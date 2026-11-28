"""All multi-UE tests on a 4-antenna burst recording, burst by burst (memory stays small), aggregated at the end.
Reference UE (A) = the UE with the most DM-RS snapshots in the burst (the downloading tablet); every other UE
with >= MINB snapshots is a second UE (B).
 T0 census           UEs, snapshots, PRBs per grant, precoder type (PRG jump size inside its own grants)
 T1 same-PRG         B vs A on one PRG, nearest A snapshot within 2 ms: share not explained by one number
 T2 bridge distance  B's number(s) from an A snapshot ~D away, checked against an A snapshot within 1 ms
                     variants W (B raw, one number), S (B self-stitched per snapshot, one number), P (one number per PRG)
 T3 edge, controlled B's own grant cut in two, known number, gap g PRGs
 T4 edge, real       A and B side by side in a slot: edge vs anchor vs combined vs oracle against a full-band A reference
 T5 coverage         PRGs per slot / per 150 ms window: A alone, A + mergeable B, A + all B
usage: gs_all.py <rec dir> <out prefix>"""
import sys, struct, json, glob, os, pickle, numpy as np, collections as C
R=sys.argv[1]; OUT=sys.argv[2]; N=272; NG=N//4; L=3*N; CH=[int(c) for c in os.environ.get("CH","0,1,2,3").split(",")]; NA=len(CH); MINB=int(os.environ.get("MINB","10"))
q=np.arange(L); rng=np.random.default_rng(7); NF=1024; DF=120e3; BIN=3e8/(NF*DF)
def load(path):
    raw=open(path,"rb").read(); i=0; out={}; base=None
    while i+12<=len(raw):
        tag=raw[i:i+4]
        if tag in (b"GSD1",b"GSD2"):
            two=tag==b"GSD2"
            if i+20>len(raw): break
            sfn,sl,ns,rn,npb,occ,coh=struct.unpack_from("<HBBHHff",raw,i+4); j=i+20; per=3*npb*8*(2 if two else 1)
            if j+npb+ns*(1+per)>len(raw): break
            m=np.frombuffer(raw,np.uint8,npb,j).astype(bool); j+=npb; a=sfn*20+sl
            if base is None: base=a
            a=a+20480 if a<base-1000 else a
            for _ in range(ns):
                l=raw[j]; j+=1
                if npb>=N and occ<4 and coh>=0.5:
                    out[(a,rn,l)]=(m[:N].reshape(NG,4).all(1),np.frombuffer(raw,np.complex64,3*npb,j)[:L])
                j+=per
            i=j
        elif tag==b"GSC1":
            npb=struct.unpack_from("<H",raw,i+10)[0]; i+=12+12*npb*8
        else: i+=1
    return out
def ends(d,a0,a1,x0):
    x=np.arange(a0,a1)-x0; return np.polyfit(x,d[a0:a1].real,1)[1]+1j*np.polyfit(x,d[a0:a1].imag,1)[1]
def deslope(Y,pm):
    m=np.repeat(pm,12); pr=m[1:]&m[:-1]; out=np.empty_like(Y)
    for r in range(Y.shape[0]): out[r]=Y[r]*np.exp(-1j*np.angle(np.sum((Y[r,1:]*np.conj(Y[r,:-1]))[pr]))*q)
    return out
def joint_gains(Tm,pm,steps=None):
    d=deslope(Tm,pm); g=np.ones(NG,complex)
    for b in range(1,NG):
        if not (pm[b-1] and pm[b]): g[b]=g[b-1]; continue
        Lv=np.array([ends(d[r],12*(b-1),12*b,12*b-0.5) for r in range(NA)]); Rv=np.array([ends(d[r],12*b,12*(b+1),12*b-0.5) for r in range(NA)])
        st=np.vdot(Lv,Rv)/np.vdot(Lv,Lv) if np.vdot(Lv,Lv)!=0 else 1
        ok=np.isfinite(st) and 0.1<abs(st)<10
        if steps is not None and ok: steps.append(st)
        g[b]=g[b-1]*(st if ok else 1)
    return g
SL=np.linspace(-0.08,0.08,81); gc=12*np.arange(NG)+5.5; ESL=np.exp(-1j*np.outer(SL,gc))
def epochs(A):
    EP=[]; TPL={}; e=-1
    for i,(t,a,pm,Y) in enumerate(A):
        new=True
        if i and t-A[i-1][0]<=0.02:
            sh=(CN>0)&pm
            if sh.sum()>=6:
                r=np.sum((np.conj(S_/np.maximum(np.repeat(CN,12),1)[None,:])*Y).reshape(NA,NG,12),axis=(0,2))[sh]
                c=np.abs(ESL[:,sh]@r)/(np.abs(r).sum()+1e-30); b=int(np.argmax(c))
                if c[b]>=0.9: new=False; s_=SL[b]; p0=np.angle(ESL[b,sh]@r)
        m12=np.repeat(pm,12)[None,:]
        if new: e+=1; S_=Y*m12; CN=pm.astype(int)
        else: S_=S_+Y*np.exp(-1j*(s_*q+p0))[None,:]*m12; CN=CN+pm
        EP.append(e); TPL[e]=(S_,CN)
    return EP,TPL
def frame(Xr,pr,X,px):
    m=np.repeat(pr&px,12); den=np.vdot(Xr[:,m],Xr[:,m]).real
    return (np.vdot(Xr[:,m],X[:,m])/den if den>0 else np.nan), m.sum()//12
def nmse(P,T): return float(np.sum(np.abs(P-T)**2)/np.sum(np.abs(T)**2))
def runs(pm):
    out=[];s=None
    for b in range(NG+1):
        on=b<NG and pm[b]
        if on and s is None: s=b
        if not on and s is not None: out.append((s,b)); s=None
    return out
def edge_number(XL,mL,XR,mR,K=24,ants=range(NA)):
    ants=list(ants); Z=np.zeros_like(XL); Z[:,mL]=XL[:,mL]; Z[:,mR]=XR[:,mR]
    iL=np.flatnonzero(mL); iR=np.flatnonzero(mR); l=iL[-K:]; r=iR[:K]; x0=(l[-1]+r[0])/2
    pr=(mL[1:]&mL[:-1])|(mR[1:]&mR[:-1]); Lv=[];Rv=[]
    for a in ants:
        s=np.angle(np.sum((Z[a,1:]*np.conj(Z[a,:-1]))[pr])); d=Z[a]*np.exp(-1j*s*q)
        for idx,out in ((l,Lv),(r,Rv)):
            x=idx-x0; out.append(np.polyfit(x,d[idx].real,1)[1]+1j*np.polyfit(x,d[idx].imag,1)[1])
    Lv=np.array(Lv); Rv=np.array(Rv); den=np.vdot(Lv,Lv)
    return np.vdot(Lv,Rv)/den if den!=0 else np.nan
def comb(f1,f2): return np.sqrt(abs(f1)*abs(f2))*np.exp(1j*np.angle(f1/abs(f1)+f2/abs(f2)))
def width(X):
    w=np.hanning(L); p=np.sum(np.abs(np.fft.ifft(X*w[None,:],NF,axis=1))**2,0); k=int(np.argmax(p)); r=np.roll(p,-k)/p[k]
    a=np.flatnonzero(np.r_[r[-40:],r[:40]]>=0.5); return (a[-1]-a[0]+1)*BIN
def gprof(X,pm):
    """delay power profile over the snapshot's own PRGs (Hann on the occupied span), peak-normalised, rolled to the peak"""
    m=np.repeat(pm,12); idx=np.flatnonzero(m); w=np.zeros(L); w[idx[0]:idx[-1]+1]=np.hanning(idx[-1]-idx[0]+1)
    p=np.sum(np.abs(np.fft.ifft(X*m[None,:]*w[None,:],NF,axis=1))**2,0); k=int(np.argmax(p)); return (np.roll(p,-k)/p[k])[:200].astype(np.float32)
Rz=C.defaultdict(list)     # all results: key tuple -> list of values
CENSUS=[]
BURSTS=sorted(glob.glob(R+"/b*/"),key=lambda s:int(s.rstrip("/").split("/b")[-1]))
for bi,bd in enumerate(BURSTS):
    fs=[bd+"replay_ch%d/ue_dmrs.bin"%c for c in CH]
    if not all(os.path.exists(f) for f in fs): print("skip",bd); continue
    Dd=[load(f) for f in fs]; keys=set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
    S=C.defaultdict(list)
    for k in sorted(keys):
        pm=np.logical_and.reduce([Dd[c][k][0] for c in range(NA)]); Y=np.array([Dd[c][k][1] for c in range(NA)]).astype(complex)
        if pm.sum() and np.all(np.isfinite(Y[:,np.repeat(pm,12)])): S[k[1]].append((k[0]*0.5e-3+k[2]/28e3,k[0],pm,Y))
    del Dd
    cnt=sorted(((len(v),u) for u,v in S.items()),reverse=True)
    if not cnt: continue
    UA=cnt[0][1]; A=sorted(S[UA],key=lambda x:x[0])
    EP,TPL=epochs(A); EG={k:joint_gains(S2/np.maximum(np.repeat(C2,12),1)[None,:],C2>0) for k,(S2,C2) in TPL.items()}
    XA=[]
    for i,(t,a,pm,Y) in enumerate(A):
        X=Y/np.repeat(EG[EP[i]],12)[None,:]; X[:,~np.repeat(pm,12)]=0; XA.append((t,a,pm,X))
    TA=np.array([x[0] for x in XA]); bya={}
    for i,x in enumerate(XA): bya.setdefault(x[1],[]).append(i)
    st=[]; [joint_gains(Y,pm,st) for t,a,pm,Y in A[::20]]
    CENSUS.append(dict(burst=bi,role="A",rnti=UA,n=len(A),prb=float(np.median([4*x[2].sum() for x in A])),jump=float(np.median(np.abs(np.degrees(np.angle(st))))) if st else np.nan,
        jump_db=float(np.median(np.abs(20*np.log10(np.abs(st))))) if st else np.nan,epochs=len(TPL)))
    Rz[("T6","A","jump_deg")].extend(np.degrees(np.angle(st)).tolist()); Rz[("T6","A","jump_db")].extend((20*np.log10(np.abs(st))).tolist())
    for i in range(0,len(A),3):
        if A[i][2].sum()>=16:
            Rz[("T6","A","prof raw")].append(gprof(A[i][3],A[i][2])); Rz[("T6","A","prof epoch")].append(gprof(XA[i][3],A[i][2]))
            g=joint_gains(A[i][3],A[i][2]); X=A[i][3]/np.repeat(g,12)[None,:]; Rz[("T6","A","prof snap")].append(gprof(X,A[i][2]))
    print("burst %d: A %#x %d snapshots, %d epochs; others %s"%(bi,UA,len(A),len(TPL),[("%#x"%u,n) for n,u in cnt[1:6]]),flush=True)
    MERGE=set()     # (a, rnti) of B snapshots mergeable (anchor within 30 ms or touching A in the slot)
    for nb,UB in cnt[1:]:
        if nb<MINB: continue
        B=sorted(S[UB],key=lambda x:x[0]); st=[]
        XB=[]
        for t,a,pb,Yb in B:
            g=joint_gains(Yb,pb,st); X=Yb/np.repeat(g,12)[None,:]; X[:,~np.repeat(pb,12)]=0; XB.append(X)
        jump=float(np.median(np.abs(np.degrees(np.angle(st))))) if st else np.nan
        typ="per-PRG" if jump>20 else "wideband"
        CENSUS.append(dict(burst=bi,role="B",rnti=UB,n=nb,prb=float(np.median([4*x[2].sum() for x in B])),jump=jump,typ=typ,
            jump_db=float(np.median(np.abs(20*np.log10(np.abs(st))))) if st else np.nan))
        Rz[("T6",typ,"jump_deg")].extend(np.degrees(np.angle(st)).tolist()); Rz[("T6",typ,"jump_db")].extend((20*np.log10(np.abs(st))).tolist())
        for ib,(t,a,pb,Yb) in enumerate(B):
            if pb.sum()>=16: Rz[("T6",typ,"prof raw")].append(gprof(Yb,pb)); Rz[("T6",typ,"prof snap")].append(gprof(XB[ib],pb))
        U=(bi,UB)
        for ib,(t,a,pb,Yb) in enumerate(B):
            # ---- T1 same-PRG similarity (raw, one PRG, 4 antennas x 12)
            j=int(np.argmin(np.abs(TA-t)))
            if abs(TA[j]-t)<=0.002:
                for b in np.flatnonzero(pb&A[j][2])[:20]:
                    x=A[j][3][:,12*b:12*b+12].ravel(); y=Yb[:,12*b:12*b+12].ravel()
                    Rz[("T1",typ,"B vs A")].append(1-abs(np.vdot(x,y))**2/(np.vdot(x,x).real*np.vdot(y,y).real))
            # ---- mergeability: anchor within 30 ms or touching in the slot
            near=[k for k in np.flatnonzero(np.abs(TA-t)<=0.03) if (XA[k][2]&pb).sum()>=4]
            touch=a in bya and any(abs(r1[0]-r2[1])==0 or abs(r2[0]-r1[1])==0 for r1 in runs(XA[bya[a][0]][2]) for r2 in runs(pb))
            if near or touch: MERGE.add((a,UB))
            # ---- T2 bridge distance
            close=[k for k in np.flatnonzero(np.abs(TA-t)<=0.001) if (XA[k][2]&pb).sum()>=4]
            for D in (0.002,0.01,0.03,0.1,0.3,1.0):
                far=[k for k in np.flatnonzero((np.abs(TA-t)>=0.8*D)&(np.abs(TA-t)<=1.25*D)) if (XA[k][2]&pb).sum()>=4]
                if not close or not far: continue
                jc=max(close,key=lambda k:(XA[k][2]&pb).sum()); jf=max(far,key=lambda k:(XA[k][2]&pb).sum())
                r,n_=frame(XA[jf][3],XA[jf][2],XA[jc][3],XA[jc][2])
                if not np.isfinite(r) or n_<8: continue
                m=np.repeat(pb&XA[jc][2]&XA[jf][2],12)
                if m.sum()<48: continue
                Tr=XA[jc][3][:,m]
                Rz[("T2",typ,D,"A own change")].append(nmse(XA[jf][3][:,m]*r,Tr))
                for var,Z in (("W",Yb),("S",XB[ib])):
                    fc,_=frame(XA[jc][3],XA[jc][2],Z,pb); ff,_=frame(XA[jf][3],XA[jf][2],Z,pb)
                    if not (np.isfinite(fc) and np.isfinite(ff)): continue
                    est=ff/r
                    Rz[("T2",typ,D,var,"phase")].append(abs(np.degrees(np.angle(est/fc)))); Rz[("T2",typ,D,var,"amp")].append(abs(20*np.log10(abs(est/fc))))
                    Rz[("T2",typ,D,var,"err")].append(nmse(Z[:,m]/est,Tr))
                    if D==0.002:
                        o=np.vdot(Z[:,m],Tr)/np.vdot(Z[:,m],Z[:,m]); Rz[("T2",typ,"oracle",var)].append(nmse(Z[:,m]*o,Tr))
                # P: one number per PRG from the far A snapshot
                P=np.zeros_like(Tr); Po=np.zeros_like(Tr); ok=True
                for ii,b in enumerate(np.flatnonzero(pb&XA[jc][2]&XA[jf][2])):
                    sl_=slice(12*b,12*b+12); xf=XA[jf][3][:,sl_]; y=Yb[:,sl_]; tr=XA[jc][3][:,sl_]
                    fb=np.vdot(xf,y)/np.vdot(xf,xf); P[:,12*ii:12*ii+12]=y/fb*r
                    Po[:,12*ii:12*ii+12]=y*np.vdot(y,tr)/np.vdot(y,y)
                Rz[("T2",typ,D,"P","err")].append(nmse(P,Tr))
                if D==0.002: Rz[("T2",typ,"oracle","P")].append(nmse(Po,Tr))
            # ---- T3 edge, controlled on B's own grant (self-stitched S and raw W), known number
            if ib%2==0:
                for g in (0,1,2,4,8):
                    for s,e_ in runs(pb):
                        if e_-s<4+g: continue
                        mid=s+(e_-s-g)//2; mL=np.zeros(L,bool); mL[12*s:12*mid]=True; mR=np.zeros(L,bool); mR[12*(mid+g):12*e_]=True
                        cc=rng.uniform(0.5,2)*np.exp(1j*rng.uniform(-np.pi,np.pi))
                        for var,Z in (("W",Yb),("S",XB[ib])):
                            f=edge_number(Z,mL,Z*cc,mR,24)
                            if np.isfinite(f): Rz[("T3",typ,g,var,"phase")].append(abs(np.degrees(np.angle(f/cc)))); Rz[("T3",typ,g,var,"amp")].append(abs(20*np.log10(abs(f/cc))))
                        break
            # ---- T4 edge, real side by side
            if a in bya:
                t0,_,pa,X0=XA[bya[a][0]]
                if (pa&pb).any() or pa.sum()<4: continue
                best=None
                for sa,ea in runs(pa):
                    for sb,eb in runs(pb):
                        cand=(sb-ea,"AB",(sa,ea),(sb,eb)) if ea<=sb else (sa-eb,"BA",(sa,ea),(sb,eb))
                        if best is None or cand[0]<best[0]: best=cand
                gap,order,(sa,ea),(sb,eb)=best
                mA=np.zeros(L,bool); mA[12*sa:12*ea]=True; mB=np.zeros(L,bool); mB[12*sb:12*eb]=True; K=min(24,12*(ea-sa),12*(eb-sb))
                full=[k for k in np.flatnonzero(np.abs(TA-t)<=0.005) if XA[k][2][(pa|pb)].all() and XA[k][1]!=a]
                ref=None
                if full:
                    k=min(full,key=lambda k:abs(TA[k]-t)); rk,_=frame(XA[k][3],XA[k][2],X0,pa); ref=XA[k][3]*rk
                for var,Z in (("W",Yb),("S",XB[ib])):
                    fe=edge_number(X0,mA,Z,mB,K) if order=="AB" else 1/edge_number(Z,mB,X0,mA,K)
                    nr=[k for k in np.flatnonzero(np.abs(TA-t)<=0.03) if XA[k][1]!=a]
                    cs=[(k,)+frame(XA[k][3],XA[k][2],Z,pb) for k in nr]; cs=[x for x in cs if x[2]>=4 and np.isfinite(x[1])]
                    fa_=np.nan; dta=np.nan; fitres=np.nan
                    if cs:
                        k2,f,n=min(cs,key=lambda x:(abs(TA[x[0]]-t)>0.002,-x[2])); r0,n0=frame(XA[k2][3],XA[k2][2],X0,pa)
                        if np.isfinite(r0) and n0>=4:
                            fa_=f/r0; dta=abs(TA[k2]-t); mm=np.repeat(XA[k2][2]&pb,12)
                            fitres=nmse(Z[:,mm]/f,XA[k2][3][:,mm])
                    row=dict(burst=bi,rnti=UB,t=t,gap=4*gap,runA=4*(ea-sa),runB=4*(eb-sb),var=var,typ=typ,dta=dta,fitres=fitres)
                    ests={"edge":fe,"anchor":fa_,"combined":comb(fe,fa_) if np.isfinite(fe) and np.isfinite(fa_) else (fe if np.isfinite(fe) else fa_)}
                    if np.isfinite(fe) and np.isfinite(fa_): row["edge_vs_anchor_deg"]=float(np.degrees(np.angle(fe/fa_))); row["edge_vs_anchor_db"]=float(20*np.log10(abs(fe/fa_)))
                    if ref is not None:
                        mb=np.repeat(pb,12); mab=np.repeat(pa|pb,12); Tr=ref[:,mb]
                        o=np.vdot(Z[:,mb],Tr)/np.vdot(Z[:,mb],Z[:,mb]); row["err_oracle"]=nmse(Z[:,mb]*o,Tr)
                        row["w_ref"]=width(ref*mab[None,:]); row["w_A"]=width(X0)
                        for nm,f in ests.items():
                            if np.isfinite(f):
                                row["err_"+nm]=nmse(Z[:,mb]/f,Tr); Xm=X0.copy(); Xm[:,mb]=Z[:,mb]/f; row["w_"+nm]=width(Xm)
                        Xr=X0.copy(); Xr[:,mb]=Z[:,mb]*np.sqrt(np.sum(np.abs(X0)**2)/max(np.sum(np.abs(Z[:,mb])**2),1e-30)*pb.sum()/pa.sum()); row["err_raw"]=nmse(Xr[:,mb],Tr)
                    Rz["T4"].append(row)
    # ---- T5 coverage
    slots=sorted({x[1] for u in S for x in S[u]})
    for sel in ("A alone","A + mergeable B","A + all B"):
        per=[]; win=C.defaultdict(lambda: np.zeros(NG,bool))
        for a in slots:
            m=np.zeros(NG,bool)
            for i in bya.get(a,[]): m|=XA[i][2]
            if sel!="A alone":
                for u in S:
                    if u==UA: continue
                    for t,a2,pb,Yb in S[u]:
                        if a2==a and (sel=="A + all B" or (a,u) in MERGE): m|=pb
            per.append(m.sum()); win[int(a*0.5e-3/0.15)]|=m
        per=np.array(per)
        Rz[("T5",sel,"PRBs per slot")].extend((4*per).tolist())
        Rz[("T5",sel,"window cover")].extend([w.mean() for w in win.values()])
pickle.dump(dict(Rz=dict(Rz),CENSUS=CENSUS),open(OUT+".pkl","wb"))
print("saved",OUT+".pkl")
