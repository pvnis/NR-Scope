"""4-antenna analysis of PDSCH DM-RS: how many basis beams a UE's precoders span (spatial-frequency rank per PRG over
snapshots), and whether another UE's channel lies in the same beam space (the condition for merging UEs coherently).
Input: <rec>/replay_ch{0..3}/ue_dmrs.bin, GSD2 records (ports 1000 and 1001) of every UE, one replay per antenna.
usage: gs_4ch.py <rec dir> <out prefix>"""
import sys, struct, json, numpy as np, collections as C
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=sys.argv[1]; OUT=sys.argv[2]; N=272; NG=N//4; L=3*N; NA=4
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
                    p0=np.frombuffer(raw,np.complex64,3*npb,j)[:L].astype(complex)
                    p1=np.frombuffer(raw,np.complex64,3*npb,j+3*npb*8)[:L].astype(complex) if two else None
                    out[(a,rn,l)]=(m[:N].reshape(NG,4).all(1),p0,p1)
                j+=per
            i=j
        elif tag==b"GSC1":
            npb=struct.unpack_from("<H",raw,i+10)[0]; i+=12+12*npb*8
        else: i+=1
    return out
import glob, os
BURSTS=sorted(glob.glob(R+"/b*/")) or [R+"/"]
Dd=[{} for _ in range(NA)]
for bi,bd in enumerate(BURSTS):
    for c in range(NA):
        f=bd+"replay_ch%d/ue_dmrs.bin"%c
        if os.path.exists(f):
            for (a,rn,l),v in load(f).items(): Dd[c][(a+bi*1000000,rn,l)]=v
print("bursts:",len(BURSTS))
keys=set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
print("DM-RS snapshots per channel: %s; on all 4: %d"%([len(d) for d in Dd],len(keys)))
ue=C.Counter(k[1] for k in keys); print("per UE (on all 4 antennas):",[("%#x"%u,n) for u,n in ue.most_common(8)])
# snapshot -> PRG mask (all antennas), Y (antennas x L) per port
S=C.defaultdict(list)
for k in sorted(keys):
    pm=np.logical_and.reduce([Dd[c][k][0] for c in range(NA)])
    for port in (0,1):
        if any(Dd[c][k][1+port] is None for c in range(NA)): continue
        Y=np.array([Dd[c][k][1+port] for c in range(NA)])
        if np.all(np.isfinite(Y[:,np.repeat(pm,12)])): S[k[1]].append((k[0],k[2],port,pm,Y))
# ===================== cross-slot anchoring: put a second UE's grant in the tablet's frame
UA,UB=0x7466,0x7176; NF=1024; DF=120e3; BIN=3e8/(NF*DF); q=np.arange(L); rng=np.random.default_rng(7)
def ends(d,a0,a1,x0):
    x=np.arange(a0,a1)-x0; return np.polyfit(x,d[a0:a1].real,1)[1]+1j*np.polyfit(x,d[a0:a1].imag,1)[1]
def deslope(Y,pm):
    m=np.repeat(pm,12); pr=m[1:]&m[:-1]; out=np.empty_like(Y)
    for r in range(NA): out[r]=Y[r]*np.exp(-1j*np.angle(np.sum((Y[r,1:]*np.conj(Y[r,:-1]))[pr]))*q)
    return out
def joint_gains(Tm,pm):
    d=deslope(Tm,pm); g=np.ones(NG,complex)
    for b in range(1,NG):
        if not (pm[b-1] and pm[b]): g[b]=g[b-1]; continue
        Lv=np.array([ends(d[r],12*(b-1),12*b,12*b-0.5) for r in range(NA)]); Rv=np.array([ends(d[r],12*b,12*(b+1),12*b-0.5) for r in range(NA)])
        st=np.vdot(Lv,Rv)/np.vdot(Lv,Lv) if np.vdot(Lv,Lv)!=0 else 1
        g[b]=g[b-1]*(st if np.isfinite(st) and 0.1<abs(st)<10 else 1)
    return g
# tablet snapshots (port 1000), time-ordered, epochs on 4 antennas, epoch gains
A=sorted([(a*0.5e-3+l/28e3,a,pm,Y) for a,l,p,pm,Y in S[UA] if p==0],key=lambda x:x[0])
SL=np.linspace(-0.08,0.08,81); gc=12*np.arange(NG)+5.5; ESL=np.exp(-1j*np.outer(SL,gc))
EP=[]; TPL={}; e=-1
for i,(t,a,pm,Y) in enumerate(A):
    new=True
    if i and t-A[i-1][0]<=0.02 and abs(a-A[i-1][1])<1000:
        sh=(CN>0)&pm
        if sh.sum()>=6:
            r=np.sum((np.conj(S_/np.maximum(np.repeat(CN,12),1)[None,:])*Y).reshape(NA,NG,12),axis=(0,2))[sh]
            c=np.abs(ESL[:,sh]@r)/(np.abs(r).sum()+1e-30); b=int(np.argmax(c))
            if c[b]>=0.9: new=False; s_=SL[b]; p0=np.angle(ESL[b,sh]@r)
    m12=np.repeat(pm,12)[None,:]
    if new: e+=1; S_=Y*m12; CN=pm.astype(int)
    else: S_=S_+Y*np.exp(-1j*(s_*q+p0))[None,:]*m12; CN=CN+pm
    EP.append(e); TPL[e]=(S_,CN)
EG={k:joint_gains(S2/np.maximum(np.repeat(C2,12),1)[None,:],C2>0) for k,(S2,C2) in TPL.items()}
XA=[]   # corrected tablet snapshots: coherent over the tablet's PRGs
for i,(t,a,pm,Y) in enumerate(A):
    X=Y/np.repeat(EG[EP[i]],12)[None,:]; X[:,~np.repeat(pm,12)]=0; XA.append((t,a,pm,X))
TA=np.array([x[0] for x in XA])
B=sorted([(a*0.5e-3+l/28e3,a,pm,Y) for a,l,p,pm,Y in S[UB] if p==0],key=lambda x:x[0])
print("tablet snapshots %d in %d epochs; second UE snapshots %d"%(len(A),e+1,len(B)))
def frame(Xr,pr,X,px):     # X ~ f * Xr on their common PRGs: f by least squares over 4 antennas
    m=np.repeat(pr&px,12)
    den=np.vdot(Xr[:,m],Xr[:,m]).real
    return (np.vdot(Xr[:,m],X[:,m])/den if den>0 else np.nan), m.sum()//12
def c_from(idxs,Yb,pb):    # second UE vs tablet snapshot(s): one complex number per tablet snapshot
    out=[]
    for j in idxs:
        f,n=frame(XA[j][3],XA[j][2],Yb,pb)
        if n>=4 and np.isfinite(f): out.append((j,f,n))
    return out
# ===================== edge matching between UEs: the number that puts grant R in grant L's frame, from the boundary only
def deslope_m(Y,m):
    pr=m[1:]&m[:-1]; out=np.empty_like(Y)
    for r in range(Y.shape[0]): out[r]=Y[r]*np.exp(-1j*np.angle(np.sum((Y[r,1:]*np.conj(Y[r,:-1]))[pr]))*q)
    return out
def edge_number(XL,mL,XR,mR,K=24,ants=range(NA)):
    """XL on subcarriers mL (left), XR on mR (right, all right of mL). Remove the delay slope (estimated inside
    each side, never across the boundary), fit a complex line on the last/first K subcarriers of each side,
    extrapolate both to the middle of the boundary, one number by LS over the antennas: XR ~ f*XL."""
    ants=list(ants); m=mL|mR; Z=np.zeros_like(XL); Z[:,mL]=XL[:,mL]; Z[:,mR]=XR[:,mR]
    iL=np.flatnonzero(mL); iR=np.flatnonzero(mR); l=iL[-K:]; r=iR[:K]; x0=(l[-1]+r[0])/2
    pr=np.zeros(L-1,bool); pr[:]=(mL[1:]&mL[:-1])|(mR[1:]&mR[:-1])
    Lv=[];Rv=[]
    for a in ants:
        s=np.angle(np.sum((Z[a,1:]*np.conj(Z[a,:-1]))[pr])); d=Z[a]*np.exp(-1j*s*q)
        for idx,out in ((l,Lv),(r,Rv)):
            x=idx-x0; c=np.polyfit(x,d[idx].real,1)[1]+1j*np.polyfit(x,d[idx].imag,1)[1]; out.append(c)
    Lv=np.array(Lv); Rv=np.array(Rv); den=np.vdot(Lv,Lv)
    return np.vdot(Lv,Rv)/den if den!=0 else np.nan
def runs(pm):
    out=[];s=None
    for b in range(NG+1):
        on=b<NG and pm[b]
        if on and s is None: s=b
        if not on and s is not None: out.append((s,b)); s=None
    return out
# ---- 1. controlled: split the second UE's own grant (one wideband precoder, so its pieces are coherent), known number
print("\n== 1. controlled test: second UE's grant cut in two, right piece times a known number, gap of g PRGs removed")
print("  | gap | fit length | antennas | cases | phase error median / p90 (deg) | amplitude error median (dB) |")
print("  |---|---|---|---|---|---|")
for g in (0,1,2,4,8):
    for K,ants in ((12,(0,)),(24,(0,)),(12,range(NA)),(24,range(NA)),(48,range(NA))):
        ph=[];am=[]
        for t,a,pb,Yb in B:
            for s,e_ in runs(pb):
                need=2*int(np.ceil(K/12))+g
                if e_-s<need: continue
                mid=s+(e_-s-g)//2
                mL=np.zeros(L,bool); mL[12*s:12*mid]=True; mR=np.zeros(L,bool); mR[12*(mid+g):12*e_]=True
                c=rng.uniform(0.5,2)*np.exp(1j*rng.uniform(-np.pi,np.pi))
                f=edge_number(Yb,mL,Yb*c,mR,K,ants)
                if np.isfinite(f): ph.append(abs(np.degrees(np.angle(f/c)))); am.append(abs(20*np.log10(abs(f/c))))
        if ph: print("  | %d PRBs | %d PRBs | %s | %d | %.1f / %.1f | %.2f |"%(4*g,K//12,"4" if len(list(ants))==4 else "1",len(ph),np.median(ph),np.percentile(ph,90),np.median(am)))
# ---- 2. real: tablet and second UE side by side in the same slot; edge number vs cross-slot anchoring number
print("\n== 2. real side-by-side slots: edge matching vs cross-slot anchoring (2 ms), and against a full-band tablet reference")
w=np.hanning(L); W=0.002
bya={}
for i,(t,a,pm,X) in enumerate(XA): bya.setdefault(a,[]).append(i)
rows=[]
for t,a,pb,Yb in B:
    if a not in bya: continue
    t0,_,pa,X0=XA[bya[a][0]]
    if (pa&pb).any() or pa.sum()<4: continue
    # nearest pair of runs (tablet run, UE run)
    best=None
    for sa,ea in runs(pa):
        for sb,eb in runs(pb):
            if ea<=sb: gap=sb-ea; cand=(gap,"T|U",(sa,ea),(sb,eb))
            else: gap=sa-eb; cand=(gap,"U|T",(sa,ea),(sb,eb))
            if best is None or gap<best[0]: best=cand
    gap,order,(sa,ea),(sb,eb)=best
    mA=np.zeros(L,bool); mA[12*sa:12*ea]=True; mB=np.zeros(L,bool); mB[12*sb:12*eb]=True
    K=min(24,12*(ea-sa),12*(eb-sb))
    if order=="T|U": fe=edge_number(X0,mA,Yb,mB,K)
    else: fe=1/edge_number(Yb,mB,X0,mA,K)
    # anchoring number in this slot's frame
    near=[j for j in np.flatnonzero(np.abs(TA-t)<=W) if XA[j][1]!=a]
    fa_=np.nan
    cs=[(j,)+frame(XA[j][3],XA[j][2],Yb,pb) for j in near]; cs=[x for x in cs if x[2]>=4 and np.isfinite(x[1])]
    if cs:
        j,f,n=max(cs,key=lambda x:x[2]); r0,n0=frame(XA[j][3],XA[j][2],X0,pa)
        if np.isfinite(r0) and n0>=4: fa_=f/r0
    # full-band reference
    err={}
    full=[k for k in np.flatnonzero(np.abs(TA-t)<=0.005) if XA[k][2][(pa|pb)].all() and XA[k][1]!=a]
    if full:
        k=min(full,key=lambda k:abs(TA[k]-t)); rk,_=frame(XA[k][3],XA[k][2],X0,pa); Xref=XA[k][3]*rk
        mb=np.repeat(pb,12)
        for nm,f in (("edge",fe),("anchor",fa_)):
            if np.isfinite(f): err[nm]=np.sum(np.abs(Yb[:,mb]/f-Xref[:,mb])**2)/np.sum(np.abs(Xref[:,mb])**2)
        o=np.vdot(Yb[:,mb],Xref[:,mb])/np.vdot(Yb[:,mb],Yb[:,mb]); err["oracle"]=np.sum(np.abs(Yb[:,mb]*o-Xref[:,mb])**2)/np.sum(np.abs(Xref[:,mb])**2)
    rows.append((t,gap,order,(ea-sa)*4,(eb-sb)*4,fe,fa_,err))
print("  side-by-side slots: %d; gaps (PRBs): %s"%(len(rows),C.Counter(4*r[1] for r in rows).most_common()))
print("  | time (s) | gap PRBs | order | tablet run | UE run | edge vs anchor: phase (deg), amp (dB) | UE-PRB error: edge / anchor / oracle |")
print("  |---|---|---|---|---|---|---|")
for t,gap,order,la,lb,fe,fa_,err in rows:
    d=fe/fa_ if np.isfinite(fa_) else np.nan
    print("  | %.3f | %d | %s | %d | %d | %s | %s |"%(t,4*gap,order,la,lb,
        ("%.1f, %.2f"%(np.degrees(np.angle(d)),20*np.log10(abs(d)))) if np.isfinite(d) else "no anchor",
        " / ".join(("%.0f%%"%(100*err[k])) if k in err else "-" for k in ("edge","anchor","oracle"))))
for gmax in (0,2,8,99):
    sel=[r for r in rows if r[1]<=gmax and np.isfinite(r[6])]
    if sel:
        d=np.array([r[5]/r[6] for r in sel]); E=[r[7] for r in sel if "edge" in r[7] and "anchor" in r[7]]
        print("  gap <= %d PRBs: n %d, edge vs anchor phase median %.1f deg, amp %.2f dB%s"%(4*gmax,len(sel),np.median(np.abs(np.degrees(np.angle(d)))),np.median(np.abs(20*np.log10(np.abs(d)))),
            ("; error edge %.0f%% anchor %.0f%% oracle %.0f%%"%(100*np.median([e["edge"] for e in E]),100*np.median([e["anchor"] for e in E]),100*np.median([e["oracle"] for e in E]))) if E else ""))
