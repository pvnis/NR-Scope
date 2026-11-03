import sys, struct, numpy as np, collections as C
F=sys.argv[1]; TAB=int(sys.argv[2],16); P=20480; N=273; L=3*N; NF=1024
raw=open(F,"rb").read(); i=0; hist=[]; snaps=C.defaultdict(list)
while i+20<=len(raw):
    sfn,sl,ns,rn,npb,occ,coh=struct.unpack_from("<HBBHHff",raw,i+4); i+=20
    if i+npb+ns*(1+3*npb*8)>len(raw): break
    m=np.frombuffer(raw,np.uint8,npb,i).astype(bool); i+=npb
    a=sfn*20+sl
    if hist: a+=round((float(np.median(hist[-64:]))-a)/P)*P
    hist.append(a)
    for _ in range(ns):
        l=raw[i]; i+=1; v=np.frombuffer(raw,np.complex64,3*npb,i).copy(); i+=3*npb*8
        if occ<4 and coh>=0.6: snaps[(a,l)].append((rn,m,v))
def segments(m):
    out=[];p=0
    while p<N:
        if m[p]:
            q=p
            while q+1<N and m[q+1]: q+=1
            out.append((p,q)); p=q+1
        else: p+=1
    return out
def linefit(d,idx,x0):
    x=idx-x0; return np.polyfit(x,d.real,1)[1]+1j*np.polyfit(x,d.imag,1)[1]
def step_at(x,a0,a1,a2,ph):
    """complex step between [a0,a1) and [a1,a2) on lattice samples x, slope ph removed"""
    q=np.arange(a0,a2); d=x[a0:a2]*np.exp(-1j*ph*q); n0=a1-0.5
    uL=linefit(d[:a1-a0],np.arange(a0,a1),n0); uR=linefit(d[a1-a0:],np.arange(a1,a2),n0)
    st=uR/uL if abs(uL)>0 else 1
    return st if np.isfinite(st) and 0.1<abs(st)<10 else 1  # a fit that failed: no step
def slope(x,segs):
    acc=0
    for s,e in segs:
        y=x[3*s:3*(e+1)]; acc+=np.sum(y[1:]*np.conj(y[:-1]))
    return np.angle(acc)
def stitch_seg(x,s,e,ph):
    """divide out one gain per PRG (CRB%4==0 boundaries) inside PRBs s..e of lattice x"""
    b=[3*c for c in range(s+1,e+1) if c%4==0]
    edges=[3*s]+b+[3*(e+1)]; g=1
    for j in range(1,len(edges)-1):
        if edges[j]-edges[j-1]<3 or edges[j+1]-edges[j]<3: continue
        g*=step_at(x,edges[j-1],edges[j],edges[j+1],ph); x[edges[j]:edges[j+1]]/=g
def profile(x,m):
    if not np.all(np.isfinite(x)): return None
    i=np.nonzero(np.repeat(m,3))[0]; a,b=i[0],i[-1]+1
    y=x[a:b]*np.hanning(b-a); r=np.abs(np.fft.ifft(y,NF))**2
    pk=int(np.argmax(r)); return np.roll(r,NF//2-pk)/r[pk]
def metrics(Rs):
    Pm=np.mean(Rs,axis=0); Pm/=Pm.max(); c=NF//2
    above=np.nonzero(Pm[c-40:c+41]>=0.5)[0]; w3=above[-1]-above[0]+1
    ghost=10*np.log10(Pm[c+82:c+89].max()); floor=10*np.log10(np.median(Pm[c+20:c+150]))
    return w3,ghost,floor
def step_ok(x,a0,a1,a2,ph):
    q=np.arange(a0,a2); d=x[a0:a2]*np.exp(-1j*ph*q); n0=a1-0.5
    uL=linefit(d[:a1-a0],np.arange(a0,a1),n0); uR=linefit(d[a1-a0:],np.arange(a1,a2),n0)
    st=uR/uL if abs(uL)>0 else np.nan
    return st if np.isfinite(st) else None
res=C.defaultdict(list); junc_inter=[]; junc_prg=[]; nmerged=0; nmerged_notab=0
T=sorted(snaps)
for key in T:
    gr=snaps[key]
    # single-UE reference: the tablet alone, >= 200 contiguous PRBs
    for rn,m,v in gr:
        if rn==TAB:
            for s,e in segments(m):
                if e-s+1>=200:
                    x=np.zeros(L,complex); mm=np.zeros(N,bool); mm[s:e+1]=True; x[3*s:3*(e+1)]=v[3*s:3*(e+1)]
                    ph=slope(x,[(s,e)]); stitch_seg(x,s,e,ph); res["R tablet alone, stitched"].append(profile(x,mm))
    rns=set(rn for rn,_,_ in gr)
    if len(rns)<2: continue
    X=np.zeros(L,complex); M=np.zeros(N,bool); owner=np.full(N,-1); segs=[]
    for rn,m,v in gr:
        new=m&~M
        for s,e in segments(new):
            X[3*s:3*(e+1)]=v[3*s:3*(e+1)]; owner[s:e+1]=rn; segs.append((s,e,rn))
        M|=new
    idx=np.nonzero(M)[0]
    if len(idx)==0 or idx[-1]-idx[0]+1<200 or M.sum()/(idx[-1]-idx[0]+1)<0.8: continue
    segs.sort(); nmerged+=1
    if TAB not in rns: nmerged_notab+=1
    ph=slope(X,[(s,e) for s,e,_ in segs])
    # raw PRG steps inside each UE's segments (before stitching): the precoder steps
    for s,e,rn in segs:
        for c in range(s+1,e+1):
            if c%4==0 and c-4>=s and c+4<=e+1:
                st=step_ok(X,3*(c-4),3*c,3*(c+4),ph)
                if st is not None: junc_prg.append(abs(np.degrees(np.angle(st))))
    res["M0 naive merge"].append(profile(X.copy(),M))
    X1=X.copy()
    for s,e,rn in segs: stitch_seg(X1,s,e,ph)
    res["M1 each UE stitched"].append(profile(X1.copy(),M))
    # inter-UE steps at touching junctions, after each UE is stitched
    for (s1,e1,r1),(s2,e2,r2) in zip(segs,segs[1:]):
        if s2==e1+1 and r1!=r2 and e1-s1>=3 and e2-s2>=3:
            st=step_ok(X1,3*(e1-3),3*s2,3*(s2+4),ph)
            if st is not None: junc_inter.append((abs(np.degrees(np.angle(st))),abs(20*np.log10(abs(st)+1e-12))))
    # M2: chain touching segments by boundary fits, gaps by the direct-path component
    X2=X1.copy(); q=np.arange(L)
    alpha=lambda s,e: np.mean(X2[3*s:3*(e+1)]*np.exp(-1j*ph*q[3*s:3*(e+1)]))
    ref=alpha(*segs[0][:2])
    for j in range(1,len(segs)):
        s1,e1,_=segs[j-1]; s2,e2,_=segs[j]
        if s2==e1+1 and e1-s1>=3 and e2-s2>=3: g=step_at(X2,3*(e1-3),3*s2,3*(s2+4),ph)
        else: g=alpha(s2,e2)/alpha(s1,e1)
        X2[3*s2:3*(e2+1)]/=g
    res["M2 + boundary fit, gaps by direct path"].append(profile(X2.copy(),M))
    X3=X1.copy()
    for s,e,_ in segs:
        a_=np.mean(X3[3*s:3*(e+1)]*np.exp(-1j*ph*q[3*s:3*(e+1)])); X3[3*s:3*(e+1)]/=a_/abs(a_)*abs(ref)/max(abs(a_),1e-12)*abs(a_)
    res["M3 + direct-path alignment of every segment"].append(profile(X3.copy(),M))
secs=(T[-1][0]-T[0][0])*0.0005
print("snapshots merged from >= 2 UEs (span >= 200 PRBs, fill >= 80%%): %d (%.0f/s); without the tablet: %d (%.0f/s)"%(nmerged,nmerged/secs,nmerged_notab,nmerged_notab/secs))
ji=np.array(junc_inter); jp=np.array(junc_prg)
print("phase step |deg| at junctions between two UEs (each stitched): n %d, median %.0f, p90 %.0f; amplitude step median %.1f dB"%(len(ji),np.median(ji[:,0]),np.percentile(ji[:,0],90),np.median(ji[:,1])) if len(ji) else "no touching junctions")
print("phase step |deg| at 4-PRB boundaries inside one UE, before stitching: n %d, median %.0f, p90 %.0f"%(len(jp),np.median(jp),np.percentile(jp,90)))
print("| variant | snapshots | -3 dB width of the mean direct path (bins of 2.44 m) | ghost at +208 m | sidelobe floor 50-370 m |")
print("|---|---|---|---|---|")
for k,v in res.items():
    v=[x for x in v if x is not None]; w3,g,fl=metrics(np.array(v))
    print("| %s | %d | %d | %.1f dB | %.1f dB |"%(k,len(v),w3,g,fl))
