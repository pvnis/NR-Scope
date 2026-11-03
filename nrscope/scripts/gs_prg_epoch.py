"""Does a UE's per-PRG precoder stay the same over successive slots?
usage: prg_epoch.py ue_dmrs.bin <rnti hex>
For two grants of the UE dt slots apart, the per-PRG ratio <h2,h1> over the shared PRGs: with the same precoders
its phase is one common value (CFO, timing drift as a linear term); new precoders scatter it uniformly.
The CSI-RS port 0 (a fixed beam) gives the channel's own change over dt."""
import sys, struct, numpy as np, collections as C
F=sys.argv[1]; TAB=int(sys.argv[2],16); P=20480; N=273
raw=open(F,"rb").read(); i=0; hist=[]; G=[]; S=[]
def unwrap(sfn,sl):
    a=sfn*20+sl
    if hist: a+=round((float(np.median(hist[-64:]))-a)/P)*P
    hist.append(a); return a
while i+8<=len(raw):
    t=raw[i:i+4]
    if t==b"GSD1":
        if i+20>len(raw): break
        sfn,sl,ns,rn,npb,occ,coh=struct.unpack_from("<HBBHHff",raw,i+4); j=i+20
        if j+npb+ns*(1+3*npb*8)>len(raw): break
        m=np.frombuffer(raw,np.uint8,npb,j).astype(bool); j+=npb; a=unwrap(sfn,sl); h=None
        for _ in range(ns):
            l=raw[j]; j+=1
            if l==2: h=np.frombuffer(raw,np.complex64,3*npb,j).astype(complex)
            j+=3*npb*8
        i=j
        if rn==TAB and h is not None and occ<4 and coh>=0.6: G.append((a,m,h))
    elif t==b"GSC1":
        sfn,sl,l,fc,npb=struct.unpack_from("<HBBHH",raw,i+4); j=i+12
        if j+12*npb*8>len(raw): break
        y=np.frombuffer(raw,np.complex64,12*npb,j).astype(complex).reshape(npb,12); i=j+12*npb*8
        S.append((unwrap(sfn,sl),y[:,0]+y[:,1]))   # port 3000 + scrambling ~ fixed per slot: use |ratio| only
    else: i+=1
G.sort(key=lambda g:g[0])
def prg_vec(m,h):
    out={}
    for g in range(0,N,4):
        if m[g:g+4].all(): out[g//4]=h[3*g:3*(g+4)]
    return out
V=[(a,prg_vec(m,h)) for a,m,h in G]
def spread(v1,v2):
    ks=sorted(set(v1)&set(v2))
    if len(ks)<8: return None
    r=np.array([np.vdot(v1[k],v2[k]) for k in ks]); k=np.array(ks,float)
    # remove the common phase and a linear term (timing drift) by a coarse search over the slope
    best=None
    for s in np.linspace(-0.6,0.6,121):
        z=r*np.exp(-1j*s*k); c=np.abs(np.sum(z/np.abs(z)))/len(z)
        if best is None or c>best: best=c
    return best   # mean resultant length: 1 = one phase on every PRG, ~1/sqrt(n) = random
by=C.defaultdict(list)
for x in range(len(V)):
    for y in range(x+1,min(len(V),x+60)):
        dt=V[y][0]-V[x][0]
        if dt>200: break
        b=[1,2,3,5,10,20,40,80,160]; k=min(b,key=lambda q:abs(q-dt))
        if abs(k-dt)>max(1,k*0.25): continue
        s=spread(V[x][1],V[y][1])
        if s is not None: by[k].append(s)
print("tablet grants %d"%len(V))
print("| dt (slots) | dt (ms) | pairs | per-PRG phase coherence median | pairs with coherence > 0.9 |\n|---|---|---|---|---|")
for k in sorted(by):
    v=np.array(by[k]); print("| %d | %.1f | %d | %.2f | %.0f%% |"%(k,k/2,len(v),np.median(v),100*np.mean(v>0.9)))
# epoch lengths: run of successive grants (<= 3 slots apart) whose coherence to the run start stays > 0.9
runs=[]; s0=0
for x in range(1,len(V)):
    s=spread(V[s0][1],V[x][1]) if V[x][0]-V[x-1][0]<=3 else None
    if s is None or s<0.9:
        runs.append(V[x-1][0]-V[s0][0]); s0=x
runs=np.array(runs)
print("precoder epochs (successive grants coherent with the first): n %d, median %.1f ms, p90 %.1f ms, longest %.1f ms"%(
    len(runs),np.median(runs)/2,np.percentile(runs,90)/2,runs.max()/2))
