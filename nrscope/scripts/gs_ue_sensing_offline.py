import sys, struct, numpy as np, collections as C
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
F=sys.argv[1]; OUT=sys.argv[2]; P=20480
raw=open(F,"rb").read(); i=0; recs=[]; hist=[]
while i+20<=len(raw):
    assert raw[i:i+4]==b"GSD1"
    sfn,sl,ns,rn,npb,occ,coh=struct.unpack_from("<HBBHHff",raw,i+4); i+=20
    if i+npb+ns*(1+3*npb*8)>len(raw): break  # the last record, cut at exit
    m=np.frombuffer(raw,np.uint8,npb,i).astype(bool); i+=npb
    syms={}
    for _ in range(ns):
        l=raw[i]; i+=1; syms[l]=np.frombuffer(raw,np.complex64,3*npb,i).copy(); i+=3*npb*8
    a=sfn*20+sl
    if hist: a+=round((float(np.median(hist[-64:]))-a)/P)*P
    hist.append(a); recs.append((a,m,occ,coh,syms))
N=len(recs[0][1])
print("grants",len(recs),"PRBs",N)
recs=[r for r in recs if r[2]<4 and r[3]>=0.6]  # port 1000 present, DM-RS clear
cov=np.sum([r[1] for r in recs],axis=0)
# the widest band (>= 160 PRBs) covered by the most grants
best=(0,0,0)
for b0 in range(0,N-160,8):
    for b1 in range(b0+160,N+1,8):
        n=sum(1 for r in recs if r[1][b0:b1].all())
        if n*(b1-b0) > best[0]*(best[2]-best[1]) if best[0] else True: best=(n,b0,b1)
n,b0,b1=best
print("band: PRBs %d-%d (%d PRBs, %.1f MHz), %d grants cover it fully"%(b0,b1-1,b1-b0,(b1-b0)*0.36,n))
T=[];X=[]
for a,m,occ,coh,syms in recs:
    if not m[b0:b1].all(): continue
    for l,v in syms.items():
        T.append((a+l/14.0)*0.5e-3); X.append(v[3*b0:3*b1])
o=np.argsort(T); T=np.array(T)[o]; X=np.array(X)[o]
K,L=X.shape; NF=1024
print("snapshots %d over %.1f s (%.0f/s); %d lattice points per snapshot (120 kHz apart)"%(K,T[-1]-T[0],K/(T[-1]-T[0]),L))
win=np.hanning(L)
def prof(x): return np.fft.ifft(x*win,NF)
crb=np.arange(b0,b1)
def stitch(x):
    """divide out one complex gain per 4-PRB group (CRB % 4 == 0 starts a group): each side's
    line, on the snapshot's own phase slope, extrapolated to the boundary"""
    y=x.copy()
    rho=np.sum(x[1:]*np.conj(x[:-1])); ph=np.angle(rho); q=np.arange(L); d=x*np.exp(-1j*ph*q)
    starts=[3*(c-b0) for c in crb if c%4==0 and c>b0]
    edges=[0]+starts+[L]
    g=np.ones(len(edges)-1,complex)
    for j in range(1,len(edges)-1):
        a0,a1,a2=edges[j-1],edges[j],edges[j+1]
        xl=np.arange(a0,a1); xr=np.arange(a1,a2); x0=a1-0.5
        cl=np.polyfit(xl-x0,d[a0:a1].real,1)+1j*np.polyfit(xl-x0,d[a0:a1].imag,1)
        cr=np.polyfit(xr-x0,d[a1:a2].real,1)+1j*np.polyfit(xr-x0,d[a1:a2].imag,1)
        step=cr[1]/cl[1] if abs(cl[1])>0 else 1
        g[j]=g[j-1]*step
    for j in range(len(edges)-1): y[edges[j]:edges[j+1]]/=g[j]
    return y
R0=np.array([prof(x) for x in X])
los0=int(np.argmax(np.mean(np.abs(R0)**2,axis=0)))
SH=20-los0   # roll every profile so the direct path sits at bin 20
los=20
_prof=prof
def prof(x): return np.roll(_prof(x),SH)
def variant(v):
    Y=X.copy()
    if v>=2: Y=np.array([stitch(x) for x in Y])
    R=np.array([prof(x) for x in Y])
    if v>=1: R=R/R[:,los:los+1]*np.abs(np.mean(R[:,los]))
    return R
def clutter(R,v,T):
    if v>=3:  # per-epoch: subtract a 20 ms moving mean in slow time
        out=np.empty_like(R); j0=0
        for k in range(len(T)):
            sel=np.abs(T-T[k])<0.010; out[k]=R[k]-R[sel].mean(axis=0)
        return out
    return R-R.mean(axis=0)
# metrics over 150 ms windows
fD=np.arange(-400,401,4.0)
names=["V0 raw","V1 + LoS norm","V2 + PRG stitching","V3 + 20 ms clutter removal"]
res=[]; maps={}
bin_m=3e8/(NF*120e3)
ghost=int(round(1/1.44e6/(1/(NF*120e3))))
for v in range(4):
    R=variant(v); Pm=np.mean(np.abs(R)**2,axis=0); Pm/=Pm.max()
    g=10*np.log10(np.max(Pm[(los+ghost-3)%NF:(los+ghost+4)%NF]) if (los+ghost+4)<NF else np.max(np.r_[Pm[los+ghost-3:],Pm[:3]]))
    # phase coherence of the direct path across consecutive snapshots
    ph=np.abs(np.mean(R[1:,los]*np.conj(R[:-1,los])/(np.abs(R[1:,los])*np.abs(R[:-1,los])+1e-20)))
    leak=[];floor=[]
    t0=T[0]; w=0.15
    while t0+w<T[-1]:
        sel=(T>=t0)&(T<t0+w)
        if sel.sum()>=60:
            Rw=clutter(R[sel],v,T[sel]); tw=T[sel]
            E=np.exp(-2j*np.pi*np.outer(fD,tw-tw.mean()))
            M=np.abs(E@Rw[:,:300])**2
            # static residue: energy left in the direct path's bins outside +-8 Hz, against its power before
            # clutter removal (coherent gain sel.sum()^2 for a line); and how far moving content stands above the floor
            stat=M[:,los-3:los+4]; zero=np.abs(fD)<=8
            leak.append(10*np.log10(stat[~zero].mean()/(np.mean(np.abs(R[sel][:,los])**2)*sel.sum()**2+1e-30)))
            far=M[:,los+8:300]; floor.append(10*np.log10(np.percentile(far,99.9)/np.median(far)))
            if t0 < T[0]+5 <= t0+w+0.15: maps[v]=(M,tw)
        t0+=w
    res.append((names[v],g,ph,np.nanmedian(leak),np.nanmedian(floor)))
print("bin %.2f m; PRG ghost at %d bins (%.0f m) from the direct path (bin %d)"%(bin_m,ghost,ghost*bin_m,los))
print("| variant | ghost at +208 m (dB vs direct path) | direct-path phase coherence, consecutive snapshots | static residue outside +-8 Hz (dB vs direct path) | brightest moving cells over the floor (dB) |")
print("|---|---|---|---|---|")
for n_,g,ph,lk,fl in res: print("| %s | %.1f | %.3f | %.1f | %.1f |"%(n_,g,ph,lk,fl))
fig,ax=plt.subplots(2,4,figsize=(20,8))
for v in range(4):
    R=variant(v); Pm=np.mean(np.abs(R)**2,axis=0); Pm/=Pm.max()
    rr=(np.arange(NF)-los)*bin_m
    ax[0,v].plot(rr[:300],10*np.log10(Pm[:300]+1e-12)); ax[0,v].set_ylim(-60,2); ax[0,v].set_xlim(-50,500)
    ax[0,v].axvline(ghost*bin_m,color="r",ls=":"); ax[0,v].set_title(names[v]+": mean range profile"); ax[0,v].set_xlabel("excess path (m)")
    if v in maps:
        M,tw=maps[v]; M=10*np.log10(M/M.max()+1e-12)
        ax[1,v].imshow(M.T[:200],aspect="auto",origin="lower",extent=[fD[0],fD[-1],(0-los)*bin_m,(200-los)*bin_m],vmin=-50,vmax=0,cmap="viridis")
        ax[1,v].set_xlabel("Doppler (Hz)"); ax[1,v].set_ylabel("excess path (m)"); ax[1,v].set_title("range-Doppler, 150 ms at t=5 s")
plt.tight_layout(); plt.savefig(OUT,dpi=90)
