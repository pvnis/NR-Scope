import sys, struct, numpy as np, collections as C
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
F=sys.argv[1]; OUT=sys.argv[2]; TAB=int(sys.argv[3],16); P=20480
raw=open(F,"rb").read(); i=0; recs=[]; hist=[]
while i+20<=len(raw):
    if raw[i:i+4]==b"GSC1":
        npb_=struct.unpack_from("<H",raw,i+10)[0]; i+=12+12*npb_*8; continue
    assert raw[i:i+4]==b"GSD1"
    sfn,sl,ns,rn,npb,occ,coh=struct.unpack_from("<HBBHHff",raw,i+4); i+=20
    if i+npb+ns*(1+3*npb*8)>len(raw): break  # the last record, cut at exit
    m=np.frombuffer(raw,np.uint8,npb,i).astype(bool); i+=npb
    syms={}
    for _ in range(ns):
        l=raw[i]; i+=1; syms[l]=np.frombuffer(raw,np.complex64,3*npb,i).copy(); i+=3*npb*8
    a=sfn*20+sl
    if hist: a+=round((float(np.median(hist[-64:]))-a)/P)*P
    hist.append(a); recs.append((a,m,occ,coh,syms,rn))
N=len(recs[0][1])
print("grants",len(recs),"PRBs",N)
recs=[r[:5] for r in recs if r[2]<4 and r[3]>=0.6 and r[5]==TAB]  # port 1000 present, DM-RS clear
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
    if not m[b0:b1].all() or rn!=TAB: continue
    for l,v in syms.items():
        T.append((a+l/14.0)*0.5e-3); X.append(v[3*b0:3*b1])
o=np.argsort(T); T=np.array(T)[o]; X=np.array(X)[o]
K,L=X.shape; NF=1024
print("snapshots %d over %.1f s (%.0f/s); %d lattice points per snapshot (120 kHz apart)"%(K,T[-1]-T[0],K/(T[-1]-T[0]),L))
win=np.hanning(L)
def prof(x): return np.fft.ifft(x*win,NF)
crb=np.arange(b0,b1)
def stitch(x,gains=False):
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
    if gains: return g,edges
    for j in range(len(edges)-1): y[edges[j]:edges[j+1]]/=g[j]
    return y
# --- epochs: successive snapshots whose per-PRG phases agree with the running epoch template
EDG=stitch(X[0],True)[1]; NG=len(EDG)-1
def prg_r(t,x): return np.array([np.vdot(t[EDG[j]:EDG[j+1]],x[EDG[j]:EDG[j+1]]) for j in range(NG)])
SL=np.linspace(-0.08,0.08,81); gidx=np.array([(EDG[j]+EDG[j+1])/2 for j in range(NG)])
def align(t,x):
    r=prg_r(t,x); best=(-1,0,0)
    for s_ in SL:
        z=r*np.exp(-1j*s_*gidx); c=np.abs(np.sum(z/(np.abs(z)+1e-30)))/NG
        if c>best[0]: best=(c,s_,np.angle(np.sum(z)))
    return best
epoch=np.zeros(len(X),int); e=0; tmpl=X[0].copy(); q=np.arange(L)
for k in range(1,len(X)):
    c,s_,p0=align(tmpl,X[k])
    if c<0.9 or T[k]-T[k-1]>0.02:
        e+=1; tmpl=X[k].copy()
    else:
        tmpl=tmpl+X[k]*np.exp(-1j*(s_*q+p0))
    epoch[k]=e
EG={}
walk=[]
for ee in range(e+1):
    idx=np.nonzero(epoch==ee)[0]; t=np.zeros(L,complex)
    for k in idx:
        if k==idx[0]: t=X[k].copy(); continue
        c,s_,p0=align(t,X[k]); t=t+X[k]*np.exp(-1j*(s_*q+p0))
    g,_=stitch(t,True); EG[ee]=g
    if len(idx)>=5:   # per-snapshot stitching against the epoch estimate: the phase walk at the band's far end
        for k in idx[::max(1,len(idx)//10)]:
            gk,_=stitch(X[k],True); d=gk/g; d=d/d[0]; walk.append(abs(np.degrees(np.angle(d[-1]))))
ln=np.bincount(epoch); dur=np.array([T[epoch==ee][-1]-T[epoch==ee][0] for ee in range(e+1)])
print("epochs %d: snapshots per epoch median %d, mean %.1f; duration median %.1f ms, p90 %.1f ms; snapshots in epochs of >= 5: %.0f%%"%(
    e+1,np.median(ln),ln.mean(),1e3*np.median(dur),1e3*np.percentile(dur,90),100*ln[ln>=5].sum()/len(X)))
print("per-snapshot stitching vs the epoch's: phase difference at the band's far end, median %.0f deg, p90 %.0f deg (n %d)"%(
    np.median(walk),np.percentile(walk,90),len(walk)))
def estitch(k):
    y=X[k].copy(); g=EG[epoch[k]]
    for j in range(NG): y[EDG[j]:EDG[j+1]]/=g[j]
    return y
R0=np.array([prof(x) for x in X])
los0=int(np.argmax(np.mean(np.abs(np.array([prof(stitch(x)) for x in X[::20]]))**2,axis=0)))  # the stitched direct path
SH=20-los0   # roll every profile so the direct path sits at bin 20
los=20
_prof=prof
def prof(x): return np.roll(_prof(x),SH)
def variant(v):
    Y=X.copy()
    if v in (2,3): Y=np.array([stitch(x) for x in Y])
    if v in (4,5,6): Y=np.array([estitch(k) for k in range(len(Y))])
    R=np.array([prof(x) for x in Y])
    if v>=1: R=R/R[:,los:los+1]*np.abs(np.mean(R[:,los]))
    return R
def clutter(R,v,T,EP=None):
    if v==6:  # subtract each precoder epoch's own mean
        out=np.empty_like(R)
        for ee in np.unique(EP):
            sel=EP==ee; out[sel]=R[sel]-R[sel].mean(axis=0)
        return out
    if v in (3,5):  # per-epoch: subtract a 20 ms moving mean in slow time
        out=np.empty_like(R); j0=0
        for k in range(len(T)):
            sel=np.abs(T-T[k])<0.010; out[k]=R[k]-R[sel].mean(axis=0)
        return out
    return R-R.mean(axis=0)
# metrics over 150 ms windows
fD=np.arange(-400,401,4.0)
names=["V0 raw","V1 + LoS norm","V2 + PRG stitching per snapshot","V3 = V2 + 20 ms clutter removal","V4 + PRG stitching per precoder epoch","V5 = V4 + 20 ms clutter removal","V6 = V4 + per-epoch clutter removal"]
res=[]; maps={}
bin_m=3e8/(NF*120e3)
ghost=int(round(1/1.44e6/(1/(NF*120e3))))
for v in range(7):
    R=variant(v); Pm=np.mean(np.abs(R)**2,axis=0); Pm/=Pm.max()
    g=10*np.log10(np.max(Pm[(los+ghost-3)%NF:(los+ghost+4)%NF]) if (los+ghost+4)<NF else np.max(np.r_[Pm[los+ghost-3:],Pm[:3]]))
    # phase coherence of the direct path across consecutive snapshots
    ph=np.abs(np.mean(R[1:,los]*np.conj(R[:-1,los])/(np.abs(R[1:,los])*np.abs(R[:-1,los])+1e-20)))
    leak=[];floor=[]
    t0=T[0]; w=0.15
    while t0+w<T[-1]:
        sel=(T>=t0)&(T<t0+w)
        if sel.sum()>=60:
            Rw=clutter(R[sel],v,T[sel],epoch[sel]); tw=T[sel]
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
fig,ax=plt.subplots(2,7,figsize=(35,8))
for v in range(7):
    R=variant(v); Pm=np.mean(np.abs(R)**2,axis=0); Pm/=Pm.max()
    rr=(np.arange(NF)-los)*bin_m
    ax[0,v].plot(rr[:300],10*np.log10(Pm[:300]+1e-12)); ax[0,v].set_ylim(-60,2); ax[0,v].set_xlim(-50,500)
    ax[0,v].axvline(ghost*bin_m,color="r",ls=":"); ax[0,v].set_title(names[v]+": mean range profile"); ax[0,v].set_xlabel("excess path (m)")
    if v in maps:
        M,tw=maps[v]; M=10*np.log10(M/M.max()+1e-12)
        ax[1,v].imshow(M.T[:200],aspect="auto",origin="lower",extent=[fD[0],fD[-1],(0-los)*bin_m,(200-los)*bin_m],vmin=-50,vmax=0,cmap="viridis")
        ax[1,v].set_xlabel("Doppler (Hz)"); ax[1,v].set_ylabel("excess path (m)"); ax[1,v].set_title("range-Doppler, 150 ms at t=5 s")
plt.tight_layout(); plt.savefig(OUT,dpi=90)
