"""Per-PRG alignment of PDSCH DM-RS against the 8-port CSI-RS channel.
usage: csi_align.py ue_dmrs.bin <rnti hex> [pci]
Reads GSD1 (grant DM-RS, port 1000 per CDM pair) and GSC1 (raw REs of the CSI-RS symbol) records."""
import sys, struct, numpy as np, collections as C
F=sys.argv[1]; TAB=int(sys.argv[2],16); PCI=int(sys.argv[3]) if len(sys.argv)>3 else 632
P=20480; N=273; NF=1024; MAXAGE=20  # slots between the CSI-RS and the grant

def gold(cinit,n):
    Nc=1600; x1=np.zeros(Nc+n+31,np.uint8); x2=np.zeros(Nc+n+31,np.uint8); x1[0]=1
    for i in range(31): x2[i]=(cinit>>i)&1
    for i in range(Nc+n):
        x1[i+31]=x1[i+3]^x1[i]; x2[i+31]=x2[i+3]^x2[i+2]^x2[i+1]^x2[i]
    return x1[Nc:Nc+n]^x2[Nc:Nc+n]
_seq={}
def csirs_seq(slot,l,nid,m):
    k=(slot,l)
    if k not in _seq:
        ci=((1<<10)*(14*slot+l+1)*(2*nid+1)+nid)%(1<<31); c=gold(ci,2*m).astype(float)
        _seq[k]=((1-2*c[0::2])+1j*(1-2*c[1::2]))/np.sqrt(2)
    return _seq[k]

raw=open(F,"rb").read(); i=0; hist=[]
def unwrap(sfn,sl):
    a=sfn*20+sl
    if hist: a+=round((float(np.median(hist[-64:]))-a)/P)*P
    hist.append(a); return a
grants=[]; csi={}   # csi[a] = H (N, 8) per CRB, nan where absent
while i+8<=len(raw):
    tag=raw[i:i+4]
    if tag==b"GSD1":
        if i+20>len(raw): break
        sfn,sl,ns,rn,npb,occ,coh=struct.unpack_from("<HBBHHff",raw,i+4); j=i+20
        if j+npb+ns*(1+3*npb*8)>len(raw): break
        m=np.frombuffer(raw,np.uint8,npb,j).astype(bool); j+=npb; a=unwrap(sfn,sl); syms={}
        for _ in range(ns):
            l=raw[j]; j+=1; syms[l]=np.frombuffer(raw,np.complex64,3*npb,j).astype(complex); j+=3*npb*8
        i=j
        if rn==TAB and occ<4 and coh>=0.6: grants.append((a,m,syms))
    elif tag==b"GSC1":
        sfn,sl,l,fc,npb=struct.unpack_from("<HBBHH",raw,i+4); j=i+12
        if j+12*npb*8>len(raw): break
        y=np.frombuffer(raw,np.complex64,12*npb,j).astype(complex).reshape(npb,12); i=j+12*npb*8
        a=unwrap(sfn,sl); r=csirs_seq(sl,l,PCI,2*(fc+npb)+4); H=np.full((N,8),np.nan,complex)
        for q in range(npb):
            n=fc+q
            if n>=N: break
            for g in range(4):
                z0=y[q,2*g]*np.conj(r[2*n]); z1=y[q,2*g+1]*np.conj(r[2*n+1])
                H[n,2*g]=(z0+z1)/2; H[n,2*g+1]=(z0-z1)/2
        csi[a]=H
    else:
        i+=1
print("tablet grants %d, CSI-RS symbols %d"%(len(grants),len(csi)))
ca=np.array(sorted(csi));
# --- the CSI-RS ports: power, PRB-to-PRB coherence, correlation between ports
Hs=np.array([csi[a] for a in ca])           # (T, N, 8)
pw=np.nanmean(np.abs(Hs)**2,axis=(0,1)); pw/=pw.max()
co=[np.abs(np.nansum(Hs[:,1:,p]*np.conj(Hs[:,:-1,p])))/np.nansum(np.abs(Hs[:,1:,p])**2) for p in range(8)]
print("CSI-RS port power (dB):", " ".join("%.1f"%(10*np.log10(x)) for x in pw))
print("CSI-RS port coherence PRB to PRB:", " ".join("%.2f"%x for x in co))
X=Hs.reshape(-1,8); X=X[np.all(np.isfinite(X),axis=1)]; Rc=X.conj().T@X/len(X); ev=np.linalg.eigvalsh(Rc)[::-1]
print("CSI-RS port covariance eigenvalues (dB re the first):", " ".join("%.1f"%(10*np.log10(e/ev[0])) for e in ev))
dd=np.diff(ca); print("CSI-RS period (slots), most common:", C.Counter(dd).most_common(4))

# CSI-RS ports at the DM-RS CDM-pair positions (subcarriers 1, 5, 9 of every CRB), complex linear interpolation
fp=np.array([12*n+2*(p//2)+0.5 for n in range(N) for p in range(0,1)])
def csi_at_dmrs(H):
    out=np.zeros((3*N,8),complex); fd=np.array([12*n+k for n in range(N) for k in (1,5,9)],float)
    for p in range(8):
        f=12*np.arange(N)+2*(p//2)+0.5; v=H[:,p]; ok=np.isfinite(v)
        out[:,p]=np.interp(fd,f[ok],v[ok].real)+1j*np.interp(fd,f[ok],v[ok].imag)
    return out
_ci={}
def nearest_csi(a):
    k=np.searchsorted(ca,a); best=None
    for j in (k-1,k):
        if 0<=j<len(ca) and abs(ca[j]-a)<=MAXAGE and (best is None or abs(ca[j]-a)<abs(best-a)): best=ca[j]
    if best is None: return None,None
    if best not in _ci: _ci[best]=csi_at_dmrs(csi[best])
    return _ci[best],a-best

def segments(m):
    out=[];p=0
    while p<N:
        if m[p]:
            q=p
            while q+1<N and m[q+1]: q+=1
            out.append((p,q)); p=q+1
        else: p+=1
    return out
def lsfit(A,h):
    w,*_=np.linalg.lstsq(A,h,rcond=None); r=h-A@w; return w,np.sum(np.abs(r)**2)/np.sum(np.abs(h)**2)

# --- test 1: one 8-port precoder for the whole grant (grants of >= 48 PRBs)
res=C.defaultdict(list); ages=[]; res_exc=[]
for a,m,syms in grants:
    A,age=nearest_csi(a)
    if A is None or 2 not in syms: continue
    ages.append(age)
    idx=np.nonzero(np.repeat(m,3))[0]
    if len(idx)<3*48: continue
    h=syms[2][idx]
    # one precoder over the grant? phase steps at CRB%4 boundaries vs inside the PRGs (linear phase removed)
    dv=h[1:]*np.conj(h[:-1]); ok=np.diff(idx)==1; ph=np.angle(np.sum(dv[ok]))
    crb=idx[1:]//3; bnd=(idx[1:]%3==0)&(crb%4==0)&ok; ins=ok&~bnd
    dev=np.abs(np.angle(dv*np.exp(-1j*ph)))*57.3
    exc=np.median(dev[bnd])-np.median(dev[ins]) if bnd.sum()>4 else 99
    res_exc.append(exc)
    if exc>8: continue
    if 11 in syms:   # noise floor: the same precoder on the second DM-RS symbol
        d=syms[11][idx]; res["noise"].append(np.sum(np.abs(h-d)**2)/2/np.sum(np.abs(h)**2))
    res["csi"].append(lsfit(A[idx],h)[1])
    sh=np.roll(A,3*24,axis=0); res["csi shifted 24 PRB"].append(lsfit(sh[idx],h)[1])
    f=idx.astype(float); B=3*(idx[-1]-idx[0]+1)/3.0; k0=ph/(2*np.pi)*B   # the direct-path delay, in bins of the band
    D=np.exp(2j*np.pi*np.outer(f,k0+np.arange(8)-2)/B)  # 8 delay columns from 2 bins before the direct path
    res["8 delay columns"].append(lsfit(D,h)[1])
    res["port 1 only"].append(lsfit(A[idx][:,:1],h)[1])
    # per PRG: an 8-port precoder per 4-CRB PRG (residual over all PRGs)
    num=den=0
    for g0 in range(0,N,4):
        k=[x for x in idx if g0*3<=x<(g0+4)*3]
        if len(k)<6: continue
        _,rr=lsfit(A[k],syms[2][k]); e=np.sum(np.abs(syms[2][k])**2); num+=rr*e; den+=e
    if den: res["csi per PRG"].append(num/den)
print("\nCSI-RS age at the grant (slots): median %d, |age| p90 %d, grants with CSI-RS within %d slots %d"%(
    np.median(ages),np.percentile(np.abs(ages),90),MAXAGE,len(ages)))
print("PRG phase-step excess per grant: median %.0f deg, grants with one precoder (excess < 8 deg) %d"%(np.median(res_exc),len(res["csi"])))
print("\nTest 1: DM-RS of a grant (>= 48 PRBs, first DM-RS symbol) fitted by the CSI-RS ports")
print("| model | grants | residual median | p25 | p75 |\n|---|---|---|---|---|")
for k,v in res.items():
    v=np.array(v); print("| %s | %d | %.1f dB | %.1f dB | %.1f dB |"%(k,len(v),*(10*np.log10(np.percentile(v,[50,25,75])))))
