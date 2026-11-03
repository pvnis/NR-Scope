import csv, sys, collections as C, numpy as np
O=sys.argv[1]; P=20480; N=273
# survey: strict DM-RS map per slot, SFN unwrapped against the median of the last slots
air={}; hist=[]
for l in list(open(O+"/survey.csv"))[1:]:
    r=l.rstrip("\n").split(",")
    if len(r)<40: continue
    a=int(r[0])*20+int(r[1])
    if hist: a+=round((float(np.median(hist[-64:]))-a)/P)*P
    hist.append(a)
    b=np.array([int(c) for x in r[33] for c in bin(int(x,16))[2:].zfill(4)][:N],bool); p=0
    while p<N:  # runs under 3 PRBs are noise hits
        if b[p]:
            q=p
            while q+1<N and b[q+1]: q+=1
            if q-p+1<3: b[p:q+1]=False
            p=q+1
        else: p+=1
    air[a]=b
s0=min(air); off=s0-s0%P
# grants: each takes the slot of its own DCI, placed by the DCI's air time
dcis=[r for r in csv.DictReader(open(O+"/golden_sniffer.csv")) if float(r["timestamp"])<1e6]
raw=np.array([int(r["sfn"])*20+int(r["slot"]) for r in dcis]); ta=np.array([float(r["timestamp"]) for r in dcis])
c=np.median((raw-ta/0.0005)[:200]); absd=raw+np.round((ta/0.0005+c-raw)/P)*P
key=C.defaultdict(list)
for r,aa in zip(dcis,absd): key[(r["rnti"],r["sfn"],r["slot"])].append(int(aa)+off)
for k in key: key[k].sort()
used=C.Counter(); gr=C.defaultdict(lambda: np.zeros(N,bool)); tb=C.defaultdict(lambda: np.zeros(N,bool))
TAB=sys.argv[2] if len(sys.argv)>2 else ""
for r in csv.DictReader(open(O+"/golden_sniffer_grants.csv")):
    k=(r["rnti"],r["sfn"],r["slot"]); v=key.get(k)
    if not v: continue
    a=v[min(used[k],len(v)-1)]; used[k]+=1
    if float(r["dmrs_coherence"])<0.4: continue
    for x in r["prb_runs"].split(";"):
        if x:
            s_,e=map(int,x.split("-")); gr[a][s_:e+1]=True
            if r["rnti"]==TAB: tb[a][s_:e+1]=True
def widest(m):
    best=cur=0
    for v in m:
        cur=cur+1 if v else 0; best=max(best,cur)
    return best
def psf(m):
    i=np.nonzero(m)[0]; x=np.repeat(m[i[0]:i[-1]+1],6).astype(float); x*=np.hanning(len(x))
    Pw=np.abs(np.fft.fft(x,4096))**2; Pw/=Pw.max(); h=int(np.ceil(2*4096/len(x)))+1
    return 10*np.log10(Pw[h:4096-h].max()+1e-12)
srcs=[("DM-RS on air only",lambda a:air[a]),("decoded grants (DM-RS seen)",lambda a:gr.get(a,np.zeros(N,bool))),("both",lambda a:air[a]|gr.get(a,np.zeros(N,bool)))]
if TAB: srcs.append(("the tablet's grants only",lambda a:tb.get(a,np.zeros(N,bool))))
n=len(air); secs=(max(air)-min(air)+1)*0.0005; dl_per_s=n/secs
print("%d DL slots over %.1f s (%.0f DL slots/s)"%(n,secs,dl_per_s))
print("| source | rule | slots | snapshots/s | median span | median fill | peak sidelobe |")
print("|---|---|---|---|---|---|---|")
for name,get in srcs:
    W=[];S=[];F=[];M=[]
    for a in air:
        m=get(a); W.append(widest(m))
        if m.sum()>=12:
            i=np.nonzero(m)[0]; sp=i[-1]-i[0]+1; S.append(sp); F.append(m.sum()/sp); M.append(m)
        else: S.append(0); F.append(0); M.append(None)
    W=np.array(W); S=np.array(S); F=np.array(F)
    for rule,sel in (("contiguous >= 137",W>=137),("contiguous >= 200",W>=200),("contiguous >= 250",W>=250),("span >= 137, any fill",S>=137),("span >= 137, fill >= 50%",(S>=137)&(F>=0.5)),("span >= 137, fill >= 70%",(S>=137)&(F>=0.7))):
        k=np.nonzero(sel)[0]
        if len(k)==0: print("| %s | %s | 0%% | 0 | - | - | - |"%(name,rule)); continue
        ps=[psf(M[j]) for j in k[::max(1,len(k)//400)]]
        print("| %s | %s | %.1f%% | %.0f | %d | %.2f | %.1f dB |"%(name,rule,100*len(k)/n,len(k)/secs,np.median(S[k]),np.median(F[k]),np.median(ps)))
