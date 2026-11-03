import csv,sys,numpy as np,collections as C
R=sys.argv[1]; O=R+"/replay1"; P=20480
thr=[]
for l in open(R+"/dl_rate.log"):
    f=l.split()
    if len(f)>=2 and f[1].isdigit():
        h,m,s=f[0].split(":"); thr.append((int(h)*3600+int(m)*60+float(s),int(f[1])))
t_end=15*3600+54*60+26.0; t_start=t_end-20
y=np.array([np.mean([v for t,v in thr if t_start+i<t<=t_start+i+1] or [0]) for i in range(20)])
# the DCI file in processing order: unwrap the SFN nearest to the running slot
rows=sorted(csv.DictReader(open(O+"/golden_sniffer.csv")),key=lambda r:float(r["timestamp"]))
T=[];A=[];cur=None
for r in rows:
    a=int(r["sfn"])*20+int(r["slot"])
    if cur is not None: a+=round((cur-a)/P)*P
    cur=a if cur is None else max(cur,a); T.append(float(r["timestamp"])); A.append(a); r["_a"]=a
T=np.array(T); A=np.array(A); a0=A.min()
def abs_of(r):
    a=int(r["sfn"])*20+int(r["slot"]); est=np.interp(float(r["timestamp"]),T,A); return a+round((est-a)/P)*P
dci=C.defaultdict(lambda: np.zeros(21)); dl=C.defaultdict(lambda: np.zeros(21)); sz=C.defaultdict(C.Counter)
for r in rows:
    if r["payload_hex"] and int(r["payload_hex"][0],16)>>3==1:
        sec=int((r["_a"]-a0)*0.0005)
        if 0<=sec<21: dci[r["rnti"]][sec]+=1
        sz[r["rnti"]][r["size"]]+=1
for r in csv.DictReader(open(O+"/golden_sniffer_grants.csv")):
    sec=int((abs_of(r)-a0)*0.0005)
    if 0<=sec<21: dl[r["rnti"]][sec]+=int(r["nof_rb"])
print("span of decoded slots: %.1f s"%((A.max()-a0)*0.0005))
print("tablet throughput per second (Mbps):", " ".join("%d"%v for v in y))
res=[]
for u,v in dci.items():
    if v.sum()<200: continue
    best=max((np.corrcoef(v[l:l+19],y[:19])[0,1] if v[l:l+19].std()>0 else 0,l) for l in (0,1,2))
    res.append((best[0],best[1],u,int(v.sum()),int(dl[u].sum())))
res.sort(reverse=True)
print("UEs by correlation of DL DCIs per second with the throughput (lag in s):")
for c,l,u,n,p in res[:6]:
    print("  %s  r=%.2f lag %d  DL DCIs %5d  granted PRBs %7d  sizes %s"%(u,c,l,n,p,dict(sz[u].most_common(2))))
    print("      DCIs/s: "+" ".join("%d"%x for x in dci[u][:20]))
