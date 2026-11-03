import sys, json, csv, collections as C, numpy as np
hist=[]; O=sys.argv[1]; TAB=sys.argv[2] if len(sys.argv)>2 else ""; P=20480; N=273
exec(open("/tmp/mkdmrs.py").read().split("a0=min(s[0] for s in slots)")[0].replace(
"""    if prev>=0 and a<prev-10000: wrap+=20480
    prev=a; a+=wrap""","""    if hist: a+=round((float(np.median(hist[-64:]))-a)/20480)*20480
    hist.append(a); prev=a"""))
# grants of the tablet: each takes the slot of its own DCI, placed by that DCI's air time
dcis=[r for r in csv.DictReader(open(O+"/golden_sniffer.csv")) if float(r["timestamp"])<1e6]  # tracked: air clock
raw=np.array([int(r["sfn"])*20+int(r["slot"]) for r in dcis]); ta=np.array([float(r["timestamp"]) for r in dcis])
c=np.median((raw-ta/0.0005)[:200])
absd=raw+np.round((ta/0.0005+c-raw)/P)*P
# the survey's first slot in the same frame of reference
s0=min(x[0] for x in slots); sraw=s0 % P
off=s0-sraw  # survey timeline starts in its own period 0; the DCI timeline too
key=C.defaultdict(list)
for r,aa in zip(dcis,absd): key[(r["rnti"],r["sfn"],r["slot"])].append(int(aa))
for k in key: key[k].sort()
used=C.Counter(); tab=C.defaultdict(lambda: np.zeros(N,bool)); tstat=C.Counter()
for r in csv.DictReader(open(O+"/golden_sniffer_grants.csv")):
    if r["rnti"]!=TAB: continue
    k=(r["rnti"],r["sfn"],r["slot"]); v=key.get(k)
    if not v: continue
    a=v[min(used[k],len(v)-1)]+off; used[k]+=1
    for x in r["prb_runs"].split(";"):
        if x: s_,e=map(int,x.split("-")); tab[a][s_:e+1]=True
    tstat["grants"]+=1; tstat["prb"]+=int(r["nof_rb"]); tstat["vis"]+=float(r["dmrs_coherence"])>=0.4
A=absd
for r in runs:
    m=tab.get(r[0]); r.append(int(m is not None and m[r[1]:r[2]+1].mean()>=0.5))
a0=min(s[0] for s in slots)
for s in slots: s[0]-=a0
for r in runs: r[0]-=a0
exec(open("/tmp/mkdmrs.py").read().split("# PRB-weighted summaries")[1].replace('json.dump(D,open(O+"/dmrs_viz.json","w"),separators=(",",":"))',''))
D["tablet"]=TAB; D["tablet_stats"]=[tstat["grants"],tstat["prb"],tstat["vis"]]
D["n_tab_runs"]=sum(r[10] for r in runs)
json.dump(D,open(O+"/dmrs_viz.json","w"),separators=(",",":"))
print("survey span %.1f s, DCI span %.1f s"%((max(x[0] for x in slots)-min(x[0] for x in slots))*0.0005,(A.max()-A.min())*0.0005)); print("tablet",TAB,dict(tstat),"runs under its grants",D["n_tab_runs"],"of",len(runs))
