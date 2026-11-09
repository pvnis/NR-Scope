import json,glob,sys,numpy as np,collections as C
D=sys.argv[1] if len(sys.argv)>1 else "usemap"
fs=sorted(glob.glob(D+"/sum_b*.json"),key=lambda s:int(s.split("_b")[1].split(".")[0]))
CATN=10; agg={sc:[0]*CATN for sc in ("tab","notab")}; W={sc:dict(d=[],ref=[],ne=[],lo=[],inc=[],trs=[]) for sc in agg}; blind={sc:[0]*CATN for sc in agg}; rows=[]; bursts=[]
oa=[]
for f in fs:
    d=json.load(open(f)); bi=d["burst"]; J=json.load(open(D+"/um_b%d.json"%bi))
    dl=C.Counter(x[0]//300 for x in J["slots"]); s7=C.Counter()
    for a,lo,hi,c,u in np.array(J["seg"]["notab"]).reshape(-1,5):
        if c==7: s7[a//300]+=hi-lo+1
    for k,n in dl.items():
        if n>=50: oa.append(s7[k]/(273*n))
    bl={i for i,r in enumerate(J["rnti"]) if r[1]=="blind"}
    for sc in agg:
        sinc=C.Counter()
        for a,lo,hi,c,u in np.array(J["seg"][sc]).reshape(-1,5):
            if u in bl: blind[sc][int(c)]+=int(hi-lo+1)
            if c in (1,5,6): sinc[a//300]+=hi-lo+1
        W[sc]["inc"]+=[sinc[k]/(273*n) for k,n in dl.items() if n>=50]
    r=dict(b=bi,ues=d["n_ues"],tab=d["tablet"])
    for sc in agg:
        bc={int(k):v for k,v in d[sc]["bycat"].items()}
        for k,v in bc.items(): agg[sc][k]+=v
        ws=[w for w in d[sc]["win"] if w[3]>=50]
        dd=[w[2]/(273*w[3]) for w in ws]
        W[sc]["trs"]+=[w[7]/(273*w[3]) for w in ws if len(w)>7]
        W[sc]["d"]+=dd; W[sc]["ne"]+=[w[4] for w in ws]; W[sc]["lo"]+=[w[5] for w in ws]
        r[sc]=dict(c=[bc.get(k,0) for k in range(CATN)],dens=float(np.median(dd)) if dd else 0)
    rows.append(r); bursts.append(dict(b=bi,span=J["span"],ues=d["n_ues"],
        others=[x for x in d["rnti"][1:]]))
def st(v):
    v=np.array(v,float); return dict(med=float(np.median(v)),p10=float(np.percentile(v,10)),p90=float(np.percentile(v,90)),ge71=float(np.mean(v>=0.071)),zero=float(np.mean(v<0.001)),n=len(v))
SUM=dict(agg=agg,blind=blind,win={sc:dict(dens=st(W[sc]["d"]),inc=st(W[sc]["inc"]),trs=st(W[sc]["trs"]) if W[sc]["trs"] else st([0]),ne=float(np.median(W[sc]["ne"])),lo=float(np.median(W[sc]["lo"]))) for sc in agg},
         onair=st(oa),rows=rows,bursts=bursts)
open(D+"/summary.json","w").write(json.dumps(SUM))
html=open("usemap_template.html").read().replace("__DATA__",json.dumps(SUM,separators=(",",":")))
open(D+"/index.html","w").write(html); print(len(html), SUM["win"], SUM["onair"])
