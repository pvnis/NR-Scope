import sys,pickle,numpy as np,collections as C
rows=[];cen=[]
for f in sys.argv[1:]:
    d=pickle.load(open(f,"rb")); rows+=[dict(r,src=f) for r in d["Rz"].get("T4",[])]; cen+=d["CENSUS"]
typ={c["rnti"]:c.get("typ") for c in cen if c["role"]=="B"}
g=C.defaultdict(list)
for r in rows:
    if r["var"]=="S" and r["gap"]==0 and "err_oracle" in r: g[r["rnti"]].append(r)
print("| UE | precoder | touching symbols | UE run PRBs | edge | anchor | combined | oracle | edge-anchor phase |")
print("|---|---|---|---|---|---|---|---|---|")
for u,rs in sorted(g.items(),key=lambda kv:-len(kv[1])):
    if len(rs)<10: continue
    m=lambda k:"%.0f%%"%(100*np.median([r[k] for r in rs if k in r])) if any(k in r for r in rs) else "-"
    ph=[abs(r["edge_vs_anchor_deg"]) for r in rs if "edge_vs_anchor_deg" in r]
    print("| %#x | %s | %d | %d | %s | %s | %s | %s | %s |"%(u,typ.get(u),len(rs),np.median([r["runB"] for r in rs]),m("err_edge"),m("err_anchor"),m("err_combined"),m("err_oracle"),("%.0f deg"%np.median(ph)) if ph else "-"))
