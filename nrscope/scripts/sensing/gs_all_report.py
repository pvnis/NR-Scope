"""Markdown tables from gs_all.py results (one or more .pkl files). usage: gs_all_report.py a.pkl [b.pkl ...]"""
import sys, pickle, numpy as np, collections as C
Rz=C.defaultdict(list); CEN=[]
for f in sys.argv[1:]:
    d=pickle.load(open(f,"rb"))
    for k,v in d["Rz"].items(): Rz[k].extend(v)
    CEN+= [dict(c,src=f) for c in d["CENSUS"]]
def med(v,p=None):
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    if not len(v): return "-"
    return ("%.1f"%np.median(v)) if p is None else ("%.1f / %.1f"%(np.median(v),np.percentile(v,p)))
def pct(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]; return ("%.0f%%"%(100*np.median(v))) if len(v) else "-"
print("## T0 census")
A=[c for c in CEN if c["role"]=="A"]; B=[c for c in CEN if c["role"]=="B"]
print("bursts %d; reference UE snapshots per burst median %d; PRBs per grant median %d; PRG jump median %.0f deg; epochs per burst median %d"%(
    len(A),np.median([c["n"] for c in A]),np.median([c["prb"] for c in A]),np.nanmedian([c["jump"] for c in A]),np.median([c["epochs"] for c in A])))
print("reference RNTIs:",C.Counter("%#x"%c["rnti"] for c in A).most_common())
print("\n| second UE | bursts | snapshots | PRBs/grant (median) | PRG jump inside its grants | precoder |\n|---|---|---|---|---|---|")
g=C.defaultdict(list)
for c in B: g[c["rnti"]].append(c)
for u,cs in sorted(g.items(),key=lambda x:-sum(c["n"] for c in x[1])):
    print("| %#x | %d | %d | %d | %.0f deg | %s |"%(u,len(cs),sum(c["n"] for c in cs),np.median([c["prb"] for c in cs]),np.nanmedian([c["jump"] for c in cs]),cs[0]["typ"]))
types=sorted({k[1] for k in Rz if k[0]=="T2"})
print("\n## T1 same PRG: second UE vs reference within 2 ms, share NOT explained by one number (4 ant x 12 sc)")
for t in types: print("- %s: median %s (n %d)"%(t,pct(Rz[("T1",t,"B vs A")]),len(Rz[("T1",t,"B vs A")])))
print("\n## T2 bridge distance (number from a reference snapshot ~D away; truth = reference within 1 ms)")
for t in types:
    print("\n### second UEs with %s precoder"%t)
    print("oracle (best possible) channel error: one number raw %s, self-stitched %s, per PRG %s"%(pct(Rz[("T2",t,"oracle","W")]),pct(Rz[("T2",t,"oracle","S")]),pct(Rz[("T2",t,"oracle","P")])))
    print("| D | cases | W phase med/p75 (deg) | W amp (dB) | W err | S phase | S err | P err (per PRG) | reference's own change |\n|---|---|---|---|---|---|---|---|---|")
    for D in (0.002,0.01,0.03,0.1,0.3,1.0):
        n=len(Rz[("T2",t,D,"W","phase")])
        if not n: continue
        print("| %g ms | %d | %s | %s | %s | %s | %s | %s | %s |"%(D*1e3,n,med(Rz[("T2",t,D,"W","phase")],75),med(Rz[("T2",t,D,"W","amp")]),pct(Rz[("T2",t,D,"W","err")]),
            med(Rz[("T2",t,D,"S","phase")],75),pct(Rz[("T2",t,D,"S","err")]),pct(Rz[("T2",t,D,"P","err")]),pct(Rz[("T2",t,D,"A own change")])))
print("\n## T3 edge matching, controlled (second UE's own grant cut in two, known number, 4 antennas, 2 PRBs fit)")
for t in types:
    print("\n### %s precoder\n| gap | cases | raw: phase med/p90, amp | self-stitched: phase med/p90, amp |\n|---|---|---|---|"%t)
    for gp in (0,1,2,4,8):
        n=len(Rz[("T3",t,gp,"W","phase")])
        if n: print("| %d PRBs | %d | %s deg, %s dB | %s deg, %s dB |"%(4*gp,n,med(Rz[("T3",t,gp,"W","phase")],90),med(Rz[("T3",t,gp,"W","amp")]),med(Rz[("T3",t,gp,"S","phase")],90),med(Rz[("T3",t,gp,"S","amp")])))
print("\n## T4 edge matching, real slots with reference + second UE side by side")
T4=Rz["T4"]
for var in ("W","S"):
    rows=[r for r in T4 if r["var"]==var]
    if not rows: continue
    print("\n### variant %s (%s)"%(var,"second UE raw" if var=="W" else "second UE self-stitched first"))
    print("DM-RS symbols %d (slots %d); gaps (PRBs): %s"%(len(rows),len({(r["burst"],round(r["t"]*2000)) for r in rows}),C.Counter(r["gap"] for r in rows).most_common(6)))
    print("| subset | n | edge vs anchor phase med/p90 | amp | err edge | err anchor | err combined | err oracle | err raw join | width ref / A alone / combined (m) |\n|---|---|---|---|---|---|---|---|---|---|")
    for nm,f in (("touching (gap 0)",lambda r:r["gap"]==0),("touching, both runs >= 32 PRBs",lambda r:r["gap"]==0 and min(r["runA"],r["runB"])>=32),
                 ("touching, a run < 32 PRBs",lambda r:r["gap"]==0 and min(r["runA"],r["runB"])<32),("gap 4-8 PRBs",lambda r:0<r["gap"]<=8),("gap > 8 PRBs",lambda r:r["gap"]>8),
                 ("anchor within 2 ms",lambda r:r["dta"]<=0.002),("anchor 2-30 ms only",lambda r:0.002<r["dta"]<=0.03),("no anchor",lambda r:not np.isfinite(r["dta"]))):
        s=[r for r in rows if f(r)]
        if not s: continue
        ea=[abs(r["edge_vs_anchor_deg"]) for r in s if "edge_vs_anchor_deg" in r]; ead=[abs(r["edge_vs_anchor_db"]) for r in s if "edge_vs_anchor_db" in r]
        G=lambda k:[r[k] for r in s if k in r]
        print("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s / %s / %s |"%(nm,len(s),med(ea,90),med(ead),pct(G("err_edge")),pct(G("err_anchor")),pct(G("err_combined")),pct(G("err_oracle")),pct(G("err_raw")),med(G("w_ref")),med(G("w_A")),med(G("w_combined"))))
    print("\nquality gate (no ground truth needed): keep a merge only if the second UE fits the reference on the anchor PRBs with residual below a threshold")
    print("| gate | symbols kept | share kept | err anchor | err combined | err oracle |\n|---|---|---|---|---|---|")
    base=[r for r in rows if "err_anchor" in r and np.isfinite(r.get("fitres",np.nan))]
    for th in (1e9,0.5,0.3,0.2,0.1):
        s=[r for r in base if r["fitres"]<th]
        if s: print("| %s | %d | %.0f%% | %s | %s | %s |"%("none" if th>1e8 else "fit residual < %d%%"%(100*th),len(s),100*len(s)/len(base),pct([r["err_anchor"] for r in s]),pct([r["err_combined"] for r in s]),pct([r["err_oracle"] for r in s])))
    fe=[r for r in rows if "err_edge" in r and "err_anchor" in r]
    if fe: print("edge better than anchor in %d of %d symbols with both + reference"%(sum(r["err_edge"]<r["err_anchor"] for r in fe),len(fe)))
print("\n## T5 coverage (PRGs with usable DM-RS)")
print("| set | PRBs per slot median | mean | slots >= 160 PRBs | 150 ms window covered (median) |\n|---|---|---|---|---|")
for sel in ("A alone","A + mergeable B","A + all B"):
    v=np.array(Rz[("T5",sel,"PRBs per slot")]); w=np.array(Rz[("T5",sel,"window cover")])
    if len(v): print("| %s | %d | %.0f | %.0f%% | %.0f%% |"%(sel,np.median(v),v.mean(),100*np.mean(v>=160),100*np.median(w)))

print("\n## T6 stitching bias and amplitude jumps (wideband second UE: true jumps are 0, so what stitching finds is bias)")
print("| UE kind | edges | jump phase |median| / p90 | jump amplitude |median| / p90 (dB) | jumps within 20 deg of +-90 |\n|---|---|---|---|---|")
for k in ["A"]+types:
    ph=np.array(Rz[("T6",k,"jump_deg")]); am=np.array(Rz[("T6",k,"jump_db")])
    if len(ph): print("| %s | %d | %.1f / %.1f deg | %.2f / %.2f | %.0f%% |"%("reference (tablet)" if k=="A" else "second UE, "+k,len(ph),np.median(abs(ph)),np.percentile(abs(ph),90),np.median(abs(am)),np.percentile(abs(am),90),100*np.mean(abs(abs(ph)-90)<20)))
print("\nPRG ghost (+208 m, 4-PRB period: bins 82-88 after the peak) and floor (bins 20-70), mean of peak-normalised profiles:")
print("| UE kind | stitching | snapshots | ghost | floor |\n|---|---|---|---|---|")
for k in ["A"]+types:
    for v,lab in (("prof raw","none"),("prof snap","per snapshot, 4 antennas"),("prof epoch","per epoch, 4 antennas")):
        P=Rz.get(("T6",k,v),[])
        if len(P):
            m=np.mean(np.array(P),0); print("| %s | %s | %d | %.1f dB | %.1f dB |"%("reference" if k=="A" else k,lab,len(P),10*np.log10(m[82:89].max()),10*np.log10(np.median(m[20:70]))))
