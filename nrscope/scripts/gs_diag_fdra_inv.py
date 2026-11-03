import sys, re, csv, collections as C, numpy as np
exec(open("/tmp/fdra_off.py").read().split('print("%-10s')[0].replace('O=sys.argv[1]; N=273','O=sys.argv[1]; N=273'))
src={}
for l in open(O+"/nrscope.log",errors="replace"):
    m=re.search(r"UE (0x[0-9a-f]+) FDRA: .*?(taken from|lift)",l)
    if m: src[m.group(1)]="cell" if m.group(2)=="taken from" else "own"
coh=C.defaultdict(list)
for r in csv.DictReader(open(O+"/golden_sniffer_grants.csv")): coh[r["rnti"]].append(float(r["dmrs_coherence"]))
inv=[u for u,v in coh.items() if len(v)>=20 and np.median(v)<0.4 and src.get(u)!="own"]
H2=[(t,o) for t in ("t1","t0","dyn") for o in range(1,8)]
ue=C.defaultdict(list)
for key,rs in per.items():
    for r in rs:
        if r["rnti"] in inv: ue[r["rnti"]].append((key,bits(r["payload_hex"],int(r["size"])),int(r["size"])))
print("UEs whose grants are mostly invisible (median coherence < 0.4), every layout tried on their own DL DCIs:")
for u in inv:
    sc=[]
    for h in H2:
        L=[]
        for key,b,sz in ue[u]:
            if key not in occ: continue
            oc=occ[key]; f=oc.mean()
            if f<0.02 or f>0.98: continue
            p=prbs(b,h)
            if p is None or len(p)==0: continue
            L.append((oc[p].mean()-f)/(1-f))
        if len(L)>=10: sc.append((np.mean(L),h,len(L)))
    sc.sort(reverse=True)
    sizes=sorted(set(x[2] for x in ue[u]))
    print("  %s sizes %s: best %s lift %.2f (n %d), next %s %.2f; t0 off3: %.2f, t0 off2: %.2f"%(u,sizes,sc[0][1],sc[0][0],sc[0][2],sc[1][1],sc[1][0],
          next((s for s,h,n in sc if h==("t0",3)),float("nan")),next((s for s,h,n in sc if h==("t0",2)),float("nan"))))
