import sys, csv, collections as C, numpy as np
O=sys.argv[1]; N=273
S={}
for l in list(open(O+"/survey.csv"))[1:]:
    r=l.rstrip("\n").split(",")
    if len(r)<41: continue
    k=(int(r[0]),int(r[1]))
    b=np.array([int(c) for h in r[33] for c in bin(int(h,16))[2:].zfill(4)][:N],bool)
    c8=list(map(int,r[-1].split(":")))
    S[k]=(b,c8)
G=C.defaultdict(list)
for r in csv.DictReader(open(O+"/golden_sniffer_grants.csv")):
    runs=[tuple(map(int,x.split("-"))) for x in r["prb_runs"].split(";") if x]
    G[(int(r["sfn"]),int(r["slot"]))].append((runs,float(r["dmrs_coherence"])))
noise=np.array([np.mean([c8[l-1] for l in (1,4,6,9,13)]) for b,c8 in S.values()])
print("noise: strict hits on symbols without DM-RS, per slot: mean %.1f, p90 %.0f of 273"%(noise.mean(),np.percentile(noise,90)))
n2=np.array([b.sum() for b,c in S.values()])
print("symbol 2 strict DM-RS PRBs per slot: p10 %d p25 %d median %d"%tuple(np.percentile(n2,[10,25,50])))
cls=C.Counter(); detail=C.defaultdict(list)
for k,(b,c8) in S.items():
    n=b.sum()
    if n==0 or n>40: continue
    # are the hits concentrated in few RBGs?
    rbg=C.Counter(int(p)//16 for p in np.nonzero(b)[0]); top=sum(v for _,v in rbg.most_common(2))/n
    nz=np.mean([c8[l-1] for l in (1,4,6,9,13)])
    m=np.zeros(N,bool); cohs=[]
    for runs,co in G.get(k,[]):
        for s,e in runs: m[s:e+1]=True
        cohs.append(co)
    inside=(b&m).sum()/n
    if n<=nz*2.5: c="noise level (no PDSCH DM-RS we can see)"
    elif G.get(k):
        c="decoded grant(s) there, only part of their PRBs detected" if (m.sum()>2*n and inside>0.5) else "decoded grant(s), plus DM-RS of UEs not decoded"
    elif top>0.7: c="clustered in 1-2 RBGs: an undecoded UE's small grant"
    else: c="scattered, no decoded grant: weak beam or undecoded UEs"
    cls[c]+=1; detail[c].append((n, m.sum(), inside, np.mean(cohs) if cohs else np.nan))
tot=sum(cls.values())
print("\nslots with 1-40 strict DM-RS PRBs: %d of %d DL slots"%(tot,len(S)))
for c,v in cls.most_common():
    d=np.array(detail[c],float)
    print("  %5.1f%%  %-62s DM-RS PRBs %4.1f, decoded grant PRBs %5.1f, share inside %.2f, grant coherence %.2f"%(100*v/tot,c,d[:,0].mean(),d[:,1].mean(),np.nanmean(d[:,2]),np.nanmean(d[:,3]) if np.isfinite(d[:,3]).any() else float("nan")))
# coherence of decoded grants vs how many of their PRBs pass the strict test
fr=[];co=[]
for k,v in G.items():
    if k not in S: continue
    b=S[k][0]
    for runs,c in v:
        m=np.zeros(N,bool)
        for s,e in runs: m[s:e+1]=True
        fr.append((b&m).sum()/m.sum()); co.append(c)
fr=np.array(fr); co=np.array(co)
print("\ndecoded grants: share of their PRBs detected strictly: p10 %.2f median %.2f p90 %.2f"%tuple(np.percentile(fr,[10,50,90])))
print("  grants with <30%% of their PRBs detected: %.0f%% (their mean DM-RS coherence %.2f, others %.2f)"%(100*np.mean(fr<0.3),co[fr<0.3].mean(),co[fr>=0.3].mean()))
