import csv, sys, collections as C, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
O=sys.argv[1]; N=273
# survey: strict DM-RS map per slot, and per-symbol counts
occ={}; sym={}
for l in list(open(O+"/survey.csv"))[1:]:
    r=l.rstrip("\n").split(",")
    if len(r)<39: continue
    k=(int(r[0]),int(r[1]))
    occ[k]=np.array([int(c) for h in r[33] for c in bin(int(h,16))[2:].zfill(4)][:N],bool)
    sym[k]=np.array(list(map(int,r[4:17])))
G=C.defaultdict(list)
for r in csv.DictReader(open(O+"/golden_sniffer_grants.csv")):
    runs=[tuple(map(int,x.split("-"))) for x in r["prb_runs"].split(";") if x]
    G[(int(r["sfn"]),int(r["slot"]))].append((r["rnti"],runs,int(r["nof_rb"])))
# absolute slot index (unwrap the SFN)
keys=sorted(set(occ)|set(G)); ab={}; last=None; off=0
for k in sorted(occ, key=lambda k:(k[0],k[1])): pass
seq=[]
for k in [k for k in occ]: seq.append(k)
# order by arrival in survey file
seq=list(occ.keys()); prev=-1; wrap=0
for k in seq:
    a=k[0]*20+k[1]
    if prev>=0 and a<prev-10000: wrap+=1024*20
    ab[k]=a+wrap; prev=a
# window: 160 consecutive slots (80 ms) with the most grants
A=sorted((ab[k],k) for k in ab)
best=(0,0)
for i in range(0,len(A)-160,20):
    n=len(set(u for j in range(i,i+160) for u,_,_ in G.get(A[j][1],[]))); best=max(best,(n,i))
win=A[best[1]:best[1]+160]
ues=C.Counter(u for _,k in win for u,_,_ in G.get(k,[]))
cmap=plt.get_cmap("tab20"); col={u:cmap(i%20) for i,(u,_) in enumerate(ues.most_common())}
fig,ax=plt.subplots(figsize=(16,7))
a0=win[0][0]
for a,k in win:
    x=a-a0
    o=occ[k]
    for p in np.nonzero(o)[0]: ax.add_patch(Rectangle((x,p),1,1,color="0.82",lw=0))
    for u,runs,_ in G.get(k,[]):
        for s,e in runs: ax.add_patch(Rectangle((x+0.1,s),0.8,e-s+1,color=col[u],lw=0))
ax.set_xlim(0,max(a for a,_ in win)-a0+1); ax.set_ylim(0,N)
ax.set_xlabel("slot (30 kHz, 0.5 ms each); slots 4 and 9 of every 5 ms are uplink"); ax.set_ylabel("PRB (CRB)")
ax.set_title("Sunrise PCI 632, DL grants of decoded UEs (colour = UE, FDRA type 0, RBG 16) over PDSCH DM-RS seen on air (grey), the 100 ms with the most UEs")
h=[Rectangle((0,0),1,1,color=col[u]) for u,_ in ues.most_common(12)]+[Rectangle((0,0),1,1,color="0.82")]
ax.legend(h,[u for u,_ in ues.most_common(12)]+["DM-RS on air"],loc="upper right",fontsize=7,ncol=2)
plt.tight_layout(); plt.savefig(O+"/fdra_map.png",dpi=110)
# stats
sz=[n for v in G.values() for _,_,n in v]
nue=[len(v) for v in G.values()]
gaps=[]; share=[]; expl=[]
for k,v in G.items():
    if k in occ and occ[k].sum()>0:
        m=np.zeros(N,bool)
        for _,runs,_ in v:
            for s,e in runs: m[s:e+1]=True
        expl.append((m&occ[k]).sum()/occ[k].sum())
    share.append(sum(n for _,_,n in v)/N)
    if len(v)>=2:
        iv=sorted((runs[0][0],runs[-1][1]) for _,runs,_ in v if runs)
        for (s1,e1),(s2,e2) in zip(iv,iv[1:]): gaps.append(s2-e1-1)
fig,axs=plt.subplots(2,2,figsize=(13,8))
axs[0,0].hist(sz,bins=np.arange(0,290,16),color="C0"); axs[0,0].set_title("PRBs per grant (RBG = 16 PRBs)"); axs[0,0].set_xlabel("PRBs")
axs[0,1].hist(nue,bins=np.arange(0.5,12.5,1),color="C1"); axs[0,1].set_title("decoded UEs with a DL grant, per slot (slots with any)")
axs[1,0].hist(gaps,bins=np.arange(-50,140,8),color="C2"); axs[1,0].set_title("gap between neighbouring UEs' allocations in a slot (PRBs; <0 overlap)")
axs[1,1].hist(expl,bins=20,color="C3"); axs[1,1].set_title("share of the DM-RS seen on air that decoded grants explain, per slot")
plt.tight_layout(); plt.savefig(O+"/fdra_stats.png",dpi=100)
# reference signals per symbol, over the 20 ms SSB/TRS period
base=np.median([s[l-1] for s in sym.values() for l in (1,6,9,13)])
M=np.zeros((13,40)); n=np.zeros(40)
for k,s in sym.items():
    i=ab[k]%40; M[:,i]+=np.maximum(s-base,0); n[i]+=1
M/=np.maximum(n,1)
fig,ax=plt.subplots(figsize=(14,4.5))
im=ax.imshow(M,aspect="auto",origin="lower",extent=[-0.5,39.5,0.5,13.5],cmap="viridis")
ax.set_xlabel("slot within 20 ms"); ax.set_ylabel("OFDM symbol"); plt.colorbar(im,label="PRBs with DM-RS-like coherence above noise")
ax.set_title("Reference signals per symbol: PDSCH DM-RS on 2 and 11 (pos1), TRS (CSI-RS for tracking) on 4 and 8; blank columns = UL slots")
plt.tight_layout(); plt.savefig(O+"/rs_symbols.png",dpi=100)
print("grants %d in %d slots; PRBs/grant p10/50/90 %s; UEs per slot mean %.1f max %d"%(len(sz),len(G),np.percentile(sz,[10,50,90]),np.mean(nue),max(nue)))
print("slot band share by decoded grants: median %.2f; gaps between neighbours: p10/50/90 %s, overlap share %.2f, adjacent (gap 0) %.2f"%(np.median(share),np.percentile(gaps,[10,50,90]),np.mean(np.array(gaps)<0),np.mean(np.array(gaps)==0)))
print("DM-RS on air explained by decoded grants: median %.2f"%np.median(expl))
print("window: %d DL slots, %d UEs, grants per UE:"%(len(win),len(ues)),ues.most_common(8))

# sharing on the actual PRB sets
ov=0; tot=0; gp=[]; adj=0; trans=0; nseg=[]; perue=[]
for k,v in G.items():
    if len(v)<2: continue
    lab=np.full(N,-1); cnt=np.zeros(N,int)
    for i,(u,runs,_) in enumerate(v):
        for s,e in runs: lab[s:e+1]=i; cnt[s:e+1]+=1
    ov+=(cnt>1).sum(); tot+=(cnt>0).sum()
    perue.append([sum(e-s+1 for s,e in runs) for _,runs,_ in v])
    segs=[]; i=0
    while i<N:
        if lab[i]>=0:
            j=i
            while j+1<N and lab[j+1]==lab[i]: j+=1
            segs.append((i,j,lab[i])); i=j+1
        else: i+=1
    nseg.append(len(segs))
    for (s1,e1,l1),(s2,e2,l2) in zip(segs,segs[1:]):
        if l1!=l2: gp.append(s2-e1-1); trans+=1
fig,ax=plt.subplots(figsize=(7,4)); ax.hist(gp,bins=np.arange(0,130,8),color="C2")
ax.set_title("gap between one UE's PRBs and the next UE's (PRBs), slots with >=2 UEs"); plt.tight_layout(); plt.savefig(O+"/fdra_gaps.png",dpi=100)
print("slots with >=2 decoded UEs: %d; PRBs claimed by 2+ UEs (MU-MIMO or error): %.1f%%"%(len(nseg),100*ov/max(tot,1)))
print("UE-to-UE transitions: %d, touching (gap 0) %.0f%%, gap p50 %d p90 %d PRBs; contiguous pieces per slot median %d"%(trans,100*np.mean(np.array(gp)==0),np.median(gp),np.percentile(gp,90),np.median(nseg)))
