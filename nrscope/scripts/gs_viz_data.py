import csv, sys, json, collections as C, numpy as np
O=sys.argv[1]; N=273
slots=[]; prev=-1; wrap=0; ab={}; symsum=np.zeros((13,40)); symn=np.zeros(40)
rows=[l.rstrip("\n").split(",") for l in list(open(O+"/survey.csv"))[1:]]
base=np.median([int(r[4+l-1]) for r in rows if len(r)>=39 for l in (1,6,9,13)])
for r in rows:
    if len(r)<39: continue
    sfn,sl=int(r[0]),int(r[1]); a=sfn*20+sl
    if prev>=0 and a<prev-10000: wrap+=1024*20
    a+=wrap; prev=a; ab[(sfn,sl)]=a
    c0=np.array(list(map(int,r[4:17])))
    slots.append([a,sfn,sl,r[33],int(c0[1]>base+40),int(c0[10]>base+40),int(c0[3]>base+40 and c0[7]>base+40)])
    symsum[:,a%40]+=np.maximum(c0-base,0); symn[a%40]+=1
a0=min(s[0] for s in slots)
for s in slots: s[0]-=a0
dci={}
for r in csv.DictReader(open(O+"/golden_sniffer.csv")):
    if r["payload_hex"] and int(r["payload_hex"][0],16)>>3==1: dci[(r["sfn"],r["slot"],r["rnti"])]=(int(r["size"]),int(r["al"]))
ue_idx={}; ues=[]; grants=[]
for r in csv.DictReader(open(O+"/golden_sniffer_grants.csv")):
    k=(int(r["sfn"]),int(r["slot"]))
    if k not in ab: continue
    u=r["rnti"]
    if u not in ue_idx: ue_idx[u]=len(ues); ues.append({"rnti":u,"n":0,"prb":0})
    i=ue_idx[u]; ues[i]["n"]+=1; ues[i]["prb"]+=int(r["nof_rb"])
    runs=[list(map(int,x.split("-"))) for x in r["prb_runs"].split(";") if x]
    sz,al=dci.get((r["sfn"],r["slot"],u),(0,0))
    grants.append([ab[k]-a0,i,runs,int(r["nof_rb"]),int(r["ports"]),round(float(r["dmrs_coherence"]),2),sz,al,int(r["dmrs_add"])])
log=open(O+"/nrscope.log",errors="replace").read()
for u in ues:
    u["own"]= ("UE %s FDRA: type"%u["rnti"] in log and "UE %s FDRA: type 0 rbg config 1, BWP 273 PRBs at CRB 0, field at bit"%u["rnti"] in log and "taken from" not in [l for l in log.splitlines() if "UE %s FDRA"%u["rnti"] in l][0]) if any("UE %s FDRA"%u["rnti"] in l for l in log.splitlines()) else None
ues_sorted=sorted(range(len(ues)),key=lambda i:-ues[i]["n"])
rank={old:new for new,old in enumerate(ues_sorted)}
ues=[ues[i] for i in ues_sorted]
for g in grants: g[1]=rank[g[1]]
# stats
sz=[g[3] for g in grants]; per=C.defaultdict(list)
for g in grants: per[g[0]].append(g)
gaps=[]; nue=[]
for t,v in per.items():
    nue.append(len(v))
    if len(v)<2: continue
    lab=np.full(N,-1)
    for g in v:
        for s,e in g[2]: lab[s:e+1]=g[1]
    segs=[];i=0
    while i<N:
        if lab[i]>=0:
            j=i
            while j+1<N and lab[j+1]==lab[i]: j+=1
            segs.append((i,j,lab[i])); i=j+1
        else: i+=1
    for (s1,e1,l1),(s2,e2,l2) in zip(segs,segs[1:]):
        if l1!=l2: gaps.append(s2-e1-1)
M=(symsum/np.maximum(symn,1)).round(1).tolist()
D={"n_prb":N,"slots":slots,"grants":grants,"ues":ues,"sym":M,
   "hist_prb":np.histogram(sz,bins=np.arange(0,289,16))[0].tolist(),
   "hist_gap":np.histogram(gaps,bins=np.arange(0,137,8))[0].tolist(),
   "hist_nue":np.bincount(nue,minlength=7)[1:7].tolist(),
   "n_slots":len(slots),"span":max(s[0] for s in slots)+1}
json.dump(D,open(O+"/viz.json","w"),separators=(",",":"))
print("slots",len(slots),"grants",len(grants),"ues",len(ues),"gaps",len(gaps),"bytes",len(json.dumps(D,separators=(",",":"))))

# CSI-RS per slot: TRS (sym 4/8 from sc 2), NZP CSI-RS on sym 13 sc 0-7, CSI-IM (empty 8-11 on 12/13)
csi_rows={}
for r in rows:
    if len(r)<40: continue
    k=(int(r[0]),int(r[1])); trs=0; csi=0; zp=0
    for e in r[39].split(";"):
        f=e.split(":")
        if len(f)<4: continue
        if f[0]=="z":
            if f[1] in ("12","13") and 8<=int(f[2])<=11 and float(f[3])<0.2: zp+=1
        elif float(f[3])>0.6:
            if f[0] in ("4","8") and f[1]=="2" and f[2]=="3": trs=1
            if f[0]=="13" and int(f[1])<=7 and f[2]=="1": csi+=1
    csi_rows[k]=(trs,int(csi>=4),int(zp>=3))
for s in D["slots"]:
    t=csi_rows.get((s[1],s[2]),(0,0,0)); s[6]=t[0]; s+= [t[1],t[2]]
D["csi_table"]=[
 ["TRS (tracking CSI-RS)","4 and 8","2, 6, 10 (every 4th)","1 port, density 3","slots 5 and 6 of every 40 (20 ms)"],
 ["NZP CSI-RS for CSI","13","0 to 7 (4 CDM pairs)","8 ports, frequency CDM2","mostly every 20 slots, offsets 1, 10, 11"],
 ["CSI-IM","12","8 to 11","4 empty REs","every 10 slots, offsets 1 and 2"],
 ["CSI-IM","13","8 to 11","4 empty REs","every 80 slots, offsets 0, 10, 20"],
 ["Empty REs beside the TRS","4 and 8","1, 5, 9","rate matching","with the TRS, slots 5 and 6 of every 40"]]
json.dump(D,open(O+"/viz.json","w"),separators=(",",":"))
print("csi flags: trs %d csi %d zp %d"%tuple(sum(s[i] for s in D["slots"]) for i in (6,7,8)))
