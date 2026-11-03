import sys, json, collections as C, numpy as np
O=sys.argv[1]
rows=[l.rstrip("\n").split(",") for l in list(open(O+"/survey.csv"))[1:]]
prev=-1; wrap=0; slots=[]; runs=[]
for r in rows:
    if len(r)<40: continue
    sfn,sl=int(r[0]),int(r[1]); a=sfn*20+sl
    if prev>=0 and a<prev-10000: wrap+=20480
    prev=a; a+=wrap
    trs=any(e.startswith(("4:2:3:","8:2:3:")) and float(e.split(":")[3])>0.6 for e in r[39].split(";"))
    slots.append([a,sfn,sl,int(trs)])
    f=[x for x in r if x.startswith("R")]
    if not f or len(f[0])<2: continue
    for e in f[0][1:].split(";"):
        q=e.split("/"); lo,hi=map(int,q[0].split("-"))
        comb={"0":0,"1":1,"b":2}[q[1]]; c0,c1,occ,p1,t2=map(float,q[2:7]); sy=int(q[7],16); ex=float(q[8])
        if trs: sy&=~((1<<4)|(1<<8))
        sy&=~(1<<13)  # symbol 13 carries the CSI-RS here
        ports=0 if occ<0.25 else (2 if occ>4 else 1)
        cdm=2 if (p1<0.3 or (comb==2)) else 1
        prec=3 if ex<-900 else (1 if ex>20 else (0 if abs(ex)<8 else 2))
        typ=1 if (max(c0,c1)-t2)>0.15 else 0
        runs.append([a,lo,hi,comb,int(round(100*max(c0,c1))),ports,cdm,sy,prec,typ])
a0=min(s[0] for s in slots)
for s in slots: s[0]-=a0
for r in runs: r[0]-=a0
# PRB-weighted summaries
W=lambda f: dict(C.Counter({k:0 for k in []}))
def summ(idx,labels):
    c=C.Counter()
    for r in runs: c[r[idx]]+=r[2]-r[1]+1
    return [[labels.get(k,str(k)),v] for k,v in sorted(c.items(),key=lambda kv:-kv[1])]
symname=lambda b:"+".join(str(l) for l in range(14) if b>>l&1)
D={"n_prb":273,"span":max(s[0] for s in slots)+1,"slots":slots,"runs":runs,
   "sum_ports":summ(5,{0:"1 layer, port 1000",1:"2 layers, ports 1000+1001",2:"port 1001 only"}),
   "sum_prec":summ(8,{0:"one precoder",1:"new precoder every 4 PRBs",2:"unclear",3:"run too short to tell"}),
   "sum_comb":summ(3,{0:"comb 0 (CDM group 0)",1:"comb 1 (CDM group 1)",2:"both combs"}),
   "sum_cdm":summ(6,{1:"1 (comb 1 carries data)",2:"2 (comb 1 empty or DM-RS)"}),
   "sum_typ":summ(9,{1:"type 1",0:"not clear"})}
c=C.Counter()
for r in runs: c[symname(r[7])]+=r[2]-r[1]+1
D["sum_sym"]=[[k or "2",v] for k,v in c.most_common(6)]
json.dump(D,open(O+"/dmrs_viz.json","w"),separators=(",",":"))
print("slots",len(slots),"runs",len(runs),"bytes",len(json.dumps(D,separators=(",",":"))))
for k in ("sum_ports","sum_prec","sum_comb","sum_cdm","sum_typ","sum_sym"): print(k,D[k])
