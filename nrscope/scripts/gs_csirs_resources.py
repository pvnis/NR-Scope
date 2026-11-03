import sys, collections as C, numpy as np
rows=[l.rstrip("\n").split(",") for l in list(open(sys.argv[1]+"/survey.csv"))[1:]]
prev=-1; wrap=0; occ=C.defaultdict(list); zp=C.defaultdict(list); slots=[]
MAP={0:"1 port",1:"CDM pair (even start)",2:"CDM pair (odd start)",3:"TRS density 3"}
for r in rows:
    if len(r)<40: continue
    a=int(r[0])*20+int(r[1])
    if prev>=0 and a<prev-10000: wrap+=20480
    prev=a; a+=wrap; slots.append(a)
    nz=set(); z=set()
    for e in r[39].split(";"):
        if not e: continue
        f=e.split(":")
        if f[0]=="z":
            if float(f[3])<0.2: z.add((int(f[1]),int(f[2])))
        elif float(f[3])>0.6: nz.add((int(f[0]),int(f[2]),int(f[1])))
    # group NZP by (symbol, mapping): the set of subcarriers
    g=C.defaultdict(set)
    for l,h,c in nz: g[(l,h)].add(c)
    for (l,h),cs in g.items(): occ[(l,h,tuple(sorted(cs)))].append(a)
    gz=C.defaultdict(set)
    for l,c in z: gz[l].add(c)
    for l,cs in gz.items():
        if len(cs)<=8: zp[(l,tuple(sorted(cs)))].append(a)
a0=min(slots); T=max(slots)-a0
def per(ts):
    d=np.diff(sorted(ts)); 
    if len(d)==0: return "-", "-"
    v=C.Counter(d).most_common(1)[0][0]
    return v, sorted(set(t%v for t in ts))[:3] if v>0 else "-"
print("survey spans %d slots (%.1f s); DL slots surveyed %d"%(T,T/2000,len(slots)))
print("\nNZP CSI-RS (sequence with n_ID = PCI found, coherence > 0.6):")
print("  %-8s %-24s %-28s %6s %8s %s"%("symbol","mapping","subcarriers in the PRB","times","period","slot offset"))
for k,v in sorted(occ.items(),key=lambda kv:-len(kv[1])):
    if len(v)<4: continue
    p,o=per(v); print("  %-8d %-24s %-28s %6d %8s %s"%(k[0],MAP[k[1]],",".join(map(str,k[2])),len(v),p,o))
print("\nREs left empty inside the PDSCH (zero-power CSI-RS / CSI-IM candidates, power < 0.2 of the PRB):")
for k,v in sorted(zp.items(),key=lambda kv:-len(kv[1]))[:12]:
    if len(v)<4: continue
    p,o=per(v); print("  symbol %-3d subcarriers %-22s %5d times, period %s, offsets %s"%(k[0],",".join(map(str,k[1])),len(v),p,o))
