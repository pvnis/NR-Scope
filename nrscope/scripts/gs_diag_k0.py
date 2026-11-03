import csv,sys,numpy as np,collections as C
O=sys.argv[1]; N=273
M={}
for l in list(open(O+"/survey.csv"))[1:]:
    r=l.rstrip("\n").split(",")
    M[(int(r[0]),int(r[1]))]=np.array([int(c) for x in r[33] for c in bin(int(x,16))[2:].zfill(4)][:N],bool)
def nxt(k,d):
    a=k[0]*20+k[1]+d; return ((a//20)%1024,a%20)
g=list(csv.DictReader(open(O+"/golden_sniffer_grants.csv")))
for name,sel in (("invisible",lambda r:float(r["dmrs_coherence"])<0.4),("visible",lambda r:float(r["dmrs_coherence"])>=0.4)):
    out=[]
    for d in range(0,5):
        L=[]
        for r in g:
            if not sel(r): continue
            k=nxt((int(r["sfn"]),int(r["slot"])),d)
            if k not in M: continue
            m=np.zeros(N,bool)
            for x in r["prb_runs"].split(";"):
                if x: s,e=map(int,x.split("-")); m[s:e+1]=True
            o=M[k]; f=o.mean()
            if f<0.02 or f>0.98: continue
            L.append((o[m].mean()-f)/(1-f))
        out.append("k0=%d: lift %.2f (n %d)"%(d,np.mean(L),len(L)))
    print(name+": "+";  ".join(out))
