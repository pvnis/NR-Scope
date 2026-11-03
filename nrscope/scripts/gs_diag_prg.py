import csv,sys,numpy as np
g=list(csv.DictReader(open(sys.argv[1]+"/golden_sniffer_grants.csv")))
def st(r,k):
    d,n=r[k].split(":"); return float(d),int(n)
vis=[r for r in g if float(r["dmrs_coherence"])>=0.6 and int(r["nof_rb"])>=8]
print("visible decoded grants (DM-RS coherence >= 0.6, >= 8 PRBs): %d"%len(vis))
names=[("step_inside","inside a PRB"),("step_crb4","PRB boundary at CRB%4==0"),("step_crb2","PRB boundary at CRB%4==2"),("step_odd","PRB boundary at odd CRB")]
for k,nm in names:
    num=sum(st(r,k)[0]*st(r,k)[1] for r in vis if st(r,k)[1]>0); den=sum(st(r,k)[1] for r in vis if st(r,k)[1]>0)
    print("  %-26s mean |phase step| %5.1f deg over %d steps"%(nm,num/den,den))
# per grant: does it step at the 4-PRB boundaries?
ratio=[]
for r in vis:
    a,na=st(r,"step_crb4"); b,nb=st(r,"step_crb2"); c,nc=st(r,"step_odd")
    if na>=2 and nb+nc>=4: ratio.append(a-(b*nb+c*nc)/(nb+nc))
ratio=np.array(ratio)
print("per grant: excess step at CRB%%4 boundaries over the others: p10 %.0f / median %.0f / p90 %.0f deg"%tuple(np.percentile(ratio,[10,50,90])))
print("  grants with > 20 deg excess (a new precoder every 4 PRBs): %.0f%%; within +-8 deg (one precoder across the grant): %.0f%%"%(100*np.mean(ratio>20),100*np.mean(np.abs(ratio)<8)))
import collections as C
per=C.defaultdict(list)
for r in vis:
    a,na=st(r,"step_crb4"); b,nb=st(r,"step_crb2"); c,nc=st(r,"step_odd")
    if na>=2 and nb+nc>=4: per[r["rnti"]].append(a-(b*nb+c*nc)/(nb+nc)>20)
u=[(k,np.mean(v),len(v)) for k,v in per.items() if len(v)>=10]
print("UEs with >= 10 such grants: %d; share of their grants with a precoder every 4 PRBs, per UE:"%len(u), " ".join("%.0f%%"%(100*m) for k,m,n in sorted(u,key=lambda x:-x[1])))
nb=np.array([int(r["nof_rb"]) for r in vis]); ex=[]
for r in vis:
    a,na=st(r,"step_crb4"); b,n2=st(r,"step_crb2"); c,nc=st(r,"step_odd")
    ex.append(a-(b*n2+c*nc)/max(n2+nc,1) if na>=2 and n2+nc>=4 else np.nan)
ex=np.array(ex)
for lo,hi in ((8,32),(32,128),(128,274)):
    m=(nb>=lo)&(nb<hi)&np.isfinite(ex); print("  grants of %3d-%3d PRBs: %4d, with a precoder every 4 PRBs %.0f%%"%(lo,hi-1,m.sum(),100*np.mean(ex[m]>20)))
