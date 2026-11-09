import sys, os, numpy as np, collections as C
REC=sys.argv[1]; saved=list(sys.argv)
src=open(os.path.expanduser("~/field/gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv=["gs_all.py",REC,"/dev/null"]; exec(src); sys.argv=saved
bd=REC+"/b0/"; Dd=[load(bd+"replay_ch%d/ue_dmrs.bin"%c) for c in CH]; keys=set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
S=C.defaultdict(list)
for k in sorted(keys):
    pm=np.logical_and.reduce([Dd[c][k][0] for c in range(NA)]); Y=np.array([Dd[c][k][1] for c in range(NA)]).astype(complex)
    if pm.sum() and np.all(np.isfinite(Y[:,np.repeat(pm,12)])): S[k[1]].append((k[0],k[2],pm,Y))
UA=max(S,key=lambda u:len(S[u])); A=sorted(S[UA],key=lambda x:(x[0],x[1]))
res=C.defaultdict(list)
for (a1,l1,p1,Y1),(a2,l2,p2,Y2) in zip(A[:-1],A[1:]):
    m=np.repeat(p1&p2,12)
    if m.sum()<120 or a2-a1>2: continue
    r=np.vdot(Y1[:,m],Y2[:,m]); res[(l1,l2,a1%2,a2-a1)].append(np.degrees(np.angle(r)))
for k,v in sorted(res.items(),key=lambda x:-len(x[1]))[:10]:
    v=np.array(v); print("l %d->%d slot parity %d dslot %d: n %d, mean phase %.1f, circ std %.1f"%(k[0],k[1],k[2],k[3],len(v),np.degrees(np.angle(np.mean(np.exp(1j*np.radians(v))))),np.degrees(np.sqrt(-2*np.log(abs(np.mean(np.exp(1j*np.radians(v)))))))))
