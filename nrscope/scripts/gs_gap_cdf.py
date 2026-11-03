import sys, numpy as np
exec(open("/tmp/spanfill.py").read().split("def widest(m):")[0])
def gaps(m):
    i=np.nonzero(m)[0]
    if len(i)<2: return []
    g=[]; cur=0
    for v in m[i[0]:i[-1]+1]:
        if not v: cur+=1
        elif cur: g.append(cur); cur=0
    return g
srcs=[("DM-RS on air only",lambda a:air[a]),("decoded grants (DM-RS seen)",lambda a:gr.get(a,np.zeros(N,bool))),("both",lambda a:air[a]|gr.get(a,np.zeros(N,bool)))]
if TAB: srcs.append(("tablet grants only",lambda a:tb.get(a,np.zeros(N,bool))))
X=[1,2,3,4,6,8,12,16,24,32,48,64]
print("share of gaps (and of missing PRBs) with length <= x PRBs; slots with a span >= 137 PRBs")
print("| source | gaps | "+" | ".join("<=%d"%x for x in X)+" |")
print("|---|---|"+"---|"*len(X))
for name,get in srcs:
    G=[]
    for a in air:
        m=get(a); i=np.nonzero(m)[0]
        if len(i)<2 or i[-1]-i[0]+1<137: continue
        G+=gaps(m)
    G=np.array(G)
    if len(G)==0: print("| %s | 0 |"%name); continue
    cnt=" | ".join("%.0f%%"%(100*np.mean(G<=x)) for x in X)
    prb=" | ".join("%.0f%%"%(100*G[G<=x].sum()/G.sum()) for x in X)
    print("| %s | %d | %s |"%(name,len(G),cnt)); print("| (by missing PRBs) | %d PRBs | %s |"%(G.sum(),prb))
