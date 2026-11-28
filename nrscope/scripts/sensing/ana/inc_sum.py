import json,glob,numpy as np,collections as C
R=C.defaultdict(list)
for f in glob.glob('/home/sens/isac_data/ana/inc_[0-9]*.json'):
    for k,v in json.load(open(f)).items(): R[k]+=v
q=lambda k:("%5.2f"%np.median(R[k])) if R.get(k) else "  -  "
print("windows",len(R['truth|width']),"truth width",q('truth|width'))
for c in ("split2","split4","gap"):
  for s in ("0.0","0.15","0.5"):
    for m in ("oracle","edge","raw","inc","single"):
      if R.get(f"ctl|{c}|{s}|{m}|corr"): print(f"{c:7s}{s:5s}{m:7s}"," ".join(q(f"ctl|{c}|{s}|{m}|{kk}") for kk in ("corr","dyn","recall","false","width")))
print("real windows",len(R['real|segments']),"segments",q('real|segments'))
for m in ("inc_real","inc_tab"): print(m," ".join(q(f"real|{m}|{kk}") for kk in ("corr","dyn","recall","false","width")))
