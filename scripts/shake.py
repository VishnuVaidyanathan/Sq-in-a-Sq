import sys,pickle,numpy as onp
from opt import *
f,seed,amp=sys.argv[1],int(sys.argv[2]),float(sys.argv[3])
v,s,n=pickle.load(open(f,"rb")); rng=onp.random.default_rng(seed)
q=v.copy(); q[:2*n]+=rng.normal(0,amp,2*n); q[2*n:]+=rng.normal(0,amp/2,n)
best=squeeze(q,n,s+0.05,log=lambda *a:None)
print(f,seed,amp,"->",round(best[1],5),flush=True)
pickle.dump((best[0],best[1],n),open(f"{f}.s{seed}","wb"))
