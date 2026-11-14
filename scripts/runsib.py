import numpy as onp, pickle, sys, time
from seg import seg
from opt import *
f,n,S=sys.argv[1],int(sys.argv[2]),float(sys.argv[3])
T=seg(f"/mnt/user-data/uploads/{f}-image.jpg",S)
assert len(T)==n,len(T)
v=onp.concatenate([T[:,0],T[:,1],T[:,2]])
t0=time.time(); best=squeeze(v,n,S+0.06,log=lambda *a:None)
print(n,"RESULT side",round(best[1],6),"record",S,"time",round(time.time()-t0))
pickle.dump((best[0],best[1],n),open(f"p{n}.pkl","wb"))
