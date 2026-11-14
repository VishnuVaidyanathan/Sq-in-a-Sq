import pickle,numpy as onp,math,sys
from opt import *
src,dl=sys.argv[1],float(sys.argv[2])
v,s,n=pickle.load(open(src,"rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
m=(onp.degrees(T[:,2])>12)&(onp.degrees(T[:,2])<24)
c=T[m,:2].mean(0); r=math.radians(dl); R=onp.array([[math.cos(r),-math.sin(r)],[math.sin(r),math.cos(r)]])
T[m,:2]=(T[m,:2]-c)@R.T+c; T[m,2]+=r
q=onp.concatenate([T[:,0],T[:,1],T[:,2]])
best=squeeze(q,n,13.4,log=lambda *a:None)
print(src,"main block turned",dl,"-> RESULT",round(best[1],5),"(tiles in block",int(m.sum()),")",flush=True)
pickle.dump((best[0],best[1],n),open(f"v4_{src}_{dl}.pkl","wb"))
