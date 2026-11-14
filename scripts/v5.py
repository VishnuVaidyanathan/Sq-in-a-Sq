import pickle,numpy as onp,math,sys,json
from opt import *
name=sys.argv[1]; spec=json.loads(sys.argv[2])   # {"set":{id:deg}, "add":{id:deg}}
v,s,n=pickle.load(open("v4_v3K2.pkl_-4.0.pkl","rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
order=sorted(range(n),key=lambda i:(math.floor(round(T[i,1],4)),round(T[i,0],4)))
ix={k+1:i for k,i in enumerate(order)}
for k,d in spec.get("set",{}).items(): T[ix[int(k)],2]=math.radians(d)
for k,d in spec.get("add",{}).items(): T[ix[int(k)],2]+=math.radians(d)
q=onp.concatenate([T[:,0],T[:,1],T[:,2]])
best=squeeze(q,n,s+0.15,log=lambda *a:None)
w=best[0]; ang=lambda k: round(math.degrees(w[2*n+ix[k]]),1)
print(name,"RESULT",round(best[1],5),"| final angles 78-81:",[ang(k) for k in (78,79,80,81)],"92-94:",[ang(k) for k in (92,93,94)],flush=True)
pickle.dump((best[0],best[1],n),open(f"v5_{name}.pkl","wb"))
