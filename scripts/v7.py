import pickle,numpy as onp,math,sys,json
from opt import *
name=sys.argv[1]; spec=json.loads(sys.argv[2])
v,s,n=pickle.load(open("v4_v3K2.pkl_-4.0.pkl","rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
order=sorted(range(n),key=lambda i:(math.floor(round(T[i,1],4)),round(T[i,0],4)))
ix={k+1:i for k,i in enumerate(order)}
d=onp.degrees(T[:,2])
if "corner" in spec:
    th=spec["corner"]; m=onp.where(d<-3)[0]; srt=m[onp.argsort(d[m])]; groups=[[srt[0]]]
    for i in srt[1:]:
        if d[i]-d[groups[-1][-1]]<0.6: groups[-1].append(i)
        else: groups.append([i])
    if spec.get("whole"): groups=[list(m)]
    for g in groups:
        g=onp.array(g); r=math.radians(th-d[g].mean()); c=T[g,:2].mean(0)
        R=onp.array([[math.cos(r),-math.sin(r)],[math.sin(r),math.cos(r)]]); T[g,:2]=(T[g,:2]-c)@R.T+c; T[g,2]+=r
for k,a in spec.get("set",{}).items(): T[ix[int(k)],2]=math.radians(a)
q=onp.concatenate([T[:,0],T[:,1],T[:,2]])
best=squeeze(q,n,s+0.12,log=lambda *a:None)
w=best[0]; ang=lambda k: round(math.degrees(w[2*n+ix[k]]),1)
print(name,"RESULT",round(best[1],5),"| 15,26,37:",[ang(k) for k in (15,26,37)],"13,14:",[ang(k) for k in (13,14)],"27,39,38,50:",[ang(k) for k in (27,39,38,50)],"41,53:",[ang(41),ang(53)],"tail 75-77:",[ang(k) for k in (75,76,77)],flush=True)
pickle.dump((best[0],best[1],n),open(f"v7_{name}.pkl","wb"))
