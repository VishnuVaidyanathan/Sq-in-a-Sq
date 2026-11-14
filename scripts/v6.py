import pickle,numpy as onp,math,sys
from opt import *
th=float(sys.argv[1])
v,s,n=pickle.load(open("v4_v3K2.pkl_-4.0.pkl","rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
d=onp.degrees(T[:,2]); m=onp.where((d>13)&(d<24))[0]
# cluster by angle
srt=m[onp.argsort(d[m])]; groups=[[srt[0]]]
for i in srt[1:]:
    if d[i]-d[groups[-1][-1]]<0.25: groups[-1].append(i)
    else: groups.append([i])
for g in groups:
    g=onp.array(g); r=math.radians(th-d[g].mean()); c=T[g,:2].mean(0)
    R=onp.array([[math.cos(r),-math.sin(r)],[math.sin(r),math.cos(r)]])
    T[g,:2]=(T[g,:2]-c)@R.T+c; T[g,2]=math.radians(th)
q=onp.concatenate([T[:,0],T[:,1],T[:,2]])
best=squeeze(q,n,s+0.12,log=lambda *a:None)
w=best[0]; dd=onp.degrees(w[2*n:]); mm=dd[(dd>13)&(dd<24)]
print("all main tiles set to",th,"groups",[len(g) for g in groups],"-> RESULT",round(best[1],5),"final main angles min/max",round(mm.min(),1),round(mm.max(),1),"distinct",sorted(set(onp.round(mm,1).tolist())),flush=True)
pickle.dump((best[0],best[1],n),open(f"v6_{th}.pkl","wb"))
