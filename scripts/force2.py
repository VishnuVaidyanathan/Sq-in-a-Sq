import pickle,numpy as onp,math,sys
from opt import *
v,s,n=pickle.load(open(sys.argv[1],"rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
order=sorted(range(n),key=lambda i:(math.floor(round(T[i,1],4)),round(T[i,0],4)))
idof={i:k+1 for k,i in enumerate(order)}
mu=1e5
q,sj=joint(v,n,s,mu); q,sj=joint(q,n,sj,mu)
x,y,t=q[:n],q[n:2*n],q[2*n:]
I,J=pairs_of(x,y); dep,w=pen_terms(x,y,t,sj,I,J); dep=onp.array(dep); w=onp.array(w).reshape(4,n)
F=2*mu*dep; Fw=2*mu*w   # forces; sum of right-wall forces should be ~0.5 each side (ds objective 1)
print("side",round(sj,5),"wall force totals L,R,T,B:",[round(float(Fw[k].sum()),3) for k in range(4)])
thr=0.01
load=set()
for a,b,f in zip(I,J,F):
    if f>thr: load.add(idof[a]); load.add(idof[b])
for k in range(4):
    for i in range(n):
        if Fw[k,i]>thr: load.add(idof[i])
d=onp.degrees(t)
corner=[idof[i] for i in range(n) if d[i]<-3]; main=[idof[i] for i in range(n) if d[i]>3]
print("load-bearing tiles:",len(load),"of",n)
print("corner block tiles carrying load:",sorted(set(corner)&load),"of",len(corner))
print("main block tiles carrying load:",len(set(main)&load),"of",len(main))
print("idle tiles (no load):",sorted(set(range(1,n+1))-load))
for nm,k in (("left",0),("right",1),("top",2),("bottom",3)):
    print(nm,"wall pressed by:",[(idof[i],round(float(Fw[k,i]),2)) for i in range(n) if Fw[k,i]>thr])
pickle.dump((sorted(load),[(idof[a],idof[b],float(f)) for a,b,f in zip(I,J,F) if f>thr]),open("force_"+sys.argv[2]+".pkl","wb"))
