import pickle,numpy as onp,math,sys
from opt import *
mode=sys.argv[1]
v1,s1,_=pickle.load(open("p131.pkl","rb")); n1=131; A=onp.c_[v1[:n1],v1[n1:2*n1],v1[2*n1:]]
v5,s5,_=pickle.load(open("p155.pkl.s1","rb")); n5=155; B=onp.c_[v5[:n5],v5[n5:2*n5],v5[2*n5:]]
def idx(T,x,y):
    d=onp.hypot(T[:,0]-x,T[:,1]-y); i=int(d.argmin()); assert d[i]<0.15,(x,y,d[i]); return i
extras=[A[idx(A,*p)] for p in ((3.27,6.16),(1.68,6.25),(2.55,6.90))]
wedge=A[idx(A,2.50,8.09)].copy()
rem=[idx(B,*p) for p in ((1.51,6.41),(3.27,6.40),(2.92,7.34),(2.63,8.29))]
B2=B.copy()
if mode=='K2':
    i=idx(B,4.23,6.67); B2[i]=A[idx(A,4.46,6.50)]
keep=onp.ones(n5,bool); keep[rem]=False
N=onp.vstack([B2[keep]]+[e[None] for e in extras]+[wedge[None]]); n=len(N); assert n==155
v=onp.concatenate([N[:,0],N[:,1],N[:,2]])
best=squeeze(v,n,13.10,log=lambda *a:None)
print(mode,"RESULT",round(best[1],5),flush=True)
pickle.dump((best[0],best[1],n),open(f"v3{mode}.pkl","wb"))
