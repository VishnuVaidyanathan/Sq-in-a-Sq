import pickle,numpy as onp,math,sys,time
from opt import *
mode=sys.argv[1]
v1,s1,_=pickle.load(open("p131.pkl","rb")); n1=131; A=onp.c_[v1[:n1],v1[n1:2*n1],v1[2*n1:]]
v5,s5,_=pickle.load(open("p155.pkl","rb")); n5=155; B=onp.c_[v5[:n5],v5[n5:2*n5],v5[2*n5:]]
c1=A[onp.degrees(A[:,2])<-1]; c5i=onp.where(onp.degrees(B[:,2])<-1)[0]; c5=B[c5i]
# match
extra=[]
for t in c1:
    d=onp.hypot(c5[:,0]-t[0],c5[:,1]-t[1]).min()
    if d>0.3: extra.append(t)
extra=onp.array(extra); print("131 corner",len(c1),"155 corner",len(c5),"extra tiles",len(extra),extra.round(2).tolist())
# overlap depth of each 155 tile with extras
def depth(a,b):
    I=onp.array([0]);J=onp.array([1])
    dep,_=pen_terms(onp.array([a[0],b[0]]),onp.array([a[1],b[1]]),onp.array([a[2],b[2]]),99.,I,J); return float(dep[0])
ov=onp.array([max(depth(b,e) for e in extra) for b in B])
main=onp.degrees(B[:,2])>1
cen=extra[:,:2].mean(0)
if mode=='a': rem=onp.argsort(-ov)[:len(extra)]
elif mode=='b':
    idx=onp.where(main)[0]; rem=idx[onp.argsort(onp.hypot(B[idx,0]-cen[0],B[idx,1]-cen[1]))[:len(extra)]]
elif mode=='c':
    idx=onp.where(~main & (onp.degrees(B[:,2])>-1))[0]; rem=idx[onp.argsort(-ov[idx])[:len(extra)]]
print("mode",mode,"removing",[(round(B[i,0],2),round(B[i,1],2),round(math.degrees(B[i,2]),1),round(ov[i],2)) for i in rem])
keep=onp.ones(n5,bool); keep[rem]=False
N=onp.vstack([B[keep],extra]); n=len(N); assert n==155
v=onp.concatenate([N[:,0],N[:,1],N[:,2]])
t0=time.time()
best=squeeze(v,n,13.25,log=lambda *a:print(mode,*a,flush=True))
print(mode,"RESULT side",round(best[1],6),"vs my baseline",round(s5,6),"record 12.958447","time",round(time.time()-t0))
pickle.dump((best[0],best[1],n),open(f"v1{mode}.pkl","wb"))
