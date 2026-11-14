import pickle,numpy as onp,math,sys,time
from opt import *
v1,s1,_=pickle.load(open("p131.pkl","rb")); n1=131; A=onp.c_[v1[:n1],v1[n1:2*n1],v1[2*n1:]]
v5,s5,_=pickle.load(open("p155.pkl","rb")); n5=155; B=onp.c_[v5[:n5],v5[n5:2*n5],v5[2*n5:]]
c1=A[onp.degrees(A[:,2])<-1]; c5=B[onp.degrees(B[:,2])<-1]
extra=onp.array([t for t in c1 if onp.hypot(c5[:,0]-t[0],c5[:,1]-t[1]).min()>0.3])
deg=onp.degrees(B[:,2]); mi=onp.where((deg>14)&(deg<21))[0]
th=B[mi,2].mean(); u=onp.array([math.cos(th),math.sin(th)]); w=onp.array([-math.sin(th),math.cos(th)])
pu=B[mi,:2]@u; pw=B[mi,:2]@w
def rows(p,q):
    r=onp.round(p-p.min()).astype(int); out=[]
    for k in range(r.max()+1):
        ids=mi[r==k]; out.append(ids[onp.argsort(q[r==k])])
    return out
Ru=rows(pu,pw); Rw=rows(pw,pu)
if __name__=="__main__":
    print("main tiles",len(mi),"angle",round(math.degrees(th),2))
    print("rows by u:",[len(r) for r in Ru]); print("rows by w:",[len(r) for r in Rw])
    for nm,R in (("u",Ru),("w",Rw)):
        print(nm,[(round(B[r[0],0],1),round(B[r[0],1],1),round(B[r[-1],0],1),round(B[r[-1],1],1)) for r in R[:3]])
    print("extra",extra.round(2).tolist())
