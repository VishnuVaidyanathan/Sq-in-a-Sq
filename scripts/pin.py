# free packing, one condition: right tip of tile 27 sits  delta  past the flat wall (x = s-9+delta)
import autograd.numpy as np, numpy as onp, math, sys, pickle, json
from autograd import value_and_grad
from scipy.optimize import minimize
from opt import pairs_of, pen_terms
src,name,delta=sys.argv[1],sys.argv[2],float(sys.argv[3]); ex=json.loads(sys.argv[4]) if len(sys.argv)>4 else {}
v,s0,n=pickle.load(open(src,"rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
R0=pickle.load(open("v4_v3K2.pkl_-4.0.pkl","rb")); T0=onp.c_[R0[0][:n],R0[0][n:2*n],R0[0][2*n:]]
order=sorted(range(n),key=lambda i:(math.floor(round(T0[i,1],4)),round(T0[i,0],4))); ix={k+1:i for k,i in enumerate(order)}
d=onp.degrees(T[:,2]); cb=onp.where(d<-13)[0]
T[cb,1]+=ex.get("dy",0.0); T[cb,0]+=ex.get("dx",0.0)
if "turn" in ex:
    r=math.radians(ex["turn"]); c=T[cb,:2].mean(0); Rm=onp.array([[math.cos(r),-math.sin(r)],[math.sin(r),math.cos(r)]]); T[cb,:2]=(T[cb,:2]-c)@Rm.T+c; T[cb,2]+=r
for k,a in ex.get("set",{}).items(): T[ix[int(k)],2]=math.radians(a)
I28=ix[28]; I41=ix[41]
p27=ix[27]; mask=onp.ones(n,bool); mask[p27]=False; oth=onp.where(mask)[0]
i28=None;i41=None
def build(p,s):
    xo=p[:n-1]; yo=p[n-1:2*n-2]; t=p[2*n-2:3*n-2]
    a=t[p27]
    # x of every tile except 27, y of every tile except 27
    def at(arr,i): return arr[i if i<p27 else i-1]
    a28=t[I28]; a41=t[I41]
    ybl28=at(yo,I28)-0.5*np.sin(a28)+0.5*np.cos(a28)     # bottom-left corner of 28
    ytl41=at(yo,I41)-0.5*np.sin(a41)-0.5*np.cos(a41)     # top-left corner of 41
    ymid=0.5*(ybl28+ytl41)
    x27=s-9.0+delta-0.5*(np.cos(a)-np.sin(a))            # right tip of 27 is delta past the wall line
    y27=ymid-0.5*(np.sin(a)+np.cos(a))                   # and sits exactly in the 28/41 joint
    x=np.concatenate([xo[:p27],np.reshape(x27,(1,)),xo[p27:]])
    y=np.concatenate([yo[:p27],np.reshape(y27,(1,)),yo[p27:]])
    return x,y,t
def P(p,s,I,J):
    x,y,t=build(p,s); dep,wl=pen_terms(x,y,t,s,I,J); return np.sum(dep**2)+np.sum(wl**2)
def prs(p,s):
    x,y,t=build(p,s); return pairs_of(onp.array(x),onp.array(y))
p=onp.concatenate([T[oth,0],T[oth,1],T[:,2]]); s=s0+ex.get("room",0.12)
I,J=prs(p,s); f=value_and_grad(lambda q:P(q,s,I,J)); p=minimize(f,p,jac=True,method='L-BFGS-B',options=dict(maxiter=3000,ftol=1e-30,gtol=1e-14)).x
for mu in (1e2,1e3,1e4,1e5,1e6):
    for rep in range(2):
        I,J=prs(p,s); f=value_and_grad(lambda q:q[-1]+mu*P(q[:-1],q[-1],I,J))
        r_=minimize(f,onp.append(p,s),jac=True,method='L-BFGS-B',options=dict(maxiter=4000,ftol=1e-30,gtol=1e-13,maxcor=30)).x; p,s=r_[:-1],r_[-1]
x,y,t=build(p,s); x=onp.array(x); I,J=prs(p,s); dep,wl=pen_terms(x,y,t,s,I,J); viol=float(max(onp.max(dep),onp.max(wl)))
A=onp.degrees(onp.array(t)); ang=lambda q: round(float(A[ix[q]]),1)
cbA=A[cb]
print(name,"| 27 pushed in by",delta,ex,"-> SIDE",round(float(s),5),"viol %.0e"%viol,"| gap 28-41 joint %.3f"%float((onp.array(y)[I41]-0.5*math.sin(t[I41])-0.5*math.cos(t[I41]))-(onp.array(y)[I28]-0.5*math.sin(t[I28])+0.5*math.cos(t[I28]))),"| 27:",ang(27),"28:",ang(28),"41:",ang(41),"53:",ang(53),"65:",ang(65),"| corner block angles",round(float(cbA.min()),1),"to",round(float(cbA.max()),1),flush=True)
pickle.dump((onp.concatenate([x,onp.array(y),onp.array(t)]),float(s),n),open(f"pin_{name}.pkl","wb"))
