# rigid (straight-row) block squeeze. usage: rb.py src name lo hi frame slide drot dx dy [extra json]
import autograd.numpy as np, numpy as onp, math, sys, pickle, json
from autograd import value_and_grad
from scipy.optimize import minimize
from opt import pairs_of, pen_terms
src,name=sys.argv[1],sys.argv[2]; lo,hi,frame=float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])
slide=int(sys.argv[6]); drot,dx,dy=map(float,sys.argv[7:10]); extra=json.loads(sys.argv[10]) if len(sys.argv)>10 else {}
v,s0,n=pickle.load(open(src,"rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
R0=pickle.load(open("v4_v3K2.pkl_-4.0.pkl","rb")); T0=onp.c_[R0[0][:n],R0[0][n:2*n],R0[0][2*n:]]
order=sorted(range(n),key=lambda i:(math.floor(round(T0[i,1],4)),round(T0[i,0],4))); ix={k+1:i for k,i in enumerate(order)}
for k,a in extra.get("set",{}).items(): T[ix[int(k)],2]=math.radians(a)
for k,(ax,ay) in extra.get("move",{}).items(): T[ix[int(k)],0]+=ax; T[ix[int(k)],1]+=ay
d=onp.degrees(T[:,2]); bm=(d>lo)&(d<hi)
for q in extra.get("free",[]): bm[ix[q]]=False
bi=onp.where(bm)[0]; fi=onp.where(~bm)[0]
th=T[bi,2].mean()+math.radians(frame); u=onp.array([math.cos(th),math.sin(th)]); w=onp.array([-math.sin(th),math.cos(th)])
pu=T[bi,:2]@u; pw=T[bi,:2]@w
c=onp.round(pu-pu.min()).astype(int); C=c.max()+1; k=onp.zeros(len(bi)); sig=onp.zeros(C)
for cc in range(C):
    m=c==cc; kk=onp.round(pw[m]-pw[m].min()); k[m]=kk; sig[cc]=(pw[m]-kk).mean()
o=((pu-c).mean())*u
fit=o+onp.outer(c,u)+onp.outer(sig[c]+k,w); print("fit max dev",round(float(onp.abs(fit-T[bi,:2]).max()),3),"rows",[int((c==q).sum()) for q in range(C)],file=sys.stderr)
# initial turn/shift of the block about its centroid
cen=T[bi,:2].mean(0); r=math.radians(drot); Rm=onp.array([[math.cos(r),-math.sin(r)],[math.sin(r),math.cos(r)]])
o=(o-cen)@Rm.T+cen+onp.array([dx,dy]); th0=th+r
nf=len(fi)
FIX=(math.radians(extra['fixang']+frame) if 'fixang' in extra else None)
if FIX is not None:
    r=FIX-th; Rm=onp.array([[math.cos(r),-math.sin(r)],[math.sin(r),math.cos(r)]]); o=(((pu-c).mean())*u-cen)@Rm.T+cen+onp.array([dx,dy]); th0=FIX
def build(p):
    fx,fy,ft=p[:nf],p[nf:2*nf],p[2*nf:3*nf]; ox,oy=p[3*nf],p[3*nf+1]; ph=(FIX if FIX is not None else p[3*nf+2])
    sg=p[3*nf+3:3*nf+3+C] if slide else sig
    a=sg[c]+k
    return np.concatenate([fx,ox+c*np.cos(ph)-a*np.sin(ph)]),np.concatenate([fy,oy+c*np.sin(ph)+a*np.cos(ph)]),np.concatenate([ft,ph+0*a])
p=onp.concatenate([T[fi,0],T[fi,1],T[fi,2],[o[0],o[1],th0],sig if slide else []])
def P(p,s,I,J):
    x,y,t=build(p); dep,wl=pen_terms(x,y,t,s,I,J); return np.sum(dep**2)+np.sum(wl**2)
def prs(p):
    x,y,t=build(p); return pairs_of(onp.array(x),onp.array(y))
def mv(p,s):
    I,J=prs(p); x,y,t=build(p); dep,wl=pen_terms(x,y,t,s,I,J); return float(max(onp.max(dep),onp.max(wl)))
s=s0+float(extra.get("room",0.25))
I,J=prs(p); f=value_and_grad(lambda q:P(q,s,I,J)); p=minimize(f,p,jac=True,method='L-BFGS-B',options=dict(maxiter=3000,ftol=1e-30,gtol=1e-14)).x
for mu in (1e2,1e3,1e4,1e5,1e6):
    for rep in range(2):
        I,J=prs(p); f=value_and_grad(lambda q:q[-1]+mu*P(q[:-1],q[-1],I,J))
        r_=minimize(f,onp.append(p,s),jac=True,method='L-BFGS-B',options=dict(maxiter=4000,ftol=1e-30,gtol=1e-13,maxcor=30)).x; p,s=r_[:-1],r_[-1]
x,y,t=build(p); idx=onp.concatenate([fi,bi]); X=onp.zeros(n);Y=onp.zeros(n);A=onp.zeros(n); X[idx]=x;Y[idx]=y;A[idx]=t-math.radians(frame)*bm[idx]
ang=lambda q: round(math.degrees(A[ix[q]]),1)
print(name,"| block",len(bi),"tiles, start turn",drot,"shift",(dx,dy),"->","SIDE",round(float(s),5),"viol %.0e"%mv(p,s),"| block angle",round(math.degrees(FIX if FIX is not None else p[3*nf+2])-frame,2),"| 28:",ang(28),"41:",ang(41),"16:",ang(16),"53:",ang(53),flush=True)
pickle.dump((onp.concatenate([X,Y,A]),float(s),n),open(f"rb_{name}.pkl","wb"))
