import autograd.numpy as np, numpy as onp, math, sys, pickle, time
from autograd import value_and_grad
from scipy.optimize import minimize
from opt import pairs_of, pen_terms
src=sys.argv[1]; slide=int(sys.argv[2])   # slide=1: columns may slide; 0: fully rigid
v,s0,n=pickle.load(open(src,"rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
d=onp.degrees(T[:,2]); bm=(d>13)&(d<24); bi=onp.where(bm)[0]; fi=onp.where(~bm)[0]
th=T[bi,2].mean(); u=onp.array([math.cos(th),math.sin(th)]); w=onp.array([-math.sin(th),math.cos(th)])
pu=T[bi,:2]@u; pw=T[bi,:2]@w
c=onp.round(pu-pu.min()).astype(int); C=c.max()+1
k=onp.zeros(len(bi)); sig=onp.zeros(C)
for cc in range(C):
    m=c==cc; kk=onp.round(pw[m]-pw[m].min()); k[m]=kk; sig[cc]=(pw[m]-kk).mean()
pu0=(pu-c).mean()
o=pu0*u  # origin so that pos = o + c*u + (sig+k)*w
nf=len(fi); nb=len(bi)
def build(p):
    fx,fy,ft=p[:nf],p[nf:2*nf],p[2*nf:3*nf]
    ox,oy,ph=p[3*nf],p[3*nf+1],p[3*nf+2]
    sg=p[3*nf+3:3*nf+3+C] if slide else sig
    a=sg[c]+k
    bx=ox+c*np.cos(ph)-a*np.sin(ph); by=oy+c*np.sin(ph)+a*np.cos(ph)
    return np.concatenate([fx,bx]),np.concatenate([fy,by]),np.concatenate([ft,ph+0*bx])
p0=onp.concatenate([T[fi,0],T[fi,1],T[fi,2],[o[0],o[1],th],sig if slide else []])
def P(p,s,I,J):
    x,y,t=build(p); dep,wl=pen_terms(x,y,t,s,I,J); return np.sum(dep**2)+np.sum(wl**2)
def prs(p):
    x,y,t=build(p); return pairs_of(onp.array(x),onp.array(y))
def mv(p,s):
    I,J=prs(p); x,y,t=build(p); dep,wl=pen_terms(x,y,t,s,I,J); return float(max(onp.max(dep),onp.max(wl)))
def relax(p,s,iters=4000):
    I,J=prs(p); f=value_and_grad(lambda q:P(q,s,I,J))
    return minimize(f,p,jac=True,method='L-BFGS-B',options=dict(maxiter=iters,ftol=1e-30,gtol=1e-14,maxcor=30)).x
def joint(p,s,mu):
    I,J=prs(p); f=value_and_grad(lambda q:q[-1]+mu*P(q[:-1],q[-1],I,J))
    r=minimize(f,onp.append(p,s),jac=True,method='L-BFGS-B',options=dict(maxiter=4000,ftol=1e-30,gtol=1e-13,maxcor=30)).x
    return r[:-1],r[-1]
s=s0+0.15; p=relax(p0,s)
for mu in (1e2,1e3,1e4,1e5,1e6): p,s=joint(p,s,mu)
lo=s; hi=s+0.02; best=None
while best is None:
    q=relax(p,hi)
    if mv(q,hi)<1e-7: best=(q,hi)
    else: hi+=0.02
while hi-lo>1e-6:
    mid=(lo+hi)/2; q=relax(best[0],mid,6000)
    if mv(q,mid)<1e-7: hi=mid; best=(q,mid)
    else: lo=mid
q,sb=best; x,y,t=build(q)
print(src,"block tiles",nb,"columns",C,"slide" if slide else "fully rigid","-> RESULT",round(sb,5),"(bent version",round(s0,5),") block angle",round(math.degrees(q[3*nf+2]),2),flush=True)
# save in original tile order
X=onp.zeros(n);Y=onp.zeros(n);A=onp.zeros(n); idx=onp.concatenate([fi,bi]); X[idx]=x;Y[idx]=y;A[idx]=t
pickle.dump((onp.concatenate([X,Y,A]),sb,n),open(f"rig_{slide}_{src}","wb"))
