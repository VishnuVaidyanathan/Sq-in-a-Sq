import autograd.numpy as np, numpy as onp, math, time
from autograd import value_and_grad
from scipy.optimize import minimize
def pairs_of(x,y,r=1.75):
    n=len(x); d=onp.hypot(x[:,None]-x[None,:],y[:,None]-y[None,:])
    I,J=onp.where(onp.triu(d<r,1)); return I,J
def pen_terms(x,y,t,s,I,J):
    dx=x[J]-x[I]; dy=y[J]-y[I]
    ci,si=np.cos(t[I]),np.sin(t[I]); cj,sj=np.cos(t[J]),np.sin(t[J])
    dt=t[I]-t[J]; hh=0.5*(np.abs(np.cos(dt))+np.abs(np.sin(dt)))+0.5
    o1=hh-np.abs(dx*ci+dy*si); o2=hh-np.abs(-dx*si+dy*ci)
    o3=hh-np.abs(dx*cj+dy*sj); o4=hh-np.abs(-dx*sj+dy*cj)
    dep=np.minimum(np.minimum(o1,o2),np.minimum(o3,o4))
    dep=np.maximum(dep,0.0)
    h=0.5*(np.abs(np.cos(t))+np.abs(np.sin(t)))
    w=np.concatenate([np.maximum(h-x,0),np.maximum(x+h-s,0),np.maximum(h-y,0),np.maximum(y+h-s,0)])
    return dep,w
def P(v,n,s,I,J):
    x,y,t=v[:n],v[n:2*n],v[2*n:3*n]
    dep,w=pen_terms(x,y,t,s,I,J)
    return np.sum(dep**2)+np.sum(w**2)
def maxviol(v,n,s):
    x,y,t=v[:n],v[n:2*n],v[2*n:3*n]
    I,J=pairs_of(onp.array(x),onp.array(y))
    dep,w=pen_terms(x,y,t,s,I,J)
    return float(max(onp.max(dep) if len(dep) else 0,onp.max(w)))
def relax(v,n,s,iters=3000,gtol=1e-14):
    I,J=pairs_of(v[:n],v[n:2*n])
    f=value_and_grad(lambda q:P(q,n,s,I,J))
    r=minimize(f,v,jac=True,method='L-BFGS-B',options=dict(maxiter=iters,ftol=1e-30,gtol=gtol,maxcor=30))
    return r.x,r.fun
def joint(v,n,s,mu,iters=4000):
    I,J=pairs_of(v[:n],v[n:2*n])
    def F(q): return q[-1]+mu*P(q[:-1],n,q[-1],I,J)
    f=value_and_grad(F)
    r=minimize(f,onp.append(v,s),jac=True,method='L-BFGS-B',options=dict(maxiter=iters,ftol=1e-30,gtol=1e-13,maxcor=30))
    return r.x[:-1],r.x[-1]
def squeeze(v,n,s,log=print,tol=1e-6,feas=1e-7):
    v,p=relax(v,n,s); log(f"start s={s:.5f} pen={p:.3e} maxviol={maxviol(v,n,s):.2e}")
    for mu in (1e2,1e3,1e4,1e5,1e6):
        v,s=joint(v,n,s,mu); log(f"mu={mu:.0e} s={s:.6f} maxviol={maxviol(v,n,s):.2e}")
    # feasibility bisection
    lo=s; hi=s+0.02; best=None
    while True:
        vv,p=relax(v,n,hi)
        if maxviol(vv,n,hi)<feas: best=(vv,hi); break
        hi+=0.02
    while hi-lo>tol:
        mid=(lo+hi)/2; vv,p=relax(best[0],n,mid,iters=6000)
        if maxviol(vv,n,mid)<feas: hi=mid; best=(vv,mid)
        else: lo=mid
    return best
