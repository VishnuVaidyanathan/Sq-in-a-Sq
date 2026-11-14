# all tiles free, except: 13 and 14 held at a given angle, and 14's upper-left edge must touch tile 1's bottom-right corner.
# optional: tip of 27 held in the 28/41 joint, delta past the wall line.
import autograd.numpy as np, numpy as onp, math, sys, pickle, json
from autograd import value_and_grad
from scipy.optimize import minimize
from opt import pairs_of, pen_terms
src,name,ang1314=sys.argv[1],sys.argv[2],float(sys.argv[3]); ex=json.loads(sys.argv[4]) if len(sys.argv)>4 else {}
delta=ex.get("pin27",None)
v,s0,n=pickle.load(open(src,"rb")); T=onp.c_[v[:n],v[n:2*n],v[2*n:]]
R0=pickle.load(open("v4_v3K2.pkl_-4.0.pkl","rb")); T0=onp.c_[R0[0][:n],R0[0][n:2*n],R0[0][2*n:]]
order=sorted(range(n),key=lambda i:(math.floor(round(T0[i,1],4)),round(T0[i,0],4))); ix={k+1:i for k,i in enumerate(order)}
I1,I13,I14,I27,I28,I41=[ix[k] for k in (1,13,14,27,28,41)]
A=math.radians(ang1314)
# start: turn 13 and 14 together about their common centre, lift them
pr=[I13,I14]; c=T[pr,:2].mean(0); r=A-T[pr,2].mean(); Rm=onp.array([[math.cos(r),-math.sin(r)],[math.sin(r),math.cos(r)]])
T[pr,:2]=(T[pr,:2]-c)@Rm.T+c; T[pr,2]=A; T[pr,1]-=ex.get("lift",0.08)
fixed_t=onp.zeros(n,bool); fixed_t[pr]=True
def build(p,s):
    x=p[:n]; y=p[n:2*n]; tf=p[2*n:3*n]
    t=np.where(fixed_t,A,tf)
    # 14 touches the bottom-right corner of 1 with its upper-left edge
    px=x[I1]+0.5; py=y[I1]+0.5
    y14=py-(-0.5+(px-x[I14])*np.sin(A))/np.cos(A)
    sel14=(onp.arange(n)==I14)
    y=np.where(sel14,y14,y)
    if delta is not None:
        a=t[I27]
        ybl28=y[I28]-0.5*np.sin(t[I28])+0.5*np.cos(t[I28]); ytl41=y[I41]-0.5*np.sin(t[I41])-0.5*np.cos(t[I41])
        x27=s-9.0+delta-0.5*(np.cos(a)-np.sin(a)); y27=0.5*(ybl28+ytl41)-0.5*(np.sin(a)+np.cos(a))
        sel27=(onp.arange(n)==I27); x=np.where(sel27,x27,x); y=np.where(sel27,y27,y)
    return x,y,t
def P(p,s,I,J):
    x,y,t=build(p,s); dep,wl=pen_terms(x,y,t,s,I,J); return np.sum(dep**2)+np.sum(wl**2)
def prs(p,s):
    x,y,t=build(p,s); return pairs_of(onp.array(x),onp.array(y))
p=onp.concatenate([T[:,0],T[:,1],T[:,2]]); s=s0+ex.get("room",0.10)
I,J=prs(p,s); f=value_and_grad(lambda q:P(q,s,I,J)); p=minimize(f,p,jac=True,method='L-BFGS-B',options=dict(maxiter=3000,ftol=1e-30,gtol=1e-14)).x
for mu in (1e2,1e3,1e4,1e5,1e6):
    for rep in range(2):
        I,J=prs(p,s); f=value_and_grad(lambda q:q[-1]+mu*P(q[:-1],q[-1],I,J))
        r_=minimize(f,onp.append(p,s),jac=True,method='L-BFGS-B',options=dict(maxiter=4000,ftol=1e-30,gtol=1e-13,maxcor=30)).x; p,s=r_[:-1],r_[-1]
x,y,t=[onp.array(q) for q in build(p,s)]; I,J=prs(p,s); dep,wl=pen_terms(x,y,t,s,I,J); viol=float(max(onp.max(dep),onp.max(wl)))
D=onp.degrees(t); g=lambda q: round(float(D[ix[q]]),1)
top14=y[I14]-0.5*(abs(math.sin(A))+abs(math.cos(A)))
print(name,"| 13,14 held at",ang1314,("| 27 in joint by %s"%delta if delta is not None else "| 27 free"),"-> SIDE",round(float(s),5),"viol %.0e"%viol,"| top of 14 at y=%.3f"%top14,"| 15,26,37:",[g(15),g(26),g(37)],"27:",g(27),"41:",g(41),"2:",g(2),"16:",g(16),flush=True)
pickle.dump((onp.concatenate([x,y,t]),float(s),n),open(f"touch_{name}.pkl","wb"))
