import sys,time,pickle
from v2 import *
mode=sys.argv[1]
if mode=='A': rem=[Ru[0][0],Ru[1][0],Ru[2][0]]; add=extra
elif mode=='B': rem=[Ru[0][-1],Ru[1][-1],Ru[2][-1]]; add=extra
elif mode=='C': rem=[Ru[0][0]]; add=extra[[0]]
elif mode=='E': rem=list(Ru[0]); add=extra
elif mode=='G': rem=[Ru[0][0],Ru[1][0]]; add=extra[[0,2]]
elif mode=='H': rem=[Ru[0][0],Ru[0][1],Ru[1][0]]; add=extra
keep=onp.ones(n5,bool); keep[rem]=False
N=onp.vstack([B[keep],add]); n=len(N); assert n==155
v=onp.concatenate([N[:,0],N[:,1],N[:,2]])
best=squeeze(v,n,13.25,log=lambda *a:None)
print(mode,"removed",[(round(B[i,0],1),round(B[i,1],1)) for i in rem],"added",len(add),"RESULT",round(best[1],5),flush=True)
pickle.dump((best[0],best[1],n),open(f"v2{mode}.pkl","wb"))
