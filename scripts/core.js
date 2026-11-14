// ---- squeeze core (shared by page and test) ----
function dep2(p,q){
  const dx=q.x-p.x,dy=q.y-p.y;if(dx*dx+dy*dy>2.01)return 0;
  const a=p.a*0.017453292519943295,b=q.a*0.017453292519943295,d=a-b,hh=.5+(Math.abs(Math.cos(d))+Math.abs(Math.sin(d)))/2;
  const ca=Math.cos(a),sa=Math.sin(a),cb=Math.cos(b),sb=Math.sin(b);
  const o=Math.min(hh-Math.abs(dx*ca+dy*sa),hh-Math.abs(-dx*sa+dy*ca),hh-Math.abs(dx*cb+dy*sb),hh-Math.abs(-dx*sb+dy*cb));
  return o>0?o:0;
}
function wallv(t,S){const r=t.a*0.017453292519943295,h=(Math.abs(Math.cos(r))+Math.abs(Math.sin(r)))/2;
  return [Math.max(0,h-t.x),Math.max(0,t.x+h-S),Math.max(0,h-t.y),Math.max(0,t.y+h-S)];}
function neigh(T){const N=T.map(()=>[]);for(let i=0;i<T.length;i++)for(let j=i+1;j<T.length;j++){const dx=T[i].x-T[j].x,dy=T[i].y-T[j].y;if(dx*dx+dy*dy<2.6){N[i].push(j);N[j].push(i);}}return N;}
function localE(T,N,i,S){const t=T[i];let e=0;for(const j of N[i]){const d=dep2(t,T[j]);e+=d*d;}const w=wallv(t,S);return e+w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3];}
function maxViol(T,N,S){let m=0;for(let i=0;i<T.length;i++){for(const j of N[i])if(j>i){const d=dep2(T[i],T[j]);if(d>m)m=d;}const w=wallv(T[i],S);for(const v of w)if(v>m)m=v;}return m;}
function relax(T,S,iters,tol){
  const h=1e-6,n=T.length;let N=neigh(T),v=1;
  const G=new Float64Array(3*n);
  for(let it=0;it<iters;it++){
    if(it%25===0){N=neigh(T);v=maxViol(T,N,S);if(v<tol)return v;}
    for(let i=0;i<n;i++){const t=T[i];if(localE(T,N,i,S)===0){G[3*i]=G[3*i+1]=G[3*i+2]=0;continue;}
      let o=t.x;t.x=o+h;let ep=localE(T,N,i,S);t.x=o-h;let em=localE(T,N,i,S);t.x=o;G[3*i]=(ep-em)/(2*h);
      o=t.y;t.y=o+h;ep=localE(T,N,i,S);t.y=o-h;em=localE(T,N,i,S);t.y=o;G[3*i+1]=(ep-em)/(2*h);
      o=t.a;const ha=h*57.29577951308232;t.a=o+ha;ep=localE(T,N,i,S);t.a=o-ha;em=localE(T,N,i,S);t.a=o;G[3*i+2]=(ep-em)/(2*h);}
    for(let i=0;i<n;i++){const t=T[i];t.x-=.22*G[3*i];t.y-=.22*G[3*i+1];t.a-=.35*G[3*i+2]*57.29577951308232;}
  }
  N=neigh(T);return maxViol(T,N,S);
}
