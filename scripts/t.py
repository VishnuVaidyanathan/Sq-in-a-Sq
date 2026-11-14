import math,csv
r2=math.sqrt(2)
up={5:2+1/r2,10:3+1/r2,11:3.8771,17:4.6756,18:3.5+math.sqrt(7)/2,19:3+4*r2/3,26:3.5+1.5*r2,27:5+1/r2,28:3+2*r2,29:5.9344,
37:6.5987,38:6+1/r2,39:6.8189,40:4+2*r2,41:6.9473,50:53/7,51:7.7044,52:7+1/r2,53:7.8231,54:7.8488,55:7.9871,
65:5+5/r2,66:3+4*r2,67:8+1/r2,68:7.5+math.sqrt(7)/2,69:8.8562,70:8.9121,71:8.9633,
82:6+5/r2,83:4+4*r2,84:9+1/r2,85:5.5+3*r2,86:8.5+math.sqrt(7)/2,87:9.8520,88:9.9018,89:5+7/r2}
lo={2:2,5:2+1/r2,6:3,10:3+1/r2,11:2+4/math.sqrt(5),13:4,14:4,17:(40*r2+19)/17,19:6*r2-4,21:4.7438,22:5,23:5,
26:2*r2+(27+2*math.sqrt(10))/13,28:2*r2+6/math.sqrt(5),31:5.6415,33:6,34:6,37:2*r2+(113+10*math.sqrt(3))/37,40:2*r2+8/math.sqrt(5),
46:7,47:7,50:2*r2+(101+3*math.sqrt(14))/25,62:8,65:2*r2+71/13,79:9,82:9.2667,98:10}
rows=[];L=0
U={}
m=1e9
for n in range(100,0,-1):
    u=min(up.get(n,math.ceil(math.sqrt(n))),m); m=u; U[n]=u
for n in range(1,101):
    k=math.ceil(math.sqrt(n)); r=k*k-n
    f=math.sqrt(n)
    if n>=8 and r>0: f=max(f,0.5+math.sqrt(n-math.isqrt(n)+0.25))
    L=max(L,f,lo.get(n,0))
    rows.append((n,k,r,round(L,4),round(U[n],4),round(U[n]-L,4),round(k-U[n],4),round(U[n]**2-n,3)))
csv.writer(open('t.csv','w')).writerows([("n","k","r","floor","best","gap","saving_vs_grid","waste_area")]+rows)
for x in rows: print(*x)
