import math
from fractions import Fraction
D="""5 2.70711;10 3.70711;11 3.87708;17 4.67553;18 4.82288;19 4.88562;26 5.62132;27 5.70711;28 5.82444;29 5.93383;37 6.59862;38 6.70711;39 6.81072;40 6.82843;41 6.92669;50 7.57143;51 7.70080;52 7.70711;53 7.82288;54 7.84667;55 7.94577;65 8.53553;66 8.65685;67 8.70711;68 8.79880;69 8.82719;70 8.88167;71 8.94407;82 9.53553;83 9.63476;84 9.70711;85 9.74264;86 9.82288;87 9.83882;88 9.88815;89 9.94975;
101 10.53553;102 10.61139;103 10.70378;104 10.70711;105 10.80762;106 10.82298;107 10.84667;108 10.92592;109 10.94975;110 10.99679;
122 11.53553;123 11.60140;124 11.65685;125 11.70711;126 11.77474;127 11.82288;128 11.82509;129 11.88131;130 11.91119;131 11.95653;132 11.99137;
145 12.53553;146 12.60091;147 12.65685;148 12.65685;149 12.70711;150 12.77817;151 12.82288;152 12.83096;153 12.88167;154 12.93167;155 12.95845;156 12.98208;
170 13.53553;171 13.57143;172 13.61899;173 13.65685;174 13.70711;175 13.77817;176 13.82288;177 13.82303;178 13.84667;179 13.89534;180 13.93509;181 13.95691;182 13.97419;
197 14.53553;198 14.57143;199 14.61899;200 14.65685;201 14.70711;202 14.72792;203 14.77817;204 14.82288;205 14.82445;206 14.87222;207 14.89395;208 14.93777;209 14.95862;210 14.97413;
226 15.53553;227 15.57107;228 15.60902;229 15.65685;230 15.68293;231 15.70711;232 15.77817;233 15.77817;234 15.82288;235 15.82661;236 15.87608;237 15.91421;238 15.93966;239 15.95635;240 15.97556;241 15.99081;
257 16.53553;258 16.57107;259 16.60257;260 16.65685;261 16.68293;262 16.70711;263 16.74264;264 16.77817;265 16.77817;266 16.82306;267 16.84667;268 16.87931;269 16.90597;270 16.94062;271 16.95500;272 16.96972;273 16.98821;
290 17.53553;291 17.53553;292 17.60257;293 17.63415;294 17.65685;295 17.70711;296 17.70711;297 17.74106;298 17.77817;299 17.82288;300 17.82412;301 17.86889;302 17.88675;303 17.93126;304 17.94911;305 17.96066;306 17.96914;307 17.98272"""
rec={int(a):float(b) for a,b in (x.split() for x in D.replace('\n','').split(';'))}
B={}
m=1e9
for n in range(324,0,-1):
    m=min(m,rec.get(n,math.ceil(math.sqrt(n)))); B[n]=m
eps=2e-5
def form(s):
    f=s-math.floor(s)
    # rational small denom
    for d in range(2,60):
        p=round(f*d)
        if abs(f-p/d)<1.5e-5 and math.gcd(p,d)==1: return f"rational /{d}"
    r2=math.sqrt(2)
    for d in (1,2,3,4):
        for j in range(1,40):
            v=j*r2/d; 
            if abs((v-f)-round(v-f))<1.5e-5: return f"45deg: {j}/{d}*sqrt2"
    for c,nm in ((math.sqrt(7)/2,'sqrt7/2'),):
        if abs((c+.5-f)-round(c+.5-f))<1.5e-5: return nm
    return "odd-angle (poly root)"
kinds={}
for n,s in sorted(rec.items()):
    k=math.ceil(s); r=k*k-n
    mm=n-2*k+3
    if n+1 in B and abs(B[n+1]-s)<eps and (n+1) in rec: kind='drop'  # same as n+1
    elif mm>=1 and B[mm]<=s-1+eps: kind='L-child'
    else: kind='SEED'
    kinds[n]=(kind,k,r,form(s))
import collections
for k in range(3,19):
    row=[(n,)+kinds[n] for n in kinds if kinds[n][1]==k]
    c=collections.Counter(x[1] for x in row)
    print(f"k={k:2d} tilted={len(row):2d} rmin={min(x[3] for x in row):2d} rmax={max(x[3] for x in row):2d} seeds={c['SEED']} L={c['L-child']} drop={c['drop']}")
print()
for n in sorted(kinds):
    kd,k,r,fm=kinds[n]
    if kd=='SEED': print(n,rec[n],'k',k,'r',r,fm)
print(collections.Counter(v[0] for v in kinds.values()))
print(collections.Counter(v[3].split(':')[0] for v in kinds.values() if v[0]=='SEED'))
# frac families recurring
fr=collections.defaultdict(list)
for n,s in rec.items(): fr[round(s-math.floor(s),4)].append(n)
for f,l in sorted(fr.items()):
    if len(l)>=3: print(f,l)
