# 第9章　前の本の平面全体を、30面体を並べた立体の一枚の切り口として持ち上げる
import json, itertools, time
from collections import Counter, defaultdict
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b): return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
def cmp_pair(p,q):
    a=p[0]-q[0]; b=p[1]-q[1]
    if a==0: return (b>0)-(b<0)
    if b==0: return (a>0)-(a<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    big=5*a*a-b*b
    return (1 if a>0 else -1) if big>0 else ((1 if b>0 else -1) if big<0 else 0)
def zt(c): return (-c[3], c[0]+c[3], c[1]-c[3], c[2]+c[3])
def zpow(n):
    c=(1,0,0,0)
    for _ in range(n%10): c=zt(c)
    return c
add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b)); sc=lambda k,a:tuple(k*x for x in a)
def shadow(a):
    z=(0,0,0,0)
    for k in range(1,6): z=add(z,sc(a[k],zpow(2*k)))
    return z
index=lambda z:(z[0]-z[1]+z[2]-z[3])%5
# 第6章の取り戻し
def det(M):
    if len(M)==1: return M[0][0]
    return sum((-1)**j*M[0][j]*det([r[:j]+r[j+1:] for r in M[1:]]) for j in range(len(M)))
B=[zpow(2*k) for k in range(1,5)]; Mt=[[B[k][r] for k in range(4)] for r in range(4)]; dM=det(Mt)
def adjm(M):
    n=len(M); return [[(-1)**(i+j)*det([r[:i]+r[i+1:] for k,r in enumerate(M) if k!=j]) for j in range(n)] for i in range(n)]
A=adjm(Mt)
def lift(z,h):
    b=tuple(sum(A[i][j]*z[j] for j in range(4))//dM for i in range(4)); a0,S=h; r=S-sum(b)
    if r%5: return None
    t=r//5; return (a0,)+tuple(x+t for x in b)+(t,)
SIG=[tuple(v) for v in json.load(open("ch01_rt.json"))["V"]]
rings=[tuple(r) for r in json.load(open("expanded.json"))["rings"]]
cells={tuple(map(int,k.split(","))):v for k,v in json.load(open("carrier_1245.json"))["cells"].items()}
print("円環",len(rings),"個、五角形",len(cells),"枚")
# P1
ri=Counter(index(r) for r in rings)
p1=len(ri)==1
print("P1 円環の中心の指数：",dict(ri),"→","OK" if p1 else "NG")
# 五角形の頂点（2倍）
def pverts(p,k): return [sc(2,add(p,zpow(k+2*j))) for j in range(5)]
allv=set()
for p,k in cells.items(): allv.update(pverts(p,k))
res=Counter()
for v in allv:
    b=tuple(sum(A[i][j]*v[j] for j in range(4))//dM for i in range(4)); res[sum(b)%5]+=1
p2=len(res)==1
print("P2 五角形の頂点",len(allv),"個の「b の和 mod 5」：",dict(res),"→","OK（一枚の平面に乗る）" if p2 else "NG")
s0=next(iter(res)); H0=(0,s0 if s0<=2 else s0-5)
print("   平面の高さの組：",H0)
# P3：30面体を上か下に置いて切る
pidx=Counter(index(p) for p in cells)
print("   五角形の中心の指数：",dict(pidx))
def rt_center(p,updown):
    return lift(sc(2,p),(H0[0]+updown, H0[1]+3*updown))
# 上下：ch7 の円環（中心 0、指数 0 の担体）では指数 4 の五角形が上の段の下。ここでは指数の値で分け、二つの値の一方を上にする
vals=sorted(pidx)
def choose(up_val):
    UD={p:(1 if index(p)==up_val else -1) for p in cells}
    ok=0; RTc={}
    for p,k in cells.items():
        c=rt_center(p,UD[p]); RTc[p]=c
        if c is None: continue
        V=[add(c,sg) for sg in SIG]; on=[v for v in V if (v[0],sum(v[1:]))==H0]
        if len(on)==5 and {shadow(v) for v in on}==set(pverts(p,k)): ok+=1
    return ok,RTc,UD
best=None
for uv in vals:
    ok,RTc,UD=choose(uv)
    print("   指数",uv,"の五角形を上に置くと、切り口が五角形に一致：",ok,"/",len(cells))
    if ok==len(cells): best=(uv,RTc,UD)
p3=best is not None
print("P3 →","OK" if p3 else "NG")
json.dump({"H0":H0,"up_index":best[0]},open("ch09_plane.json","w"))
# 負の対照（P2）：円環一つを番地 1 だけずらして平面に混ぜると、頂点が一枚の平面に乗らない
phi2=(2,0,1,-1)
def zmul(x,y):
    r=(0,0,0,0)
    for i in range(4): r=add(r,sc(x[i],y)); y=zt(y)
    return r
r0=add(rings[0],(1,0,0,0))
extra=set()
for k in range(10):
    p=add(r0,zmul(phi2,zpow(k))); extra.update(pverts(p,k))
res2=Counter()
for v in allv|extra:
    b=tuple(sum(A[i][j]*v[j] for j in range(4))//dM for i in range(4)); res2[sum(b)%5]+=1
print("負の対照：ずらした円環を混ぜたときの余り：",dict(res2),"→ P2 は","NG" if len(res2)>1 else "OK(!)")
