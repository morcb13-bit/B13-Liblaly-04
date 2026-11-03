# 第7章　真ん中で切る（K1・K2）と、様々な角度から切る（記録）
import json, itertools
from collections import Counter
from functools import cmp_to_key
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
add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sc=lambda k,a:tuple(k*x for x in a)
def shadow(a):
    z=(0,0,0,0)
    for k in range(1,6): z=add(z,sc(a[k],zpow(2*k)))
    return z
def zmul(x,y):
    r=(0,0,0,0)
    for i in range(4):
        r=add(r,sc(x[i],y)); y=zt(y)
    return r
R=json.load(open("ch01_rt.json")); SIG=[tuple(v) for v in R["V"]]; E=R["E"]
C=[tuple(c) for c in json.load(open("ch02_ring.json"))["C"]]
# ---- 前の本の円環を、前の本の規則だけで作り直す ----
phi2=(2,0,1,-1)
def book3_ring(off):
    ring=[]
    for k in range(10):
        cen=zmul(phi2,zpow(k))
        ring.append([add(cen,zpow(k+off+2*j)) for j in range(5)])   # 外接半径 1 の五角形
    return ring
def shared_edges(ring):
    cnt=0
    for k in range(10):
        A=set(ring[k]); B=set(ring[(k+1)%10]); cnt+= len(A&B)==2
    return cnt
best=None
for off in range(2):
    se=shared_edges(book3_ring(off)); print("前の本の円環の候補 off =",off,"：隣と辺を共有",se,"/ 10")
    if se==10: best=off
ring3=book3_ring(best); other=book3_ring(1-best)
V3={v for P in ring3 for v in P}; V3o={v for P in other for v in P}
# ---- K1 ----
cut=set()
for c in C:
    for sg in SIG:
        a=add(c,sg)
        if a[0]==0 and sum(a[1:])==0: cut.add(a)
cs={shadow(a) for a in cut}
target={sc(2,v) for v in V3}
k1= len(cut)==30 and cs==target
print("K1 真ん中の平面に乗る頂点：",len(cut),"個／影＝前の本の円環の頂点30個の2倍と一致：",cs==target,"→","OK" if k1 else "NG")
print("負の対照：36°回した輪と照合：",cs=={sc(2,v) for v in V3o},"→ K1 は","NG" if cs!={sc(2,v) for v in V3o} else "OK(!)")
# ---- K2 ----
def index(z): return (z[0]-z[1]+z[2]-z[3])%5
idx=[index(zmul(phi2,zpow(k))) for k in range(10)]
up=[idx[k] for k in range(0,10,2)]; dn=[idx[k] for k in range(1,10,2)]
k2= len(set(up))==1 and len(set(dn))==1 and up[0]!=dn[0]
print("K2 円環の五角形の中心の指数：上の段の下",set(up),"下の段の上",set(dn),"→","OK" if k2 else "NG")
# ---- 後半：30面体1個を三つの軸で切る（記録） ----
def slices(n):
    h=[dot5(v,n) for v in SIG]
    lv=sorted(set(h),key=cmp_to_key(cmp_pair))
    out=[]
    for i,L in enumerate(lv):
        on=sum(1 for x in h if x==L)
        cross=sum(1 for a,b,_ in E if cmp_pair(h[a],L)*cmp_pair(h[b],L)<0)
        out.append(("面上" if i in (0,len(lv)-1) else "高さ",L,on,on+cross))
        if i<len(lv)-1:
            cross2=sum(1 for a,b,_ in E if (cmp_pair(h[a],L)<=0 and cmp_pair(h[b],lv[i+1])>=0) or (cmp_pair(h[b],L)<=0 and cmp_pair(h[a],lv[i+1])>=0))
            out.append(("あいだ",None,0,cross2))
    return lv,out
v3=next(v for v,k in zip(SIG,R["K"]) if k==3)            # 三回軸の向き（S=3 の頂点）
f=R["F"][0]; n2=tuple(sum(SIG[i][k] for i in f) for k in range(6))   # 二回軸の向き（面の中心）
AXES=[("五回軸 e0",(1,0,0,0,0,0)),("三回軸",v3),("二回軸",n2)]
rec={}
for name,n in AXES:
    lv,out=slices(n)
    seq=[o[3] for o in out if o[0]=="あいだ"]
    print(f"[記録] {name}：頂点の高さ {len(lv)} 段。段のあいだで切った切り口の角の数（下から）：",seq)
    rec[name]=seq
# ---- 輪全体を水平に切る（記録） ----
ringV=[(k,add(c,sg)) for k,c in enumerate(C) for sg in SIG]
H=[dot5(a,(1,0,0,0,0,0)) for _,a in ringV]
lv=sorted(set(H),key=cmp_to_key(cmp_pair))
print("[記録] 輪全体：頂点の高さ",len(lv),"段")
ring_rec=[]
for i,L in enumerate(lv):
    pieces=Counter()
    for k,c in enumerate(C):
        h=[dot5(add(c,sg),(1,0,0,0,0,0)) for sg in SIG]
        on=sum(1 for x in h if x==L)
        cross=sum(1 for a,b,_ in E if cmp_pair(h[a],L)*cmp_pair(h[b],L)<0)
        if on+cross>=3: pieces[on+cross]+=1
    ring_rec.append((L,dict(pieces)))
    print("   高さの組",L,"：切り口",dict(pieces) if pieces else "なし（点か辺だけ）")
json.dump({"ring3":ring3,"rec":rec},open("ch07_cut.json","w"))
