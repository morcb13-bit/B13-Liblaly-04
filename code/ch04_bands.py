# 第4章　30面体を五回軸で段に割る。判定は整数だけ。
import json, itertools
from collections import Counter, defaultdict
from functools import cmp_to_key
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b):
    return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
R=json.load(open("ch01_rt.json")); V=[tuple(v) for v in R["V"]]; F=R["F"]
cen=[tuple(sum(V[i][k] for i in f) for k in range(6)) for f in F]   # 4頂点の和（単位は1/4）
# 値 x√5+y の大小を整数だけで
def cmp_pair(p,q):
    a=p[0]-q[0]; b=p[1]-q[1]   # a√5+b の符号
    if a==0: return (b>0)-(b<0)
    if b==0: return (a>0)-(a<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    # 符号が逆：5a² と b² を比べる
    big=5*a*a-b*b
    return (1 if a>0 else -1) if big>0 else ((1 if b>0 else -1) if big<0 else 0)
def bands(axis):
    h=[dot5(c,axis) for c in cen]          # √5·(c·axis)
    levels=sorted(set(h),key=cmp_to_key(cmp_pair),reverse=True)
    band=[levels.index(x) for x in h]
    return levels,band
E0=(1,0,0,0,0,0)
levels,band=bands(E0)
cnt=[band.count(k) for k in range(len(levels))]
print("B1 高さの組（上から）：",levels)
print("   面の数：",cnt,"→","OK" if cnt==[5,5,10,5,5] else "NG")
# 辺を共有する面
def edges(f): return {frozenset((f[k],f[(k+1)%4])) for k in range(4)}
adj=defaultdict(list)
for a in range(30):
    for b in range(a+1,30):
        sh=edges(F[a])&edges(F[b])
        if sh: adj[a].append((b,next(iter(sh)))); adj[b].append((a,next(iter(sh))))
pairs=Counter(tuple(sorted((band[a],band[b]))) for a in range(30) for b,_ in adj[a] if a<b)
skip=sum(n for (x,y),n in pairs.items() if y-x>1)
print("B2 段のつながり：",dict(sorted(pairs.items())),"／飛び越え",skip,"→","OK" if skip==0 else "NG")
# 段の中のつながり
def ring_of(k):
    fs=[f for f in range(30) if band[f]==k]
    e=[(a,b) for a in fs for b,_ in adj[a] if band[b]==k and a<b]
    deg=Counter(); [deg.update([a,b]) for a,b in e]
    if not e: return "なし",None
    if all(deg[f]==2 for f in fs) and len(e)==len(fs): return f"{len(fs)}枚の輪",fs
    return "その他",fs
inner=[ring_of(k)[0] for k in range(5)]
print("B3 段の中：",inner,"→","OK" if inner==["5枚の輪","なし","10枚の輪","なし","5枚の輪"] else "NG")
# 輪の通り方：輪の隣2枚と共有する2辺が、頂点を共有する（折れる）か、離れている（貫く）か
def how(k):
    fs=[f for f in range(30) if band[f]==k]; out=Counter()
    for f in fs:
        es=[e for b,e in adj[f] if band[b]==k]
        out["折れる" if es[0]&es[1] else "貫く"]+=1
    return dict(out)
hw=[how(k) for k in (0,2,4)]
print("B4 段0・段2・段4 の通り方：",hw,"→","OK" if hw==[{"折れる":5},{"貫く":10},{"折れる":5}] else "NG")
# 影
def zt(c): return (-c[3], c[0]+c[3], c[1]-c[3], c[2]+c[3])
def zpow(n):
    c=(1,0,0,0)
    for _ in range(n%10): c=zt(c)
    return c
def shadow(a):
    z=(0,0,0,0)
    for k in range(1,6): z=tuple(x+a[k]*y for x,y in zip(z,zpow(2*k)))
    return z
def closed(k,step):
    sh={shadow(cen[f]) for f in range(30) if band[f]==k}
    one=next(iter(sh)); orb=set(); z=one
    for _ in range(10):
        orb.add(z)
        for _ in range(step): z=zt(z)
    return sh==orb and len(sh)==len([f for f in range(30) if band[f]==k])
b5=closed(2,1) and closed(0,2) and closed(4,2)
print("B5 段2は ζ 倍で10個、段0・段4は ζ² 倍で5個 →","OK" if b5 else "NG")
s0={shadow(cen[f]) for f in range(30) if band[f]==0}; s4={shadow(cen[f]) for f in range(30) if band[f]==4}
print("   段0と段4の影の向き：段0に ζ を一回掛けると段4 ？",{zt(z) for z in s0}==s4)
# 負の対照：二回軸 e0+e1
l2,b2=bands((1,1,0,0,0,0))
c2=[b2.count(k) for k in range(len(l2))]
print("負の対照（二回軸 e0+e1）：段",len(l2),"通り、面の数",c2,"→ B1 は","NG" if c2!=[5,5,10,5,5] else "OK(!)")
l3,b3=bands((1,1,1,0,0,0)) if False else (None,None)
json.dump({"band":band},open("ch04_bands.json","w"))
