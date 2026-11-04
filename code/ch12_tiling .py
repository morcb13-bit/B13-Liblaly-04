# 第12章の担体：黄金のひし形の六面体による立体のペンローズを、6次元の整数の点から作る
# 窓の判定は (√5の係数, 整数) の組で、整数だけで行う
import itertools, json, math
from fractions import Fraction
import numpy as np
r5=math.sqrt(5)
# 実の軸と共役の軸（共役は表示と係数の読み取りにだけ使う）
t=math.atan(2)
E=[np.array([0,0,1.])]+[np.array([math.sin(t)*math.cos(2*math.pi*k/5),math.sin(t)*math.sin(2*math.pi*k/5),math.cos(t)]) for k in range(1,6)]
Es=[np.array([0,0,1.])]+[np.array([2/r5*math.cos(4*math.pi*k/5),2/r5*math.sin(4*math.pi*k/5),-1/r5]) for k in range(1,6)]
# 三重積 D_ijk = det(e*_i,e*_j,e*_k) を c·(a+b√5)/2 の整数の組で読む
c=abs(np.linalg.det(np.array([Es[0],Es[1],Es[2]])))
def recog(x):
    for b in range(-12,13):
        a=2*x/c-b*r5
        if abs(a-round(a))<1e-7: return (round(a),b)
    raise ValueError(x)
PAIRS=list(itertools.combinations(range(6),2))
D={}
for i,j in PAIRS:
    D[(i,j)]=[recog(np.linalg.det(np.array([Es[i],Es[j],Es[k]]))) for k in range(6)]
def sgn(a,b):   # a+b√5 の符号
    if b==0: return (a>0)-(a<0)
    if a==0: return (b>0)-(b<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    big=a*a-5*b*b
    return (1 if a>0 else -1) if big>0 else -(1 if a>0 else -1)
def absq(p): return p if sgn(*p)>=0 else (-p[0],-p[1])
HALF={}
for ij,row in D.items():
    s=(0,0)
    for p in row: q=absq(p); s=(s[0]+q[0],s[1]+q[1])
    HALF[ij]=s           # Σ|D|（2倍のスケールでの h の2倍）
GAM=[3,7,11,17,23,29]; DEN=97   # 一般の位置へずらす量 γ = GAM/DEN
def accept(x):
    for ij,row in D.items():
        A=sum((DEN*x[k]+GAM[k])*row[k][0] for k in range(6)); B=sum((DEN*x[k]+GAM[k])*row[k][1] for k in range(6))
        H=(DEN*HALF[ij][0],DEN*HALF[ij][1])
        if sgn(H[0]-A,H[1]-B)<0 or sgn(H[0]+A,H[1]+B)<0: return False   # |A+B√5| ≤ H
    return True
SIGN=1
STRICT=False   # True：係数がちょうど 0 になる境目の六面体を入れない（鏡で対称な γ のとき）
def lin(x,row):
    A=sum((DEN*x[k]+SIGN*GAM[k])*row[k][0] for k in range(6)); B=sum((DEN*x[k]+SIGN*GAM[k])*row[k][1] for k in range(6))
    return (SIGN*A,SIGN*B)
def neg(p): return (-p[0],-p[1])
def tile_ok(x,tri):
    # 六面体 (x; i,j,k) が入る ⟺ π*(x+γ) が、残りの3本 l,m,n の張る平行六面体の中（係数が 0 以上 1 未満）
    l,m,n=[q for q in range(6) if q not in tri]
    Dl=D[(l,m)][n]; sD=sgn(*Dl)
    rows=[D[(m,n)],[neg(p) for p in D[(l,n)]],D[(l,m)]]
    for row in rows:
        N=lin(x,row)
        if sgn(*N)*sD<0 or (STRICT and sgn(*N)==0): return False
        rest=(DEN*Dl[0]-N[0],DEN*Dl[1]-N[1])
        if sgn(*rest)*sD<=0: return False
    return True
def phys(x): return sum(x[k]*E[k] for k in range(6))
def build(R):
    start=(0,)*6
    # 原点が窓に入らないこともあるので、近くから入る点を探す
    if not accept(start):
        for d in itertools.product([-1,0,1],repeat=6):
            if accept(d): start=d; break
    seen={start}; stack=[start]
    while stack:
        x=stack.pop()
        for k in range(6):
            for s in (1,-1):
                y=tuple(x[m]+(s if m==k else 0) for m in range(6))
                if y in seen: continue
                if np.linalg.norm(phys(y))>R: continue
                if accept(y): seen.add(y); stack.append(y)
    P=seen; cells=[]
    for x in P:
        for tri in itertools.combinations(range(6),3):
            if tile_ok(x,tri): cells.append((x,tri))
    return P,cells
def faces_of(cell):
    x,(i,j,k)=cell
    def sh(v,m): return tuple(v[n]+(1 if n==m else 0) for n in range(6))
    return [(x,(i,j)),(sh(x,k),(i,j)),(x,(i,k)),(sh(x,j),(i,k)),(x,(j,k)),(sh(x,i),(j,k))]
if __name__=="__main__":
    import sys
    R=float(sys.argv[1]) if len(sys.argv)>1 else 6.0
    P,cells=build(R)
    from collections import Counter, defaultdict
    F=defaultdict(list)
    for n,cl in enumerate(cells):
        for f in faces_of(cl): F[f].append(n)
    mult=Counter(len(v) for v in F.values())
    print("半径",R,"：格子点",len(P),"、六面体",len(cells),"、面の持ち主の数ごと",dict(mult))
    # 太と細：三重積の大きさ（物理側）で2種類
    vol=Counter(round(abs(np.linalg.det(np.array([E[i],E[j],E[k]]))),6) for _,(i,j,k) in cells)
    print("体積の種類",dict(vol))
    json.dump({"cells":[[list(x),list(t)] for x,t in cells]},open("ch12_tiling.json","w"))
