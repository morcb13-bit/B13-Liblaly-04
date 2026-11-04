# 第17章：二つの中心 A, B のまわりの同心球（半径 m·d、d = AB/N）で、楕円体と双曲面の帯を作る
# 点は立体のペンローズの頂点（6 個の整数）。半径の比べ方は (√5 の係数, 整数) の組で、整数だけ。
import numpy as np, json
from collections import Counter
import ch12_tiling as T
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b): return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
def sgn(a,b):   # a√5 + b の符号
    if a==0: return (b>0)-(b<0)
    if b==0: return (a>0)-(a<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    big=5*a*a-b*b
    return (1 if a>0 else -1) if big>0 else ((1 if b>0 else -1) if big<0 else 0)
N=4     # AB の長さ = N·d。同心球の半径は d, 2d, 3d, …
def shell(x,c,ab):   # x と c の距離 r の殻の番号 m（m·d ≤ r < (m+1)·d）。ab = √5|AB|² の組
    d=tuple(p-q for p,q in zip(x,c)); P,Q=dot5(d,d); m=0
    # (m+1)² |AB|² ≤ N² r²  ⇔  N²(P√5+Q) − (m+1)²(ab) ≥ 0
    while sgn(N*N*P-(m+1)**2*ab[0], N*N*Q-(m+1)**2*ab[1])>=0: m+=1
    return m
PERM=[0,4,3,2,1,5]
mirror=lambda x:tuple(x[PERM[k]] for k in range(6))
def points(gam,R):
    T.GAM=gam; P,_=T.build(R); return P
if __name__=="__main__":
    Pts=points([3,7,11,11,7,29],10.0)
    sym=all(mirror(x) in Pts for x in Pts)
    print("頂点",len(Pts),"個、鏡で自分に重なる：",sym)
    ph=lambda x:T.phys(x)
    A=min((x for x in Pts if ph(x)[1]>0.5 and mirror(x)!=x),key=lambda x:np.linalg.norm(ph(x)-np.array([0,2.6,0])))
    B=mirror(A); AB=dot5(tuple(p-q for p,q in zip(B,A)),tuple(p-q for p,q in zip(B,A)))
    sh=lambda x,c:shell(x,c,AB)
    R0=max(np.linalg.norm(ph(x)) for x in Pts); inner=[x for x in Pts if np.linalg.norm(ph(x))<R0-1.0]
    KA={x:sh(x,A) for x in inner}; KB={x:sh(x,B) for x in inner}
    print("√5|AB|² の組：",AB,"、B は A から見て番号",sh(B,A),"の殻（N =",N,"）")
    c1=all(KA[x]+KB[x]>=N-1 and abs(KA[x]-KB[x])<=N for x in inner)
    print("C1 すべての点で m + n ≥ N − 1、|m − n| ≤ N（三角形の関係）：",c1,"→","OK" if c1 else "NG")
    fixed=[x for x in inner if mirror(x)==x]
    c2=all(KA[x]==KB[x] for x in fixed)
    print("C2 鏡の面の上の点",len(fixed),"個がすべて m = n：",c2,"→","OK" if c2 else "NG")
    c3=all(KA[mirror(x)]==KB[x] for x in inner if mirror(x) in KA)
    print("C3 鏡に映すと m と n が入れ替わる：",c3,"→","OK" if c3 else "NG")
    ell=Counter(KA[x]+KB[x] for x in inner); hyp=Counter(KA[x]-KB[x] for x in inner)
    print("[記録] 点の数",len(inner))
    print("[記録] 楕円体の帯 j = m+n ごとの点の数：",dict(sorted(ell.items())))
    print("[記録] 双曲面の帯 k = m−n ごとの点の数：",dict(sorted(hyp.items())))
    B2=min((x for x in Pts if x!=B),key=lambda x:np.linalg.norm(ph(x)-ph(B)))
    c2n=all(sh(x,A)==sh(x,B2) for x in fixed)
    print("負の対照：B を隣の頂点にずらす → 鏡の面の点がすべて m = n：",c2n,"→","NG" if not c2n else "OK(!)")
    json.dump({"N":N,"AB":AB,"ell":dict(sorted(ell.items())),"hyp":dict(sorted(hyp.items())),"mirror_points":len(fixed)},
              open("ch17_confocal.json","w"),ensure_ascii=False)
