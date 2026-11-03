# 第2章　30面体1個（最初の1）と、10個の輪（次の1）
# 長さは半分を単位にした6次元の整数の組。判定は整数だけ。
import itertools
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b):  # √5·(a·b) を (√5の係数, 整数) で
    return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
def height(a): return (a[0], sum(a[1:]))   # 高さ = a0 + (a1+…+a5)/√5 の整数の組
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
# 30面体の頂点（中心からの符号の組）：第1章の S=5, 3
SIG=[p for p in itertools.product([-1,1],repeat=6)
     if sum(p[i]*p[j]*s[i][j] for i in range(6) for j in range(i+1,6)) in (5,3)]
def verts(c): return {add(c,p) for p in SIG}
# ζ の掛け算（Z[ζ10] の番地）
def zt(c): return (-c[3], c[0]+c[3], c[1]-c[3], c[2]+c[3])
def zpow(n):
    c=(1,0,0,0)
    for _ in range(n%10): c=zt(c)
    return c
def shadow(a):
    z=(0,0,0,0)
    for k in range(1,6): z=tuple(x+a[k]*y for x,y in zip(z,zpow(2*k)))
    return z
def M(a):  # e1..e5 を3つずらし、全体を反転
    b=[0]*6; b[0]=a[0]
    for k in range(1,6): b[1+(k-1+3)%5]=a[k]
    return tuple(-x for x in b)

c0=(1,1,-1,-1,1,3)
C=[c0]
for _ in range(10): C.append(M(C[-1]))
r0= C[10]==c0 and len(set(C[:10]))==10
print("R0 M を10回で戻る：",C[10]==c0,"／途中の異なる数：",len(set(C[:10])),"→","OK" if r0 else "NG")
for k in range(10): print("   c%d ="%k,C[k])
D=[sub(C[k+1],C[k]) for k in range(10)]
tot=(0,)*6
for d in D: tot=add(tot,d)
L=set(dot5(d,d) for d in D)
r1= tot==(0,)*6 and len(L)==1
print("R1 一歩の和：",tot,"／長さの組：",L,"→","OK" if r1 else "NG")
VS=[verts(c) for c in C[:10]]
sh=[len(VS[k]&VS[(k+1)%10]) for k in range(10)]
r2= sh==[4]*10
print("R2 隣どうしの共有頂点：",sh,"→","OK" if r2 else "NG")
far=[len(VS[a]&VS[b]) for a in range(10) for b in range(a+2,10) if not (a==0 and b==9)]
r3= len(far)==35 and all(x==0 for x in far)
print("R3 隣でない",len(far),"組の共有頂点の最大：",max(far),"→","OK" if r3 else "NG")
H=[height(c) for c in C[:10]]
r4= all(H[k]==H[0] for k in range(0,10,2)) and all(H[k]==(-H[0][0],-H[0][1]) for k in range(1,10,2)) and H[0][0]>0
print("R4 高さの組（a0, Σ）：",H[:2],"…","→","OK" if r4 else "NG")
SH=[shadow(c) for c in C[:10]]
r5= all(SH[(k+1)%10]==zt(SH[k]) for k in range(10))
print("R5 影が一つ進むごとに ζ 倍：",sum(SH[(k+1)%10]==zt(SH[k]) for k in range(10)),"/ 10 →","OK" if r5 else "NG")
print("   影の番地 c0：",SH[0])
# 負の対照
d=D[0]; X=c0
for _ in range(10): X=add(X,d)
n1= sub(X,c0)!=(0,)*6
print("N1 まっすぐ10歩の和：",sub(X,c0),"→ R1 は","NG" if n1 else "OK(!)")
w=(2,2,2,2,2,2)
n2=len(verts(c0)&verts(add(c0,w)))
print("N2 5回軸の向きに一歩：共有頂点",n2,"→ R2 は","NG" if n2!=4 else "OK(!)")

# R3 が NG だったので原因を調べる（事前登録の外。後から立てた確認）
print("--- R3 の原因 ---")
pairs={}
for a in range(10):
    for b in range(a+1,10):
        pairs[(a,b)]=len(VS[a]&VS[b])
gap=lambda a,b:min((b-a)%10,(a-b)%10)
from collections import Counter
byg=Counter()
for (a,b),n in pairs.items(): byg[(gap(a,b),n)]+=1
print("番号の離れ（1〜5）ごとの共有頂点：",sorted(byg.items()))
print("R3′ 番号が2離れた組（同じ段の隣）：10組すべてで共有頂点4 →",
      "OK" if byg[(2,4)]==10 else "NG")
print("R3″ 番号が3以上離れた25組：共有頂点0 →",
      "OK" if byg[(3,0)]+byg[(4,0)]+byg[(5,0)]==25 else "NG")
deg=[sum(1 for b in range(10) if b!=a and len(VS[a]&VS[b])==4) for a in range(10)]
print("各30面体が面で接する相手の数：",deg)
up=[C[k] for k in range(0,10,2)]
print("上の段5個の一歩の長さの組：",set(dot5(sub(up[(i+1)%5],up[i]),sub(up[(i+1)%5],up[i])) for i in range(5)),
      "（ジグザグの一歩と同じ (16, 32) なら同じ面の接し方）")
