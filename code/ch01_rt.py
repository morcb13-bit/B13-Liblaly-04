# 第1章　黄金のひし形30面体を、5回軸6本の整数の符号の組で書く
# 判定は整数だけ。最後の凸包の照合だけが観察者側の実数。
import itertools
from collections import Counter

def zt(c):  # ζ（10乗で1）を掛ける：(c0,c1,c2,c3) -> (-c3, c0+c3, c1-c3, c2+c3)
    return (-c[3], c[0]+c[3], c[1]-c[3], c[2]+c[3])
def zpow(n):
    c=(1,0,0,0)
    for _ in range(n%10): c=zt(c)
    return c
def addv(a,b): return tuple(x+y for x,y in zip(a,b))
def scal(k,a): return tuple(k*x for x in a)

def solid(n, s):
    """n 本の軸、内積の符号 s[i][j]（i<j）。64通りの点を S で分ける。"""
    pts=list(itertools.product([-1,1],repeat=n))
    S={p:sum(p[i]*p[j]*s[i][j] for i in range(n) for j in range(i+1,n)) for p in pts}
    cls=Counter(S.values())
    top=sorted(cls,reverse=True)
    return pts,S,cls,top

def edges_faces(V):
    Vs=set(V); n=len(V[0])
    def flip(p,ids): return tuple(-x if i in ids else x for i,x in enumerate(p))
    E=sum(1 for p in V for i in range(n) if flip(p,{i}) in Vs)//2
    F=0
    for p in V:
        for i in range(n):
            for j in range(i+1,n):
                q1,q2,q3=flip(p,{i}),flip(p,{j}),flip(p,{i,j})
                if q1 in Vs and q2 in Vs and q3 in Vs: F+=1
    return E,F//4

def report(name,n,s,unit):
    print("=="+name+"==")
    pts,S,cls,top=solid(n,s)
    print("S の分布（大きい順）：",[(k,cls[k]) for k in top])
    a,b=cls[top[0]],cls[top[1]]
    V=[p for p in pts if S[p] in (top[0],top[1])]
    t1=sorted([a,b])==[12,20]
    print("T1 上の二組：",a,b,"→","OK" if t1 else "NG")
    E,F=edges_faces(V)
    t2=(E,F)==(60,30)
    print("T2 頂点",len(V),"辺",E,"面",F,"オイラー",len(V)-E+F,"→","OK" if t2 else "NG")
    return V,t1,t2

# 30面体：e0 は真上、e1..e5 は 72° ごと。内積 ±1/√5 の符号
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(i+1,6):
        if i==0: s[i][j]=1
        else: s[i][j]=1 if (j-i)%5 in (1,4) else -1
V,t1,t2=report("黄金のひし形30面体",6,s,5)
# T3：対角線 e_i±e_j の長さ二乗 ×√5 = 2√5 ± 2s → (√5の係数, 整数) の組
# 比が φ² = (3+√5)/2 か：2·(長い) = (3+√5)·(短い) を Z[√5] で照合
def mul(x,y): return (x[0]*y[1]+x[1]*y[0], x[0]*y[0]*5+x[1]*y[1])  # (a√5+b)(c√5+d), 組は(√5の係数,整数)
ok3=0;tot3=0
for i in range(6):
    for j in range(i+1,6):
        L=(2,2*s[i][j]); Sh=(2,-2*s[i][j])
        if s[i][j]<0: L,Sh=Sh,L
        tot3+=1; ok3+= (scal(2,L)==mul((1,3),Sh))
print("T3 面15組×2のうち、対角線の比が φ²：",ok3,"/",tot3,"→","OK" if ok3==tot3 else "NG")
# T4：真上からの影  Σ σ_k ζ^{2k}（k=1..5）、e0 の影は 0
sh=Counter()
for p in V:
    z=(0,0,0,0)
    for k in range(1,6): z=addv(z,scal(p[k],zpow(2*k)))
    sh[z]+=1
mult=Counter(sh.values())
centre=sh.get((0,0,0,0),0)
print("T4 影の番地の数：",len(sh),"重なりの数ごと：",dict(mult),"中心に乗る頂点：",centre)
# 10番地ずつの輪：ζ² 回し（＝72°）と符号反転（＝180°）で閉じる輪に分ける
def orbit(z):
    out=set(); w=z
    for _ in range(5):
        out.add(w); out.add(scal(-1,w)); w=zt(zt(w))
    return frozenset(out)
rings=Counter()
for z,m in sh.items():
    if z!=(0,0,0,0): rings[(orbit(z),m)]+=1
print("   輪（番地数, 重なり）：",sorted((len(o),m) for (o,m) in rings))
t4= len(sh)==21 and centre==2 and sorted((len(o),m) for (o,m) in rings)==[(10,1),(10,2)]
print("T4 →","OK" if t4 else "NG")

# 負の対照：立方体の対角線4本（内積 −1/3）→ 菱形十二面体
sc=[[0]*4 for _ in range(4)]
for i in range(4):
    for j in range(i+1,4): sc[i][j]=-1
Vc,c1,c2=report("負の対照：菱形十二面体",4,sc,3)
# T3 の対照：対角線の長さ二乗 ×3 = 2·3 ± 2 → 8 と 4、比 2。φ² は無理数なので整数比 2 とは一致しない
print("T3 対照：長い/短い = 8/4（整数の比）。φ²=(3+√5)/2 は √5 を含むので一致しない → NG")
print("負の対照が NG を返したか：T1",("NG" if not c1 else "OK(!)"),"T2",("NG" if not c2 else "OK(!)"),"T3 NG")

# 観察者側：凸包の頂点と一致するか（実数）
import numpy as np
from scipy.spatial import ConvexHull
t=np.arctan(2)
ax=[np.array([0,0,1.])]+[np.array([np.sin(t)*np.cos(2*np.pi*k/5),np.sin(t)*np.sin(2*np.pi*k/5),np.cos(t)]) for k in range(1,6)]
pts=list(itertools.product([-1,1],repeat=6))
P=np.array([sum(p[i]*ax[i]/2 for i in range(6)) for p in pts])
H=set(ConvexHull(P).vertices)
hv=set(pts[i] for i in H)
print("観察者側 凸包の頂点",len(hv),"個、T1 の32個と一致：",hv==set(V))
