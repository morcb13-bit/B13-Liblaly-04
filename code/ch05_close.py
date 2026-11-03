# 第5章　30面体は平行移動で隣へ、正十二面体と正四面体は反射で隣へ。閉じるかを整数で見る。
import itertools, json
from math import gcd
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b):  # √5·(a·b) = P√5 + Q を (P, Q) で
    return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
def E(i): return tuple(1 if k==i else 0 for k in range(6))
add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
sc=lambda k,a:tuple(k*x for x in a)
# Q0：r_i = Σ_j s_ij e_j − √5 e_i が、どの軸とも内積 0（(√5の係数,整数) の組で）
# √5·(r_i·e_m) = √5·(Σ s_ij e_j · e_m) − 5 (e_i·e_m)
ok0=True
for i in range(6):
    w=tuple(s[i][j] if j!=i else 0 for j in range(6))
    for m in range(6):
        P,Q=dot5(w,E(m)); Pi,Qi=dot5(E(i),E(m))   # √5(e_i·e_m)=Pi√5+Qi → 5(e_i·e_m)=√5(Pi√5+Qi)=5Pi+Qi√5
        # √5(r·e_m) = (P√5+Q) − (5Pi + Qi√5) → (P−Qi)√5 + (Q−5Pi)
        if (P-Qi, Q-5*Pi)!=(0,0): ok0=False
print("Q0 √5·e_i = Σ s_ij e_j（6本すべて）→","OK" if ok0 else "NG")
# Q1：30面体の面の一歩（平行移動）
SIG=[p for p in itertools.product([-1,1],repeat=6) if sum(p[i]*p[j]*s[i][j] for i in range(6) for j in range(i+1,6)) in (5,3)]
F=json.load(open("ch01_rt.json"))["F"]; V=[tuple(v) for v in json.load(open("ch01_rt.json"))["V"]]
def verts(c,S=SIG): return {add(c,p) for p in S}
ok1=0
for f in F:
    cen=tuple(sum(V[i][k] for i in f) for k in range(6))   # 面の中心×4（反転した軸は0、残りは±4）
    step=tuple(x//2 for x in cen)                           # 中心から面までの2倍＝隣の中心（半分単位で ±2 が4つ、0 が2つ）
    ok1+= len(verts((0,)*6)&verts(step))==4
print("Q1 面の一歩（±2 が4つ・0 が2つの整数の組）で共有頂点4：",ok1,"/ 30 →","OK" if ok1==30 else "NG")
# Q2：正十二面体の反射（面に垂直な軸 e_i、原点を通る面で：回転の部分だけを見る）
# x' = x − 2(x·e_i)e_i、 x·e_i = P + Q/√5、 e_i/√5 = (Σ s_ij e_j)/5  → 5x' = 5x − 10P e_i − 2Q Σ s_ij e_j
def refl_dod(x,den,i):
    P,Q=dot5(x,E(i))      # √5(x·e_i) = P√5+Q → x·e_i = (P + Q/√5)/den
    w=tuple(s[i][j] if j!=i else 0 for j in range(6))
    y=sub(sub(sc(5,x),sc(10*P,E(i))),sc(2*Q,w))
    return y,den*5
def reduce(x,den):
    g=den
    for v in x: g=gcd(g,abs(v))
    return tuple(v//g for v in x),den//g
x0=E(2); x,den=x0,1; back=None; dens=[]
for n in range(1,61):
    x,den=refl_dod(x,den,0); x,den=refl_dod(x,den,1); x,den=reduce(x,den); dens.append(den)
    if den==1 and x==x0 and back is None: back=n
okden=all(dens[k]==5**(2*(k+1)) for k in range(60))
print("Q2 正十二面体：二つの面で交互に反射、60回までに戻った回：",back,"／分母が毎回 25 倍のまま：",okden,"→","OK" if back is None and okden else "NG")
print("   参考：最初の分母",dens[:4])
# Q3：正四面体の一本道（ボールディク＝コクセターのらせん）
T=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
# 頂点を 3^n 倍の整数で持つ。新しい頂点 v4 = (2/3)(v1+v2+v3) − v0
pts=[(sc(1,p),1) for p in T]   # (整数ベクトル, 分母)
def lift(p,den): return sc(den//p[1],p[0])
for n in range(60):
    a,b,c,d=pts[-4:]
    D=3*max(x[1] for x in (a,b,c,d))
    sv=add(add(lift(b,D),lift(c,D)),lift(d,D))
    nv=sub(sc(2,tuple(v//3 for v in sv)) if all(v%3==0 for v in sv) else None, lift(a,D)) if False else None
    # 分母 D で：v4·D = (2/3)(b+c+d)·D − a·D  → 3 で割り切れるように D を 3 倍してある
    nv=sub(tuple(2*v//3 for v in sv), lift(a,D))
    g=D
    for v in nv: g=gcd(g,abs(v))
    pts.append((tuple(v//g for v in nv),D//g))
def d2(p,q):
    L=max(p[1],q[1]); u=sub(lift(p,L),lift(q,L)); return sum(v*v for v in u),L*L
regular=True
for k in range(len(pts)-3):
    q=pts[k:k+4]; L=1
    for p in q: L=max(L,p[1])
    vals=set()
    for i in range(4):
        for j in range(i+1,4):
            a2,l2=d2(q[i],q[j]); vals.add(a2*(L*L)//l2)   # 共通分母 L² にそろえた長さの二乗
    regular&= len(vals)==1
keyset={(p[0],p[1]) for p in pts}
denoms=[p[1] for p in pts[4:]]
grow=all(denoms[k]==3**(k+1) for k in range(len(denoms)))
print("Q3 正四面体の一本道：",len(pts),"頂点、どの4頂点も正四面体：",regular,"／頂点の重なり：",len(pts)-len(keyset),"／分母が毎回 3 倍：",grow,
      "→","OK" if regular and len(keyset)==len(pts) and grow else "NG")
# 負の対照1：立方体を x と y の面で交互に反射
def refl_cube(v,ax): return tuple(-x if k==ax else x for k,x in enumerate(v))
v0=(1,2,3); v=v0; back=None
for n in range(1,61):
    v=refl_cube(refl_cube(v,0),1)
    if v==v0: back=n; break
print("負の対照：立方体（x と y の面で交互に反射）は",back,"回で戻る → 検定は「閉じる」を返す")
# 負の対照2：30面体の輪の操作 M
def M(a):
    b=[0]*6; b[0]=a[0]
    for k in range(1,6): b[1+(k-1+3)%5]=a[k]
    return tuple(-x for x in b)
c0=(1,1,-1,-1,1,3); c=c0
for n in range(1,61):
    c=M(c)
    if c==c0: print("負の対照：30面体の輪の操作 M は",n,"回で戻る（分母なし）"); break

# Q2 の分母が NG だったので原因を調べる（事前登録の外。後から立てた確認）
print("--- Q2 の分母 ---")
print("二回の反射ごとの分母（n=1〜8）：",dens[:8])
exact5=all(dens[k]==5**(k+1) for k in range(60))
print("分母はちょうど 5^n（二回の反射で 25 を掛け、5 が一つ約分で消える）：",exact5)
print("分母は増え続け 1 に戻らない（＝元の整数の組に戻らない）：",all(dens[k]<dens[k+1] for k in range(59)))
