# 第6章　番地をふる。6個の整数 ＝ 平面の番地（Z[ζ10] の4個）＋ 高さの組（2個）
import json, itertools
from collections import Counter
def zt(c): return (-c[3], c[0]+c[3], c[1]-c[3], c[2]+c[3])
def zpow(n):
    c=(1,0,0,0)
    for _ in range(n%10): c=zt(c)
    return c
def addv(a,b): return tuple(x+y for x,y in zip(a,b))
def sc(k,a): return tuple(k*x for x in a)
Z2=[zpow(2*k) for k in range(1,6)]    # ζ², ζ⁴, ζ⁶, ζ⁸, ζ¹⁰=1
def shadow(a):
    z=(0,0,0,0)
    for k in range(5): z=addv(z,sc(a[k+1],Z2[k]))
    return z
def height(a): return (a[0], sum(a[1:]))
# 平面の番地 z を ζ²,ζ⁴,ζ⁶,ζ⁸ で書き直す（整数の基底の取り替え）：4×4 の整数行列を解く
B=[Z2[k] for k in range(4)]     # 列ベクトル
def det(M):
    if len(M)==1: return M[0][0]
    return sum((-1)**j*M[0][j]*det([r[:j]+r[j+1:] for r in M[1:]]) for j in range(len(M)))
Mt=[[B[k][r] for k in range(4)] for r in range(4)]   # z = Mt · b
dM=det(Mt)
print("基底 ζ²,ζ⁴,ζ⁶,ζ⁸ の行列式：",dM,"（±1 なら整数で取り戻せる）")
def adj(M):
    n=len(M)
    return [[(-1)**(i+j)*det([r[:i]+r[i+1:] for k,r in enumerate(M) if k!=j]) for j in range(n)] for i in range(n)]
A=adj(Mt)
def to_b(z):   # b = Mt^{-1} z（整数）
    return tuple(sum(A[i][j]*z[j] for j in range(4))//dM for i in range(4))
def back(z,h,mod=5):
    b=to_b(z); a0,S=h; r=S-sum(b)
    if r%mod: return None
    t=r//5 if mod==5 else r//mod
    return (a0,)+tuple(x+t for x in b)+(t,)
def roundtrip(a,mod=5): return back(shadow(a),height(a),mod)==a
# 点を集める
SIG=[tuple(v) for v in json.load(open("ch01_rt.json"))["V"]]
C2=[tuple(c) for c in json.load(open("ch02_ring.json"))["C"]]
ROWS=json.load(open("ch03_stack.json"))["rows"]
pts=set(SIG)|set(C2)|{addv(c,s) for c in C2 for s in SIG}|{tuple(c) for r in ROWS for c in r}
print("A1 取り戻しの式：b = (ζ²,ζ⁴,ζ⁶,ζ⁸ での係数)、t = (高さの和 − b の和)/5、a = (a0, b1+t, b2+t, b3+t, b4+t, t)")
okA2=sum(roundtrip(a) for a in pts)
print("A2 行って戻る：",okA2,"/",len(pts),"→","OK" if okA2==len(pts) else "NG")
# A3：範囲の中で、立体の点になれる組の割合
cnt=Counter()
for z in itertools.product(range(-2,3),repeat=4):
    for S in range(-10,10):
        cnt[back(z,(0,S)) is not None]+=1
tot=cnt[True]+cnt[False]
print("A3 平面の番地 5⁴ 通り × 高さの和 20 通り：立体の点になれる",cnt[True],"/",tot,"（5 個に 1 個なら",tot//5,"）→","OK" if cnt[True]*5==tot else "NG")
# A4：輪の中心の影 ＝ 2φ²ζ^k
phi2=(2,0,1,-1)
targ=[sc(2,phi2)]
for k in range(9): targ.append(zt(targ[-1]))
sh=[shadow(c) for c in C2]
m=[targ.index(sh[0])] if sh[0] in targ else []
ok4= bool(m) and all(sh[k]==targ[(k+m[0])%10] for k in range(10))
print("A4 輪の中心の影：",sh[0],"＝ 2φ²·ζ^"+(str(m[0]) if m else "?"),"から一つずつ ζ 倍：",ok4,"→","OK" if ok4 else "NG")
# A5：面の番地 {i,j} と表裏
F=json.load(open("ch01_rt.json"))["F"]; V=[tuple(v) for v in json.load(open("ch01_rt.json"))["V"]]
band=json.load(open("ch04_bands.json"))["band"]
names=[]
for f in F:
    cen=tuple(sum(V[i][k] for i in f) for k in range(6))
    ij=tuple(k for k in range(6) if cen[k]==0)
    # 表裏：残る4軸の符号の並び（中心の向き）
    names.append((ij,tuple(1 if x>0 else -1 for x in cen if x!=0)))
uniq=len(set(names))==30 and all(len(n[0])==2 for n in names)
pairs=Counter(n[0] for n in names)
# 段と矛盾しない：表裏で段が鏡写し（段 b ↔ 4−b）
fl={n:b for n,b in zip(names,band)}
mirror=all(fl[(ij,tuple(-x for x in sg))]==4-b for (ij,sg),b in fl.items())
print("A5 面の番地 30 個が重ならない：",uniq,"／軸の組15通りに2枚ずつ：",set(pairs.values())=={2},"／表裏で段が b ↔ 4−b：",mirror,"→","OK" if uniq and mirror else "NG")
# 負の対照
shc=Counter(shadow(a) for a in SIG)
print("負の対照1：平面の番地だけでは、第1章の頂点32が",len(shc),"番地に重なる →","NG" if len(shc)<32 else "OK(!)")
ok3=sum(roundtrip(a,3) for a in pts)
print("負の対照2：5 ではなく 3 で割った余りにすると、行って戻るのは",ok3,"/",len(pts),"→","NG" if ok3<len(pts) else "OK(!)")
