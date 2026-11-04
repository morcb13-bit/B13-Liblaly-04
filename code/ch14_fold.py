# 第14章　立体の点を φ で伸ばし、13 の番地で読む
import numpy as np
from collections import Counter
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
S=np.array(s)
print("S² ＝ 5I：",np.array_equal(S@S,5*np.eye(6,dtype=int)))
# すべて同じ偶奇の点の基底：b0=(1,…,1), b_k=2e_k (k=1..5)
B=np.zeros((6,6),dtype=int); B[:,0]=1
for k in range(1,6): B[k,k]=2
def phi_vec(a):   # φ·a = (a + S a)/2
    v=a+S@a; assert all(v%2==0); return v//2
Binv_num=np.round(np.linalg.inv(B)*2).astype(int)   # B^{-1} = Binv_num/2
def coords(a):    # a = B c を整数 c で
    c=Binv_num@a; assert all(c%2==0); return c//2
Phi=np.column_stack([coords(phi_vec(B[:,k])) for k in range(6)])
I=np.eye(6,dtype=int)
k1= np.array_equal(Phi@Phi,Phi+I)
print("K1 Φ（基底 b0, 2e1…2e5 での整数の行列）：\n",Phi,"\n   Φ² ＝ Φ + I：",k1,"→","OK" if k1 else "NG")
def cycles(M,p,dim=6,sub=None):
    n=p**dim; idx=np.arange(n,dtype=np.int64)
    X=np.zeros((n,dim),dtype=np.int64); r=idx.copy()
    for k in range(dim): X[:,k]=r%p; r//=p
    if sub is not None: X=X@sub.T % p   # 部分格子：座標を埋め込む
    Y=(X@M.T)%p
    if sub is not None:
        return None
    nxt=np.zeros(n,dtype=np.int64)
    for k in range(dim-1,-1,-1): nxt=nxt*p+Y[:,k]
    seen=np.zeros(n,dtype=bool); lens=Counter()
    for st in range(n):
        if seen[st]: continue
        L=0; x=st
        while not seen[x]: seen[x]=True; x=nxt[x]; L+=1
        lens[L]+=1
    return lens
P2=(Phi@Phi)
l2=cycles(P2,13); l1=cycles(Phi,13)
print("K2 13⁶ 通り：Φ² の輪",dict(l2),"／Φ の輪",dict(l1))
k2= l2=={1:1,14:344772} and l1=={1:1,28:172386}
print("   →","OK" if k2 else "NG")
# K3：真ん中の平面 a0=0, a1+…+a5=0（同じ偶奇の点）の部分格子
# 部分格子の基底を Φ の座標で：c で書いた点のうち a0=0, Σa=0 のもの
def a_of(c): return B@c
cand=[]
import itertools
for c in itertools.product(range(-2,3),repeat=6):
    a=a_of(np.array(c))
    if a[0]==0 and a[1:].sum()==0 and any(c): cand.append(np.array(c))
# 整数の基底を選ぶ（ヘルミート標準形の代わりに、小さいものから一次独立に4本）
basis=[]
for c in sorted(cand,key=lambda v:abs(v).sum()):
    M=np.array(basis+[c])
    if np.linalg.matrix_rank(M)>len(basis): basis.append(c)
    if len(basis)==4: break
Bp=np.array(basis).T   # 6×4
# 部分格子の上の Φ（4×4）：Φ Bp = Bp R を解き、整数か確かめる
R=np.linalg.lstsq(Bp,Phi@Bp,rcond=None)[0]; Rint=np.round(R).astype(int)
ok_int=np.array_equal(Bp@Rint,Phi@Bp)
# 基底が部分格子全体を張るか（指数1）：同じ範囲の候補がすべて整数係数で書けるか
spans=all(np.allclose(np.round(np.linalg.lstsq(Bp,c,rcond=None)[0]),np.linalg.lstsq(Bp,c,rcond=None)[0]) for c in cand)
l3=cycles(Rint@Rint,13,dim=4)
k3= ok_int and spans and l3=={1:1,14:2040}
print("K3 真ん中の平面：Φ が平面に留まる",ok_int,"、基底が平面全体を張る",spans,"、13⁴ 通りの Φ² の輪",dict(l3),"→","OK" if k3 else "NG")
# K4：初期値の違い
d=np.array([0,0,0,0,0,1])   # c の座標で 1 だけ違う
seq=[d]; 
for t in range(30): seq.append(Phi@seq[-1])
fib=all(np.array_equal(seq[t+2],seq[t+1]+seq[t]) for t in range(28))
back=[t for t in range(1,30) if np.array_equal(seq[t]%13,d%13)]
print("K4 差に Φ を掛け続けると、差は前の二つの和（φ 倍ずつ伸びる）：",fib,"／13 で読んで差が戻る歩：",back[:2],"→","OK" if fib and back and back[0]==28 else "NG")
print("   差の大きさ（係数の最大）：",[int(abs(seq[t]).max()) for t in range(0,30,4)])
# 負の対照：11 で読む
l11=cycles(P2,11)
print("負の対照：11 で読んだ Φ² の輪の長さ",dict(l11),"→","NG（長さがそろわない）" if len([k for k in l11 if k>1])>1 else "OK(!)")
# 事前登録の負の対照（11 で Φ²）は、長さが 5 にそろってしまい NG を返さなかった。
# 11 では x²−x−1 が二つの根 4, 8 に分かれ、どちらも φ² の位数が 5 になるため。後から立てた確認：
l11b=cycles(Phi,11)
print("後から立てた確認：11 で読んだ Φ の輪の長さ",dict(l11b),"（13 では長さ 28 だけ）")
