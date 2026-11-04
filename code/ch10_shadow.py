# 第10章　影と切り口。輪の外側の五角形を真上から見た影と、真ん中の切り口を比べる
import json
from collections import Counter
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
def zmul(x,y):
    r=(0,0,0,0)
    for i in range(4): r=add(r,sc(x[i],y)); y=zt(y)
    return r
height=lambda a:(a[0],sum(a[1:]))
SIG=[tuple(v) for v in json.load(open("ch01_rt.json"))["V"]]
C=[tuple(c) for c in json.load(open("ch02_ring.json"))["C"]]
TOP=(1,1,1,1,1,1); BOT=sc(-1,TOP)
assert TOP in SIG and BOT in SIG
def around(apex):   # 五つ角の頂点と辺でつながる5頂点（符号が一つだけ違う）
    return [v for v in SIG if sum(1 for i in range(6) if v[i]!=apex[i])==1]
up=lambda k:C[k][0]>0
cut=[[add(C[k],v) for v in around(BOT if up(k) else TOP)] for k in range(10)]
far=[[add(C[k],v) for v in around(TOP if up(k) else BOT)] for k in range(10)]
# 前の本の円環の二つの候補（第7章）
phi2=(2,0,1,-1)
def book3_ring(off): return [[add(zmul(phi2,zpow(k)),zpow(k+off+2*j)) for j in range(5)] for k in range(10)]
cand={off:{sc(2,v) for P in book3_ring(off) for v in P} for off in (0,1)}
shF={shadow(v) for P in far for v in P}; shC={shadow(v) for P in cut for v in P}
g1= len(shF)==30 and shF==cand[1]
print("G1 外側の五角形の影の頂点",len(shF),"個 ＝ もう一方の候補（off=1）の2倍：",shF==cand[1],"→","OK" if g1 else "NG")
# 事前登録の「30個」は見込み違い：外側の五角形は隣と頂点を共有しないので 10×5=50 個になる（G2 と同じ事情）。後から直した判定：
g1b= len(shF)==50 and shF==cand[1]
print("G1′（個数を 50 に直した判定）：",len(shF),"個、off=1 の2倍と一致 →","OK" if g1b else "NG")
print("負の対照：切り口の円環に同じ照合：",shC==cand[1],"→","NG" if shC!=cand[1] else "OK(!)","（切り口は off=0 と一致：",shC==cand[0],"）")
def share(R): return [len({shadow(v) for v in R[k]}&{shadow(v) for v in R[(k+1)%10]}) for k in range(10)]
sF=share(far); sC=share(cut)
g2= sF==[0]*10 and sC==[2]*10
print("G2 隣どうしの共有頂点：影",sF,"切り口",sC,"→","OK" if g2 else "NG")
def csum(P):
    z=(0,0,0,0)
    for v in P: z=add(z,shadow(v))
    return z
g3= all(csum(far[k])==sc(5,shadow(C[k])) and csum(cut[k])==sc(5,shadow(C[k])) for k in range(10))
print("G3 5頂点の影の和 ＝ 中心の影の5倍（影・切り口とも）：",g3,"→","OK" if g3 else "NG")
def rel(P,k): return {sub(shadow(v),shadow(C[k])) for v in P}
z5=zpow(5)
g4=sum(rel(far[k],k)=={zmul(z5,d) for d in rel(cut[k],k)} for k in range(10))
print("G4 影の五角形 ＝ 切り口の五角形を中心で ζ⁵ 倍：",g4,"/ 10 →","OK" if g4==10 else "NG")
hF=Counter((up(k),height(v)) for k in range(10) for v in far[k]); hC={height(v) for P in cut for v in P}
hu={h for (u,h) in hF if u}; hd={h for (u,h) in hF if not u}
g5= len(hu)==1 and len(hd)==1 and next(iter(hd))==sc(-1,next(iter(hu))) and hC=={(0,0)}
print("G5 影の五角形の高さの組：上の段",hu,"下の段",hd,"／切り口",hC,"→","OK" if g5 else "NG")
json.dump({"far":far,"cut":cut,"C":C},open("ch10_shadow.json","w"))
