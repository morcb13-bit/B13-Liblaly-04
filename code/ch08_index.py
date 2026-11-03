# 第8章　前の本の指数は、立体の高さの和を5で割った余り
import json, itertools
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
index=lambda z:(z[0]-z[1]+z[2]-z[3])%5
wrong=lambda z:(z[0]+z[1]+z[2]+z[3])%5
# L1
l1=[index(zpow(k)) for k in range(1,11)]
ok1= l1==[(-1)**k%5 for k in range(1,11)] and all(index(zpow(2*k))==1 for k in range(1,6))
print("L1 ζ^1〜ζ^10 の指数：",l1,"→","OK" if ok1 else "NG")
# L2
box=list(itertools.product(range(-2,3),repeat=6))
ok2=sum(index(shadow(a))==sum(a[1:])%5 for a in box)
print("L2 範囲内",len(box),"通りで 指数＝高さの和 mod 5：",ok2,"→","OK" if ok2==len(box) else "NG")
# L3
SIG=[tuple(v) for v in json.load(open("ch01_rt.json"))["V"]]
C=[tuple(c) for c in json.load(open("ch02_ring.json"))["C"]]
ROWS=json.load(open("ch03_stack.json"))["rows"]
pts={add(c,s) for c in C for s in SIG}|set(C)|{add(tuple(c),s) for r in ROWS for c in r for s in SIG}
ok3=sum(index(shadow(a))==sum(a[1:])%5 for a in pts)
print("L3 これまでの点",len(pts),"個：",ok3,"→","OK" if ok3==len(pts) else "NG")
# L4
up={index(shadow(c)) for c in C[0::2]}; dn={index(shadow(c)) for c in C[1::2]}
mid={index(shadow(a)) for a in pts if a[0]==0 and sum(a[1:])==0 and a in {add(c,s) for c in C for s in SIG}}
ok4= len(up)==1 and len(dn)==1 and up!=dn and mid=={0}
print("L4 上の段の中心の指数",up,"下の段",dn,"真ん中の平面の頂点",mid,"→","OK" if ok4 else "NG")
print("   上の段の中心の高さの和",{sum(c[1:]) for c in C[0::2]},"下の段",{sum(c[1:]) for c in C[1::2]})
# 30面体1個の頂点：高さの和ごとの指数
from collections import Counter
cnt=Counter((sum(s[1:]),index(shadow(s))) for s in SIG)
print("[記録] 30面体の頂点（中心が原点）の (高さの和, 指数) と個数：",sorted(cnt.items()))
# 負の対照
n1=sum(wrong(shadow(a))==sum(a[1:])%5 for a in box)
n2=sum(index(shadow(a))==a[0]%5 for a in box)
print("負の対照1：a0+a1+a2+a3 を使う：",n1,"/",len(box),"→ L2 は","NG" if n1<len(box) else "OK(!)")
print("負の対照2：高さの和の代わりに a0：",n2,"/",len(box),"→ L2 は","NG" if n2<len(box) else "OK(!)")
