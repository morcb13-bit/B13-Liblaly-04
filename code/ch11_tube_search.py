# 輪の上に輪を重ねる一歩を、整数の組で探す（探索の記録）
import json, itertools
from collections import Counter
exec(open("ch05_close.py").read().split("# Q0")[0])   # s, dot5, add, sub
def cmp_pair(p,q):
    a=p[0]-q[0]; b=p[1]-q[1]
    if a==0: return (b>0)-(b<0)
    if b==0: return (a>0)-(a<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    big=5*a*a-b*b
    return (1 if a>0 else -1) if big>0 else ((1 if b>0 else -1) if big<0 else 0)
SIG=[tuple(v) for v in json.load(open("ch01_rt.json"))["V"]]
F=json.load(open("ch01_rt.json"))["F"]
C=[tuple(c) for c in json.load(open("ch02_ring.json"))["C"]]
STEPS=set()
for f in F:
    cen=tuple(sum(SIG[i][k] for i in f) for k in range(6)); STEPS.add(tuple(x//2 for x in cen))
TOUCH=(16,32)
VS=lambda c:{add(c,s) for s in SIG}
def check(A,B):
    sh=Counter(); bad=0
    for a in A:
        for b in B:
            d=sub(a,b); L=dot5(d,d); n=len(VS(a)&VS(b)); sh[n]+=1
            if n not in (0,4) or cmp_pair(L,TOUCH)<0: bad+=1
    return bad,sh
cands=set()
for a in C:
    for b in C:
        for f in STEPS:
            v=sub(add(a,f),b)
            if v[0]>0 or (v[0]==0 and sum(v[1:])>0): cands.add(v)
res=[]
for v in cands:
    B=[add(c,v) for c in C]
    bad,sh=check(C,B)
    if bad==0 and sh[4]>0: res.append((sh[4],v,(v[0],sum(v[1:]))))
res.sort(key=lambda r:(-r[0],r[2]))
for r in res[:12]: print(r)
print(len(res))

def shift(a,j):
    b=[0]*6; b[0]=a[0]
    for k in range(1,6): b[1+(k-1+j)%5]=a[k]
    return tuple(b)
v0=res[0][1]
def stack(seq,N):
    rings=[C]; off=(0,)*6
    for n in range(1,N):
        off=add(off,seq(n-1)); rings.append([add(c,off) for c in C])
    return rings
def tube_check(rings):
    bad=0; cont=Counter()
    for i in range(len(rings)):
        for j in range(i+1,len(rings)):
            b,sh=check(rings[i],rings[j]); bad+=b; cont[j-i]+=sh[4]
    return bad,dict(cont)
for name,seq in [("同じ一歩",lambda n:v0),("右回り",lambda n:shift(v0,n)),("左回り",lambda n:shift(v0,-n)),("一つおき",lambda n:shift(v0,2*n))]:
    print(name,tube_check(stack(seq,8)))
