# 終章　600胞体（イコシアン 120 個）を w で切り、二枚の切り口の影を重ねると黄金のひし形 30 面体になるか
# 数はすべて a + b√5（a, b は整数）の組で持ち、座標は 4 倍しておく。足し算・掛け算・符号だけで判定する。
import itertools, json
from collections import Counter
mul=lambda x,y:(x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
add=lambda x,y:(x[0]+y[0],x[1]+y[1]); neg=lambda x:(-x[0],-x[1]); Z=(0,0)
def sgn(x):   # a + b√5 の符号
    a,b=x
    if b==0: return (a>0)-(a<0)
    if a==0: return (b>0)-(b<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    t=a*a-5*b*b
    return (1 if a>0 else -1) if t>0 else (1 if b>0 else -1)
cmp=lambda x,y:sgn(add(x,neg(y)))
def icosians():
    P=set()
    for i in range(4):
        for s in (4,-4):
            v=[Z]*4; v[i]=(s,0); P.add(tuple(v))
    for sg in itertools.product((2,-2),repeat=4): P.add(tuple((s,0) for s in sg))
    base=[Z,(2,0),(-1,1),(1,1)]          # 0, 1, 1/φ, φ（×2 したものを、さらに ×2 して 4 倍）
    even=[p for p in itertools.permutations(range(4)) if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])%2==0]
    for p in even:
        for sg in itertools.product((1,-1),repeat=3):
            v=[Z]*4; vals=[base[0]]+[x if s>0 else neg(x) for x,s in zip(base[1:],sg)]
            for k in range(4): v[p[k]]=vals[k]
            P.add(tuple(v))
    return sorted(P)
d2=lambda a,b:(lambda d:add(add(mul(d[0],d[0]),mul(d[1],d[1])),mul(d[2],d[2])))([add(x,neg(y)) for x,y in zip(a,b)])
n2=lambda a:d2(a,(Z,Z,Z))
# 第1章の 30 面体：√5|x|² = sq·√5 + cross  →  a + b√5 の組では (cross, sq)
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j); s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def rt_len(v): return (sum(v[i]*v[j]*s[i][j] for i in range(6) for j in range(6) if i!=j), sum(x*x for x in v))
def check_rt(I,Dd,label):
    U=I+Dd; dists=[(d2(a,b),ia,ib) for (ia,a),(ib,b) in itertools.combinations(list(enumerate(U)),2)]
    m=dists[0][0]
    for d,_,_ in dists:
        if cmp(d,m)<0: m=d
    E=[(a,b) for d,a,b in dists if d==m]
    deg=Counter(x for e in E for x in e)
    ni=len(I)
    bip=all((a<ni)!=(b<ni) for a,b in E)
    degok=all(deg[k]==5 for k in range(ni)) and all(deg[k]==3 for k in range(ni,len(U)))
    nb={k:{b if a==k else a for a,b in E if k in (a,b)} for k in range(len(U))}
    faces=sum(1 for a,b in itertools.combinations(range(ni),2) if len(nb[a]&nb[b])==2)
    return len(E),bip,degok,faces
if __name__=="__main__":
    print("基準（走らせる前に決めたもの）")
    print(" R1 イコシアン 120 個を実部 w で切ると、切り口の点の数は 1・12・20・12・30・12・20・12・1")
    print(" R2 w = 1/2 の切り口（20 点）と w = 1/(2φ) の切り口（12 点）を、w を落として重ねる（影）。いちばん短い距離の組はちょうど 60 本で、")
    print("    どれも 12 点の側と 20 点の側を結び、12 点の側はどれも 5 本、20 点の側はどれも 3 本。12 点の組で共通の相手をちょうど 2 つ持つものが 30 組（ひし形の面）")
    print(" R3 12 点の側と 20 点の側の、中心からの長さの二乗の比が、第 1 章の 30 面体の S＝5 の頂点と S＝3 の頂点の比と、整数の組でちょうど一致する")
    print(" 負の対照：12 点の切り口を w = φ/2 のほうに替えると、R3 は NG を返す")
    print()
    Q=icosians(); print("[記録] イコシアン",len(Q),"個、長さの二乗がすべて 16（4 倍の座標で 1）：",all(n2(q[1:])==add((16,0),neg(mul(q[0],q[0]))) for q in Q))
    # w を大きい順に並べる（比較は整数の組で）
    import functools
    ws=sorted({q[0] for q in Q},key=functools.cmp_to_key(cmp),reverse=True)
    sl={w:[q[1:] for q in Q if q[0]==w] for w in ws}
    cnt=[len(sl[w]) for w in ws]
    r1=cnt==[1,12,20,12,30,12,20,12,1]
    print("R1 w（4 倍、a+b√5 の組）",ws,"\n   切り口の点の数",cnt,"→","OK" if r1 else "NG")
    Dd=sl[(2,0)]; I=sl[(-1,1)]; I2=sl[(1,1)]
    ne,bip,degok,faces=check_rt(I,Dd,"")
    r2=ne==60 and bip and degok and faces==30
    print("R2 影を重ねた 32 点：最短の組",ne,"本、12 点と 20 点を結ぶ",bip,"、本数 5/3",degok,"、ひし形の面",faces,"→","OK" if r2 else "NG")
    R=json.load(open("ch01_rt.json")); V=[tuple(v) for v in R["V"]]; K=R["K"]
    i5={rt_len(v) for v,k in zip(V,K) if k==5}; i3={rt_len(v) for v,k in zip(V,K) if k==3}
    a,b=list(i5)[0],list(i3)[0]; ci={n2(p) for p in I}; cd={n2(p) for p in Dd}; c,d=list(ci)[0],list(cd)[0]
    r3=len(i5)==len(i3)==len(ci)==len(cd)==1 and mul(c,b)==mul(d,a)
    print("R3 600胞体の影：12 点",c,"／20 点",d,"　第1章の30面体：S=5",a,"／S=3",b,"　比が一致",mul(c,b)==mul(d,a),"→","OK" if r3 else "NG")
    c2=n2(I2[0]); ok_n=mul(c2,b)==mul(d,a)
    print("負の対照：12 点を w = φ/2 の切り口",c2,"に替える → 比が一致",ok_n,"→","OK(!)" if ok_n else "NG")
    f=lambda x:x[0]+x[1]*5**0.5
    json.dump({"ws":[f(w)/4 for w in ws],"slices":[[[round(f(x)/4,6) for x in p] for p in sl[w]] for w in ws]},open("ch19_shadow.json","w"))
