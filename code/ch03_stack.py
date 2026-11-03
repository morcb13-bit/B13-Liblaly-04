# 第3章　1から13まで積む。数 n を n 個の30面体の鎖で表す（奇数はジグザグ、偶数はまっすぐ）
# 長さは半分を単位にした6次元の整数の組。判定は整数だけ。
import itertools, json
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b):  # √5·(a·b) を (√5の係数, 整数) で
    return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
def height(a): return (a[0], sum(a[1:]))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)
SIG=[p for p in itertools.product([-1,1],repeat=6)
     if sum(p[i]*p[j]*s[i][j] for i in range(6) for j in range(i+1,6)) in (5,3)]
def verts(c): return {add(c,p) for p in SIG}
def shift(a,j):
    b=[0]*6; b[0]=a[0]
    for k in range(1,6): b[1+(k-1+j)%5]=a[k]
    return tuple(b)
def gt(x,y):  # x√5+X > y√5+Y を整数だけで：(p√5 > q) の比べ方
    p=x[0]-y[0]; q=y[1]-x[1]          # p√5 > q ?
    if p>=0 and q<0: return True
    if p<=0 and q>=0: return False
    if p>0: return 5*p*p>q*q           # 両辺正
    return 5*p*p<q*q                   # p<0, q<0：-p√5 < -q
TOUCH=(16,32)
u0=(0,2,2,0,-2,-2); u1=shift(u0,1); u2=shift(u0,2); w=neg(shift(u0,3))  # w：u0 と u1 のまん中の向き（36°）
def chain(n,a=u0,b=u1):
    # 奇数：a, b を交互（ジグザグ）。偶数：w だけ（まっすぐ。ジグザグの進む向きと同じ向き）
    c=[(0,)*6]
    for k in range(n-1):
        c.append(add(c[-1], w if n%2==0 else (a if k%2==0 else b)))
    return c
def check(c):
    n=len(c)
    adj=[len(verts(c[k])&verts(c[k+1])) for k in range(n-1)]
    far=[(len(verts(c[a])&verts(c[b])), dot5(sub(c[b],c[a]),sub(c[b],c[a]))) for a in range(n) for b in range(a+2,n)]
    s1=all(x==4 for x in adj)
    s2=all(sh==0 and gt(L,TOUCH) for sh,L in far)
    s3=all(height(x)==(0,0) for x in c)
    steps=[sub(c[k+1],c[k]) for k in range(n-1)]
    folds=set(dot5(neg(steps[k]),steps[k+1]) for k in range(n-2))
    return s1,s2,s3,folds
print("一歩 u0 =",u0,"高さ",height(u0),"長さの組",dot5(u0,u0))
print("u1 =",u1,"（e1〜e5 を一つずらす：72°）")
ANG108=dot5(neg(u0),u1); ANG180=dot5(neg(w),w)
print("w =",w,"（u0 と u1 のまん中の向き）高さ",height(w),"長さの組",dot5(w,w))
print("折れ角108°の組（−u0・u1）：",ANG108,"／180°の組（−w・w）：",ANG180)
ALL=[];tot=0;rows=[]
for n in range(1,14):
    c=chain(n); s1,s2,s3,f=check(c); tot+=n
    kind="ジグザグ" if n%2 and n>2 else ("まっすぐ" if n>1 else "1個")
    s4 = f==set() or (n%2==1 and f=={ANG108}) or (n%2==0 and f=={ANG180})
    rows.append((n,s1,s2,s3,s4))
    print(f"n={n:2d} {kind:5s} S1 {'OK' if s1 else 'NG'}  S2 {'OK' if s2 else 'NG'}  S3 {'OK' if s3 else 'NG'}  S4 {'OK' if s4 else 'NG'}  折れ角の組 {sorted(f)}")
    ALL.append(c)
print("S1〜S4 すべての n で OK：",all(all(r[1:]) for r in rows))
print("S5 30面体の数：",tot,"（1+2+…+13 =",13*14//2,"）→","OK" if tot==91 else "NG")
fib=[1,1]
while fib[-1]<13: fib.append(fib[-1]+fib[-2])
ok5=all(len(ALL[fib[i]-1])==len(ALL[fib[i-1]-1])+len(ALL[fib[i-2]-1]) for i in range(2,len(fib)))
print("   段",fib,"：鎖の個数が前の二つの和 →","OK" if ok5 else "NG")
# 閉じる：同じ向きに曲がり続けると5歩で閉じる（第2章の上の段の輪）
c=[(0,)*6]
for k in range(5): c.append(add(c[-1],shift(u0,k)))
print("u0,u1,u2,u3,u4 の5歩で戻る：",c[5]==(0,)*6)
# 負の対照：u0 と u2 を交互（折れ角36°）
nc=chain(5,u0,u2)
s1,s2,s3,f=check(nc)
d=dot5(sub(nc[2],nc[0]),sub(nc[2],nc[0]))
print("負の対照（u0,u2 交互）：S2",'OK(!)' if s2 else 'NG',"／一つ飛びの中心の長さの組",d,"接する距離",TOUCH,"より近い：",gt(TOUCH,d))
json.dump({"rows":ALL},open("ch03_stack.json","w"))
