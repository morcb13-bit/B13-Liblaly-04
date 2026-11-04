# 第12章　立体のペンローズ（黄金のひし形の六面体）の上の整数の波
import json, math
from collections import Counter, defaultdict, deque
from ch12_tiling import faces_of
cells=[(tuple(x),tuple(t)) for x,t in json.load(open("ch12_tiling.json"))["cells"]]
F=defaultdict(list)
for n,cl in enumerate(cells):
    for f in faces_of(cl): F[f].append(n)
mult=Counter(len(v) for v in F.values())
print("C0 六面体",len(cells),"個、面の持ち主の数ごと",dict(mult),"→","OK" if max(mult)<=2 else "NG")
fat=sum(1 for _,t in cells if t in {(0,1,2),(0,2,3),(0,3,4),(0,4,5),(0,1,5)} )
nb=defaultdict(list)
for f,v in F.items():
    if len(v)==2: a,b=v; nb[a].append(b); nb[b].append(a)
# いちばん大きいつながり
comp=[-1]*len(cells); best=None
for s in range(len(cells)):
    if comp[s]>=0: continue
    q=[s]; comp[s]=s; members=[s]
    while q:
        a=q.pop()
        for b in nb[a]:
            if comp[b]<0: comp[b]=s; q.append(b); members.append(b)
    if best is None or len(members)>len(best): best=members
keep=set(best)
deg={a:len(nb[a]) for a in keep}
L=1
for d in set(deg.values()): L=L*d//math.gcd(L,d)
print("   つながった六面体",len(keep),"個、隣の数ごと",dict(sorted(Counter(deg.values()).items())),"、L =",L)
# 有向辺：(a,b) は a から b へ。値は「b に入ってくる値」
def step(inc,Lc,exact=True):
    out={}
    for b in keep:
        d=deg[b]; s=sum(inc.get((a,b),0) for a in nb[b])
        for a in nb[b]:
            if exact: v=(2*Lc//d)*s-Lc*inc.get((a,b),0)
            else: v=((2*Lc)//d)*s-Lc*inc.get((a,b),0)
            if v: out[(b,a)]=v
    return out
def back(inc,Lc):
    # 逆向き：辺を入れ替えて同じ散らし方を掛け、L² で割る
    prev={}
    for b in keep:
        d=deg[b]; outs=[inc.get((b,a),0) for a in nb[b]]; s=sum(outs)
        for a,o in zip(nb[b],outs):
            num=(2*Lc//d)*s-Lc*o
            if num%(Lc*Lc): return None
            if num: prev[(a,b)]=num//(Lc*Lc)
    return prev
# 出発：中心にいちばん近い、隣の数 6 の六面体に、すべての辺から 1 を入れる
import numpy as np
from ch12_tiling import phys
def ctr(n):
    x,(i,j,k)=cells[n]; return phys(x)+(phys([1 if m==i else 0 for m in range(6)])+phys([1 if m==j else 0 for m in range(6)])+phys([1 if m==k else 0 for m in range(6)]))/2
start=min((n for n in keep if deg[n]==6),key=lambda n:np.linalg.norm(ctr(n)))
inc0={(a,start):1 for a in nb[start]}
dist={start:0}; q=deque([start])
while q:
    a=q.popleft()
    for b in nb[a]:
        if b not in dist: dist[b]=dist[a]+1; q.append(b)
T=12
hist=[inc0]; inc=inc0
for t in range(T): inc=step(inc,L); hist.append(inc)
S=[sum(h.values()) for h in hist]; Q=[sum(v*v for v in h.values()) for h in hist]
w1=all(S[t+1]==L*S[t] for t in range(T)); w2=all(Q[t+1]==L*L*Q[t] for t in range(T))
print("W1 和が毎歩 L 倍（",T,"歩）：",w1,"→","OK" if w1 else "NG")
print("W2 二乗和が毎歩 L² 倍：",w2,"→","OK" if w2 else "NG")
x=hist[-1]; ok3=True
for t in range(T):
    x=back(x,L)
    if x is None: ok3=False; break
ok3= ok3 and {k:v for k,v in x.items() if v}==inc0
print("W3 逆に",T,"歩走らせて最初に戻る：",ok3,"→","OK" if ok3 else "NG")
w4=all(dist.get(b,10**9)<=t for t,h in enumerate(hist) for (a,b) in h)
print("W4 t 歩目の値は出発から t 歩以内の六面体だけ：",w4,"→","OK" if w4 else "NG")
# 記録：2色に塗れるか
col={start:0}; q=deque([start]); bip=True
while q:
    a=q.popleft()
    for b in nb[a]:
        if b not in col: col[b]=1-col[a]; q.append(b)
        elif col[b]==col[a]: bip=False
print("[記録] 2色に塗り分けられる：",bip)
if bip:
    alt=all(col[b]==t%2 for t,h in enumerate(hist) for (a,b) in h)
    print("[記録] 値は毎歩反対の色の六面体へ移る：",alt)
reach=[len({b for (a,b) in h}) for h in hist]
print("[記録] 各歩で値のある六面体の数：",reach)
# 負の対照：L=12 のまま切り捨て
h=inc0; Qn=[sum(v*v for v in h.values())]
for t in range(T): h=step(h,12,exact=False); Qn.append(sum(v*v for v in h.values()))
n2=all(Qn[t+1]==144*Qn[t] for t in range(T))
y=h; okb=True
for t in range(T):
    y=back(y,12)
    if y is None: okb=False; break
print("負の対照（L=12、2·12/d を切り捨て）：二乗和が毎歩144倍",n2,"／逆に戻れる",okb and y=={k:v for k,v in inc0.items()},"→","NG" if not n2 or not okb else "OK(!)")
json.dump({"keep":sorted(keep),"start":start,"L":L,"T":T,"col":{str(k):v for k,v in col.items()},
           "hist":[{f"{a},{b}":v for (a,b),v in h.items()} for h in hist[:9]]},open("ch12_wave.json","w"))
