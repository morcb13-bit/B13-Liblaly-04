# 第16章　正十二面体（30面体の S=3 の頂点 20 個）の中に、正四面体と立方体を整数だけで探す
import json, itertools
from collections import Counter
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b): return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
def cmp_pair(p,q):
    a=p[0]-q[0]; b=p[1]-q[1]
    if a==0: return (b>0)-(b<0)
    if b==0: return (a>0)-(a<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    big=5*a*a-b*b
    return (1 if a>0 else -1) if big>0 else ((1 if b>0 else -1) if big<0 else 0)
sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
R=json.load(open("ch01_rt.json")); V=[tuple(v) for v in R["V"]]; K=R["K"]
D20=[v for v,k in zip(V,K) if k==3]; I12=[v for v,k in zip(V,K) if k==5]
d2=lambda a,b:dot5(sub(a,b),sub(a,b))
def tetras(P):
    out=[]
    for q in itertools.combinations(range(len(P)),4):
        L={d2(P[a],P[b]) for a,b in itertools.combinations(q,2)}
        if len(L)==1: out.append(frozenset(q))
    return out
T=tetras(D20)
print("A1 正十二面体の20頂点の中の正四面体：",len(T),"組 →","OK" if len(T)==10 else "NG")
print("   正四面体の辺の長さの組：",{d2(D20[min(t)],D20[sorted(t)[1]]) for t in T})
def is_cube(idx):
    P=[D20[i] for i in idx]
    L=list({d2(a,b) for a,b in itertools.combinations(P,2)})
    m=L[0]
    for x in L[1:]:
        if cmp_pair(x,m)<0: m=x
    nb=[sum(1 for b in P if b!=a and d2(a,b)==m) for a in P]
    edges=sum(1 for a,b in itertools.combinations(P,2) if d2(a,b)==m)
    return all(x==3 for x in nb) and edges==12
neg=lambda v:tuple(-x for x in v)
cubes=[]; mirror_ok=0
for a,b in itertools.combinations(T,2):
    if a&b: continue
    if is_cube(a|b):
        cubes.append(a|b)
        if {neg(D20[i]) for i in a}=={D20[i] for i in b}: mirror_ok+=1
print("A2 二つの正四面体を合わせた立方体：",len(cubes),"組、うち互いの鏡写し（符号の反転）",mirror_ok,"→","OK" if len(cubes)==5 and mirror_ok==5 else "NG")
ct=Counter(i for t in T for i in t); cc=Counter(i for c in cubes for i in c)
a3= set(ct.values())=={2} and set(cc.values())=={2} and len(ct)==20 and len(cc)==20
print("A3 各頂点が正四面体に",set(ct.values()),"回、立方体に",set(cc.values()),"回 →","OK" if a3 else "NG")
sh=[len(a&b) for a,b in itertools.combinations(cubes,2)]
anti=all(all(neg(D20[i]) in {D20[j] for j in (a&b)} for i in (a&b)) for a,b in itertools.combinations(cubes,2))
print("A4 二つの立方体が共有する頂点の数：",Counter(sh),"／共有する2点は向かい合う：",anti,"→","OK" if set(sh)=={2} and anti else "NG")
TN=tetras(I12)
print("負の対照：正二十面体の12頂点の中の正四面体：",len(TN),"組 →","NG" if not TN else "OK(!)")
json.dump({"D20":D20,"T":[sorted(t) for t in T],"cubes":[sorted(c) for c in cubes]},open("ch16_atoms.json","w"))
