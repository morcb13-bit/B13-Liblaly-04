# 第13章の見立ての確かめ：入ってきた面の向かいの面から抜ける（まっすぐ抜ける）と、経路はどこを進むか
import json, numpy as np, itertools
from collections import defaultdict, Counter
import ch12_tiling as T
P,cells=T.build(10.0)
F=defaultdict(list)
for n,cl in enumerate(cells):
    for f in T.faces_of(cl): F[f].append(n)
E=lambda m:T.phys([1 if r==m else 0 for r in range(6)])
def ctr(n):
    x,(i,j,k)=cells[n]; return T.phys(x)+(E(i)+E(j)+E(k))/2
def opposite(n,f):
    x,(i,j,k)=cells[n]; fs=T.faces_of(cells[n]); idx=fs.index(f); return fs[idx^1]
def rail(n,f):   # 面 f から入って、向かいの面から抜け続ける
    seq=[n]; seen={(n,f)}
    while True:
        g=opposite(n,f); own=[m for m in F[g] if m!=n]
        if not own: break
        n=own[0]; f=g; seq.append(n)
        if (n,f) in seen: break
        seen.add((n,f))
    return seq,f[1]
# 面の種類 {i,j} はレールの上で変わらない
start=min(range(len(cells)),key=lambda n:np.linalg.norm(ctr(n)))
same=0;tot=0;dirs=Counter();lens=[]
rails=[]
for f in T.faces_of(cells[start]):
    s,ij=rail(start,f); tot+=1
    rails.append((s,f[1]))
# 多くの六面体から：まっすぐのずれ（直線からの最大のはなれ）と向き
dev=[];ang=[]
for n in range(0,len(cells),25):
    for f in T.faces_of(cells[n])[:6:2]:
        s,ij=rail(n,f)
        if len(s)<12: continue
        pts=np.array([ctr(m) for m in s]); d=pts[-1]-pts[0]; d/=np.linalg.norm(d)
        perp=pts-pts[0]-np.outer((pts-pts[0])@d,d); dev.append(np.abs(np.linalg.norm(perp,axis=1)).max())
        i,j=f[1]; ax=np.cross(E(i),E(j)); ax/=np.linalg.norm(ax); ang.append(np.degrees(np.arccos(min(1,abs(d@ax)))))
        lens.append(len(s))
print("レールの数",len(dev),"長さ（中央値）",int(np.median(lens)))
print("直線からの最大のはなれ：中央値",round(float(np.median(dev)),3),"最大",round(float(max(dev)),3),"（六面体の辺の長さ 1）")
print("進む向きと e_i×e_j の向きの角度：中央値",round(float(np.median(ang)),2),"度、最大",round(float(max(ang)),2),"度")
