# 第15章　立体の風洞：鏡で対称な立体のペンローズに風を通し、物体のまわりの横向きの力を見る
import json, math, numpy as np
from collections import defaultdict, Counter
import ch12_tiling as T
PERM=[0,4,3,2,1,5]          # 鏡：e1↔e4, e2↔e3
def mirror_cell(c):
    x,t=c; y=tuple(x[PERM[k]] for k in range(6)); return (y,tuple(sorted(PERM[k] for k in t)))
def build(gam,R=8.0,strict=False):
    T.GAM=gam; T.STRICT=strict; P,cells=T.build(R); return cells
def carrier(cells):
    F=defaultdict(list)
    for n,cl in enumerate(cells):
        for f in T.faces_of(cl): F[f].append(n)
    return F
E=lambda m:T.phys([1 if r==m else 0 for r in range(6)])
def ctr(c):
    x,(i,j,k)=c; return T.phys(x)+(E(i)+E(j)+E(k))/2
def run(cells,F,removed,sources,steps,check_back=False):
    keep=[n for n in range(len(cells)) if n not in removed]; ks=set(keep)
    nb=defaultdict(list)
    for f,v in F.items():
        if len(v)==2 and v[0] in ks and v[1] in ks: a,b=v; nb[a].append(b); nb[b].append(a)
    deg={a:len(nb[a]) for a in nb}
    L=1
    for d in set(deg.values()): L=L*d//math.gcd(L,d)
    def step(inc):
        out={}
        for b in nb:
            d=deg[b]; s=sum(inc.get((a,b),0) for a in nb[b])
            if s==0 and not any(inc.get((a,b),0) for a in nb[b]): continue
            for a in nb[b]:
                v=(2*L//d)*s-L*inc.get((a,b),0)
                if v: out[(b,a)]=v
        return out
    def back(inc):
        prev={}
        for b in nb:
            outs=[inc.get((b,a),0) for a in nb[b]]; s=sum(outs)
            if not any(outs): continue
            d=deg[b]
            for a,o in zip(nb[b],outs):
                num=(2*L//d)*s-L*o
                if num%(L*L): return None
                if num: prev[(a,b)]=num//(L*L)
        return prev
    inc={}
    for b,a in sources:
        if b in nb and a in nb[b]: inc[(a,b)]=inc.get((a,b),0)+1
    hist=[inc]
    for t in range(steps): hist.append(step(hist[-1]))
    ok_back=None
    if check_back:
        x=hist[-1]
        for t in range(steps):
            x=back(x)
            if x is None: break
        ok_back= x is not None and {k:v for k,v in x.items() if v}==hist[0]
    return hist,L,ok_back,nb
def cell_vals(h):
    s=defaultdict(int)
    for (a,b),v in h.items(): s[b]+=v
    return s
if __name__=="__main__":
    GAM_SYM=[3,7,11,11,7,29]   # 鏡で変わらない：γ1=γ4, γ2=γ3
    # 境目にちょうど乗る六面体は入れない（鏡で対称なまま、重なりを出さない）
    cells=build(GAM_SYM,strict=True); F=carrier(cells)
    idx={c:n for n,c in enumerate(cells)}
    mir={n:idx.get(mirror_cell(c)) for n,c in enumerate(cells)}
    mult=Counter(len(v) for v in F.values())
    r0= all(m is not None for m in mir.values()) and max(mult)<=2
    print("R0 六面体",len(cells),"個：鏡で自分に重なる",all(m is not None for m in mir.values()),"、面の持ち主",dict(mult),"→","OK" if r0 else "NG")
    C={n:ctr(c) for n,c in enumerate(cells)}
    xs=[p[0] for p in C.values()]; xmin=min(xs)
    # 風：左端の厚さ1の層に、左側の隣から入れる（鏡で対称）
    def wind_sources(cells,F,removed=set()):
        Fm=F; src=[]
        for f,v in Fm.items():
            if len(v)==2:
                a,b=v
                for p,q in ((a,b),(b,a)):
                    if C[q][0]<xmin+2.0 and C[p][0]<C[q][0]-1e-9: src.append((q,p))   # x が同じ面（横向きの面）は風の入口にしない
        return src
    src=wind_sources(cells,F)
    print("   風の入口：",len(src),"面、鏡で対称",set((mir[q],mir[p]) for q,p in src)==set(src))
    # 物体：中心 O（鏡の面の上）の球
    O=np.array([0.0,0.0,0.0])
    ball={n for n,p in C.items() if np.linalg.norm(p-O)<1.6}
    side={n for n,p in C.items() if np.linalg.norm(p-(O+np.array([0.3,1.5,0])))<0.9}-ball
    sidem={mir[n] for n in side}
    pos={n for n,p in C.items() if p[1]>1e-9}
    def A(h):
        s=cell_vals(h); return sum(s.get(n,0)**2-s.get(mir[n],0)**2 for n in pos)
    T_STEPS=24
    res={}
    for name,rem in [("対称な球",ball),("出っ張り",ball|side),("鏡写し",ball|sidem)]:
        hist,L,okb,nb=run(cells,F,rem,src,T_STEPS,check_back=True)
        S=[sum(h.values()) for h in hist]; Q=[sum(v*v for v in h.values()) for h in hist]
        r1=all(S[t+1]==L*S[t] for t in range(T_STEPS)) and all(Q[t+1]==L*L*Q[t] for t in range(T_STEPS))
        res[name]=[A(h) for h in hist]
        print(f"[{name}] 取り除いた六面体 {len(rem)}、L={L}：R1 {'OK' if r1 else 'NG'}、R2 逆に戻る {'OK' if okb else 'NG'}")
    r3=all(a==0 for a in res["対称な球"])
    print("R3 対称な球の A：すべての歩で 0 →",r3,"→","OK" if r3 else "NG")
    r4=any(a!=0 for a in res["出っ張り"]) and all(a==-b for a,b in zip(res["出っ張り"],res["鏡写し"]))
    print("R4 出っ張りの A が 0 でない",any(a!=0 for a in res["出っ張り"]),"、鏡写しでちょうど −A",all(a==-b for a,b in zip(res["出っ張り"],res["鏡写し"])),"→","OK" if r4 else "NG")
    sgn=[(1 if a>0 else -1 if a<0 else 0) for a in res["出っ張り"]]
    print("   出っ張りの A の符号（歩ごと）：","".join({1:"+",-1:"−",0:"0"}[x] for x in sgn))
    # R5：重ね合わせ
    half1=src[:len(src)//2]; half2=src[len(src)//2:]
    h1,_,_,_=run(cells,F,ball,half1,12); h2,_,_,_=run(cells,F,ball,half2,12); h12,_,_,_=run(cells,F,ball,half1+half2,12)
    sm=defaultdict(int)
    for k,v in h1[-1].items(): sm[k]+=v
    for k,v in h2[-1].items(): sm[k]+=v
    r5={k:v for k,v in sm.items() if v}==h12[-1]
    print("R5 二か所から入れた風 ＝ それぞれの和（12 歩目、すべての辺）：",r5,"→","OK" if r5 else "NG")
    # 負の対照：鏡で変わる γ
    cells2=build([3,7,11,17,23,29]); F2=carrier(cells2); idx2={c:n for n,c in enumerate(cells2)}
    sym2=all(mirror_cell(c) in idx2 for c in cells2)
    C2={n:ctr(c) for n,c in enumerate(cells2)}
    mir2={n:idx2.get(mirror_cell(c)) for n,c in enumerate(cells2)}
    xmin2=min(p[0] for p in C2.values())
    src2=[]
    for f,v in F2.items():
        if len(v)==2:
            a,b=v
            for p,q in ((a,b),(b,a)):
                if C2[q][0]<xmin2+2.0 and C2[p][0]<C2[q][0]-1e-9: src2.append((q,p))
    ball2={n for n,p in C2.items() if np.linalg.norm(p-O)<1.6}
    h,_,_,_=run(cells2,F2,ball2,src2,T_STEPS)
    pos2={n for n,p in C2.items() if p[1]>1e-9}
    def A2(hh):
        s=cell_vals(hh); return sum(s.get(n,0)**2-(s.get(mir2[n],0)**2 if mir2[n] is not None else 0) for n in pos2)
    a2=[A2(hh) for hh in h]
    print("負の対照（鏡で変わる γ）：担体が鏡で重なる",sym2,"、対称な球の A が 0 でない歩",sum(1 for a in a2 if a!=0),"/",len(a2),"→","NG" if any(a!=0 for a in a2) else "OK(!)")
    json.dump({"A":{k:[str(a) for a in v] for k,v in res.items()}},open("ch15_wind.json","w"))
