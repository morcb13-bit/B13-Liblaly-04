# 第18章　銅 ── 第16章の正四面体を、6 個の整数のままシェルピンスキーのギャスケットにして並べる
# 判定は整数だけ。長さは √5 倍した二乗を (√5 の係数, 整数) の組で比べる。
import json, itertools
s=[[0]*6 for _ in range(6)]
for i in range(6):
    for j in range(6):
        if i!=j:
            a,b=min(i,j),max(i,j)
            s[i][j]=1 if a==0 else (1 if (b-a)%5 in (1,4) else -1)
def dot5(a,b): return (sum(a[i]*b[i] for i in range(6)), sum(a[i]*b[j]*s[i][j] for i in range(6) for j in range(6) if i!=j))
sub=lambda a,b:tuple(x-y for x,y in zip(a,b)); add=lambda a,b:tuple(x+y for x,y in zip(a,b))
mul=lambda a,k:tuple(x*k for x in a)
d2=lambda a,b:dot5(sub(a,b),sub(a,b))
R=json.load(open("ch16_atoms.json")); D20=[tuple(v) for v in R["D20"]]
T0=[D20[i] for i in R["T"][0]]                       # 第16章の正四面体（水素）
def gasket(tet,m):          # m 段のギャスケット。大きい正四面体を 2^(m−1) 倍して、辺の中点もすべて整数にする
    S=2**(m-1); T=[[mul(v,S) for v in tet]]
    for _ in range(m-1):
        T=[[t[i]]+[tuple((x+y)//2 for x,y in zip(t[i],t[j])) for j in range(4) if j!=i] for t in T for i in range(4)]
    return T
pts=lambda T:{p for t in T for p in t}
refl=lambda T,u,v:[[sub(add(u,v),p) for p in t] for t in T]   # 辺 uv の中点で折り返す：p → u + v − p
if __name__=="__main__":
    print("基準（走らせる前に決めたもの）")
    print(" G1 3 段のギャスケット：小さな正四面体 16 個、どれも 6 本の辺が等しく、大きな正四面体の辺の 1/4（長さの二乗で 16 倍）。点は 34 個")
    print(" G2 辺の中点で折り返して並べる：隣どうしは点を 5 個だけ分け合い、2 つ目からの原子は 1 個ごとに新しい点を 29 個持ち込む。隣でない原子は点を分け合わない")
    print(" G3 2 回折り返すと平行移動になり、原子の列はまっすぐ伸びる（1 つおきの原子が同じ整数の差でずれる）")
    print(" 負の対照：辺の中点ではなく頂点で折り返すと、分け合う点は 1 個、新しい点は 33 個になり、G2 は NG を返す")
    print()
    rows=[]
    for m in (1,2,3,4):
        T=gasket(T0,m); P=pts(T); big=d2(mul(T0[0],2**(m-1)),mul(T0[1],2**(m-1)))
        ok=all(len({d2(a,b) for a,b in itertools.combinations(t,2)})==1 for t in T)
        sm=d2(T[0][0],T[0][1]); A,B=T0[0],T0[1]; S=2**(m-1)
        N2=refl(T,mul(A,S),mul(B,S)); P2=pts(N2)
        rows.append((m,len(T),len(P),len(P&P2),len(P2-P)))
        print(f"[記録] {m} 段：小さな正四面体 {len(T):3d}、点 {len(P):3d}、正四面体がそろう {ok}、隣と分け合う点 {len(P&P2)}、新しい点 {len(P2-P)}")
    m=3; S=4; T=gasket(T0,m); P=pts(T)
    big=d2(mul(T0[0],S),mul(T0[1],S)); sm={d2(a,b) for t in T for a,b in itertools.combinations(t,2)}
    g1=len(T)==16 and len(P)==34 and len(sm)==1 and mul(list(sm)[0],16)==big
    print("G1 小さな正四面体",len(T),"個、点",len(P),"個、辺の長さの組",sm,"（大きな辺",big,"）→","OK" if g1 else "NG")
    # 銅の導線：辺 (0,1) の中点、次は向かいの辺 (2,3) の像の中点、と交互に折り返す
    K=6; atoms=[T]; V=[mul(v,S) for v in T0]; cur=V[:]
    for k in range(1,K):
        if k%2==1: u,v=cur[0],cur[1]
        else:      u,v=cur[2],cur[3]
        nxt=[sub(add(u,v),p) for p in cur]
        atoms.append(refl(atoms[-1],u,v)); cur=nxt
    PS=[pts(a) for a in atoms]
    sh=[len(PS[k]&PS[k+1]) for k in range(K-1)]; far=[len(PS[i]&PS[j]) for i in range(K) for j in range(i+2,K)]
    new=[len(PS[k]-set().union(*PS[:k])) for k in range(1,K)]
    g2=all(x==5 for x in sh) and all(x==0 for x in far) and all(x==29 for x in new)
    print("G2 隣と分け合う点",sh,"、隣でない組",set(far),"、新しく持ち込む点",new,"→","OK" if g2 else "NG")
    tot=len(set().union(*PS)); print("[記録] 原子",K,"個の点の合計",tot,"＝ 34 ＋ 29 ×",K-1)
    shift=[sub(min(PS[k+2]),min(PS[k])) for k in range(K-2)]
    tr=[all(add(p,shift[0]) in PS[k+2] for p in PS[k]) for k in range(K-2)]
    g3=all(tr) and len(set(shift))==1
    print("G3 1 つおきの原子のずれ",shift[0],"（すべての k で同じ、点の集まりが重なる）：",g3,"→","OK" if g3 else "NG")
    # 負の対照：頂点で折り返す
    A=mul(T0[0],S); N2=[[sub(mul(A,2),p) for p in t] for t in T]; P2=pts(N2)
    ng=len(P&P2)==5 and len(P2-P)==29
    print("負の対照：頂点で折り返す → 分け合う点",len(P&P2),"、新しい点",len(P2-P),"→","OK(!)" if ng else "NG")
    json.dump({"atoms":[[list(map(list,t)) for t in a] for a in atoms],"rows":rows,"shift":shift[0]},open("ch18_gasket.json","w"))
