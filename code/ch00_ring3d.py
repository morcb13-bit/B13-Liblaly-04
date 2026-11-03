# 序章の図（fig00_ring3d.html）の座標の確認。図に入れた頂点・面・中心をそのまま読む。
# 実数の座標の確認で、整数の検定ではない（整数の照合は第7章）。
import json, sys
D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "ring3d_D.json"))
V, C = D["V"], D["C"]
def cut(k):
    c = C[k]; P = D["PB"] if c[2] > 0 else D["PT"]
    return [tuple(round(c[i] + V[j][i], 4) for i in range(3)) for j in P]
pts = [cut(k) for k in range(10)]
on = sum(all(abs(p[2]) < 1e-4 for p in pts[k]) for k in range(10))
print("C1 切り口の五角形が真ん中の平面に乗る：", on, "/ 10")
def same(p, q): return all(abs(p[i] - q[i]) < 1e-3 for i in range(3))
sh = [sum(any(same(p, q) for q in pts[(k + 1) % 10]) for p in pts[k]) for k in range(10)]
print("C2 隣どうしが共有する頂点の数：", sh, "→ 2 が", sh.count(2), "/ 10")
allp = []
for k in range(10):
    for p in pts[k]:
        if not any(same(p, q) for q in allp): allp.append(p)
print("C3 円環の頂点の数（重なりを除く）：", len(allp))
up = sum(1 for c in C if c[2] > 0)
print("上段の30面体：", up, "個（下面で切る）、下段：", 10 - up, "個（上面で切る）")
