# 第7章の静止図：30面体1個を三つの軸に垂直な平面で切った切り口を並べる（表示用、実数）
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib import font_manager as fm
t=np.arctan(2)
AX=np.array([[0,0,1.]]+[[np.sin(t)*np.cos(2*np.pi*k/5),np.sin(t)*np.sin(2*np.pi*k/5),np.cos(t)] for k in range(1,6)])
R=json.load(open("ch01_rt.json")); P=np.array(R["V"])@AX/2; E=[e[:2] for e in R["E"]]
v3=np.array(next(v for v,k in zip(R["V"],R["K"]) if k==3))@AX/2
f=R["F"][0]; n2=np.array([sum(R["V"][i][k] for i in f) for k in range(6)])@AX/2
axes=[("五回軸",np.array([0,0,1.])),("三回軸",v3/np.linalg.norm(v3)),("二回軸",n2/np.linalg.norm(n2))]
def section(n,h):
    pts=[]
    for a,b in E:
        za,zb=P[a]@n-h,P[b]@n-h
        if za*zb<0: l=za/(za-zb); pts.append(P[a]+l*(P[b]-P[a]))
    pts=np.array(pts)
    u=np.cross(n,[0,0,1] if abs(n[2])<0.9 else [1,0,0]); u/=np.linalg.norm(u); w=np.cross(n,u)
    q=np.c_[pts@u,pts@w]; c=q.mean(0); q=q[np.argsort(np.arctan2(q[:,1]-c[1],q[:,0]-c[0]))]
    return q
jp=fm.FontProperties(family="Noto Sans CJK JP")
ink="#1e2724"; fill="#d9826a"; bg="#f8f7f2"
fig=plt.figure(figsize=(13,5.6),dpi=150); fig.patch.set_facecolor(bg)
ncol=9; cw=1/(ncol+1.2)
for r,(name,n) in enumerate(axes):
    lv=np.unique(np.round(P@n,9))
    mids=(lv[:-1]+lv[1:])/2
    y0=1-(r+1)*0.31
    fig.text(0.012,y0+0.13,name,fontproperties=jp,fontsize=15,color=ink,weight="bold")
    for j,h in enumerate(mids):
        q=section(n,h)
        ax=fig.add_axes([0.11+j*cw*0.98,y0,cw*0.9,0.26]); ax.set_xlim(-1.7,1.7); ax.set_ylim(-2.15,1.6); ax.set_aspect("equal"); ax.axis("off")
        ax.add_patch(Polygon(q,closed=True,fc=fill,ec=ink,lw=1.4))
        ax.text(0,-2.12,f"{len(q)} 角形",ha="center",va="bottom",fontproperties=jp,fontsize=11,color=ink)
fig.text(0.5,0.015,"黄金のひし形30面体を、それぞれの軸に垂直な平面で下から上へ切った切り口（頂点の高さで区切った段のあいだごとに一枚）",
         ha="center",fontproperties=jp,fontsize=10.5,color="#5d6660")
fig.savefig("fig07b_sections.png",facecolor=bg)
print("ok")
