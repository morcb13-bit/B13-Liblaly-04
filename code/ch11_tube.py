# 第11章　立体からチューブへ。30面体の輪を整数の一歩で重ねていく
import json
from collections import Counter
exec(open("ch11_tube_search.py").read().split("cands=set()")[0])   # s, dot5, add, sub, cmp_pair, SIG, C, check, TOUCH
def zt(c): return (-c[3], c[0]+c[3], c[1]-c[3], c[2]+c[3])
def zpow(n):
    c=(1,0,0,0)
    for _ in range(n%10): c=zt(c)
    return c
sc=lambda k,a:tuple(k*x for x in a)
def shadow(a):
    z=(0,0,0,0)
    for k in range(1,6): z=add(z,sc(a[k],zpow(2*k)))
    return z
def shift(a,j):
    b=[0]*6; b[0]=a[0]
    for k in range(1,6): b[1+(k-1+j)%5]=a[k]
    return tuple(b)
v0=(4,2,4,4,2,0)
B=[add(c,v0) for c in C]
bad,sh=check(C,B)
z1= bad==0 and sh[4]==6 and (v0[0],sum(v0[1:]))==(4,12)
print("Z1 一歩 v =",v0,"高さの組",(v0[0],sum(v0[1:])),"：面で接する組",sh[4],"、食い込み",bad,"→","OK" if z1 else "NG")
def stack(seq,N):
    rings=[C]; off=(0,)*6; offs=[off]
    for n in range(1,N):
        off=add(off,seq(n-1)); offs.append(off); rings.append([add(c,off) for c in C])
    return rings,offs
def tube_check(rings):
    bad=0; cont=Counter()
    for i in range(len(rings)):
        for j in range(i+1,len(rings)):
            b,sh=check(rings[i],rings[j]); bad+=b; cont[j-i]+=sh[4]
    return bad,cont
N=8
rings,offs=stack(lambda n:shift(v0,n),N)
bad,cont=tube_check(rings)
z2= bad==0 and cont[1]==6*(N-1) and all(cont[d]==0 for d in range(2,N))
print("Z2 回しながら",N,"段：食い込み",bad,"、面で接する組（段の差ごと）",dict(cont),"→","OK" if z2 else "NG")
s5=(0,)*6
for k in range(5): s5=add(s5,shift(v0,k))
z3= shadow(s5)==(0,0,0,0)
print("Z3 5歩の和",s5,"の影",shadow(s5),"→","OK（5段で真上に戻る）" if z3 else "NG")
print("   5段ごとの高さの組：",(s5[0],sum(s5[1:])))
# Z4：各段の真ん中の平面（高さの組＝その段のずらし分）の切り口
ok4=0
for r,off in zip(rings,offs):
    H=(off[0],sum(off[1:])); cut=set()
    for c in r:
        for sg in SIG:
            a=add(c,sg)
            if (a[0],sum(a[1:]))==H: cut.add(a)
    sh_cut={sub(shadow(a),shadow(off)) for a in cut}
    base={shadow(add(c,sg)) for c in C for sg in SIG if (add(c,sg)[0],sum(add(c,sg)[1:]))==(0,0)}
    ok4+= len(cut)==30 and sh_cut==base
print("Z4 各段の真ん中の切り口が、最初の円環をずらしたもの（頂点30個）：",ok4,"/",N,"→","OK" if ok4==N else "NG")
# 記録：回さずに重ねる
r2,_=stack(lambda n:v0,N); bad2,cont2=tube_check(r2)
print("[記録] 回さずに",N,"段：食い込み",bad2,"、面で接する組",dict(cont2))
print("[記録] 一歩の影",shadow(v0),"は、輪の30面体の中心の影と同じ長さ：",shadow(v0) in {shadow(c) for c in C} or shadow(sc(-1,v0)) in {shadow(c) for c in C})
# 負の対照
vh=(2,1,2,2,1,0); bad3,sh3=check(C,[add(c,vh) for c in C])
print("負の対照：一歩を半分",vh,"→ 食い込み",bad3,"→ Z1 は","NG" if bad3>0 else "OK(!)")
json.dump({"offs":offs,"v0":v0},open("ch11_tube.json","w"))
