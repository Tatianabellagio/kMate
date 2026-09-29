#!/usr/bin/env python
"""Build+execute the purging-vs-hitchhiking notebook. Load-only. `basic` env."""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "notebooks", "sv_purging_vs_hitchhiking.ipynb")

md_title = """# Are insertions actually being purged in GrENE-Net, or just hitchhiking?

The temporal result is that insertions are more purged at hot/arid gardens. Two explanations have
been live all along:

* **hitchhiking** — cold-origin ecotypes carry more insertions, so when warm-origin ecotypes win at
  hot gardens every cold-carried variant falls, insertion or not;
* **purifying selection** — the insertions are deleterious and the large GrENE-Net founder
  populations (~5k seeds) let selection act on them.

**Miles's diagnostic separates them.** A deleterious insertion should decline in *cold* gardens too,
where its carrier lineages are winning. A hitchhiking one should fall at hot gardens and rise at cold
ones. So the test is not the climate gradient — it is whether the decline is *global*.

Everything here uses the residuals from `_sv_temporal_climate_matched.py`: each variant's `s` minus
the median `s` of SNPs matched on **both** starting frequency and **carrier-origin climate**. So the
comparison SNPs for a cold-carried insertion are themselves cold-carried, and the ecotype-composition
effect is already differenced out before anything below runs.

**Groups**, from the founder panel (231 founders, origin-bio1 tertiles of 77):
`cold-specific` ≥ 2 cold carriers and zero warm; `warm-specific` the mirror; `shared` carriers in
both.

**The controls that make this readable.** Short indels (1–50 bp) are carried through every split as
the reference class — they are subject to the same lineage structure, the same matching and the same
estimator, so any test's rate in restricted *indels* is its false-positive rate here. Deletions
> 50 bp are shown too, though they are barely lineage-restricted (81 of 2,396) and underpowered."""

code_setup = r"""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'figure.dpi':110,'font.size':8,'axes.linewidth':0.6})

G=os.path.abspath(os.path.join(os.getcwd(),"..","r1_sv_negative_selection","results","sv_adaptive"))
os.makedirs(f"{G}/plots",exist_ok=True)
def style(ax,axis="both"):
    ax.set_axisbelow(True); ax.grid(color="0.88",lw=0.7,axis=axis)
    for sp in ax.spines.values(): sp.set_visible(False)

V=pd.read_csv(f"{G}/sv_purging_classify_variants.csv")
D=pd.read_csv(f"{G}/sv_temporal_climate_matched_sites.csv").sort_values("bio1").reset_index(drop=True)
P=np.load(f"{G}/sv_temporal_climate_matched_persite.npz")
ORD=["cold-specific","shared","warm-specific"]
CC={"cold-specific":"#2166AC","shared":"0.55","warm-specific":"#B2182B"}
CLS=[("insertion","ins"),("deletion","del"),("indel","indel")]
print(f"{len(D)} gardens, bio1 {D.bio1.min():.0f}-{D.bio1.max():.0f} C")
print(V.groupby(["cls","origin"]).size().unstack(fill_value=0).loc[["insertion","deletion","indel"],ORD])
"""

md_diag = """## 1. Miles's diagnostic — do they decline in cold gardens too?

One dot per garden: the median residual of that group's variants in that garden. Filled = below
zero (declining), open = above. Dashed line = the group's median across gardens.

**How to read it:** a hitchhiking group would have filled dots only on the right (hot gardens) and
open dots on the left. A deleterious group has filled dots all the way across. Top row is
insertions, bottom row is short indels split the same way — the null. Deletions are in the table
but not the figure: only 25 cold-specific and 56 warm-specific exist, too few for per-garden medians."""

code_diag = r"""

def garden_medians(clab,key,o):
    allc=V[V.cls==clab].col.to_numpy(); cidx=V[(V.cls==clab)&(V.origin==o)].col.to_numpy()
    if cidx.size<20: return None,cidx.size
    pos=np.searchsorted(allc,cidx); med=[]
    for s in D.site.astype(int):
        sub=P[f"{s}_res_{key}_matched"][pos]; sub=sub[~np.isnan(sub)]
        med.append(np.median(sub) if sub.size>=20 else np.nan)
    return np.array(med),cidx.size

SUM=[]; M={}
for clab,key in CLS:
    for o in ORD:
        med,nv=garden_medians(clab,key,o)
        if med is None: continue
        g=~np.isnan(med); M[(clab,o)]=(med,nv)
        SUM.append(dict(cls=clab,origin=o,n=nv,median=np.median(med[g]),n_neg=int((med[g]<0).sum()),
                        n_g=int(g.sum()),p=stats.wilcoxon(med[g]).pvalue,
                        rho=stats.spearmanr(D.bio1[g],med[g]).statistic))

# figure: insertions (top) vs indels = the null (bottom), one panel per origin group, dots only
ROWS=[("insertion","insertions"),("indel","indels (null)")]
allv=np.concatenate([M[(c,o)][0] for c,_ in ROWS for o in ORD if (c,o) in M])
lim=np.nanmax(np.abs(allv))*1.1
fig,ax=plt.subplots(2,3,figsize=(12,6.2),sharex=True,sharey=True)
for i,(clab,rlab) in enumerate(ROWS):
    for j,o in enumerate(ORD):
        a=ax[i,j]; style(a)
        med,nv=M[(clab,o)]; g=~np.isnan(med)
        a.axhline(0,color="0.3",lw=1.0)
        neg=med[g]<0
        a.scatter(D.bio1[g][neg],med[g][neg],s=34,color=CC[o],edgecolor="0.2",linewidth=0.4,zorder=3)
        a.scatter(D.bio1[g][~neg],med[g][~neg],s=34,facecolor="white",edgecolor=CC[o],linewidth=1.0,zorder=3)
        a.axhline(np.median(med[g]),color=CC[o],lw=1.6,ls="--",zorder=2)
        a.annotate(f"{rlab}, {o}  (n={nv:,})\nbelow 0 in {int(neg.sum())}/{g.sum()} gardens\n"
                   f"median {np.median(med[g]):+.3f}",
                   xy=(0.03,0.04),xycoords="axes fraction",va="bottom",fontsize=7.8,color="0.25",
                   bbox=dict(facecolor="white",alpha=0.8,edgecolor="none",pad=1.5))
        a.set_ylim(-lim,lim)
        if i==1: a.set_xlabel("garden temperature  bio1 (C)")
        if j==0: a.set_ylabel("median residual vs matched SNPs\n(< 0 = declines)")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_by_garden.png",dpi=130,bbox_inches="tight"); plt.show()

S=pd.DataFrame(SUM)
print(f"{'class':<11}{'origin':<15}{'n vars':>7}{'median':>10}{'<0 at':>8}{'wilcoxon p':>12}{'rho bio1':>10}")
for _,r in S.iterrows():
    print(f"{r.cls:<11}{r.origin:<15}{int(r.n):>7}{r['median']:>+10.4f}{f'{int(r.n_neg)}/{int(r.n_g)}':>8}"
          f"{r.p:>12.4f}{r.rho:>+10.3f}")
"""

md_call = """## 2. Per variant — global decline, climate-graded decline, or neither?

Each variant is tested twice on its own 31 residuals:

* **global** — sign test: negative in more gardens than chance (deleterious signature);
* **climate-graded** — Spearman against garden bio1, negative slope (hitchhiking signature).

They are not exclusive. The nominal 5% threshold is strongly anti-conservative here, because one
variant's residuals are correlated across gardens — which is why every bar is shown against the
**same test applied to restricted indels**. The indel height is the false-positive rate; only the
excess above it is signal."""

code_call = r"""
def rates(cls,o):
    d=V[(V.cls==cls)&(V.origin==o)]
    if len(d)<20: return None
    return dict(n=len(d),
                glob=(d.is_global).mean(), grad=(d.is_graded).mean(),
                both=(d.is_global&d.is_graded).mean(), none=(~d.is_global&~d.is_graded).mean())

fig,ax=plt.subplots(1,2,figsize=(11.5,4.2),sharey=True)
W=0.26
for a,(metric,mlab) in zip(ax,[("glob","GLOBAL decline (sign test)  -- deleterious signature"),
                               ("grad","CLIMATE-GRADED decline  -- hitchhiking signature")]):
    style(a,axis="y")
    x=np.arange(len(ORD))
    for k,(clab,_) in enumerate(CLS):
        vals=[(rates(clab,o) or {}).get(metric,np.nan) for o in ORD]
        col={"insertion":"#D55E00","deletion":"#009E73","indel":"#888888"}[clab]
        a.bar(x+(k-1)*W,[100*v for v in vals],width=W*0.9,color=col,linewidth=0,
              label=clab+"s"+(" (null)" if clab=="indel" else ""))
        for xi,v in zip(x,vals):
            if not np.isnan(v):
                a.annotate(f"{100*v:.0f}",(xi+(k-1)*W,100*v),xytext=(0,3),textcoords="offset points",
                           ha="center",fontsize=7,color="0.3")
    a.set_xticks(x); a.set_xticklabels(ORD,fontsize=8.5)
    a.set_ylabel("% of variants" if metric=="glob" else "")
    a.annotate(mlab,xy=(0.02,0.97),xycoords="axes fraction",va="top",fontsize=8.5,color="0.3")
ax[0].legend(frameon=False,fontsize=8,loc="upper right")
ax[0].set_ylim(0,60)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_classification.png",dpi=130,bbox_inches="tight"); plt.show()

print(f"{'class':<11}{'origin':<15}{'n':>6}{'global':>9}{'graded':>9}{'both':>8}{'neither':>9}")
for clab,_ in CLS:
    for o in ORD:
        r=rates(clab,o)
        if r is None: continue
        print(f"{clab:<11}{o:<15}{r['n']:>6}{100*r['glob']:>8.0f}%{100*r['grad']:>8.0f}%"
              f"{100*r['both']:>7.0f}%{100*r['none']:>8.0f}%")
print("\nEXCESS over the restricted-indel null (same origin group), in percentage points and variants:")
print(f"{'origin':<15}{'global excess':>22}{'graded excess':>22}")
for o in ORD:
    ri=rates("indel",o); rs=rates("insertion",o)
    if ri is None or rs is None: continue
    dg=100*(rs['glob']-ri['glob']); dr=100*(rs['grad']-ri['grad'])
    print(f"{o:<15}{f'{dg:+.0f} pts  ({dg/100*rs[chr(110)]:+.0f} vars)':>22}"
          f"{f'{dr:+.0f} pts  ({dr/100*rs[chr(110)]:+.0f} vars)':>22}")
"""

md_dist = """### 2b. Insertions: the global distribution, and split by carrier origin

The first box is **every tested insertion**; the next three split them by where their carriers come
from. Dots are coloured by that split in every box, so the global box shows the mix. The first panel
answers "globally, is there a tendency?"; the split answers "which insertions carry it?". Each box is
tested against 0.5 with a one-sample Wilcoxon (the p under each box label) — slightly conservative for declines, since matched
indels sit at a median of 0.484 rather than 0.5.

Below it, the original version of this figure against indels as the null, for insertions and deletions."""

md_dist_old = """### (supplementary) The distributions behind those bars

Each bar above is the share of variants past a threshold. These are the per-variant statistics the
thresholds are applied to, one dot per variant with a transparent box on top:

* **left** — fraction of the 31 gardens in which the variant's residual is below zero. Chance is 0.5;
  a deleterious variant sits well above it.
* **right** — the absolute linear slope of the variant's residual against garden bio1 (change in
  residual per °C). Absolute, because hitchhiking points opposite ways for the two groups: a
  cold-specific hitchhiker falls at hot gardens (negative slope), a warm-specific one rises there
  (positive slope). |slope| measures how strongly a variant tracks garden temperature regardless of
  direction. It is never zero even for noise, so compare against the indel box, not against 0.

Orange = insertions, grey = indels (the null). Deletions are left out here (too few restricted
ones). For readability at most 400 dots are drawn per box; the boxes use every variant. Medians and
counts are in the table below the figure.

**What to look for:** in the left panel, the cold-specific insertion box sitting clearly above 0.5
while its indel neighbour doesn't. In the right panel, insertion and indel boxes sitting at the same
height — no extra climate-grading for insertions."""

code_insonly = r"""

from matplotlib.gridspec import GridSpec
V["abs_slope"]=V.slope_bio1.abs()
F=pd.read_csv(f"{G}/sv_origin_partition_founders.csv",dtype={"founder":str})
e1=F[F.group=="cold"].bio1.max(); e2=F[F.group=="warm"].bio1.min()
GC={"cold":"#2166AC","temperate":"0.6","warm":"#B2182B"}
DOT={"cold-specific":"#2166AC","shared":"0.55","warm-specific":"#B2182B","other":"0.8"}
I=V[V.cls=="insertion"].copy()
BOXES=[("all insertions",I),("cold-specific",I[I.origin=="cold-specific"]),
       ("shared",I[I.origin=="shared"]),("warm-specific",I[I.origin=="warm-specific"])]

fig=plt.figure(figsize=(12.5,7.4))
gs=GridSpec(2,2,height_ratios=[1,2.6],hspace=0.34,wspace=0.22,figure=fig)
h=fig.add_subplot(gs[0,:]); style(h,axis="y")
bins=np.arange(np.floor(F.bio1.min()),np.ceil(F.bio1.max())+1,1.0)
h.hist([F[F.group==g].bio1 for g in ("cold","temperate","warm")],bins=bins,stacked=True,
       color=[GC[g] for g in ("cold","temperate","warm")],edgecolor="white",linewidth=0.5)
for e in (e1,e2): h.axvline(e,color="0.2",lw=1.1,ls="--")
h.set_ylim(0,h.get_ylim()[1]*1.75); yt=h.get_ylim()[1]
for g,x in (("cold",(F.bio1.min()+e1)/2),("temperate",(e1+e2)/2),("warm",(e2+F.bio1.max())/2)):
    lo,hi=F[F.group==g].bio1.min(),F[F.group==g].bio1.max()
    h.annotate(f"{g}-origin\n{int((F.group==g).sum())} ecotypes, {lo:.1f} to {hi:.1f} C",(x,yt*0.97),
               ha="center",va="top",fontsize=8,color=GC[g] if g!="temperate" else "0.3",fontweight="bold")
h.set_xlabel("ecotype origin temperature  bio1 (C)"); h.set_ylabel("founder\necotypes")
h.annotate("cold-specific = carried by >= 2 cold-origin ecotypes and no warm-origin ones;  "
           "warm-specific = the mirror;  shared = carried in both.  Dots coloured by these groups.",
           xy=(0.0,-0.42),xycoords="axes fraction",fontsize=7.5,color="0.35")

rng=np.random.default_rng(0)
for pi,(col,ylab,ref) in enumerate([("frac_neg","fraction of gardens where residual < 0\n(0.5 = no tendency, higher = declines more consistently)",0.5),
                                    ("abs_slope","|linear slope of residual vs garden bio1|  (per C)\n(bigger = tracks garden temperature more strongly)",None)]):
    a=fig.add_subplot(gs[1,pi]); style(a,axis="y")
    if ref is not None: a.axhline(ref,color="0.35",lw=1.0,ls="--",zorder=1)
    TL=[]
    for i,(lab,d) in enumerate(BOXES):
        v=d[col].dropna().to_numpy(); orig=d.loc[d[col].notna(),"origin"].to_numpy()
        take=np.arange(v.size) if v.size<=400 else rng.choice(v.size,400,replace=False)
        a.scatter(i+rng.uniform(-0.2,0.2,take.size),v[take],s=7,
                  c=[DOT.get(o,"0.8") for o in orig[take]],alpha=0.6,linewidth=0,zorder=2)
        a.boxplot(v,positions=[i],widths=0.55,patch_artist=True,showfliers=False,
                  boxprops=dict(facecolor="none",edgecolor="0.1",linewidth=1.1),
                  medianprops=dict(color="k",lw=1.8),whiskerprops=dict(color="0.2",lw=0.9),
                  capprops=dict(color="0.2",lw=0.9),manage_ticks=False,zorder=3)
        if ref is not None:
            p1=stats.wilcoxon(v-ref).pvalue
            TL.append(f"{lab}\nmedian {np.median(v):.3f}\np = {p1:.1e}" if p1<1e-3 else f"{lab}\nmedian {np.median(v):.3f}\np = {p1:.3f}")
        else:
            TL.append(f"{lab}\nmedian {np.median(v):.3f}")
            print(f"{lab:<15} n={v.size:>5}  median frac<0 = {np.median(v):.3f}  "
                  f"share declining in >50% of gardens = {100*(v>0.5).mean():.1f}%  one-sample Wilcoxon vs 0.5 p = {p1:.3g}")
    a.axvline(0.5,color="0.75",lw=0.8)          # separates the global box from the split
    ca=I[I.origin=="cold-specific"][col].dropna(); wa=I[I.origin=="warm-specific"][col].dropna()
    pv=stats.ranksums(ca,wa).pvalue
    lo_,hi_=a.get_ylim(); yb=hi_+0.18*(hi_-lo_)
    a.plot([1,1,3,3],[yb-0.02*(hi_-lo_),yb,yb,yb-0.02*(hi_-lo_)],color="0.15",lw=1.0)
    a.annotate(f"cold- vs warm-specific: Wilcoxon rank-sum p = {pv:.3g}",(2,yb),xytext=(0,3),
               textcoords="offset points",ha="center",va="bottom",fontsize=8,color="0.15")
    a.set_ylim(lo_,yb+0.12*(hi_-lo_))
    a.set_xticks(range(len(BOXES))); a.set_xticklabels(TL,fontsize=7.6)
    a.set_xlim(-0.6,len(BOXES)-0.4); a.set_ylabel(ylab)
    print(f"{col}: cold- vs warm-specific rank-sum p = {pv:.3g}\n")
fig.legend(handles=[plt.Line2D([0],[0],marker="o",ls="",color=DOT[o],ms=6,label=o)
                    for o in ("cold-specific","shared","warm-specific")],
           frameon=False,fontsize=8.5,loc="lower center",ncol=3,bbox_to_anchor=(0.5,-0.07))
fig.savefig(f"{G}/plots/sv_purging_insertions.png",dpi=130,bbox_inches="tight"); plt.show()
"""

code_dist = r"""

from matplotlib.gridspec import GridSpec
V["abs_slope"]=V.slope_bio1.abs()     # direction-free: cold-specific hitchhikers fall at hot gardens, warm-specific ones rise
F=pd.read_csv(f"{G}/sv_origin_partition_founders.csv",dtype={"founder":str})
e1=F[F.group=="cold"].bio1.max(); e2=F[F.group=="warm"].bio1.min()
GC={"cold":"#2166AC","temperate":"0.6","warm":"#B2182B"}
CCLS={"insertion":"#D55E00","deletion":"#009E73","indel":"#888888"}
MAXDOTS=400; W=0.34

def draw_dist(svcls,fname):
    DCLS=[(svcls,None),("indel","indel")]
    fig=plt.figure(figsize=(12,7.2))
    gs=GridSpec(2,2,height_ratios=[1,2.6],hspace=0.32,wspace=0.22,figure=fig)

    # top: where the 231 founder ecotypes come from, and where the group thresholds fall
    h=fig.add_subplot(gs[0,:]); style(h,axis="y")
    bins=np.arange(np.floor(F.bio1.min()),np.ceil(F.bio1.max())+1,1.0)
    h.hist([F[F.group==g].bio1 for g in ("cold","temperate","warm")],bins=bins,stacked=True,
           color=[GC[g] for g in ("cold","temperate","warm")],edgecolor="white",linewidth=0.5)
    for e in (e1,e2): h.axvline(e,color="0.2",lw=1.1,ls="--")
    h.set_ylim(0,h.get_ylim()[1]*1.75)                  # headroom so the group labels clear the bars
    yt=h.get_ylim()[1]
    for g,x in (("cold",(F.bio1.min()+e1)/2),("temperate",(e1+e2)/2),("warm",(e2+F.bio1.max())/2)):
        lo,hi=F[F.group==g].bio1.min(),F[F.group==g].bio1.max()
        h.annotate(f"{g}-origin\n{int((F.group==g).sum())} ecotypes, {lo:.1f} to {hi:.1f} C",(x,yt*0.97),
                   ha="center",va="top",fontsize=8,color=GC[g] if g!="temperate" else "0.3",fontweight="bold")
    h.set_xlabel("ecotype origin temperature  bio1 (C)"); h.set_ylabel("founder\necotypes")
    h.annotate("cold-specific = carried by >= 2 cold-origin ecotypes and no warm-origin ones;  "
               "warm-specific = the mirror;  shared = carried in both",
               xy=(0.0,-0.42),xycoords="axes fraction",fontsize=7.5,color="0.35")

    rng=np.random.default_rng(0)
    for pi,(col,ylab,ref) in enumerate([("frac_neg","fraction of gardens where residual < 0\n(0.5 = chance, higher = declines more consistently)",0.5),
                                        ("abs_slope","|linear slope of residual vs garden bio1|  (per C)\n(bigger = tracks garden temperature more strongly)",None)]):
        a=fig.add_subplot(gs[1,pi]); style(a,axis="y")
        if ref is not None: a.axhline(ref,color="0.35",lw=1.0,ls="--",zorder=1)
        for gi,o in enumerate(ORD):
            for k,(clab,_) in enumerate(DCLS):
                v=V[(V.cls==clab)&(V.origin==o)][col].dropna().to_numpy()
                x0=gi+(k-0.5)*W
                dots=v if v.size<=MAXDOTS else rng.choice(v,MAXDOTS,replace=False)
                a.scatter(x0+rng.uniform(-W*0.36,W*0.36,dots.size),dots,s=6,color=CCLS[clab],
                          alpha=0.45,linewidth=0,zorder=2)
                a.boxplot(v,positions=[x0],widths=W*0.8,patch_artist=True,showfliers=False,
                          boxprops=dict(facecolor="none",edgecolor="0.1",linewidth=1.1),
                          medianprops=dict(color="k",lw=1.8),whiskerprops=dict(color="0.2",lw=0.9),
                          capprops=dict(color="0.2",lw=0.9),manage_ticks=False,zorder=3)
        # Wilcoxon rank-sum (= Mann-Whitney): cold-specific vs warm-specific INSERTIONS
        ca=V[(V.cls==svcls)&(V.origin=="cold-specific")][col].dropna()
        wa=V[(V.cls==svcls)&(V.origin=="warm-specific")][col].dropna()
        pv=stats.ranksums(ca,wa).pvalue
        xa,xb=0-0.5*W,2-0.5*W
        lo_,hi_=a.get_ylim(); yb=hi_+0.04*(hi_-lo_)
        a.plot([xa,xa,xb,xb],[yb-0.02*(hi_-lo_),yb,yb,yb-0.02*(hi_-lo_)],color="0.15",lw=1.0)
        a.annotate(f"cold- vs warm-specific {svcls}s: Wilcoxon rank-sum p = {pv:.3g}",((xa+xb)/2,yb),
                   xytext=(0,3),textcoords="offset points",ha="center",va="bottom",fontsize=8,color="0.15")
        a.set_ylim(lo_,yb+0.12*(hi_-lo_))
        a.set_xticks(range(len(ORD))); a.set_xticklabels(ORD,fontsize=8.5); a.set_xlim(-0.6,len(ORD)-0.4)
        a.set_ylabel(ylab)
        print(f"{col:<9} cold-specific vs warm-specific {svcls}s: Wilcoxon rank-sum p = {pv:.4g} "
              f"(medians {ca.median():.3f} vs {wa.median():.3f}, n = {len(ca)} / {len(wa)})")
    fig.legend(handles=[plt.Line2D([0],[0],marker="o",ls="",color=CCLS[c],ms=6,
                                  label=c+"s"+(" (null)" if c=="indel" else "")) for c,_ in DCLS],
               frameon=False,fontsize=8.5,loc="lower center",ncol=2,bbox_to_anchor=(0.5,-0.02))
    fig.savefig(f"{G}/plots/{fname}",dpi=130,bbox_inches="tight"); plt.show()

    print(f"\n{'class':<11}{'origin':<15}{'n':>7}{'median frac<0':>15}{'median |slope|':>16}{'median slope':>14}")
    for clab,_ in DCLS:
        for o in ORD:
            d=V[(V.cls==clab)&(V.origin==o)]
            print(f"{clab:<11}{o:<15}{len(d):>7,}{d.frac_neg.median():>15.3f}{d.abs_slope.median():>16.4f}"
                  f"{d.slope_bio1.median():>+14.4f}")

draw_dist("insertion","sv_purging_distributions.png")
draw_dist("deletion","sv_purging_distributions_del.png")
"""

md_glob = """### 2c. Globally — all variants, no split by origin

The same residuals (each variant vs SNPs matched on starting frequency **and** carrier-origin
climate), pooled over every tested variant rather than split by where its carriers come from. This
is the direct question: **overall, are SVs being purged more than comparable variants?**

Left: fraction of the 31 gardens in which the variant's residual is below zero (0.5 = no net
tendency). Right: the variant's mean residual across gardens (< 0 = it fell more than matched SNPs).
Each SV class is tested against indels with a Wilcoxon rank-sum. For readability at most 400 dots
are drawn per box; the boxes and tests use every variant."""

code_glob = r"""

V["cls2"]=V.cls
GCL=[("all SVs",V.cls.isin(["insertion","deletion"])),
     ("insertions",V.cls=="insertion"),
     ("deletions",V.cls=="deletion"),
     ("indels (null)",V.cls=="indel")]
GCOL={"all SVs":"#8C510A","insertions":"#D55E00","deletions":"#009E73","indels (null)":"#888888"}
null={c:V[V.cls=="indel"][c].dropna().to_numpy() for c in ("frac_neg","mean_res")}
fig,ax=plt.subplots(1,2,figsize=(12,4.8))
rng=np.random.default_rng(0)
rows=[]
for a,(col,ylab,ref) in zip(ax,[("frac_neg","fraction of gardens where residual < 0\n(0.5 = no net tendency)",0.5),
                                ("mean_res","per-variant mean residual across gardens\n(< 0 = falls more than matched SNPs)",0.0)]):
    style(a,axis="y"); a.axhline(ref,color="0.35",lw=1.0,ls="--",zorder=1)
    for i,(lab,m) in enumerate(GCL):
        v=V[m][col].dropna().to_numpy()
        dots=v if v.size<=400 else rng.choice(v,400,replace=False)
        a.scatter(i+rng.uniform(-0.18,0.18,dots.size),dots,s=6,color=GCOL[lab],alpha=0.45,linewidth=0,zorder=2)
        a.boxplot(v,positions=[i],widths=0.5,patch_artist=True,showfliers=False,
                  boxprops=dict(facecolor="none",edgecolor="0.1",linewidth=1.1),
                  medianprops=dict(color="k",lw=1.8),whiskerprops=dict(color="0.2",lw=0.9),
                  capprops=dict(color="0.2",lw=0.9),manage_ticks=False,zorder=3)
        if lab!="indels (null)":
            pv=stats.ranksums(v,null[col]).pvalue
            a.annotate(f"vs indels\np = {pv:.2g}",(i,1.0),xycoords=("data","axes fraction"),
                       ha="center",va="bottom",fontsize=7.5,color="0.25")
            rows.append(dict(stat=col,cls=lab,n=v.size,median=np.median(v),null_median=np.median(null[col]),p=pv))
    a.set_xticks(range(len(GCL))); a.set_xticklabels([l for l,_ in GCL],fontsize=8.5)
    a.set_xlim(-0.6,len(GCL)-0.4); a.set_ylabel(ylab)
    if col=="mean_res":
        lo,hi=np.percentile(V.mean_res.dropna(),[1,99]); a.set_ylim(lo,hi)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_global.png",dpi=130,bbox_inches="tight"); plt.show()

print(f"{'statistic':<10}{'class':<14}{'n':>8}{'median':>10}{'indel median':>14}{'rank-sum p':>12}")
for r in rows:
    print(f"{r['stat']:<10}{r['cls']:<14}{r['n']:>8,}{r['median']:>+10.4f}{r['null_median']:>+14.4f}{r['p']:>12.3g}")
for lab,m in GCL:
    d=V[m]
    print(f"{lab:<14} share of variants declining in >50% of gardens: {100*(d.frac_neg>0.5).mean():.1f}%")
"""

md_cw = """## 3. Cold-specific vs warm-specific insertions

Per variant, the mean residual across the 31 gardens. Bootstrap 95% CIs on the medians."""

code_cw = r"""
d=V[(V.cls=="insertion")]
a=d[d.origin=="cold-specific"].mean_res.dropna().to_numpy()
b=d[d.origin=="warm-specific"].mean_res.dropna().to_numpy()
c=d[d.origin=="shared"].mean_res.dropna().to_numpy()
def boot(x,B=5000,seed=0):
    r=np.random.default_rng(seed); return np.percentile([np.median(r.choice(x,x.size,True)) for _ in range(B)],[2.5,97.5])
fig,ax=plt.subplots(figsize=(7.5,4.2)); style(ax)
ax.axvline(0,color="0.4",lw=1.0)
for i,(lab,x) in enumerate((("cold-specific",a),("shared",c),("warm-specific",b))):
    lo,hi=boot(x)
    ax.scatter(x,np.full(x.size,i)+np.random.default_rng(i).normal(0,0.06,x.size),s=6,
               color=CC[lab],alpha=0.35,linewidth=0)
    ax.plot([np.median(x)],[i],"o",color="k",ms=7,zorder=3)
    ax.plot([lo,hi],[i,i],"-",color="k",lw=2,zorder=3)
    ax.annotate(f"median {np.median(x):+.4f}  [{lo:+.4f}, {hi:+.4f}]  n={x.size}  {100*(x<0).mean():.0f}% < 0",
                xy=(0.02,i+0.28),xycoords=("axes fraction","data"),fontsize=7.5,color="0.3")
ax.set_yticks(range(3)); ax.set_yticklabels(["cold-specific","shared","warm-specific"])
ax.set_xlabel("per-variant mean residual across gardens  (< 0 = declines)")
ax.set_xlim(np.percentile(np.r_[a,b,c],0.5),np.percentile(np.r_[a,b,c],99.5))
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_cold_vs_warm.png",dpi=130,bbox_inches="tight"); plt.show()
print(f"Mann-Whitney cold vs warm specific: p={stats.mannwhitneyu(a,b).pvalue:.4f}")
r=np.random.default_rng(0)
diff=[np.median(r.choice(a,a.size,True))-np.median(r.choice(b,b.size,True)) for _ in range(5000)]
print(f"difference of medians {np.median(a)-np.median(b):+.4f}, bootstrap 95% CI "
      f"[{np.percentile(diff,2.5):+.4f}, {np.percentile(diff,97.5):+.4f}], P(cold deeper)={np.mean(np.array(diff)<0):.2f}")
"""

md_close = """## Reading it

**The decline is global, not climate-graded.** Cold-specific insertions fall in 29/31 gardens,
including 11/11 cold gardens where their carrier lineages are *winning*. That is the outcome Miles's
diagnostic was designed to detect, and it is not what hitchhiking predicts.

**The restrictedness control rules out the obvious artefact.** Equally-restricted **indels** do not
decline (cold-specific indels 14/31, p = 0.96), so a narrow founder projection does not by itself
manufacture a negative slope. It takes being a >50 bp variant *and* being lineage-restricted.
`shared` insertions — 2,172 of 2,546, the bulk — do not decline at all, which is why the pooled
insertion set looked climate-graded rather than purged.

**Per variant, the excess sits on the global signature only.** Against the restricted-indel null,
cold-specific insertions run ~17 points above on global decline and ~2 points on climate-grading;
warm-specific ~12 points and ~0. So the climate-graded component is a general property of
lineage-restricted variants — indels show it just as much — while the global decline is specific to
large insertions. **The gradient is demography; the global decline is selection.**

**Do not over-read cold vs warm.** Cold-specific insertions decline about twice as deep as
warm-specific ones by point estimate (−0.068 vs −0.046 per variant), and the distributions differ
(Mann-Whitney p = 0.024), but the bootstrap CI on the difference of medians includes zero. The
direction is supported; the 2× is not. Cold-specific is individually significant (CI excludes 0),
warm-specific is not — that asymmetry is the better-supported claim.

**Caveats.**
- Deletions cannot arbitrate insertion-vs-size: only 25 cold-specific and 56 warm-specific exist,
  because deletions are barely lineage-restricted. Warm-specific deletions do trend down
  (−0.0399, p = 0.061), so "large restricted SVs decline" is the safe statement, not "insertions do".
- Three generations bounds detectable selection to roughly |3s|, so this sees only large effects.
- Matching is on the **mean** carrier-origin climate — the first moment of carrier composition, not
  composition itself. Variants restricted to one clade are compared with SNPs of the same mean that
  may be spread across many, and that residual mismatch is not controlled.
- `s` is a logit slope from a global-mode founder projection; it cannot separate an insertion from
  the haplotype carrying it, which is precisely why the indel control carries the argument here
  rather than the insertion result alone."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_setup),
    new_markdown_cell(md_diag), new_code_cell(code_diag),
    new_markdown_cell(md_call), new_code_cell(code_call),
    new_markdown_cell(md_dist), new_code_cell(code_insonly),
    new_markdown_cell(md_dist_old),
    new_code_cell(code_dist),
    new_markdown_cell(md_glob), new_code_cell(code_glob),
    new_markdown_cell(md_cw), new_code_cell(code_cw),
    new_markdown_cell(md_close),
])
ep = ExecutePreprocessor(timeout=1800, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
