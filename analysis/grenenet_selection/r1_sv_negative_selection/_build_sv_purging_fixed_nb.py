#!/usr/bin/env python
"""Build+execute `sv_purging_fixed.ipynb`: the purging-vs-hitchhiking analysis redone at the ALLELE
unit, consistently, after fixing the position-level carrier filter.
Load-only from the `_rec2` outputs of `_sv_temporal_climate_matched.py`, `_sv_purging_classify.py`,
`_sv_purging_annotate.py` (SV_FILTER=record SV_MIN_MAC=2). `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "notebooks", "sv_purging_fixed.ipynb")

md_title = """# Are insertions purged in GrENE-Net? — Truvari event unit, filter fixed

This is the corrected version of `sv_purging_vs_hitchhiking.ipynb`, which is kept unchanged for the
record.

## What was wrong

The temporal analyses filtered variants with `lib.founder_panel_keep()`, which checks the founder
carrier count **per genomic position**, not per variant. At a multi-allelic position every allele
passed if *any* allele there cleared the floor. Example, Chr2:7158959 — one ~2 kb insertion present
as 15 alleles (different inserted sequences):

| allele carriers | 21 | 5 | 4 | 3 | 2 | 2 | 2 | 2 | 2 | 1 | 1 | 1 | 1 | 1 | 1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

The 21-carrier allele passed a "≥ 12 carriers" filter, and so did all 14 others. One of them, with
3 carriers, was the top "cold-specific purged insertion" of the old analysis.

That **mixes units**: an allele was admitted on its *bubble's* frequency but then analysed as an
*allele* — labelled from its own carriers, matched on its own frequency. It was concentrated exactly
where the old result lived: 78% of the old "cold-specific" insertions had fewer than 12 carriers of
their own, and those carried the decline.

## The rule here

**One unit throughout — the SV event.** arch3 splits one insertion into one panel record per
assembly path, so the allele (record) is too fine: one insertion appears as several "alleles" with
small sequence differences. The bubble is too coarse in principle. Here the alleles of each pangenome
bubble are collapsed with **`truvari collapse` defaults** (size similarity ≥ 0.95, sequence similarity
≥ 0.95 by edit distance, ≤ 500 bp apart; greedy, no chaining, earliest allele as representative —
`_build_event_map.py`). A founder carries an event if it carries any of its alleles; the event's pool
frequency is the sum of its alleles' frequencies (founders almost never carry two alleles of one
event: 0.6% of carrier pairs). Short indels and SNPs stay per record.

Every variant is filtered on its **own** carriers and call rate (≥ 2 carriers, called in ≥ 90% of 231
founders), labelled cold-specific / shared / warm-specific from its own carriers, matched on its own
starting frequency and carrier-origin climate, and scored on its own trajectory
(`_build_event_units.py` → `_sv_temporal_climate_matched.py` → `_sv_purging_classify.py` →
`_sv_purging_annotate.py`, `SV_UNIT=event SV_MIN_MAC=2`, outputs `_event2`). The earlier allele-level
run (`_rec2`) is kept for comparison.

**Twin notebook.** `sv_purging_fixed_allele.ipynb` is the same analysis, cell for cell, with every
insertion / deletion record as its own variant (no grouping) — open the two side by side.

The statistic throughout: each variant's residual `s` against SNPs matched on starting frequency and
carrier-origin climate, summarised as the **fraction of the 31 gardens in which it declines**
(0.5 = no tendency)."""

code_setup = r"""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
from matplotlib.gridspec import GridSpec
import statsmodels.formula.api as smf
plt.rcParams.update({'figure.dpi':110,'font.size':8,'axes.linewidth':0.6})
G=os.path.abspath(os.path.join(os.getcwd(),"..","r1_sv_negative_selection","results","sv_adaptive"))
os.makedirs(f"{G}/plots",exist_ok=True)
TAG="_event2"      # Truvari events (SV_UNIT=event, SV_MIN_MAC=2); "_rec2" = the earlier allele run
def style(ax,axis="both"):
    ax.set_axisbelow(True); ax.grid(color="0.88",lw=0.7,axis=axis)
    for sp in ax.spines.values(): sp.set_visible(False)
def dotbox(a,x,v,color,rng,w=0.55,maxdots=400,colors=None):
    take=np.arange(v.size) if v.size<=maxdots else rng.choice(v.size,maxdots,replace=False)
    c=color if colors is None else [colors[i] for i in take]
    a.scatter(x+rng.uniform(-w*0.36,w*0.36,take.size),v[take],s=7,c=c,alpha=0.55,linewidth=0,zorder=2)
    a.boxplot(v,positions=[x],widths=w,patch_artist=True,showfliers=False,
              boxprops=dict(facecolor="none",edgecolor="0.1",linewidth=1.1),medianprops=dict(color="k",lw=1.8),
              whiskerprops=dict(color="0.2",lw=0.9),capprops=dict(color="0.2",lw=0.9),manage_ticks=False,zorder=3)
V=pd.read_csv(f"{G}/sv_purging_annotated{TAG}.csv")
V["abs_slope"]=V.slope_bio1.abs()
ORD=["cold-specific","shared","warm-specific"]
DOT={"cold-specific":"#2166AC","shared":"0.55","warm-specific":"#B2182B","other":"0.8"}
F=pd.read_csv(f"{G}/sv_origin_partition_founders.csv",dtype={"founder":str})
D=pd.read_csv(f"{G}/sv_temporal_climate_matched_sites{TAG}.csv").sort_values("bio1").reset_index(drop=True)
P=np.load(f"{G}/sv_temporal_climate_matched_persite{TAG}.npz")
"""

md_check = """## 0. Check — is every tested variant judged on its own carriers?

The filter floor is 2 carriers; every tested variant (SV event, or indel record) must clear it on its
**own** carrier count — for an event, the union of its alleles' carriers. The old analysis would fail
this check. Also shown: how the tested set compares with the old one."""

code_check = r"""
bad=int((V.n_carriers<2).sum())
print(f"tested variants: {len(V):,}   with < 2 carriers of their own: {bad}   -> {'OK' if bad==0 else 'FAIL'}")
print(f"call rate: min {V.call_rate.min():.3f}  (floor 0.90)")
print("\ntested set by class and carrier origin:")
print(pd.crosstab(V.cls,V.origin).loc[["insertion","deletion","indel"],ORD+["other"]].to_string())
old=pd.read_csv(f"{G}/sv_purging_annotated.csv")   # the old, position-filtered MAC>=12 run
oi=old[(old.cls=="insertion")&(old.origin=="cold-specific")]
print(f"\nold 'cold-specific' insertions: {len(oi)}; with >= 12 carriers of their own: {int((oi.n_carriers>=12).sum())}")
print("own-carrier distribution of cold-specific insertions now:",
      V[(V.cls=="insertion")&(V.origin=="cold-specific")].n_carriers.describe()[["count","50%","max"]].to_dict())
"""

md_main = """## 1. Insertions: globally, and split by carrier origin

Top: where the 231 founder ecotypes come from, and the tertile boundaries that define the groups.
Bottom left: fraction of gardens in which each insertion declines — every tested insertion first,
then split; dots coloured by the split. Each box is tested against 0.5 (one-sample Wilcoxon);
cold- vs warm-specific with a Wilcoxon rank-sum. Bottom right: how strongly each insertion tracks
garden temperature: the linear slope of its residual against garden bio1 (> 0 = it does better in
warmer gardens, < 0 = in colder ones), each group tested against 0."""

code_main = r"""
e1=F[F.group=="cold"].bio1.max(); e2=F[F.group=="warm"].bio1.min()
GC={"cold":"#2166AC","temperate":"0.6","warm":"#B2182B"}
I=V[V.cls=="insertion"].copy()
BOXES=[("all insertions",I)]+[(o,I[I.origin==o]) for o in ORD]
fig=plt.figure(figsize=(12.5,7.6))
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
h.annotate("cold-specific = variant carried by >= 2 cold-origin ecotypes and no warm-origin ones;  "
           "warm-specific = the mirror;  shared = carried in both.  Groups use each variant's OWN carriers.",
           xy=(0.0,-0.42),xycoords="axes fraction",fontsize=7.5,color="0.35")
rng=np.random.default_rng(0)
for pi,(col,ylab,ref) in enumerate([("frac_neg","fraction of gardens where residual < 0\n(0.5 = no tendency, higher = declines more consistently)",0.5),
                                    ("slope_bio1","linear slope of residual vs garden bio1  (per C)\n(> 0 = does better in warmer gardens, < 0 = in colder)",0.0)]):
    a=fig.add_subplot(gs[1,pi]); style(a,axis="y")
    if ref is not None: a.axhline(ref,color="0.35",lw=1.0,ls="--",zorder=1)
    TL=[]
    for i,(lab,d) in enumerate(BOXES):
        d=d[d[col].notna()]; v=d[col].to_numpy()
        dotbox(a,i,v,None,rng,colors=[DOT.get(o,"0.8") for o in d.origin])
        if ref is not None:
            p1=stats.wilcoxon(v-ref).pvalue
            fm="+.4f" if col=="slope_bio1" else ".3f"
            TL.append(f"{lab}\n(n={v.size:,})\nmedian {np.median(v):{fm}}\np = {p1:.1e}" if p1<1e-3 else
                      f"{lab}\n(n={v.size:,})\nmedian {np.median(v):{fm}}\np = {p1:.3f}")
        else:
            TL.append(f"{lab}\nmedian {np.median(v):.3f}")
    a.axvline(0.5,color="0.75",lw=0.8)
    ca=I[I.origin=="cold-specific"][col].dropna(); wa=I[I.origin=="warm-specific"][col].dropna()
    pv=stats.ranksums(ca,wa).pvalue
    lo_,hi_=a.get_ylim(); yb=hi_+0.18*(hi_-lo_)
    a.plot([1,1,3,3],[yb-0.02*(hi_-lo_),yb,yb,yb-0.02*(hi_-lo_)],color="0.15",lw=1.0)
    a.annotate(f"cold- vs warm-specific: Wilcoxon rank-sum p = {pv:.3g}",(2,yb),xytext=(0,3),
               textcoords="offset points",ha="center",va="bottom",fontsize=8,color="0.15")
    a.set_ylim(lo_,yb+0.12*(hi_-lo_))
    a.set_xticks(range(len(BOXES))); a.set_xticklabels(TL,fontsize=7.4); a.set_xlim(-0.6,len(BOXES)-0.4)
    a.set_ylabel(ylab)
    print(f"{col}: cold- vs warm-specific rank-sum p = {pv:.3g}")
    if col=="slope_bio1":
        for lab,d in BOXES:
            v=d[col].dropna(); print(f"   slope {lab:<15} median {v.median():+.4f}  >0 in {100*(v>0).mean():.1f}%  vs 0 p={stats.wilcoxon(v).pvalue:.2g}")
fig.legend(handles=[plt.Line2D([0],[0],marker="o",ls="",color=DOT[o],ms=6,label=o) for o in ORD],
           frameon=False,fontsize=8.5,loc="lower center",ncol=3,bbox_to_anchor=(0.5,-0.07))
fig.savefig(f"{G}/plots/sv_purging_fixed_insertions.png",dpi=130,bbox_inches="tight"); plt.show()
for lab,d in BOXES:
    v=d.frac_neg.dropna()
    print(f"{lab:<15} n={len(v):>5,}  median {v.median():.3f}  declining in >50% of gardens: {100*(v>0.5).mean():.1f}%  "
          f"vs 0.5 p={stats.wilcoxon(v-0.5).pvalue:.2g}")
"""

md_main_ctx = """## 1b. The same figure, insertions and deletions separately, by genomic context

Exactly §1's figure — founder-origin histogram, fraction of gardens declining, slope vs garden
bio1 — repeated for each (SV class × context) combination: insertions and deletions, each in CDS,
intron, and intergenic sequence. The histogram is identical every time (shown for reference); only
the variant set going into the two boxplots changes. Sample sizes shrink a lot in CDS (as few as
~35–55 per origin group, see the panel titles) — read those for direction, not precision."""

code_main_ctx = r"""
CTX3=("CDS","intron","intergenic")
def make_ctx_fig(cls_,cx):
    Ic=V[(V.cls==cls_)&(V.context==cx)].copy()
    BOXES=[(f"all {cls_}s",Ic)]+[(o,Ic[Ic.origin==o]) for o in ORD]
    fig=plt.figure(figsize=(12.5,7.6))
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
                   ha="center",va="top",fontsize=8,color=GC[g] if g!="temperate" else "0.3")
    h.set_xlabel("ecotype origin temperature  bio1 (C)"); h.set_ylabel("founder\necotypes")
    h.annotate(f"{cls_}s in {cx} only.  cold-specific = variant carried by >= 2 cold-origin ecotypes and no warm-origin ones;  "
               "warm-specific = the mirror;  shared = carried in both.  Groups use each variant's OWN carriers.",
               xy=(0.0,-0.42),xycoords="axes fraction",fontsize=7.5,color="0.35")
    rng=np.random.default_rng(0)
    for pi,(col,ylab,ref) in enumerate([("frac_neg","fraction of gardens where residual < 0\n(0.5 = no tendency, higher = declines more consistently)",0.5),
                                        ("slope_bio1","linear slope of residual vs garden bio1  (per C)\n(> 0 = does better in warmer gardens, < 0 = in colder)",0.0)]):
        a=fig.add_subplot(gs[1,pi]); style(a,axis="y")
        a.axhline(ref,color="0.35",lw=1.0,ls="--",zorder=1); TL=[]
        for i,(lab,d) in enumerate(BOXES):
            dd=d[d[col].notna()]; v=dd[col].to_numpy()
            if v.size==0:
                TL.append(f"{lab}\n(n=0)"); continue
            dotbox(a,i,v,None,rng,colors=[DOT.get(o,"0.8") for o in dd.origin])
            fm="+.4f" if col=="slope_bio1" else ".3f"
            if v.size>=2:
                p1=stats.wilcoxon(v-ref).pvalue
                TL.append(f"{lab}\n(n={v.size:,})\nmedian {np.median(v):{fm}}\np = {p1:.1e}" if p1<1e-3 else
                          f"{lab}\n(n={v.size:,})\nmedian {np.median(v):{fm}}\np = {p1:.3f}")
            else:
                TL.append(f"{lab}\n(n={v.size:,})\nmedian {np.median(v):{fm}}")
        a.axvline(0.5,color="0.75",lw=0.8)
        ca=Ic[Ic.origin=="cold-specific"][col].dropna(); wa=Ic[Ic.origin=="warm-specific"][col].dropna()
        pv=stats.ranksums(ca,wa).pvalue if len(ca)>=2 and len(wa)>=2 else np.nan
        lo_,hi_=a.get_ylim(); yb=hi_+0.18*(hi_-lo_)
        a.plot([1,1,3,3],[yb-0.02*(hi_-lo_),yb,yb,yb-0.02*(hi_-lo_)],color="0.15",lw=1.0)
        lab_p=f"cold- vs warm-specific: Wilcoxon rank-sum p = {pv:.3g}" if pv==pv else "cold- vs warm-specific: n too small"
        a.annotate(lab_p,(2,yb),xytext=(0,3),textcoords="offset points",ha="center",va="bottom",fontsize=8,color="0.15")
        a.set_ylim(lo_,yb+0.12*(hi_-lo_))
        a.set_xticks(range(len(BOXES))); a.set_xticklabels(TL,fontsize=7.4); a.set_xlim(-0.6,len(BOXES)-0.4)
        a.set_ylabel(ylab)
        print(f"{cls_} {cx} {col}: {lab_p}")
    fig.legend(handles=[plt.Line2D([0],[0],marker="o",ls="",color=DOT[o],ms=6,label=o) for o in ORD],
               frameon=False,fontsize=8.5,loc="lower center",ncol=3,bbox_to_anchor=(0.5,-0.07))
    fig.savefig(f"{G}/plots/sv_purging_fixed_{cls_}_{cx}.png",dpi=130,bbox_inches="tight"); plt.show()
    for lab,d in BOXES:
        v=d.frac_neg.dropna()
        if len(v)<2: print(f"{lab:<15} n={len(v)}"); continue
        print(f"{lab:<15} n={len(v):>5,}  median {v.median():.3f}  declining in >50% of gardens: {100*(v>0.5).mean():.1f}%  "
              f"vs 0.5 p={stats.wilcoxon(v-0.5).pvalue:.2g}")

for cls_ in ("insertion","deletion"):
    for cx in CTX3:
        print(f"\n=== {cls_} x {cx} ===")
        make_ctx_fig(cls_,cx)
"""

md_null = """## 2. Against the null — insertions vs indels at the same own-carrier count

Rarer variants decline more consistently in this statistic whatever they are, so insertions are
compared with **indels carried by the same number of founders** and from the same origin group.
Indels are the null: same estimator, same matching, same filter — and essentially biallelic, so
allele and event coincide for them."""

code_null = r"""
BANDS=[(2,5),(6,11),(12,10**6)]
blab=lambda lo,hi: f"{lo}-{hi}" if hi<10**6 else f"{lo}+"
rng=np.random.default_rng(0)
fig,ax=plt.subplots(1,3,figsize=(13,4.4),sharey=True)
rows=[]
for a,o in zip(ax,ORD):
    style(a,axis="y"); a.axhline(0.5,color="0.35",lw=1.0,ls="--",zorder=1); TL=[]
    for i,(lo,hi) in enumerate(BANDS):
        sel=lambda c: V[(V.cls==c)&(V.origin==o)&(V.n_carriers>=lo)&(V.n_carriers<=hi)].frac_neg.dropna().to_numpy()
        vi,vd=sel("insertion"),sel("indel")
        if vi.size>=5: dotbox(a,i-0.2,vi,DOT[o],rng,w=0.34)
        if vd.size>=5: dotbox(a,i+0.2,vd,"0.7",rng,w=0.34)
        p=stats.ranksums(vi,vd).pvalue if vi.size>=10 and vd.size>=10 else np.nan
        TL.append(f"{blab(lo,hi)} carriers\nins n={vi.size}  indel n={vd.size}\np = {p:.2g}")
        rows.append(dict(origin=o,carriers=blab(lo,hi),n_ins=vi.size,ins=np.median(vi) if vi.size else np.nan,
                         n_indel=vd.size,indel=np.median(vd) if vd.size else np.nan,p=p))
    a.annotate(f"{o}   (coloured = insertions, grey = indels)",xy=(0.02,1.01),xycoords="axes fraction",
               fontsize=8,color=DOT[o] if o!="shared" else "0.3")
    a.set_xticks(range(len(BANDS))); a.set_xticklabels(TL,fontsize=7.3)
    if o=="cold-specific": a.set_ylabel("fraction of gardens where residual < 0\n(0.5 = no tendency)")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_fixed_vs_indels.png",dpi=130,bbox_inches="tight"); plt.show()
print(pd.DataFrame(rows).round(3).to_string(index=False))
"""

md_garden = """## 3. Per garden — do they decline in cold gardens too?

One dot per garden: the median residual of that group's variants there. Filled = the group declines
in that garden. Insertions on top, indels (the null) below; one shared y-axis. A purged group sits
below zero across the whole temperature axis; a hitchhiking one only on the side where its carrier
lineages lose."""

code_garden = r"""
def garden_medians(cls,o):
    key={"insertion":"ins","indel":"indel","deletion":"del"}[cls]
    d=V[(V.cls==cls)&(V.origin==o)]
    if len(d)<20: return None
    ri=d.ridx.to_numpy(); med=[]
    for s in D.site.astype(int):
        r=P[f"{s}_res_{key}_matched"][ri]; r=r[~np.isnan(r)]
        med.append(np.median(r) if r.size>=20 else np.nan)
    return np.array(med)
ROWS=[("insertion","insertions"),("indel","indels (null)")]
M={(c,o):garden_medians(c,o) for c,_ in ROWS for o in ORD}
allv=np.concatenate([m for m in M.values() if m is not None]); lim=np.nanmax(np.abs(allv))*1.1
fig,ax=plt.subplots(2,3,figsize=(12,6.2),sharex=True,sharey=True)
for i,(c,rl) in enumerate(ROWS):
    for j,o in enumerate(ORD):
        a=ax[i,j]; style(a); med=M[(c,o)]; g=~np.isnan(med); neg=med[g]<0
        a.axhline(0,color="0.3",lw=1.0)
        a.scatter(D.bio1[g][neg],med[g][neg],s=34,color=DOT[o],edgecolor="0.2",linewidth=0.4,zorder=3)
        a.scatter(D.bio1[g][~neg],med[g][~neg],s=34,facecolor="white",edgecolor=DOT[o],linewidth=1.0,zorder=3)
        a.axhline(np.median(med[g]),color=DOT[o],lw=1.6,ls="--",zorder=2)
        pw=stats.wilcoxon(med[g]).pvalue
        a.annotate(f"{rl}, {o}  (n={int(((V.cls==c)&(V.origin==o)).sum()):,})\nbelow 0 in {int(neg.sum())}/{g.sum()} gardens"
                   f"   median {np.median(med[g]):+.3f}   p = {pw:.2g}",
                   xy=(0.03,0.04),xycoords="axes fraction",va="bottom",fontsize=7.6,color="0.25",
                   bbox=dict(facecolor="white",alpha=0.8,edgecolor="none",pad=1.5))
        a.set_ylim(-lim,lim)
        if i==1: a.set_xlabel("garden temperature  bio1 (C)")
        if j==0: a.set_ylabel("median residual vs matched SNPs\n(< 0 = declines)")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_fixed_by_garden.png",dpi=130,bbox_inches="tight"); plt.show()
"""

md_garden_insdel = """## 3b. Insertions vs deletions, before and after carrier-climate matching

Same figure as §3, but insertions (top) and **deletions** (bottom) instead of the indel null, shown
twice: first against SNPs matched on starting frequency **only** (`p0only` -- no carrier-climate
control), then against SNPs matched on starting frequency **and carrier-origin climate**
(`matched`, as in §3). Comparing the two versions is the carrier-climate-matching check from
`sv_temporal_climate_matched.ipynb` at the per-class, per-garden level."""

code_garden_insdel = r"""
def per_garden_fig(tag,ROWS,fname,ylab):
    def garden_medians(cls,o):
        key={"insertion":"ins","indel":"indel","deletion":"del"}[cls]
        d=V[(V.cls==cls)&(V.origin==o)]
        if len(d)<20: return None
        ri=d.ridx.to_numpy(); med=[]
        for s in D.site.astype(int):
            r=P[f"{s}_res_{key}_{tag}"][ri]; r=r[~np.isnan(r)]
            med.append(np.median(r) if r.size>=20 else np.nan)
        return np.array(med)
    M={(c,o):garden_medians(c,o) for c,_ in ROWS for o in ORD}
    allv=np.concatenate([m for m in M.values() if m is not None]); lim=np.nanmax(np.abs(allv))*1.1
    fig,ax=plt.subplots(2,3,figsize=(12,6.2),sharex=True,sharey=True)
    for i,(c,rl) in enumerate(ROWS):
        for j,o in enumerate(ORD):
            a=ax[i,j]; style(a); med=M[(c,o)]; g=~np.isnan(med); neg=med[g]<0
            a.axhline(0,color="0.3",lw=1.0)
            a.scatter(D.bio1[g][neg],med[g][neg],s=34,color=DOT[o],edgecolor="0.2",linewidth=0.4,zorder=3)
            a.scatter(D.bio1[g][~neg],med[g][~neg],s=34,facecolor="white",edgecolor=DOT[o],linewidth=1.0,zorder=3)
            a.axhline(np.median(med[g]),color=DOT[o],lw=1.6,ls="--",zorder=2)
            pw=stats.wilcoxon(med[g]).pvalue
            a.annotate(f"{rl}, {o}  (n={int(((V.cls==c)&(V.origin==o)).sum()):,})\nbelow 0 in {int(neg.sum())}/{g.sum()} gardens"
                       f"   median {np.median(med[g]):+.3f}   p = {pw:.2g}",
                       xy=(0.03,0.04),xycoords="axes fraction",va="bottom",fontsize=7.6,color="0.25",
                       bbox=dict(facecolor="white",alpha=0.8,edgecolor="none",pad=1.5))
            a.set_ylim(-lim,lim)
            if i==1: a.set_xlabel("garden temperature  bio1 (C)")
            if j==0: a.set_ylabel(ylab)
    fig.tight_layout(); fig.savefig(f"{G}/plots/{fname}",dpi=130,bbox_inches="tight"); plt.show()

ROWS_ID=[("insertion","insertions"),("deletion","deletions")]
print("=== p0-only matching (no carrier-climate control) ===")
per_garden_fig("p0only",ROWS_ID,"sv_purging_fixed_by_garden_ins_del_p0only.png",
               "median residual vs matched SNPs\n(p0 only; < 0 = declines)")
print("=== p0 + carrier-climate matching ===")
per_garden_fig("matched",ROWS_ID,"sv_purging_fixed_by_garden_ins_del_matched.png",
               "median residual vs matched SNPs\n(p0 + carrier climate; < 0 = declines)")
"""

md_size = """## 4. Size and inserted sequence (TE vs other)

Insertion size from the panel; inserted-sequence class from `insertion_tair10_class.csv` (joined on the
arch3 panel record index). Size is tested **within** each origin group, because restricted insertions
are larger than shared ones, and jointly with TE content and own carrier count, because size and TE
content are correlated."""

code_size = r"""

S=pd.read_csv(f"{G}/sv_purging_insertions{TAG}_seqclass.csv")
bins=[50,100,200,500,2000,1e7]; lab=["51-100","101-200","201-500","501-2k",">2k"]
S["sbin"]=pd.cut(S.size_bp,bins,labels=lab)
rng=np.random.default_rng(0)
def mean_ci(v,rng,B=2000):
    bs=rng.choice(v,(B,v.size)).mean(1); return v.mean(),np.percentile(bs,2.5),np.percentile(bs,97.5)
fig,ax=plt.subplots(1,3,figsize=(13,4.3),sharey=True)
for a,o in zip(ax,ORD):
    style(a,axis="y"); a.axhline(0.5,color="0.35",lw=1.0,ls="--",zorder=1); d=S[S.origin==o]; TL=[]
    for i,bb in enumerate(lab):
        v=d[d.sbin==bb].frac_neg.dropna().to_numpy()
        if v.size>=5: dotbox(a,i,v,DOT[o],rng)
        TL.append(f"{bb}\nn={v.size}")
    mm=smf.ols("frac_neg ~ lsize + te + lcar",data=d).fit()
    a.annotate(f"{o}   size {mm.params.lsize:+.3f}/log10 bp (p={mm.pvalues.lsize:.2g}), TE {mm.params.te:+.3f} (p={mm.pvalues.te:.2g})",
               xy=(0.02,1.01),xycoords="axes fraction",fontsize=7.6,color=DOT[o] if o!="shared" else "0.3")
    a.set_xticks(range(len(lab))); a.set_xticklabels(TL,fontsize=7.5); a.set_xlabel("insertion size (bp)")
    if o=="cold-specific": a.set_ylabel("fraction of gardens where residual < 0\n(0.5 = no tendency)")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_fixed_size.png",dpi=130,bbox_inches="tight"); plt.show()
print("size vs TE content: Spearman rho =",round(S[["lsize","te"]].corr(method="spearman").iloc[0,1],3))
m=smf.ols("frac_neg ~ (lsize + te) * C(origin, Treatment('shared')) + lcar",data=S).fit()
print("\nfrac_neg ~ (log size + TE) x origin + log own carriers")
print(m.summary2().tables[1][["Coef.","P>|t|"]].round(4).to_string())
"""

md_size_read = """**Reading it (events).** Size and TE content barely matter. Across all insertions, larger ones and
TE-derived ones decline very slightly *less* (size −0.013 per log10 bp, p = 0.03; TE −0.016, p = 0.02),
and neither interacts with origin (all p ≥ 0.06). The allele-level pattern — larger cold-specific
insertions falling more and larger warm-specific ones rising more — does not hold on events."""

md_ctx = """## 5. Genomic context, rechecked on the fixed set

Context as in `sv_purging_followups.ipynb` (CDS > UTR > intron > reference TE > promoter ≤ 2 kb >
intergenic). Insertions as dots and boxes, with the **mean ± 95% bootstrap CI** as a white dot; the
grey diamond is the mean ± 95% CI of **indels in the same context and origin group** — the null.
Group sizes differ a lot (cold-specific: CDS 47, UTR 29 vs promoter 240): that doesn't bias a mean or median, but the
small contexts have much wider CIs, which the bars show. Tested jointly with size and own carrier count, since size
now matters and differs between contexts."""

code_ctx = r"""

CTX=["CDS","UTR","intron","TE","promoter","intergenic"]
A=S.copy()
rng=np.random.default_rng(0)
def mean_ci(v,rng,B=2000):
    bs=rng.choice(v,(B,v.size)).mean(1); return v.mean(),np.percentile(bs,2.5),np.percentile(bs,97.5)
fig,ax=plt.subplots(1,3,figsize=(13,4.3),sharey=True)
for a,o in zip(ax,ORD):
    style(a,axis="y"); a.axhline(0.5,color="0.35",lw=1.0,ls="--",zorder=1); d=A[A.origin==o]; TL=[]
    for i,c in enumerate(CTX):
        v=d[d.context==c].frac_neg.dropna().to_numpy()
        if v.size>=5:
            dotbox(a,i,v,DOT[o],rng)
            mu,lo,hi=mean_ci(v,rng)
            a.errorbar([i],[mu],yerr=[[mu-lo],[hi-mu]],fmt="o",ms=6,mfc="white",mec="k",mew=1.3,ecolor="k",elinewidth=1.6,capsize=3,zorder=5)
        nd=V[(V.cls=="indel")&(V.origin==o)&(V.context==c)].frac_neg.dropna().to_numpy()
        if nd.size>=10:
            mu,lo,hi=mean_ci(nd,rng)
            a.errorbar([i+0.38],[mu],yerr=[[mu-lo],[hi-mu]],fmt="D",ms=5.5,color="0.45",elinewidth=1.3,capsize=2.5,zorder=4)
        TL.append(f"{c}\nn={v.size}")
    a.annotate(o,xy=(0.02,1.01),xycoords="axes fraction",fontsize=8.5,color=DOT[o] if o!="shared" else "0.3")
    a.set_xticks(range(len(CTX))); a.set_xticklabels(TL,fontsize=7.5)
    if o=="cold-specific": a.set_ylabel("fraction of gardens where residual < 0\n(0.5 = no tendency)")
fig.tight_layout()
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],marker="o",ls="",mfc="white",mec="k",mew=1.3,ms=6,label="insertion mean ± 95% CI"),
                    Line2D([],[],marker="D",ls="",color="0.45",ms=5.5,label="indel mean ± 95% CI (same context)")],
           frameon=False,fontsize=8,loc="lower center",bbox_to_anchor=(0.5,0.98),ncol=2)
fig.savefig(f"{G}/plots/sv_purging_fixed_context.png",dpi=130,bbox_inches="tight"); plt.show()
t=A.groupby(["context","origin"]).frac_neg.agg(["median","size"])
print("insertions: median fraction of gardens declining (n)")
print(t["median"].unstack().reindex(CTX)[ORD].round(3).astype(str).add(" (").add(t["size"].unstack().reindex(CTX)[ORD].astype(str)).add(")").to_string())
A["ctx3"]=A.context.map({"CDS":"coding","UTR":"coding","intron":"intron","TE":"outside","promoter":"outside","intergenic":"outside"})
for f,lab_ in (("frac_neg ~ C(ctx3, Treatment('outside')) * C(origin, Treatment('shared')) + lcar","context x origin + carriers"),
               ("frac_neg ~ C(ctx3, Treatment('outside')) * C(origin, Treatment('shared')) + lsize * C(origin, Treatment('shared')) + lcar","... + size x origin")):
    m=smf.ols(f,data=A).fit(); tb=m.summary2().tables[1][["Coef.","P>|t|"]]
    print(f"\n{lab_}:"); print(tb[tb.index.str.contains("ctx3")].round(4).to_string())
# genic (CDS+UTR+intron) vs strict intergenic, cold-specific
print("\ncold-specific, genic (CDS+UTR+intron) vs intergenic: mean frac_neg [95% CI]")
for cls_ in ("insertion","indel"):
    d=V[(V.cls==cls_)&(V.origin=="cold-specific")]
    g=d[d.context.isin(["CDS","UTR","intron"])].frac_neg.dropna().to_numpy(); n_=d[d.context=="intergenic"].frac_neg.dropna().to_numpy()
    mg,lg,hg=mean_ci(g,rng); mn,ln,hn=mean_ci(n_,rng)
    print(f"  {cls_:<9} genic {mg:.3f} [{lg:.3f},{hg:.3f}] n={g.size}   intergenic {mn:.3f} [{ln:.3f},{hn:.3f}] n={n_.size}   diff {mg-mn:+.3f}  rank-sum p={stats.ranksums(g,n_).pvalue:.2g}")
# insertions vs indels within context, cold-specific only (the group that declines)
print("\ncold-specific: insertions vs indels in the same context (rank-sum)")
for c in CTX:
    a_=A[(A.origin=="cold-specific")&(A.context==c)].frac_neg; b_=V[(V.cls=="indel")&(V.origin=="cold-specific")&(V.context==c)].frac_neg
    if len(a_)>=10 and len(b_)>=10:
        print(f"  {c:<11} ins {a_.median():.3f} (n={len(a_)})  indel {b_.median():.3f} (n={len(b_)})  p={stats.ranksums(a_,b_).pvalue:.2g}")
"""

md_freq = """## 6. Genic vs intergenic, by frequency — are the purged genic insertions the rare ones?

If genic cold-specific insertions were deleterious, the classic expectation is that they sit at
**rarer** natural frequencies than intergenic ones, and that the genic-vs-intergenic difference in
decline depends on frequency. Cold-specific insertions (and cold-specific indels as the reference),
genic = CDS + UTR + intron, intergenic = strict intergenic (promoter and TE excluded).

- **Left / middle:** own-carrier spectrum, genic vs intergenic, among insertions that decline in
  ≥ 65% of gardens ("globally purged") and among the rest.
- **Right:** mean fraction of gardens declining per carrier bin (bootstrap 95% CI); solid =
  insertions, dotted = indels. Below: `frac_neg ~ genic × log10 carriers`.

Note the temporal set needs p0 ≥ 0.02, which removes most 2–3-carrier variants — the rarest end is
not testable here."""

code_freq = r"""
V6=pd.read_csv(f"{G}/sv_purging_annotated{TAG}.csv")
V6["reg"]=np.where(V6.context.isin(["CDS","UTR","intron"]),"genic",np.where(V6.context=="intergenic","intergenic","other"))
V6["lcar"]=np.log10(V6.n_carriers); V6["purged"]=V6.frac_neg>=0.65
CI=V6[(V6.cls=="insertion")&(V6.origin=="cold-specific")&V6.reg.isin(["genic","intergenic"])].copy()
N=V6[(V6.cls=="indel")&(V6.origin=="cold-specific")&V6.reg.isin(["genic","intergenic"])].copy()
print("cold-specific insertions by own-carrier count:",CI.n_carriers.value_counts().sort_index().to_dict())
BINS=[(2,4),(5,5),(6,6),(7,8),(9,11),(12,40)]
lab=lambda lo,hi: f"{lo}" if lo==hi else f"{lo}-{hi}"
COL={"genic":"#B2182B","intergenic":"#2166AC"}
fig,ax=plt.subplots(1,3,figsize=(14,4.3))
# A: spectrum, purged vs not, genic vs intergenic
for a,(pl,pv) in zip(ax[:2],(("globally purged (>=65% of gardens)",True),("not purged (<65%)",False))):
    a.set_axisbelow(True); a.grid(color="0.88",lw=0.7,axis="y")
    for sp in a.spines.values(): sp.set_visible(False)
    w=0.38
    for k,r in enumerate(("genic","intergenic")):
        x=CI[(CI.reg==r)&(CI.purged==pv)].n_carriers
        share=[((x>=lo)&(x<=hi)).mean() for lo,hi in BINS]
        a.bar(np.arange(len(BINS))+(k-0.5)*w,share,width=w*0.95,color=COL[r],alpha=0.8,label=f"{r} (n={len(x)})")
    g=CI[(CI.reg=="genic")&(CI.purged==pv)].n_carriers; n=CI[(CI.reg=="intergenic")&(CI.purged==pv)].n_carriers
    a.annotate(f"cold-specific insertions, {pl}\ngenic vs intergenic spectrum: rank-sum p = {stats.ranksums(g,n).pvalue:.2g}",
               xy=(0.98,0.97),xycoords="axes fraction",ha="right",va="top",fontsize=7.8,color="0.25")
    a.set_xticks(range(len(BINS))); a.set_xticklabels([lab(*b) for b in BINS])
    a.set_xlabel("own carriers among the 231 founders"); a.set_ylabel("share of variants"); a.legend(frameon=False,fontsize=7.5,loc="center right")
# B: decline as a function of frequency
a=ax[2]; a.set_axisbelow(True); a.grid(color="0.88",lw=0.7)
for sp in a.spines.values(): sp.set_visible(False)
a.axhline(0.5,color="0.35",lw=1.0,ls="--")
rng=np.random.default_rng(0); rows=[]
for data,ls,tag in ((CI,"-","insertions"),(N,":","indels")):
    for r in ("genic","intergenic"):
        mu,lo_,hi_=[],[],[]
        for lo,hi in BINS:
            v=data[(data.reg==r)&(data.n_carriers>=lo)&(data.n_carriers<=hi)].frac_neg.to_numpy()
            if v.size<8: mu.append(np.nan); lo_.append(np.nan); hi_.append(np.nan); continue
            bs=[rng.choice(v,v.size,True).mean() for _ in range(1000)]
            mu.append(v.mean()); lo_.append(np.percentile(bs,2.5)); hi_.append(np.percentile(bs,97.5))
            rows.append(dict(cls=tag,region=r,carriers=lab(lo,hi),n=v.size,mean=v.mean()))
        x=np.arange(len(BINS))+(0.08 if r=="intergenic" else -0.08)
        a.errorbar(x,mu,yerr=[np.array(mu)-np.array(lo_),np.array(hi_)-np.array(mu)],fmt="o",ls=ls,color=COL[r],
                   alpha=1 if tag=="insertions" else 0.45,capsize=2,label=f"{tag}, {r}")
a.set_xticks(range(len(BINS))); a.set_xticklabels([lab(*b) for b in BINS])
a.set_xlabel("own carriers among the 231 founders"); a.set_ylabel("mean fraction of gardens declining")
a.legend(frameon=False,fontsize=7,loc="lower left")
a.annotate("cold-specific: decline vs frequency\nsolid = insertions, dotted = indels",xy=(0.98,0.97),xycoords="axes fraction",
           ha="right",va="top",fontsize=7.8,color="0.25")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_cold_genic_by_frequency.png",dpi=130,bbox_inches="tight"); plt.show()
print(pd.DataFrame(rows).pivot_table(index="carriers",columns=["cls","region"],values=["mean","n"]).round(3).to_string())
for data,tag in ((CI,"insertions"),(N,"indels")):
    d=data.assign(genic=(data.reg=="genic").astype(int))
    m=smf.ols("frac_neg ~ genic * lcar",data=d).fit()
    print(f"\n{tag}: frac_neg ~ genic x log10 carriers")
    print(m.summary2().tables[1][["Coef.","P>|t|"]].round(4).to_string())
"""

md_freq_read = """**Reading it (events).** Numbers per bin are small (14–29), CIs overlap throughout, and the
genic × frequency interaction is not significant (insertions p = 0.29, indels p = 0.99). Among
globally purged cold-specific insertions, genic ones are not rarer than intergenic ones. No
frequency-dependent genic signal in the trajectory statistic."""

md_sfs = """## 7. The site-frequency spectra — by context, by carrier origin, purged vs not

The natural-panel spectra from `_sv_sfs_context_origin.py` (all 231 founders; called in ≥ 90%,
≥ 1 carrier; SNPs and short indels a random 5% of records; insertions and deletions (> 50 bp) as **events** —
the alleles of a pangenome bubble collapsed with `truvari collapse` defaults, size and sequence
similarity ≥ 0.95, `_build_event_map.py`; a founder carries an event if it carries any of its alleles).
Events, not alleles: arch3 splits one insertion into one record per assembly path, and counting those
records separately pushes insertions artificially rare. Built directly from the panel, no MAC floor.
x = carrier founders (log).

- **7a** all variants pooled, per context (the view in `sv_purging_followups.ipynb`).
- **7b** the same split by carrier origin (cold-specific ≥ 2 cold & 0 warm carriers; warm-specific
  the mirror; shared = carriers in both). Table: share with ≤ 3 carriers, fold over SNPs of the same
  origin and context.
- **7c** cold-specific variants in the temporal set (events / indel records), globally purged (declining in ≥ 65% of
  gardens) vs not, insertions vs indels. Only ~3–25 carriers are testable (p0 ≥ 0.02)."""

code_sfs = r"""
S7=np.load(f"{G}/sv_sfs_context_origin_event.npz")   # "*_bubble" classes here = Truvari EVENTS
CTX7=["CDS","UTR","intron","TE","promoter","intergenic"]
OC={"cold-specific":"#2166AC","shared":"0.3","warm-specific":"#B2182B"}
LN=[("snp","SNPs","0.3","-",1.6),("indel","short indels","0.65","-",1.6),("ins_bubble","insertions (events)","#D95F02","-",2.0),
    ("del_bubble","deletions (events)","#1B9E77","-",1.6)]
def osel(c,o):
    nc,nw=S7[f"{c}_nc"],S7[f"{c}_nw"]
    return {"all":np.ones(nc.size,bool),"cold-specific":(nc>=2)&(nw==0),"warm-specific":(nw>=2)&(nc==0),"shared":(nc>0)&(nw>0)}[o]
def ecdf(a,v,**k):
    v=np.sort(np.asarray(v,float)); a.step(np.r_[1,v],np.r_[0,np.arange(1,v.size+1)/v.size],where="post",**k)
def logx(a,ticks,lim):
    a.set_xscale("log"); a.set_xticks(ticks); a.set_xticklabels(ticks); a.set_xlim(*lim); a.minorticks_off()
def sfs_grid(origins,fname):
    fig,ax=plt.subplots(len(origins),6,figsize=(17,2.9*len(origins)+0.4),sharex=True,sharey=True,squeeze=False); rows=[]
    for r,o in enumerate(origins):
        for k,x in enumerate(CTX7):
            a=ax[r,k]; style(a)
            for c,lab,col,ls,lw in LN:
                n=S7[f"{c}_n"][osel(c,o)&(S7[f"{c}_ctx"]==k)]
                if n.size<20: continue
                ecdf(a,n,color=col,ls=ls,lw=lw,label=lab); rows.append(dict(origin=o,context=x,cls=c,N=n.size,le3=(n<=3).mean()))
            n_ins=(osel("ins_bubble",o)&(S7["ins_bubble_ctx"]==k)).sum()
            a.annotate(f"{x}\ninsertions n={n_ins:,}",xy=(0.97,0.05),xycoords="axes fraction",ha="right",fontsize=7.5,color="0.25")
            if k==0: a.set_ylabel(f"{o}\ncumulative distribution",color=OC.get(o,"0.2"))
            if r==len(origins)-1: a.set_xlabel("carriers among the 231 founders")
            logx(a,[1,2,5,10,20,50,100,200] if o=="all" or len(origins)>1 else [2,5,10,20],(0.9,235))
    ax[0,0].set_ylim(0,1.01)
    h,l=ax[0,0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",bbox_to_anchor=(0.5,0.995),ncol=5,frameon=False,fontsize=8.5)
    fig.tight_layout(); fig.savefig(f"{G}/plots/{fname}",dpi=120,bbox_inches="tight"); plt.show()
    return pd.DataFrame(rows)
print("7a — all variants pooled")
T0=sfs_grid(["all"],"sv_sfs_context_pooled.png")
print("7b — split by carrier origin")
T=sfs_grid(ORD,"sv_sfs_context_by_origin.png")
T=pd.concat([T0,T]); snp=T[T.cls=="snp"].set_index(["origin","context"]).le3
T["fold"]=T.le3.values/snp.reindex(pd.MultiIndex.from_frame(T[["origin","context"]])).values
for c in ("ins_bubble","indel","del_bubble"):
    print(f"\n{c}: share with <= 3 carriers, fold over SNPs (same origin & context)")
    print(T[T.cls==c].pivot(index="origin",columns="context",values="fold")[CTX7].reindex(["all"]+ORD).round(2).to_string())
"""

code_sfs_purged = r"""
V7=pd.read_csv(f"{G}/sv_purging_annotated{TAG}.csv"); V7["purged"]=V7.frac_neg>=0.65
fig,ax=plt.subplots(1,6,figsize=(17,3.4),sharey=True)
print("cold-specific, temporal set: median own carriers purged vs not (n), rank-sum p")
for k,x in enumerate(CTX7):
    a=ax[k]; style(a); line=f"  {x:<10}"; txt=[x]
    for cls_,col in (("insertion","#D95F02"),("indel","0.55")):
        d=V7[(V7.cls==cls_)&(V7.origin=="cold-specific")&(V7.context==x)]
        p_,q_=d[d.purged].n_carriers,d[~d.purged].n_carriers
        ecdf(a,p_,color=col,lw=2.0,label=f"{cls_}, purged (≥ 65% of gardens)")
        ecdf(a,q_,color=col,lw=1.3,ls="--",label=f"{cls_}, not purged")
        pv=stats.ranksums(p_,q_).pvalue
        line+=f"  {cls_[:3]} {p_.median():.0f} vs {q_.median():.0f} (n={len(p_)}/{len(q_)}, p={pv:.2g})"
        txt.append(f"{cls_[:3]} n={len(p_)}/{len(q_)}, p={pv:.2g}")
    print(line)
    a.annotate("\n".join(txt),xy=(0.03,0.97),xycoords="axes fraction",ha="left",va="top",fontsize=7.3,color="0.25")
    logx(a,[2,3,5,10,20],(1.8,30)); a.set_xlabel("own carriers among the 231 founders")
ax[0].set_ylabel("cold-specific\ncumulative distribution",color=OC["cold-specific"]); ax[0].set_ylim(0,1.01)
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",bbox_to_anchor=(0.5,0.99),ncol=4,frameon=False,fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_sfs_cold_purged_vs_not_by_context.png",dpi=120,bbox_inches="tight"); plt.show()
"""

md_sfs_cds = """**7b′ — CDS only: insertions vs SNPs, per origin.** Insertion events solid, SNPs dashed, coloured
by carrier origin; intergenic alongside as the reference context. Compare each solid line with the
dashed line of **its own colour** — the three origin classes cover different frequency ranges by
definition (restricted variants can't exceed their group's 77 founders), so cold vs shared vs warm
lines are not comparable to each other directly."""

code_sfs_cds = r"""
fig,ax=plt.subplots(1,2,figsize=(10,4.2),sharey=True)
print("share with <= 3 carriers (insertion events vs SNPs), and KS test ins vs SNP within origin")
for a,x in zip(ax,["CDS","intergenic"]):
    style(a); k=CTX7.index(x)
    for o in ORD:
        ni=S7["ins_bubble_n"][osel("ins_bubble",o)&(S7["ins_bubble_ctx"]==k)]
        ns=S7["snp_n"][osel("snp",o)&(S7["snp_ctx"]==k)]
        ecdf(a,ni,color=OC[o],lw=2.2,label=f"{o}: insertions (n={ni.size:,})")
        ecdf(a,ns,color=OC[o],lw=1.3,ls="--",label=f"{o}: SNPs (n={ns.size:,})")
        ks=stats.ks_2samp(ni,ns).pvalue
        print(f"  {x:<10} {o:<14} ins {(ni<=3).mean():.2f}  SNP {(ns<=3).mean():.2f}  fold {(ni<=3).mean()/(ns<=3).mean():.2f}  "
              f"median {np.median(ni):.0f} vs {np.median(ns):.0f}  KS p={ks:.2g}")
    a.annotate(x,xy=(0.97,0.05),xycoords="axes fraction",ha="right",fontsize=9,color="0.25")
    logx(a,[2,5,10,20,50,100,200],(1.8,235)); a.set_xlabel("carriers among the 231 founders")
ax[0].set_ylabel("cumulative distribution"); ax[0].set_ylim(0,1.01)
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",bbox_to_anchor=(0.5,0.99),ncol=3,frameon=False,fontsize=8)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_sfs_cds_ins_vs_snp_by_origin.png",dpi=130,bbox_inches="tight"); plt.show()
"""

md_sfs_read = """**Reading it (events).**
- **Pooled (7a):** insertion events are rarer than SNPs in CDS (1.47× at ≤ 3 carriers) and UTR
  (1.32×), barely elsewhere (TE 1.02×, intergenic 1.08×).
- **By origin (7b):** that genic rare-excess belongs to **shared** insertions — CDS 2.65×, UTR 1.70×,
  promoter 1.41×, but intergenic 1.02× and TE 1.05×: the classic purifying pattern. Shared deletions
  show a weaker version (CDS 1.75×). Cold- and warm-specific insertions sit close to SNPs (≤ 1.27×).
- **7b′:** shared CDS insertion events have a median of 13 carriers vs 25 for shared CDS SNPs
  (KS p = 6e-15); in intergenic sequence 13 vs 17.
- **Purged vs not (7c):** among cold-specific insertions, purged ones are slightly rarer in promoter
  (median 6 vs 7, p = 0.015) and intergenic sequence (p = 0.05), not in CDS, UTR, intron or TE.
- So the insertions with a purified natural spectrum are the **shared genic** ones — the class whose
  founder burden and population load fall in GrENE-Net (`sv_burden_lmm.ipynb`) — not the
  cold-specific insertions."""

md_close = """## Reading it (events, `truvari collapse` defaults)

**The check passes.** Every tested variant clears the carrier floor on its own (union) carriers
(0 of 36,532 below it).

**1. Per variant, insertions show no purging.** Across all 7,800 insertion events the median fraction
of gardens declining is 0.484 — if anything slightly *below* 0.5 (p = 0.01).

**2. No lineage-direction effect in that statistic.** Cold-specific insertion events: median 0.516
(p = 0.19); warm-specific: 0.484 (p = 0.27); cold vs warm p = 0.09.

**2b. But they follow their home climate.** The signed slope of each insertion's residual against
garden temperature separates the groups: cold-specific insertions do relatively better in **colder**
gardens (median −0.0068 per °C, 39% positive, p = 6e-13), warm-specific ones in **warmer** gardens
(+0.0063, 58% positive, p = 8e-8; cold vs warm p = 1e-19); shared ones lean slightly warm (+0.0021).
This is beyond SNPs matched on the *mean* climate of their carriers — a climate-graded (home-garden)
pattern, not a global decline. Selection on the insertions and tighter linkage to climate-adapted
haplotypes are not separable here.

**3. Size, TE content and genomic context add little.** Larger and TE-derived insertions decline very
slightly less; coding vs outside +0.017 (p = 0.08); no origin interaction survives.

**4. The natural spectrum points at shared genic insertions** (§7: shared CDS 2.65× rarer than SNPs,
intergenic 1.02×) — not at the restricted classes.

**What changed from the allele-level version.** With alleles, cold-specific insertions declined in 27/31
gardens (median 0.581), warm-specific ones rose in 24/31, and size amplified both. On events none of
that holds: 70% of the allele-level "cold-specific" insertions are fragments of events that are
actually shared (median 3 own carriers vs 18 for their event), so those patterns tracked the cold
haplotype of a shared insertion, not the insertion itself.

**Where the selection signal is.** The per-variant trajectory statistic (one noisy fraction per
variant) does not detect it. The founder-level burden model and the population load test do: shared
CDS insertion load predicts lower founder fitness and falls beyond matched SNPs in 27/31 gardens —
see `sv_burden_lmm.ipynb`.

**Relation to the old notebook.** `sv_purging_vs_hitchhiking.ipynb` reported a strong cold-specific
decline (median 0.645, p = 9e−12). That came from the position-level filter bug plus allele counting;
it survives neither fix."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_setup),
    new_markdown_cell(md_check), new_code_cell(code_check),
    new_markdown_cell(md_main), new_code_cell(code_main),
    new_markdown_cell(md_main_ctx), new_code_cell(code_main_ctx),
    new_markdown_cell(md_null), new_code_cell(code_null),
    new_markdown_cell(md_garden), new_code_cell(code_garden),
    new_markdown_cell(md_garden_insdel), new_code_cell(code_garden_insdel),
    new_markdown_cell(md_size), new_code_cell(code_size), new_markdown_cell(md_size_read),
    new_markdown_cell(md_ctx), new_code_cell(code_ctx),
    new_markdown_cell(md_freq), new_code_cell(code_freq), new_markdown_cell(md_freq_read),
    new_markdown_cell(md_sfs), new_code_cell(code_sfs), new_markdown_cell(md_sfs_cds), new_code_cell(code_sfs_cds), new_code_cell(code_sfs_purged), new_markdown_cell(md_sfs_read),
    new_markdown_cell(md_close),
])
ep = ExecutePreprocessor(timeout=1800, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
