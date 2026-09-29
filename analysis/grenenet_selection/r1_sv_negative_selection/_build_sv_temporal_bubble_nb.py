#!/usr/bin/env python
"""Build+execute the bubble-vs-allele temporal notebook.

Load-only from `_sv_temporal_bubble.py` (per-site shifts) and `_sv_temporal_bubble_units.py`
(per-bubble table). `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "notebooks", "sv_temporal_bubble.ipynb")

md_title = """# Does the temporal SV-purging climate gradient survive counting events instead of alleles?

`temporal_s_consolidated.ipynb` reports two things about the per-site de-trended shift in `s`
(median `s` minus the same-$p_0$-decile SNP baseline, averaged over deciles):

1. the **net** shift is ~0 — no genome-wide SV purifying excess;
2. but it is **climate-graded** — negative (more purged) at hot/arid sites, positive at cold/wet
   ones, ρ ≈ −0.48 with bio1 and ≈ +0.65 with bio18, and **carried by insertions**.

Those are computed on biallelic panel records, i.e. **alleles**. arch3 is the biallelic
representation of a multi-allelic graph, so one insertion event appears once per distinct inserted
sequence. This notebook recomputes the same statistic on **bubbles** (events).

**`s` cannot be re-aggregated from per-allele `s`** — it is a logit slope, and the slope of a sum is
not the sum of slopes. So the trajectories are aggregated first and the estimator re-run:

$$\\text{AF}_\\text{bubble}(\\text{pool}) = \\sum_{a \\in \\text{bubble}} \\text{AF}_a , \\qquad
p_0(\\text{bubble}) = \\sum_{a \\in \\text{bubble}} p_{0,a}$$

then the same plot-replicate logit-slope estimator. Both units are computed in **one pass** over the
pool matrices, so they see identical data.

**The control.** The allele-level result is recomputed here too. If it does not land on the
published ρ, nothing else in this notebook means anything.

**The validity constraint.** Summing allele frequencies is only legitimate if no haplotype carries
two of the bubble's alleles. arch3 records are *atomic* variants nested in a bubble, and a single
path can carry several of them — so this is **not** automatic. It is tested directly on the founder
panel (does any founder carry more than one of the bubble's records?) and non-exclusive bubbles are
**dropped, not clipped**: their summed frequency is not a frequency. Without that filter, summed
$p_0$ exceeds 1 for 2.6% of bubbles, reaching 17.3."""

code_setup = r"""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
from scipy import stats
plt.rcParams.update({'figure.dpi':110,'font.size':8,'axes.linewidth':0.6})

G=os.path.abspath(os.path.join(os.getcwd(),"..","r1_sv_negative_selection","results","sv_adaptive"))
os.makedirs(f"{G}/plots",exist_ok=True)
def style(ax,axis="both"):
    ax.set_axisbelow(True); ax.grid(color="0.88",lw=0.7,axis=axis)
    for sp in ax.spines.values(): sp.set_visible(False)

D=pd.read_csv(f"{G}/sv_temporal_bubble_sites.csv")
S=pd.read_csv(f"{G}/sv_temporal_bubble_summary.csv").set_index("measure")
U=np.load(f"{G}/sv_temporal_bubble_units.npz",allow_pickle=True)
k=U["k"]; excl=U["excl"]; p0b=U["p0b"]; bcls=U["bcls"].astype(str); keep_b=U["keep_b"]; keep_ne=U["keep_noexcl"]
CU={"allele":"#D55E00","bubble":"#0072B2"}
print(f"{len(D)} sites; bio1 {D.bio1.min():.1f}-{D.bio1.max():.1f} C, bio18 {D.bio18.min():.0f}-{D.bio18.max():.0f} mm")
print(f"{k.size:,} bubbles from {int(k.sum()):,} SV alleles; exclusive {excl.sum():,} ({100*excl.mean():.1f}%)")
print(f"tested bubbles: {keep_b.sum():,} (exclusivity required) vs {keep_ne.sum():,} (not required)")
"""

md_control = """## 1. Control — does the allele-level recomputation reproduce the published result?

Same estimator, same 31 sites, same filters as `temporal_s_consolidated`. The published values are
median shift +0.0004, negative at 15/31 sites, ρ = −0.48 (bio1) and +0.65 (bio18)."""

code_control = r"""
pub={"median_shift":0.0004,"n_neg":15,"rho_bio1":-0.48,"rho_bio18":0.65}
r=S.loc["allele_sv"]
print(f"{'':<16}{'published':>12}{'recomputed':>13}{'delta':>10}")
for kk,lab in [("median_shift","median shift"),("n_neg","sites < 0"),("rho_bio1","rho bio1"),("rho_bio18","rho bio18")]:
    print(f"{lab:<16}{pub[kk]:>12.4f}{r[kk]:>13.4f}{r[kk]-pub[kk]:>+10.4f}")
print(f"\np-values: bio1 {r.p_bio1:.4f}, bio18 {r.p_bio18:.6f}  ({int(r.n_sites)} sites)")
print("reproduced" if abs(r.rho_bio1-pub['rho_bio1'])<0.02 and abs(r.rho_bio18-pub['rho_bio18'])<0.02
      else "DOES NOT REPRODUCE -- stop here")
"""

md_main = """## 2. The gradient, allele vs bubble

Each point is a site: x = its climate, y = its de-trended shift in `s` (below 0 = SVs more purged
than same-frequency SNPs). Left column is the allele unit, right is the bubble unit; top row is
temperature (bio1), bottom is precipitation of the warmest quarter (bio18, low = arid summer).

**Insertions first, then deletions.** The published gradient is attributed to insertions; deletions
are the structure control, and they do not behave like a null here — see the closing notes.
The two figures share a y-axis scale per row within a figure but not between figures."""

code_main = r"""

def draw_gradient(cls,lab,fname):
    fig,ax=plt.subplots(2,2,figsize=(10.5,7),sharey="row")
    for i,(cv,clab,cmap) in enumerate([("bio1","site mean annual temperature  bio1 (C)","coolwarm"),
                                       ("bio18","precipitation of warmest quarter  bio18 (mm, low = arid)","BrBG_r")]):
        for j,unit in enumerate(("allele","bubble")):
            a=ax[i,j]; style(a); col=f"{unit}_{cls}"
            d=D[[cv,col]].dropna()
            a.axhline(0,color="0.6",lw=0.8,ls=":")
            a.scatter(d[cv],d[col],c=d[cv],cmap=cmap,s=52,zorder=3,edgecolor="0.3",linewidth=0.4)
            b_,a0=np.polyfit(d[cv],d[col],1); xs=np.array([d[cv].min(),d[cv].max()])
            a.plot(xs,a0+b_*xs,color=CU[unit],lw=1.6,ls="--")
            r=stats.spearmanr(d[cv],d[col])
            a.annotate(f"{lab} {unit}s\n"+r"$\rho$"+f" = {r.statistic:+.3f}   p = {r.pvalue:.3g}\n"
                       f"{int((d[col]<0).sum())}/{len(d)} sites below 0",
                       xy=(0.03,0.97),xycoords="axes fraction",va="top",fontsize=8.5,
                       color=CU[unit],bbox=dict(facecolor="white",alpha=0.75,edgecolor="none",pad=1.5))
            a.set_xlabel(clab)
            if j==0: a.set_ylabel("de-trended shift in s\n(below 0 = more purged than same-freq SNPs)")
    fig.tight_layout(); fig.savefig(f"{G}/plots/{fname}",dpi=130,bbox_inches="tight"); plt.show()

draw_gradient("ins","insertion","sv_temporal_bubble_gradient.png")
draw_gradient("del","deletion","sv_temporal_bubble_gradient_del.png")
"""

md_grid = """## 3. Per-site grid — the same figure as `temporal_s_consolidated`, both units

The figure from the source notebook, rebuilt here twice: once with **SV = alleles** (which should
reproduce it) and once with **SV = bubbles**. Same estimator, same sliding-window median over raw
$(s, p_0)$ pairs, same SNP-only baseline, sites ordered cold → hot.

**The two grids share one y-axis**, computed over both, so the attenuation is visible rather than
hidden by autoscaling. Reading each panel: the orange line below 0 in the rare (left) part of the
axis is the purging signal; the source notebook's claim is that this deepens at the warm sites.

Indels stay at the allele unit in both grids — the bubble map only covers records with |dlen| > 50 bp,
so there is no bubble aggregation for them. They are the reference class, not a second test."""

code_grid = r"""

from matplotlib.lines import Line2D
P=np.load(f"{G}/sv_temporal_bubble_persite.npz")
COL={"SNP":"#888888","indel":"#0072B2","SV":"#D55E00"}
XTICKS=[0.03,0.05,0.1,0.2,0.5]; XLIM=(0.026,0.62)
meta=D.sort_values("bio1").reset_index(drop=True); sites=meta.site.tolist()

def rolling_curve(s,p0,thin=15):
    o=np.argsort(p0); p0s,ss=p0[o],s[o]
    win=max(51,len(ss)//8); win+=1-win%2
    med=pd.Series(ss).rolling(win,center=True,min_periods=win//3).median().to_numpy()
    return p0s[::thin],med[::thin]

def site_curves(site,svkey):
    p0b_,base=rolling_curve(P[f"{site}_s_SNP"],P[f"{site}_p0_SNP"]); lb=np.log(p0b_)
    out={}
    for c,key in (("SNP","SNP"),("indel","indel"),("SV",svkey)):
        p0c,medc=rolling_curve(P[f"{site}_s_{key}"],P[f"{site}_p0_{key}"])
        out[c]=(p0c,medc-np.interp(np.log(p0c),lb,base))
    return out

def _vis(p0c,rc):
    m=(p0c>=XLIM[0])&(p0c<=XLIM[1]); return p0c[m],rc[m]

CUR={u:{s_:site_curves(s_,key) for s_ in sites} for u,key in (("allele","SVallele"),("bubble","SV"))}
allv=np.concatenate([_vis(*CUR[u][s_][cl])[1] for u in CUR for s_ in sites for cl in ("SV","indel")])
pad=0.05*(allv.max()-allv.min()); YLIM=(allv.min()-pad,allv.max()+pad)

def draw_grid(unit,fname):
    ncol=7; nrow=int(np.ceil(len(sites)/ncol))
    norm=mpl.colors.Normalize(vmin=meta.bio1.min(),vmax=meta.bio1.max()); cmap=mpl.cm.coolwarm
    fig,axes=plt.subplots(nrow,ncol,figsize=(18,2.4*nrow),sharex=True,sharey=True)
    axes=axes.ravel()
    for ax in axes[len(sites):]: ax.axis("off")
    for i,s_ in enumerate(sites):
        ax=axes[i]; cv=CUR[unit][s_]
        p0sv,rsv=_vis(*cv["SV"]); p0in,rin=_vis(*cv["indel"]); p0sn,rsn=_vis(*cv["SNP"])
        ax.axhline(0,color="k",lw=0.7,ls=":")
        ax.fill_between(p0sv,0,rsv,color=COL["SV"],alpha=0.20,lw=0)
        ax.plot(p0sn,rsn,"-",color=COL["SNP"],lw=0.9)
        ax.plot(p0in,rin,"-",color=COL["indel"],lw=1.1)
        ax.plot(p0sv,rsv,"-",color=COL["SV"],lw=1.7)
        b1=meta.bio1.iloc[i]; c=cmap(norm(b1))
        ax.annotate(f"site {int(meta.site.iloc[i])}   bio1={b1:.0f}"+u"°C",xy=(0.03,0.96),
                    xycoords="axes fraction",va="top",ha="left",fontsize=6.8,color=c,fontweight="bold")
        for sp in ax.spines.values(): sp.set_color(c); sp.set_linewidth(1.4)
        ax.set_xscale("log"); ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
        ax.set_xticklabels([str(t) for t in XTICKS],fontsize=6); ax.tick_params(labelsize=6)
    axes[0].set_ylim(*YLIM)
    axes[min(len(sites),ncol)-1].annotate(f"SV = {unit}s",xy=(0.97,0.04),xycoords="axes fraction",
                                          ha="right",va="bottom",fontsize=9,color=COL["SV"],fontweight="bold")
    fig.supxlabel("initial frequency $p_0$ (log)",y=0.005,fontsize=10)
    fig.supylabel(f"s - same-$p_0$ SNP baseline   (below 0 = more purged)   [SV = {unit}s]",x=0.006,fontsize=10)
    fig.legend(handles=[Line2D([0],[0],color=COL[c],lw=2,label=c if c!="SV" else f"SV ({unit}s)")
                        for c in ("SNP","indel","SV")],
               loc="upper right",ncol=3,fontsize=9,frameon=False,bbox_to_anchor=(0.995,1.004))
    fig.tight_layout(rect=[0.02,0.02,1,0.99])
    fig.savefig(f"{G}/plots/{fname}",dpi=130,bbox_inches="tight"); plt.show()

draw_grid("allele","sv_temporal_bubble_grid_allele.png")
draw_grid("bubble","sv_temporal_bubble_grid_bubble.png")
print(f"shared ylim {YLIM[0]:+.3f} to {YLIM[1]:+.3f} (true min/max over both grids -- no clipping)")
for unit in ("allele","bubble"):
    dp=[np.trapz(np.clip(_vis(*CUR[unit][s_]["SV"])[1],None,0),np.log(_vis(*CUR[unit][s_]["SV"])[0])) for s_ in sites]
    r=stats.spearmanr(meta.bio1,dp)
    print(f"SV={unit:<7} area below 0 vs bio1: rho={r.statistic:+.3f} p={r.pvalue:.4f}")
"""

md_paired = """## 4. The same sites, both units

Each point is one site, x = its allele-level shift, y = its bubble-level shift, coloured by
temperature. On the dotted 1:1 line the two units agree. Points in the upper-left / lower-right
quadrants changed sign between units."""

code_paired = r"""
fig,ax=plt.subplots(1,2,figsize=(10.5,4.4))
for a,cls,lab in zip(ax,("ins","del"),("insertions","deletions")):
    style(a); d=D[["bio1",f"allele_{cls}",f"bubble_{cls}"]].dropna()
    lim=np.abs(np.r_[d[f"allele_{cls}"],d[f"bubble_{cls}"]]).max()*1.15
    a.axhline(0,color="0.8",lw=0.8); a.axvline(0,color="0.8",lw=0.8)
    a.plot([-lim,lim],[-lim,lim],color="0.6",lw=0.9,ls=":")
    sc=a.scatter(d[f"allele_{cls}"],d[f"bubble_{cls}"],c=d.bio1,cmap="coolwarm",s=55,zorder=3,
                 edgecolor="0.3",linewidth=0.4)
    a.set_xlim(-lim,lim); a.set_ylim(-lim,lim)
    flip=int(((d[f"allele_{cls}"]<0)!=(d[f"bubble_{cls}"]<0)).sum())
    a.annotate(f"{lab}\nsign changes at {flip}/{len(d)} sites\nslope(bubble~allele) = "
               f"{np.polyfit(d[f'allele_{cls}'],d[f'bubble_{cls}'],1)[0]:+.2f}",
               xy=(0.03,0.97),xycoords="axes fraction",va="top",fontsize=8.5,color="0.25",
               bbox=dict(facecolor="white",alpha=0.75,edgecolor="none",pad=1.5))
    a.set_xlabel("allele-level shift"); a.set_ylabel("bubble-level shift")
    plt.colorbar(sc,ax=a,label="bio1 (C)",fraction=0.046,pad=0.03)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_temporal_bubble_paired.png",dpi=130,bbox_inches="tight"); plt.show()
"""

md_all = """## 5. Every measure side by side

Spearman ρ of the per-site shift against each climate axis, for all measures computed. Filled
markers are p < 0.05. `*_ins_single` / `*_ins_multi` are the cheap non-independence diagnostic:
allele-level insertions split by whether their bubble holds one allele or several — if the gradient
lived only in the multi-allele stratum, allele-splitting would be the explanation."""

code_all = r"""
ORD=["allele_sv","bubble_sv","allele_ins","bubble_ins","allele_del","bubble_del",
     "allele_ins_k1","bubble_ins_k1","allele_ins_multi","bubble_ins_multi",
     "allele_del_k1","bubble_del_k1","allele_del_multi","bubble_del_multi"]
ORD=[m for m in ORD if m in S.index]
fig,ax=plt.subplots(1,2,figsize=(10.5,3.8),sharey=True)
yy=np.arange(len(ORD))[::-1]
for a,cv in zip(ax,("bio1","bio18")):
    style(a,axis="x"); a.axvline(0,color="0.6",lw=0.9)
    rho=S.loc[ORD,f"rho_{cv}"].to_numpy(); pv=S.loc[ORD,f"p_{cv}"].to_numpy()
    for y,r_,p_,m in zip(yy,rho,pv,ORD):
        c=CU["bubble"] if m.startswith("bubble") else CU["allele"]
        a.plot([0,r_],[y,y],color=c,lw=1.2,alpha=0.6)
        a.scatter([r_],[y],s=70,facecolor=c if p_<0.05 else "white",edgecolor=c,linewidth=1.4,zorder=3)
        a.annotate(f"{r_:+.2f}" + ("" if p_<0.05 else "  n.s."),(r_,y),xytext=(6 if r_>=0 else -6,0),
                   textcoords="offset points",ha="left" if r_>=0 else "right",va="center",fontsize=7,color="0.3")
    a.set_xlim(-1,1); a.set_xlabel(r"Spearman $\rho$ of per-site shift vs "+cv)
ax[0].set_yticks(yy); ax[0].set_yticklabels(ORD,fontsize=8)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_temporal_bubble_rho.png",dpi=130,bbox_inches="tight"); plt.show()

print(f"{'measure':<20}{'median shift':>14}{'<0 at':>9}{'rho bio1':>11}{'p':>9}{'rho bio18':>12}{'p':>9}")
for m in ORD:
    r=S.loc[m]
    print(f"{m:<20}{r.median_shift:>+14.4f}{f'{int(r.n_neg)}/{int(r.n_sites)}':>9}"
          f"{r.rho_bio1:>+11.3f}{r.p_bio1:>9.4f}{r.rho_bio18:>+12.3f}{r.p_bio18:>9.4f}")
"""

md_crux = """## 6. The crux — where the gradient lives

The pooled bubble row hides the result. Split insertion events by how many alleles their bubble
holds, and score each split **both ways**:

- **k = 1** (the allele *is* the bubble): the two units are the same events with the same
  trajectories, so their `s` is identical to floating-point — verified at 5.55e-17 over 25,264
  events. Any difference between these two bars would have to be the filter, not the unit.
- **k > 1**: the alleles and the summed event are genuinely different quantities. If the alleles
  carry the gradient but the summed event does not, the climate-correlated change is
  **redistribution among alleles at one locus** — one inserted sequence rising as another falls —
  rather than net gain or loss of the insertion.

That is the difference between selection on *the insertion* and selection on *which haplotype
carries it*."""

code_crux = r"""

def draw_crux(cls,lab,fname):
    CRUX=[("k = 1\n(allele = bubble)",f"allele_{cls}_k1",f"bubble_{cls}_k1"),
          ("k > 1\n(multi-allele)",f"allele_{cls}_multi",f"bubble_{cls}_multi")]
    fig,ax=plt.subplots(1,2,figsize=(10.5,4.0),sharey=True)
    for a,cv in zip(ax,("bio1","bio18")):
        style(a,axis="y"); a.axhline(0,color="0.5",lw=0.9)
        x=np.arange(len(CRUX)); w=0.34
        for off,unit,idx in ((-w/2,"allele",1),(w/2,"bubble",2)):
            vals=[S.loc[c[idx],f"rho_{cv}"] if c[idx] in S.index else np.nan for c in CRUX]
            pvs=[S.loc[c[idx],f"p_{cv}"] if c[idx] in S.index else np.nan for c in CRUX]
            a.bar(x+off,vals,width=w*0.92,color=CU[unit],linewidth=0,alpha=0.9,label=f"{unit}s")
            for xi,(v,pv) in enumerate(zip(vals,pvs)):
                t=f"{v:+.2f}\n"+("p<0.001" if pv<0.001 else f"p={pv:.3f}")
                a.annotate(t,(xi+off,v),xytext=(0,4 if v>=0 else -17),textcoords="offset points",
                           ha="center",fontsize=7.5,color="0.25")
        a.set_xticks(x); a.set_xticklabels([c[0] for c in CRUX],fontsize=8.5)
        a.set_ylim(-1,1)
        a.annotate(f"{lab}, vs {cv}",xy=(0.03,0.97),xycoords="axes fraction",va="top",
                   fontsize=8.5,color="0.3")
    ax[0].set_ylabel(r"Spearman $\rho$ of per-site shift vs climate")
    ax[0].legend(frameon=False,fontsize=8,loc="lower left")
    fig.tight_layout(); fig.savefig(f"{G}/plots/{fname}",dpi=130,bbox_inches="tight"); plt.show()
    print(f"{lab} events, gradient by allele count and unit")
    print(f"{'measure':<26}{'rho bio1':>11}{'p':>9}{'rho bio18':>12}{'p':>9}{'sites <0':>10}")
    for _,ca,cb_ in CRUX:
        for c in (ca,cb_):
            if c not in S.index: continue
            r=S.loc[c]
            print(f"{c:<26}{r.rho_bio1:>+11.3f}{r.p_bio1:>9.4f}{r.rho_bio18:>+12.3f}{r.p_bio18:>9.4f}"
                  f"{f'{int(r.n_neg)}/{int(r.n_sites)}':>10}")
    n1=int((keep_b&(bcls==cls)&(k==1)).sum()); nm=int((keep_b&(bcls==cls)&(k>1)).sum())
    print(f"  tested {lab} bubbles: {n1:,} with k=1, {nm:,} with k>1 ({100*nm/(n1+nm):.0f}% multi-allele)\n")

draw_crux("ins","insertions","sv_temporal_bubble_crux.png")
draw_crux("del","deletions","sv_temporal_bubble_crux_del.png")
"""

md_freq = """## 7. The frequency floor — are private SVs in this?

Both units require founder **MAC ≥ 12 of 231** (MAF ≳ 5%), call rate ≥ 90%, and
$p_0 \\in [0.02, 0.98]$ in the seed pool. So this is a common-variant analysis and **no private SV
survives**: zero singletons in the tested set, either unit.

The two filters are not equivalent, though. The allele side uses `lib.founder_panel_keep`, which is
**position-based** — it passes a record if *any* panel record at that position clears MAC ≥ 12. At a
multi-allelic position a rare allele therefore rides in on a common sibling's frequency. The bubble
side counts that event's own carriers, so its floor is exact.

This is worth knowing precisely because the rare, lineage-specific alleles are the ones most likely
to tag a haplotype rather than mark an event."""

code_freq = r"""

ka=U["keep_a"]; na=U["n_alt"]
fig,ax=plt.subplots(1,2,figsize=(10.5,3.8))
for a,(v,lab,c) in zip(ax,[(na[ka],"tested alleles (own founder MAC)","allele"),
                           (U["cb"][keep_b],"tested bubbles (own founder MAC)","bubble")]):
    style(a,axis="y")
    a.hist(v,bins=np.arange(0,233,4),color=CU[c],linewidth=0)
    a.axvline(12,color="#D55E00",lw=1.3,ls="--")
    below=int((v<12).sum())
    a.annotate(f"MAC >= 12 floor\n{below:,} of {v.size:,} below it ({100*below/v.size:.1f}%)\n"
               f"private (MAC = 1): {int((v==1).sum()):,}",
               xy=(0.35,0.95),xycoords="axes fraction",va="top",fontsize=8.5,color="0.3")
    a.set_xlabel("carriers among the 231 founders"); a.set_ylabel("count"); a.set_yscale("log")
    a.annotate(lab,xy=(0.35,0.55),xycoords="axes fraction",fontsize=8.5,color=CU[c])
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_temporal_bubble_macfloor.png",dpi=130,bbox_inches="tight"); plt.show()
print(f"tested alleles {int(ka.sum()):,}: own MAC < 12 for {int((na[ka]<12).sum()):,} "
      f"({100*(na[ka]<12).mean():.1f}%) -- these entered via the position-based filter; "
      f"private (MAC=1): {int((na[ka]==1).sum()):,}")
print(f"tested bubbles {int(keep_b.sum()):,}: own MAC < 12 for {int((U['cb'][keep_b]<12).sum()):,} (exact floor)")
"""

md_valid = """## 8. Validity — when is a bubble frequency even defined as a sum?

A bubble's alleles add only if no haplotype carries two of them. Tested on the founder panel: a
bubble is *exclusive* if no founder carries more than one of its records. Exclusivity is universal
at one allele and breaks down as the allele count grows — and the $p_0 \\ge 0.02$ floor on bubbles
selects preferentially for exactly the high-allele bubbles where it breaks."""

code_valid = r"""
KB=[(1,1),(2,2),(3,4),(5,9),(10,10**6)]
fig,ax=plt.subplots(1,3,figsize=(13,3.8))
labs=[f"{lo}" if lo==hi else (f"{lo}-{hi}" if hi<10**5 else f"{lo}+") for lo,hi in KB]
sel=[(k>=lo)&(k<=hi) for lo,hi in KB]
style(ax[0],axis="y")
ax[0].bar(range(len(KB)),[100*excl[s].mean() for s in sel],color="#0072B2",width=0.62,linewidth=0)
for i,s in enumerate(sel):
    ax[0].annotate(f"{100*excl[s].mean():.1f}%\nn={int(s.sum()):,}",(i,100*excl[s].mean()),
                   xytext=(0,4),textcoords="offset points",ha="center",fontsize=7.5,color="0.25")
ax[0].set_xticks(range(len(KB))); ax[0].set_xticklabels(labs); ax[0].set_ylim(0,118)
ax[0].set_xlabel("alleles in the bubble"); ax[0].set_ylabel("% of bubbles that are exclusive\n(summation valid)")
style(ax[1],axis="y")
ax[1].bar(range(len(KB)),[np.median(p0b[s]) for s in sel],color="0.55",width=0.62,linewidth=0)
ax[1].axhline(0.02,color="#D55E00",lw=1.2,ls="--")
ax[1].annotate("p0 >= 0.02 floor for a bubble to be tested\n-> the floor selects high-allele bubbles,\nwhich are the ones that are not exclusive",
               xy=(0.03,0.95),xycoords="axes fraction",va="top",fontsize=8,color="0.3")
ax[1].set_xticks(range(len(KB))); ax[1].set_xticklabels(labs); ax[1].set_yscale("log")
ax[1].set_xlabel("alleles in the bubble"); ax[1].set_ylabel("median $p_0$ of the bubble")
style(ax[2],axis="y")
tb=[int((keep_ne&(bcls==c)).sum()) for c in ("ins","del")]
te=[int((keep_b&(bcls==c)).sum()) for c in ("ins","del")]
x=np.arange(2)
ax[2].bar(x-0.17,tb,width=0.32,color="0.75",linewidth=0,label="passes the frequency/MAC filters")
ax[2].bar(x+0.17,te,width=0.32,color="#0072B2",linewidth=0,label="and is exclusive (kept)")
for xi,(u,v) in enumerate(zip(tb,te)):
    ax[2].annotate(f"{u:,}",(xi-0.17,u),xytext=(0,3),textcoords="offset points",ha="center",fontsize=7.5,color="0.3")
    ax[2].annotate(f"{v:,}\n({100*v/u:.0f}%)",(xi+0.17,v),xytext=(0,3),textcoords="offset points",ha="center",fontsize=7.5,color="0.3")
ax[2].set_xticks(x); ax[2].set_xticklabels(["insertion bubbles","deletion bubbles"])
ax[2].set_ylabel("tested bubbles"); ax[2].legend(frameon=False,fontsize=7.5,loc="upper right")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_temporal_bubble_validity.png",dpi=130,bbox_inches="tight"); plt.show()

print(f"summed p0 exceeds 1 in {int((p0b>1).sum()):,} of {p0b.size:,} bubbles ({100*(p0b>1).mean():.2f}%), max {p0b.max():.1f}")
print(f"dropping non-exclusive bubbles removes {keep_ne.sum()-keep_b.sum():,} of {keep_ne.sum():,} otherwise-testable bubbles")
"""

md_close = """## Reading it

**The control holds.** The allele-level recomputation lands on the published values exactly
(median +0.0004, negative at 15/31, ρ = −0.479 / +0.650), so both units are being compared on equal
footing, with the same estimator, sites and de-trending.

**The pooled bubble row says the gradient disappears — and that row is misleading.** Insertions go
from ρ = −0.388 / +0.731 as alleles to −0.036 / −0.038 as bubbles, both n.s. But 87% of the tested
insertion bubbles are multi-allele (3,557 of 4,088), so the pooled row is essentially the k > 1 row.

**Split by allele count, the picture inverts.** For **k = 1** insertions — where the allele *is* the
event, verified identical to 5.55e-17 — the gradient is fully present as an event: ρ = −0.421 (bio1)
and +0.703 (bio18), against −0.433 / +0.691 for the same events scored as alleles. For **k > 1**
insertions the alleles carry the gradient (−0.405 / +0.745) but the summed event does not
(+0.196, p = 0.29 / −0.306, p = 0.09).

**So the answer is not "the trend is an artefact of counting alleles".** It is:

- at **biallelic** insertion loci the climate gradient is a property of the **event** — the insertion
  itself rises or falls with climate;
- at **multi-allelic** loci the climate-correlated change is **redistribution among alleles** at one
  locus — one inserted sequence rising as another falls — which cancels in the sum.

That is the distinction between selection on the insertion and selection on *which haplotype carries
it*, and it is the concrete version of the lineage-tagging worry that motivated this recomputation.
It does **not** invalidate the published gradient; it localises it.

**How much weight the k = 1 result carries.** It rests on **531** tested insertion bubbles
against 3,557 multi-allele ones. That is enough for the per-site shift (the estimator takes
medians within $p_0$ deciles and averages over them) but it is a modest set, and the two strata
are not frequency-matched — a k = 1 event is rarer by construction, since showing k alleles
requires at least k carriers. Part of the k = 1 vs k > 1 contrast is therefore a frequency
contrast, and this design cannot separate the two.

**Deletions have their own figures (§2 and §5) and run the other way — not an aggregation artefact.** `allele_del_k1` (−0.559)
and `bubble_del_k1` (−0.587) agree as they must for identical events, so the negative bio18
correlation for deletions is a property of the data, not of summation — opposite in sign to
insertions (+0.703). Read directionally that is: at arid sites insertion ALTs are purged and
deletion ALTs are not.

**At k > 1 the deletions inverted relative to insertions, which is a caution flag.** For
insertions the alleles carried the gradient and the summed event lost it. For deletions the
opposite: the alleles show nothing (−0.169 / +0.059, both n.s.) while the summed event does
(+0.099 n.s. / −0.456, p = 0.010). Aggregation *creating* a correlation where the components
have none is not what the redistribution account predicts, and it has not been explained. Until
it is, the deletion k > 1 row should not be used.

**Treat the whole deletion contrast as an observation, not a conclusion** — see the projection caveat
below, which can produce exactly this anti-correlation without any SV-specific selection.

**Private SVs are not in this.** Zero singletons survive in either unit. But 12.8% of the tested
alleles (633 of 4,946) sit below the nominal MAC ≥ 12 floor, having passed the position-based
allele filter on a sibling allele's frequency; the bubble floor is exact. Since rare alleles at
multi-allelic positions are exactly the lineage-tagging ones, this leak works in the direction
of the allele-level signal, and is a second reason the k > 1 allele row should not be read as an
event-level result.

**Validity.** Summing allele frequencies requires that no haplotype carries two of a bubble's
alleles. arch3 records are atomic variants nested in a bubble, so this is not automatic: it holds for
93.1% of bubbles, 100% at k = 1 but only 55.9% at k ≥ 10, and summed $p_0$ exceeds 1 for 2.6% of
bubbles (max 17.3). Non-exclusive bubbles are dropped. In practice the MAC and call-rate filters had
already removed most of them — the requirement costs only 81 of 6,607 otherwise-testable bubbles —
so it changes no conclusion here, but it would matter for any analysis with a lower frequency floor.

Caveats. Pool AF is a global-mode founder projection, so `s` is a projection of founder-haplotype
slopes and a SNP on the same founders shares it: SV-specific selection cannot be separated from
riding a selected haplotype, and this applies to the insertion/deletion sign contrast above. Bubbles
sit at higher $p_0$ than their alleles by construction, so the bubble set is not a subset of the
allele set (6,526 bubbles vs 4,946 alleles) and the two units do not cover the same frequency range
even after decile de-trending. The de-trending baseline is the SNP median per $p_0$ decile, so every
number is relative to SNPs, not absolute."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_setup),
    new_markdown_cell(md_control), new_code_cell(code_control),
    new_markdown_cell(md_main), new_code_cell(code_main),
    new_markdown_cell(md_grid), new_code_cell(code_grid),
    new_markdown_cell(md_paired), new_code_cell(code_paired),
    new_markdown_cell(md_all), new_code_cell(code_all),
    new_markdown_cell(md_crux), new_code_cell(code_crux),
    new_markdown_cell(md_freq), new_code_cell(code_freq),
    new_markdown_cell(md_valid), new_code_cell(code_valid),
    new_markdown_cell(md_close),
])
ep = ExecutePreprocessor(timeout=1800, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
