#!/usr/bin/env python
"""Build+execute `sv_climate_maf.ipynb`: the garden-climate split of the matched frequency change
(the context-control summary figure, `supp_fig2c_summary.png` in `sv_main_fig2.ipynb`) as the lead
figure, then the same figure restricted to MAF >= 5% in the founders, to ask whether the warm-garden
excess is carried by common variants or by rare ones, then against generation-1 plot parallelism as a
direct measure of selection strength. Shared setup/data code in `_main_fig_common.py`.
`basic` env."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_fig_common import code_setup, code_grenenet, build, new_markdown_cell, new_code_cell

md_title = """# Is the insertion signal selection pressure? Garden climate, MAF, and parallelism

The lead figure is the context control from `sv_main_fig2.ipynb`, split by garden climate: for each
event, its allele-frequency change from the seed mix to the garden minus the mean change of its
matched SNPs, averaged over the class and expressed as a percent of where the class started, with the
31 gardens grouped into terciles of garden temperature (cold / mid / hot).

The question here is whether the warm-garden excess is a **selection-pressure** effect: if the
CDS-insertion decline is selection acting on coding insertions, it should be strongest where selection
acted hardest, and it should survive restricting to variants common enough to be estimated well
(MAF ≥ 5% in the founders) — the signal would simply be measured with less noise.

Two proxies for "where selection acted hardest" are used, and they disagree. Garden temperature is the
one this figure splits on. Generation-1 **parallelism** (last section) measures it directly from how
much the plots within a garden agree, and it turns out to be uncorrelated with garden temperature
(*r* = +0.18, p = 0.34) — so "selection was strongest in the warm gardens" does not hold for that
measure, and the two sections are testing genuinely different things.

The MAF ≥ 5% version is a diagnostic, not the preferred analysis: most SVs sit at low frequency, so
the filter discards most of the data, and the matching already pairs every event with SNPs carried by
the same number of founders. The full set is what the parallelism section uses.

**Unit and filters.** SV **event** (arch3 alleles of a pangenome bubble merged with `truvari collapse`
defaults); called in ≥ 90% of the 231 founders; ≥ 2 carrier founders. Context, first match wins:
CDS > UTR > intron > reference TE > promoter (≤ 2 kb upstream) > intergenic. Frequencies are kMate's
per-pool values (`garden_af_event.npz`): the founding frequency is the mean of the 8 seed-mix runs and
a garden's frequency is the mean over its plots at the plot's last sampled generation (352 pools)."""

md_fig = """## Lead figure — matched frequency change by garden climate, all variants

For each event *v* and garden *s*, Δ<sub>v</sub> = *p*<sub>v</sub>(garden) − *p*<sub>v</sub>(seed mix)
minus the mean Δ of its matched SNPs: all SNPs of the same context with the **same carrier count** and
the **same carrier-origin-climate quintile**, so SNPs that started at the same frequency on the same
kind of founder background (events whose exact cell holds no SNP fall back to their carrier-count bin).
The differences are summed over the class and divided by the class's summed starting frequency, so the
value is a percent of where the class began — equivalently, the relative change of the class's summed
frequency minus that of the composition-matched SNP set (`load_diff`).

Small dots = gardens; large points = the tercile mean with a 95% percentile bootstrap CI over that
tercile's gardens (2,000 resamples; cold *n* = 11, mid *n* = 10, hot *n* = 10). The annotation gives the
all-garden mean, the number of gardens below zero, and the empirical matched-SNP-set null *p* (1,000
random SNP sets with the class's composition, one SNP per event drawn from that event's matches,
two-sided). Gardens share one founder pool and are not independent, so garden counts are descriptive."""

code_zone = r"""
ZC={"cold":"#2166AC","mid":"0.55","hot":"#B2182B"}; ZOFF={"cold":-0.25,"mid":0.0,"hot":0.25}
MAC=12; COMMON=(NCAR>=MAC)&(NCAR<=231-MAC); ALLV=np.ones(NCAR.size,bool)      # MAF >= 5% among the 231 founders
SUMC=["CDS","intron","intergenic"]
def bci(v,B=2000,seed=2):
    # percentile bootstrap CI of the mean over gardens (resample gardens with replacement)
    r=np.random.default_rng(seed); bs=r.choice(v,(B,v.size)).mean(1); return np.percentile(bs,2.5),np.percentile(bs,97.5)
def zcontrast(a,b,B=2000,seed=5):
    # hot minus cold: resample each tercile's gardens independently
    r=np.random.default_rng(seed); bs=r.choice(a,(B,a.size)).mean(1)-r.choice(b,(B,b.size)).mean(1)
    return bs.mean(),np.percentile(bs,2.5),np.percentile(bs,97.5)
def zone_figure(flt,head,fname):
    fig,ax=plt.subplots(1,2,figsize=(13,3.6),sharey=True); rows=[]
    for a,k in zip(ax,("ins","del")):
        style(a,axis="x"); a.axvline(0,color="0.35",lw=1,ls="--")
        for yi,ctx in enumerate(SUMC):
            ci=CTX.index(ctx); ms=(KD==k)&(CX==ci)&flt; mn=(KD=="snp")&(CX==ci)&flt
            d=load_diff(ms,mn); pe=null_p(d.mean(),matched_null(ms,mn),two_sided=True)
            rec=dict(kind=KL[k],context=ctx,n_events=int(ms.sum()),n_snps=int(mn.sum()),mean=d.mean(),
                     gardens_below=int((d<0).sum()),null_p_two_sided=pe)
            for z in ("cold","mid","hot"):
                dz=d[ZN==z]; a.scatter(dz,yi+ZOFF[z]+rng.uniform(-0.06,0.06,dz.size),s=11,color=ZC[z],alpha=0.45,linewidth=0,zorder=2)
                lo,hi=bci(dz); a.errorbar(dz.mean(),yi+ZOFF[z],xerr=[[dz.mean()-lo],[hi-dz.mean()]],fmt="o",ms=5.5,color=ZC[z],mec="0.15",mew=0.6,capsize=2,elinewidth=1.2,zorder=3)
                rec[z]=dz.mean(); rec[f"{z}_lo"]=lo; rec[f"{z}_hi"]=hi
            lo,hi=bci(d); rec["boot_lo"]=lo; rec["boot_hi"]=hi
            hm,hl,hh=zcontrast(d[ZN=="hot"],d[ZN=="cold"]); rec["hot_minus_cold"]=hm; rec["hmc_lo"]=hl; rec["hmc_hi"]=hh
            rows.append(rec)
            a.annotate(f"{100*d.mean():+.1f}%, {(d<0).sum()}/31 gardens < 0, null {pfmt(pe)}",xy=(0.01,yi-0.44),xycoords=("axes fraction","data"),ha="left",va="top",fontsize=7,color="0.3",bbox=dict(facecolor="white",alpha=0.85,edgecolor="none",pad=1))
            if yi<2: a.axhline(yi+0.5,color="0.9",lw=0.8)
        a.set_yticks(range(3)); a.set_yticklabels(SUMC,fontsize=8); a.set_ylim(2.6,-0.6)
        a.set_xlabel(f"frequency change of {KL[k]} minus matched same-context SNPs, per garden (% of start)",fontsize=8)
        a.xaxis.set_major_formatter(PercentFormatter(1.0,decimals=0))
        a.annotate(f"{head} {KL[k]}  (n = "+", ".join(f"{c} {int(((KD==k)&(CX==CTX.index(c))&flt).sum()):,}" for c in SUMC)+")",xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5,color=KC[k])
    for z,l in (("cold","cold gardens"),("mid","mid"),("hot","hot")): ax[0].errorbar([],[],xerr=[],fmt="o",ms=5,color=ZC[z],label=l)
    ax[0].legend(frameon=False,fontsize=7.5,loc="lower right")
    fig.tight_layout(); fig.savefig(f"{G}/plots/{fname}",dpi=200); plt.show()
    return pd.DataFrame(rows)
def show_zone(T):
    C=["kind","context","n_events","mean","boot_lo","boot_hi","gardens_below","null_p_two_sided","cold","mid","hot","hot_minus_cold","hmc_lo","hmc_hi"]
    print(T[C].round(4).to_string(index=False))
ZALL=zone_figure(ALLV,"all","climate_maf_all_variants.png"); show_zone(ZALL)
"""

md_maf = """## The same figure, MAF ≥ 5% in the founders

Identical statistic, matching and test; the only change is that both the events and the SNPs they are
matched against are restricted to a founder minor-allele count between 12 and 219 of 231 (MAF ≥ 5%),
the same `COMMON` filter used for the MAF ≥ 5% versions of the lead analysis in `sv_main_fig2.ipynb`.
The filter is on the **founder panel**, not on the seed mix, so it is a property of the variant rather
than of any garden.

This costs most of the events — CDS insertions drop from 715 to 201 — so wider intervals are expected.
The question is whether the *pattern* survives, in particular the hot-garden offset in CDS insertions."""

code_maf = r"""
ZM=zone_figure(COMMON,"MAF ≥ 5%","climate_maf_common.png"); show_zone(ZM)
"""

md_cmp = """## Side by side — does the warm-garden effect survive the filter?

The hypothesis under test: selection was strongest in the warm gardens, so a genuine selection effect
on coding insertions should persist (or sharpen) once noisy rare variants are removed. The table
reports each tercile mean and the hot-minus-cold contrast (bootstrap over each tercile's gardens,
2,000 resamples) for both variant sets."""

code_cmp = r"""
A=ZALL.assign(filt="all variants"); B=ZM.assign(filt="MAF >= 5%")
K=pd.concat([A,B]).set_index(["kind","context","filt"]).sort_index()
print("tercile means and the hot - cold contrast, percent of starting frequency:\n")
for k in ("insertions","deletions"):
    print(f"--- {k} ---")
    for c in SUMC:
        a=K.loc[(k,c,"all variants")]; b=K.loc[(k,c,"MAF >= 5%")]
        print(f"  {c:<11} n {int(a.n_events):>5} -> {int(b.n_events):>4}   "
              f"all {100*a['mean']:+.2f}% [{100*a.boot_lo:+.2f},{100*a.boot_hi:+.2f}] -> {100*b['mean']:+.2f}% [{100*b.boot_lo:+.2f},{100*b.boot_hi:+.2f}]")
        print(f"  {'':<11} cold {100*a.cold:+.2f} -> {100*b.cold:+.2f}   mid {100*a.mid:+.2f} -> {100*b.mid:+.2f}   hot {100*a.hot:+.2f} -> {100*b.hot:+.2f}")
        print(f"  {'':<11} hot - cold {100*a.hot_minus_cold:+.2f}% [{100*a.hmc_lo:+.2f},{100*a.hmc_hi:+.2f}] -> {100*b.hot_minus_cold:+.2f}% [{100*b.hmc_lo:+.2f},{100*b.hmc_hi:+.2f}]")
    print()
print("full table:"); print(pd.concat([A,B]).round(4).to_string(index=False))
"""

md_read = """## What the comparison says

**The warm-garden effect in CDS insertions does not survive the filter.** All variants: hot gardens
−1.7%, hot − cold = −1.6% [−2.8, −0.3], an interval that excludes zero. MAF ≥ 5%: hot gardens −0.4%,
hot − cold = −0.2% [−1.7, +1.5]. The all-garden CDS-insertion mean also moves from −0.7% [−1.3, −0.2]
(excludes zero) to −0.5% [−1.1, +0.1] (includes zero). Some of that widening is the loss of power
(715 → 201 events), but the tercile *pattern* flattens as well — the hot-cold ordering is what
disappears, not just its precision. A selection-pressure reading predicts the opposite: common variants
are estimated best, so removing the rare tail should leave the climate ordering standing.

**The one climate contrast that does survive is intergenic.** Intergenic insertions: hot − cold =
−1.4% [−2.5, −0.4] for all variants and −1.5% [−2.8, −0.2] at MAF ≥ 5%, with cold gardens *positive*
(+0.7% → +1.0%) and hot gardens negative. Intergenic sequence is the least functionally constrained
context here, so a climate contrast that is strongest there and absent in CDS once rare variants are
removed is hard to read as selection pressure on coding insertions. It fits a genome-wide,
climate-dependent sorting of founder lineages better than a context-specific functional effect — the
same reading Figure 2c attaches to a slope on temperature.

**CDS deletions move the other way and strengthen.** +0.7% [+0.1, +1.2] → +1.0% [+0.3, +1.5] at
MAF ≥ 5%, most positive in cold gardens (+1.3%). Whatever is happening to coding deletions is not the
mirror of the insertion story and is not explained by warm-garden selection pressure.

*How to weigh this.* The filter is a noise diagnostic, not a better analysis — it removes ~70% of the
events, and the carrier-count matching already compares each insertion with SNPs that started at the
same frequency, so the rare events are not unmatched noise. The informative part is not that the CDS
interval widens but that the hot-cold *ordering* flattens, which is what the parallelism test below
probes directly.

*Caveats.* Six context × kind cells are tested at each filter with no multiple-testing correction;
terciles hold 10–11 gardens each, so tercile means are noisy; and the MAF ≥ 5% set is not a random
subsample of the full set — common insertions differ from rare ones in age, length and repeat content,
so the filter changes the biology being measured as well as the noise."""

md_par = """## A direct measure of selection strength — generation-1 parallelism

The climate tercile is only a proxy for how hard selection acted. `generation_1_parallelism.txt`
(`/global/scratch/users/tbellg/gea_grene-net/key_files/`) measures it directly: per site, how much the
plots within that garden moved in the *same* direction between founding and generation 1. Under strong
selection the replicate plots agree; under drift they do not. Two versions are given, computed from SNP
allele frequencies and from ecotype frequencies; they are near-identical here (*r* = 0.98), so both are
shown. 30 of the 31 gardens have a value (site 33 is missing).

**Test.** Across gardens, Pearson and Spearman correlation of the garden's parallelism with the
garden's matched frequency change for each kind × context (`load_diff`, all variants — the carrier-count
matching already pairs each event with SNPs that started at the same frequency, so the rare events are
kept). Because garden temperature has its own slope on this statistic, each correlation is repeated as
a **partial correlation controlling for garden bio1** (both variables linearly de-trended on bio1
first). A negative correlation means the more strongly selected gardens show more decline of the class
relative to its matched SNPs.

**The specificity question is the point of the test.** If parallelism predicts the insertion decline in
CDS but not in intergenic sequence, that is selection acting on coding insertions. If it predicts the
decline everywhere, or only outside genes, it is genome-wide sorting that happens to drag insertions
with it."""
code_par = r"""
PAR=pd.read_csv("/global/scratch/users/tbellg/gea_grene-net/key_files/generation_1_parallelism.txt",sep="\t")
PW=PAR.pivot(index="site",columns="source",values="mean")
def partial(x,y,z):
    # partial correlation of x and y controlling for z (linear de-trend on z, then Pearson)
    rx=x-np.polyval(np.polyfit(z,x,1),z); ry=y-np.polyval(np.polyfit(z,y,1),z); return stats.pearsonr(rx,ry)
DP=pd.DataFrame(dict(site=GS["sites"],bio1=BIO,zone=ZN)).join(PW,on="site")
for k_ in ("ins","del"):
    for c_ in SUMC: DP[f"{k_}_{c_}"]=load_diff((KD==k_)&(CX==CTX.index(c_)),(KD=="snp")&(CX==CTX.index(c_)))
DP=DP.dropna(subset=["snp","ecotype"]).reset_index(drop=True)
fig,axs=plt.subplots(2,3,figsize=(13,6.6)); PROWS=[]
for ri,k_ in enumerate(("ins","del")):
    for cj,c_ in enumerate(SUMC):
        a=axs[ri,cj]; style(a); a.axhline(0,color="0.35",lw=1,ls="--")
        x=DP["snp"].to_numpy(); y=DP[f"{k_}_{c_}"].to_numpy()
        a.scatter(x,y,c=DP.bio1,cmap="coolwarm",vmin=BIO.min(),vmax=BIO.max(),s=26,linewidth=0.3,edgecolor="0.3",zorder=3)
        lr=stats.linregress(x,y); xx=np.linspace(x.min(),x.max(),50)
        a.plot(xx,lr.intercept+lr.slope*xx,color=KC[k_],lw=1.6,zorder=2)
        for s_ in ("snp","ecotype"):
            xs=DP[s_].to_numpy(); r_,p_=stats.pearsonr(xs,y); rs_,ps_=stats.spearmanr(xs,y); pr_,pp_=partial(xs,y,DP.bio1.to_numpy())
            PROWS.append(dict(parallelism=s_,kind=KL[k_],context=c_,pearson_r=r_,p=p_,spearman_rho=rs_,p_spearman=ps_,
                              partial_r_given_bio1=pr_,p_partial=pp_))
        r_,p_=stats.pearsonr(x,y); pr_,pp_=partial(x,y,DP.bio1.to_numpy())
        a.annotate(f"r = {r_:+.2f} ({pfmt(p_)})\npartial r | bio1 = {pr_:+.2f} ({pfmt(pp_)})",
                   xy=(0.03,0.04),xycoords="axes fraction",fontsize=6.8,color="0.3")
        a.annotate(f"{KL[k_]}, {c_}",xy=(0.03,0.93),xycoords="axes fraction",fontsize=8,color=KC[k_])
        a.yaxis.set_major_formatter(PercentFormatter(1.0,decimals=0))
        if ri==1: a.set_xlabel("generation-1 parallelism of the garden (SNP-based)",fontsize=7.5)
        if cj==0: a.set_ylabel("frequency change minus matched\nsame-context SNPs (% of start)",fontsize=7.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/parallelism_vs_matched_change.png",dpi=200,bbox_inches="tight")
fig.savefig(f"{G}/plots/parallelism_vs_matched_change.pdf",bbox_inches="tight"); plt.show()
# is parallelism just another way of saying "warm garden"?
fig2,axs2=plt.subplots(1,2,figsize=(8.4,3.1))
for a,s_ in zip(axs2,("snp","ecotype")):
    style(a); a.scatter(DP.bio1,DP[s_],c=DP.bio1,cmap="coolwarm",vmin=BIO.min(),vmax=BIO.max(),s=26,linewidth=0.3,edgecolor="0.3",zorder=3)
    r_,p_=stats.pearsonr(DP.bio1,DP[s_]); rs_,ps_=stats.spearmanr(DP.bio1,DP[s_])
    a.annotate(f"{s_}-based parallelism\nr = {r_:+.2f} ({pfmt(p_)}),  ρ = {rs_:+.2f} ({pfmt(ps_)})",xy=(0.03,0.88),xycoords="axes fraction",fontsize=7,color="0.3")
    a.set_xlabel("garden temperature  bio1 (°C)",fontsize=7.5); a.set_ylabel("generation-1 parallelism",fontsize=7.5)
fig2.tight_layout(); fig2.savefig(f"{G}/plots/parallelism_vs_bio1.png",dpi=200,bbox_inches="tight"); plt.show()
PT=pd.DataFrame(PROWS)
print("parallelism vs the matched frequency change per garden (n = %d gardens)"%len(DP))
print("negative r = more strongly selected gardens lose more of the class relative to its matched SNPs\n")
print(PT.round(4).to_string(index=False))
r_,p_=stats.pearsonr(DP.snp,DP.ecotype)
print(f"\nSNP-based vs ecotype-based parallelism: r = {r_:+.3f} ({pfmt(p_)})")
for s_ in ("snp","ecotype"):
    r_,p_=stats.pearsonr(DP.bio1,DP[s_]); print(f"{s_}-based parallelism vs garden bio1: r = {r_:+.3f} ({pfmt(p_)})  -- parallelism is not a restatement of garden temperature")
rg=np.random.default_rng(5); hi=DP[DP.snp>=DP.snp.median()]; lo=DP[DP.snp<DP.snp.median()]
print(f"\nmedian split on SNP parallelism (high n = {len(hi)}, low n = {len(lo)}), mean matched change:")
for k_ in ("ins","del"):
    for c_ in SUMC:
        a_=hi[f"{k_}_{c_}"].to_numpy(); b_=lo[f"{k_}_{c_}"].to_numpy()
        bs=rg.choice(a_,(2000,a_.size)).mean(1)-rg.choice(b_,(2000,b_.size)).mean(1)
        print(f"  {KL[k_]:<11} {c_:<11} high {100*a_.mean():+.2f}%   low {100*b_.mean():+.2f}%   difference {100*bs.mean():+.2f}% [{100*np.percentile(bs,2.5):+.2f}, {100*np.percentile(bs,97.5):+.2f}]")
print("\nper garden, sorted by parallelism:")
print(DP.sort_values("snp")[["site","bio1","zone","snp","ecotype","ins_CDS","ins_intron","ins_intergenic"]].round(4).to_string(index=False))
"""

md_par_read = """## What parallelism says

**Parallelism is not a restatement of garden temperature.** SNP-based parallelism vs garden bio1:
*r* = +0.18 (p = 0.34); ecotype-based *r* = +0.23 (p = 0.22). The two hottest-selected gardens (sites
32 and 4, parallelism 0.78 and 0.70) are warm, but across all 30 gardens the two axes are independent.
That makes parallelism a genuinely new test rather than the climate split in different clothing — and
it means "selection was strongest in the warm gardens" does not hold for this measure of selection
strength.

**CDS insertions show no relationship with it whatsoever.** *r* = +0.00 (p = 1.0) for SNP-based
parallelism, −0.08 (p = 0.69) for ecotype-based; partial correlations controlling for bio1 are +0.08
and +0.02. The median split is flat: high-parallelism gardens −0.78%, low-parallelism gardens −0.68%,
difference −0.10% [−1.17, +1.01]. If the coding-insertion decline were purifying selection, the gardens
where selection demonstrably acted hardest should show more of it. They do not.

**Intergenic insertions do, strongly.** *r* = −0.53 (p = 0.0023) SNP-based and −0.53 (p = 0.0027)
ecotype-based, and it survives controlling for climate (partial *r* = −0.51, p = 0.0038). Median split:
high-parallelism gardens −0.49% vs low-parallelism +0.76%, difference −1.27% [−2.22, −0.37]. This is the
same context that kept its climate contrast under the MAF ≥ 5% filter. So the one insertion signal that
tracks selection strength is the one *outside* genes.

**Deletions move the opposite way in the same gardens.** Intron deletions *r* = +0.45 (p = 0.012) and
intergenic deletions *r* = +0.46 (p = 0.011), both robust to controlling for bio1 (partial *r* = +0.46,
+0.47). In the strongly selected gardens, non-coding deletions rise relative to matched SNPs while
non-coding insertions fall.

**Reading.** Insertions falling and deletions rising together, outside genes, in proportion to how
hard selection acted in that garden, is the signature of **linked selection**, not of selection on the
variants themselves: intergenic insertions are disproportionately lineage-specific and sit in
repeat-rich, low-recombination regions, so when strong selection sorts founder lineages hard they get
dragged along. The carrier-count and carrier-climate matching controls the SNPs' starting frequency and
the climate of their carriers, but it does not control linkage or which lineage a variant marks.

*Caveats.* 30 gardens, 12 correlations reported with no multiple-testing correction; parallelism is
measured at generation 1 while the frequency change runs to the last sampled generation; and the
parallelism estimate itself has a CI per site (the `lower`/`upper` columns) that this test ignores by
using the point estimate."""

md_burden = """## The same question at the founder level, where kinship can be controlled

Everything above is variant-level, and variant-level tests have a structural ceiling here. From
`DRIVE_AF_README.md`, kMate estimates the founder mixture *h* per pool and then obtains every record's
frequency by the projection

> *p*<sub>v</sub> = (*h* · var_pa<sub>v</sub>) / (*h* · var_called<sub>v</sub>)

so a variant's frequency in a pool is a deterministic function of *h* and of **which founders carry
it**. Two variants with the same carrier and called sets have identical trajectories by construction;
there is no per-variant information beyond its carrier pattern. The matching used above equalises
carrier *count* and the *mean* carrier origin-climate in five quintiles, but a mean over a 231-dimensional
carrier vector is a coarse summary — two variants can share it and still be carried by disjoint sets of
founders. So a variant-level contrast cannot separate "selection on the insertion" from "selection on
the lineages that happen to carry insertions".

The founder-level burden model can, because it puts relatedness in the model explicitly. Per garden
*s*, with *y*<sub>f</sub> = log[(*H*<sub>0,f</sub> + Δ*H*<sub>s,f</sub> + 10⁻⁵)/(*H*<sub>0,f</sub> + 10⁻⁵)]:

> *y* = β₀ + β₁·z(shared CDS insertion rate) + β₂·z(shared intergenic insertion rate) +
> β₃·z(call rate, SVs) + β₄·z(call rate, SNPs) + β₅·[sequencing route] + *u* + ε,
> *u* ~ N(0, σ²<sub>g</sub>**K**), ε ~ N(0, σ²<sub>e</sub>**I**)

β₁ is therefore already a CDS-vs-intergenic contrast with kinship controlled. These fits exist
(`_sv_burden_lmm.py` → `sv_burden_lmm_by_garden_event.csv`); this section only analyses them across
gardens.

**Why the per-garden β₁ need meta-analysis rather than a correlation.** Each β₁ carries its own SE
(0.085–0.168), and those SEs are *large* relative to the spread of the β₁ themselves (SD 0.126, and the
SD expected from sampling noise alone is 0.138). Correlating β₁ against temperature treats that noise as
real between-garden variation. The proper tools are an inverse-variance pooled effect, Cochran's *Q* for
whether the gardens differ by more than their own uncertainty, and a meta-regression of β₁ on garden
bio1 — reported here at both inverse-variance and equal weights, with leave-one-garden-out and a
permutation test on the slope."""
code_burden = r"""
BK=pd.read_csv(f"{G}/sv_burden_lmm_by_garden_kinds_event.csv")     # _sv_burden_by_garden_kinds.py
ZCB={"cold":"#2166AC","mid":"0.55","hot":"#B2182B"}
def metareg(b_,v_,x_,t2=0.0):
    ww=1/(v_+t2); X=np.column_stack([np.ones(len(b_)),x_]); XtW=X.T*ww
    bh=np.linalg.solve(XtW@X,XtW@b_); cv=np.linalg.inv(XtW@X); return bh,np.sqrt(np.diag(cv))
def burden_meta(kind):
    g=BK[BK.kind==kind].reset_index(drop=True)
    b=g.beta.to_numpy(); se=g.se.to_numpy(); x=g.bio1.to_numpy(); v=se**2; w=1/v
    pooled=(w*b).sum()/w.sum(); pse=np.sqrt(1/w.sum())
    Q=(w*(b-pooled)**2).sum(); dfQ=len(b)-1; pQ=stats.chi2.sf(Q,dfQ); I2=max(0.,(Q-dfQ)/Q)*100
    return dict(g=g,b=b,se=se,x=x,v=v,pooled=pooled,pse=pse,Q=Q,dfQ=dfQ,pQ=pQ,I2=I2)
def burden_figure(kind,fname):
    M=burden_meta(kind); g,b,se,x,v=M["g"],M["b"],M["se"],M["x"],M["v"]
    bh,bhse=metareg(b,v,x); su,suse=metareg(b,np.full_like(v,v.mean()),x)
    fig=plt.figure(figsize=(12.5,4.6)); gsb=fig.add_gridspec(1,2,width_ratios=[1.35,1],wspace=0.30)
    a=fig.add_subplot(gsb[0,0]); style(a,axis="x"); panel(a,"a"); o=np.argsort(x)
    a.axvline(0,color="0.35",lw=1,ls="--")
    a.axvspan(M["pooled"]-1.96*M["pse"],M["pooled"]+1.96*M["pse"],color=KC[kind],alpha=0.13,lw=0)
    a.axvline(M["pooled"],color=KC[kind],lw=1.4)
    for yi,i in enumerate(o):
        zc_=ZCB[g.zone[i]]
        a.errorbar(b[i],yi,xerr=1.96*se[i],fmt="o",ms=4,color=zc_,ecolor=zc_,elinewidth=1,capsize=1.8,alpha=0.95)
    a.set_yticks(range(len(x))); a.set_yticklabels([f"{int(g.site[i])}" for i in o],fontsize=5.6)
    a.set_ylim(-1,len(x)); a.set_ylabel("garden (site), coldest at the bottom",fontsize=7.5)
    a.set_xlabel(f"\u03b2\u2081: effect of CDS-{KL[kind][:-1]} burden on founder log-fitness\n(per SD, kinship controlled, \u00b1 95% CI)",fontsize=7.5)
    a.annotate(f"pooled \u03b2\u2081 = {M['pooled']:+.3f} \u00b1 {M['pse']:.3f}  (p = {2*stats.norm.sf(abs(M['pooled']/M['pse'])):.1g})\n"
               f"negative in {(b<0).sum()}/{len(b)} gardens\n"
               f"Q = {M['Q']:.1f} on {M['dfQ']} df (p = {M['pQ']:.2f}), I\u00b2 = {M['I2']:.0f}%",
               xy=(0.02,0.015),xycoords="axes fraction",fontsize=6.8,color="0.3",
               bbox=dict(facecolor="white",alpha=0.85,edgecolor="none",pad=1.5))
    for z_,l_ in (("cold","cold gardens"),("mid","mid"),("hot","hot")): a.errorbar([],[],fmt="o",ms=4,color=ZCB[z_],label=l_)
    a.legend(frameon=False,fontsize=6.8,loc="upper left")
    a=fig.add_subplot(gsb[0,1]); style(a); panel(a,"b"); a.axhline(0,color="0.35",lw=1,ls="--")
    a.axhline(M["pooled"],color=KC[kind],lw=1.2,alpha=0.7)
    a.errorbar(x,b,yerr=1.96*se,fmt="none",ecolor="0.75",elinewidth=0.8,zorder=1)
    a.scatter(x,b,s=1400/(se*1000),c=x,cmap="coolwarm",vmin=BIO.min(),vmax=BIO.max(),linewidth=0.3,edgecolor="0.3",zorder=3)
    xx=np.linspace(x.min(),x.max(),50)
    a.plot(xx,bh[0]+bh[1]*xx,color="0.2",lw=1.5,label=f"inverse-variance: {bh[1]:+.4f}/\u00b0C (p = {2*stats.norm.sf(abs(bh[1]/bhse[1])):.3f})")
    a.plot(xx,su[0]+su[1]*xx,color="0.2",lw=1.2,ls=":",label=f"equal weights: {su[1]:+.4f}/\u00b0C (p = {2*stats.norm.sf(abs(su[1]/suse[1])):.3f})")
    a.set_xlabel("garden temperature  bio1 (\u00b0C)"); a.set_ylabel("\u03b2\u2081 per garden (point size \u221d precision)",fontsize=7.5)
    a.legend(frameon=False,fontsize=6.6,loc="lower left")
    TER=[("cold <12 \u00b0C",x<12),("mid 12\u201315 \u00b0C",(x>=12)&(x<=15)),("hot >15 \u00b0C",x>15)]; PL={}
    for l_,m_ in TER:
        ww=1/v[m_]; PL[l_]=((ww*b[m_]).sum()/ww.sum(),np.sqrt(1/ww.sum()))
    a.annotate("pooled by tercile\n"+"\n".join(f"{l_}: {PL[l_][0]:+.3f} \u00b1 {PL[l_][1]:.3f}" for l_,_ in TER),
               xy=(0.97,0.96),xycoords="axes fraction",ha="right",va="top",fontsize=6.6,color="0.3")
    fig.savefig(f"{G}/plots/{fname}.png",dpi=200,bbox_inches="tight")
    fig.savefig(f"{G}/plots/{fname}.pdf",bbox_inches="tight"); plt.show()
    return M,bh,bhse,su,suse,PL,TER
def burden_report(kind,M,bh,bhse,su,suse,PL,TER):
    b,se,x,v,g=M["b"],M["se"],M["x"],M["v"],M["g"]
    print(f"pooled \u03b2\u2081 (inverse variance): {M['pooled']:+.4f} \u00b1 {M['pse']:.4f}  z = {M['pooled']/M['pse']:+.2f}  p = {2*stats.norm.sf(abs(M['pooled']/M['pse'])):.3g}")
    print(f"negative in {(b<0).sum()}/{len(b)} gardens; nominally p<0.05 in {(g.p<0.05).sum()}")
    print(f"\nheterogeneity: Q = {M['Q']:.1f} on {M['dfQ']} df (p = {M['pQ']:.3g}), I\u00b2 = {M['I2']:.1f}%")
    print(f"  SD of the \u03b2\u2081 = {b.std(ddof=1):.4f}; SD expected from sampling noise alone = {np.sqrt(v.mean()):.4f}")
    print(f"  excess SD attributable to real between-garden differences = {np.sqrt(max(0.,b.var(ddof=1)-v.mean())):.4f}")
    print(f"\nmeta-regression of \u03b2\u2081 on garden bio1:")
    print(f"  inverse variance: {bh[1]:+.5f}/\u00b0C \u00b1 {bhse[1]:.5f}  p = {2*stats.norm.sf(abs(bh[1]/bhse[1])):.4g}")
    print(f"  equal weights:    {su[1]:+.5f}/\u00b0C \u00b1 {suse[1]:.5f}  p = {2*stats.norm.sf(abs(su[1]/suse[1])):.4g}")
    rp=stats.pearsonr(x,b); rsp=stats.spearmanr(x,b)
    print(f"  unweighted: Pearson r = {rp[0]:+.3f} (p = {rp[1]:.3g}), Spearman rho = {rsp.statistic:+.3f} (p = {rsp.pvalue:.3g})")
    rng2=np.random.default_rng(0)
    nl=np.array([metareg(b,v,rng2.permutation(x))[0][1] for _ in range(5000)])
    print(f"  permutation p (shuffle bio1, 5,000 draws): {(np.abs(nl)>=abs(bh[1])).mean():.4f}")
    LOO=[]
    for i in range(len(b)):
        m_=np.ones(len(b),bool); m_[i]=False; s_,sse_=metareg(b[m_],v[m_],x[m_])
        LOO.append(dict(dropped_site=int(g.site[i]),bio1=x[i],beta=b[i],se=se[i],slope=s_[1],p=2*stats.norm.sf(abs(s_[1]/sse_[1]))))
    LO=pd.DataFrame(LOO).sort_values("p",ascending=False)
    print(f"\nleave-one-garden-out on the slope: p > 0.05 after dropping {int((LO.p>0.05).sum())} of {len(b)} single gardens")
    print(LO.head(4).round(5).to_string(index=False))
    print("\npooled \u03b2\u2081 per temperature tercile:")
    for l_,m_ in TER: print(f"  {l_:<14} n = {int(m_.sum()):>2}   {PL[l_][0]:+.4f} \u00b1 {PL[l_][1]:.4f}   [{PL[l_][0]-1.96*PL[l_][1]:+.4f}, {PL[l_][0]+1.96*PL[l_][1]:+.4f}]")
    k_=list(PL); (c_,cs_),(m2_,ms_),(h_,hs_)=PL[k_[0]],PL[k_[1]],PL[k_[2]]
    dz=(h_-c_)/np.sqrt(hs_**2+cs_**2)
    print(f"  hot \u2212 cold = {h_-c_:+.4f} \u00b1 {np.sqrt(hs_**2+cs_**2):.4f}  (z = {dz:+.2f}, p = {2*stats.norm.sf(abs(dz)):.3g})")
    print(f"  mid \u2212 cold = {m2_-c_:+.4f};  hot \u2212 mid = {h_-m2_:+.4f}")
RI=burden_figure("ins","burden_by_garden_climate"); burden_report("ins",*RI)
"""

md_burden_read = """## What the founder-level test settles

**The purging effect itself is strong, and it is there in nearly every garden.** Pooled
β₁ = **−0.20 ± 0.02** (p ≈ 7 × 10⁻¹⁷), negative in **29 of 31** gardens and nominally significant on its
own in 12 of them. This is with kinship in the model and with intergenic insertion burden conditioned
out, so it is a CDS-specific contrast with relatedness controlled — the two things the variant-level
analysis could not do. It is also the shape Figure 2c calls purifying selection: an offset below zero
across the whole temperature range, present in cold gardens too.

**But the climate gradient does not survive here.** Three results say the per-garden β₁ are effectively
one number:

- **Heterogeneity.** *Q* = 30.6 on 30 df (p = 0.43), I² = 2%. The SD of the β₁ is 0.126 and the SD
  expected from their own sampling noise alone is 0.138 — the excess attributable to real between-garden
  differences is **0.000**. There is no between-garden variation left for climate to explain.
- **The slope is weighting-dependent and one-garden-dependent.** Inverse-variance meta-regression gives
  −0.0148/°C (p = 0.012, permutation p = 0.026), but at equal weights it is −0.0100/°C (p = 0.11), and
  the unweighted correlation is Spearman ρ = −0.24 (p = 0.20). Leave-one-out: dropping site 23 alone —
  cold, β₁ = +0.054, the smallest SE in the set and therefore the highest weight — takes the slope to
  p = 0.20. One garden out of 31 carries the result.
- **The shape is not a gradient.** Pooled by tercile: cold −0.148 ± 0.039, mid −0.232 ± 0.043, hot
  −0.236 ± 0.044. Mid and hot are identical (difference −0.004); hot − cold = −0.088 ± 0.059
  (z = −1.5, p = 0.13). At most this is a cold-versus-rest step, not something monotone in temperature,
  and it is not significant.

**So the climate reading and the purging reading come apart.** The purging of coding insertions is
better supported than the variant-level panels suggested, because at the founder level it survives the
control that matters (kinship) and is uniform across gardens. The *climate gradient* is the part that
does not survive: it is strong in the variant-level statistic — which cannot separate a variant from the
lineages carrying it — and absent once relatedness is in the model. The most economical account is that
the gradient reflects which founder lineages carry insertions and how those lineages fare in warm
gardens, while the burden effect on fitness is the same everywhere."""

md_del = """## The same thing for deletions

Identical model, identical fitting code (`_sv_burden_by_garden_kinds.py`, which imports the `LMM`
class and the design from `_sv_burden_lmm.py` rather than re-implementing them), with the CDS
**deletion** rate as β₁ and the intergenic deletion rate as the within-class baseline. Short indels and
SNPs are fitted the same way in the next panel, as the scale against which an SV effect should be read.
The insertion fit above is reproduced by this script to the fourth decimal (pooled −0.2012 ± 0.0241),
which is the check that the model is being reused faithfully."""

code_del = r"""
RD=burden_figure("del","burden_by_garden_climate_del"); burden_report("del",*RD)
"""

md_kinds = """### All four variant kinds, per SD of burden and per variant

The CDS-vs-intergenic burden effect is not unique to SVs, and the comparison sets the scale for how
much weight the insertion result can carry. It has to be read on **two** scales, because they order the
classes differently and only one of them is biologically meaningful for comparing variant types.

**(a) Per SD of burden** is what the model returns. On that scale SNPs look strongest. But an SD of
burden is a wildly different number of variants in each class: the median founder carries ≈79,000 of
376,000 called CDS SNPs, so an SD of that rate is ≈10,000 SNPs, whereas it carries ≈46 of 347 called
CDS insertions and an SD is ≈9.5 insertions. Per-SD effects therefore mostly reflect how much room
each class has to vary between founders.

**(b) Per extra variant** divides that out: β₁ / SD(rate) / (variants called), the cost of carrying one
more coding variant of that class. This reverses the ordering, and it is the scale on which "an
insertion is worse than a SNP" is a claim about mutations rather than about sample sizes. Plotted as
|effect| on a log axis because the classes span nearly three orders of magnitude. Caveat: this is a
local slope of a model linear in the *rate*, so it should not be multiplied up to a total load.

**(c)** The per-garden effects against garden temperature, per kind — flat everywhere, which is the
point of the two preceding sections."""


code_kinds = r"""
BD=pd.read_csv(f"{G}/founder_burden_event.csv",dtype={"founder":str})
KO=["snp","indel","ins","del"]; KRD=[]; SC=[]
for k_ in KO:
    c_=f"{k_}_shared_CDS_all"; ig=f"{k_}_shared_intergenic_all"
    carried=BD[c_].to_numpy(float); called=BD[f"called_{c_}"].to_numpy(float); rt=carried/called
    rg=BD[ig].to_numpy(float)/BD[f"called_{ig}"].to_numpy(float)
    M=burden_meta(k_); nmed=np.median(called)
    pv=M["pooled"]/rt.std()/nmed; pv_se=M["pse"]/rt.std()/nmed
    rho=stats.spearmanr(M["x"],M["b"])
    KRD.append(dict(kind=KL[k_],pooled_beta=M["pooled"],se=M["pse"],p=2*stats.norm.sf(abs(M["pooled"]/M["pse"])),
                    negative_in=str((M["b"]<0).sum())+"/"+str(len(M["b"])),Q=M["Q"],p_Q=M["pQ"],I2=M["I2"],
                    rho_bio1=rho.statistic,p_rho=rho.pvalue))
    SC.append(dict(kind=k_,label=KL[k_],carried_median=np.median(carried),called_median=nmed,rate_mean=rt.mean(),
                   variants_per_SD=rt.std()*nmed,beta_per_SD=M["pooled"],pct_per_SD=100*(np.exp(M["pooled"])-1),
                   beta_per_variant=pv,se_per_variant=pv_se,pct_per_variant=100*(np.exp(pv)-1),
                   r_callrate_sv=np.corrcoef(rt,BD.callrate_sv)[0,1],r_own_intergenic=np.corrcoef(rt,rg)[0,1]))
S_=pd.DataFrame(SC).set_index("kind")
fig=plt.figure(figsize=(14.5,3.6)); gsk=fig.add_gridspec(1,3,width_ratios=[1,1.1,1.15],wspace=0.42)
# (a) per SD of burden -- the scale the model returns
a=fig.add_subplot(gsk[0,0]); style(a,axis="x"); panel(a,"a"); a.axvline(0,color="0.35",lw=1,ls="--")
for yi,k_ in enumerate(KO):
    M=burden_meta(k_); a.errorbar(M["pooled"],yi,xerr=1.96*M["pse"],fmt="o",ms=6,color=KC[k_],capsize=3,elinewidth=1.4)
    a.annotate(f"{S_.variants_per_SD[k_]:,.0f} variants = 1 SD",xy=(0.03,yi+0.26),xycoords=("axes fraction","data"),fontsize=6,color="0.45")
a.set_yticks(range(len(KO))); a.set_yticklabels([KL[k_] for k_ in KO],fontsize=8); a.set_ylim(-0.6,len(KO)-0.4)
a.set_xlabel("pooled \u03b2\u2081 per SD of burden\n(founder log-fitness, \u00b1 95% CI)",fontsize=7.5)
# (b) per extra variant -- the scale on which the classes are comparable
a=fig.add_subplot(gsk[0,1]); style(a,axis="x"); panel(a,"b")
for yi,k_ in enumerate(KO):
    e_=abs(S_.beta_per_variant[k_]); s_=S_.se_per_variant[k_]
    a.errorbar(e_,yi,xerr=[[min(1.96*s_,e_*0.95)],[1.96*s_]],fmt="o",ms=6,color=KC[k_],capsize=3,elinewidth=1.4)
    a.annotate(f"{S_.pct_per_variant[k_]:+.3f}% per variant",xy=(e_,yi+0.26),fontsize=6.2,color="0.4",ha="center")
a.set_xscale("log"); a.set_yticks(range(len(KO))); a.set_yticklabels([KL[k_] for k_ in KO],fontsize=8)
a.set_ylim(-0.6,len(KO)-0.4); a.set_xlim(2e-5,8e-2)
a.set_xlabel("cost of carrying ONE more coding variant\n|\u03b2\u2081| per variant, log scale (\u00b1 95% CI)",fontsize=7.5)
a.annotate("ordering reverses:\nan insertion costs \u2248280\u00d7 a SNP",xy=(0.03,0.80),xycoords="axes fraction",fontsize=6.6,color="0.4")
# (c) per garden against garden temperature
a=fig.add_subplot(gsk[0,2]); style(a); panel(a,"c"); a.axhline(0,color="0.35",lw=1,ls="--")
for k_ in KO:
    M=burden_meta(k_); bh,_=metareg(M["b"],M["v"],M["x"]); xx=np.linspace(M["x"].min(),M["x"].max(),50)
    a.scatter(M["x"],M["b"],s=12,color=KC[k_],alpha=0.55,linewidth=0,zorder=2)
    a.plot(xx,bh[0]+bh[1]*xx,color=KC[k_],lw=1.8,zorder=3,label=KL[k_]+f": {bh[1]:+.4f}/\u00b0C")
a.set_xlabel("garden temperature  bio1 (\u00b0C)"); a.set_ylabel("\u03b2\u2081 per garden (per SD)",fontsize=7.5)
a.legend(frameon=False,fontsize=6.6,loc="lower left",ncol=2,columnspacing=1.0)
fig.savefig(f"{G}/plots/burden_by_kind.png",dpi=200,bbox_inches="tight")
fig.savefig(f"{G}/plots/burden_by_kind.pdf",bbox_inches="tight"); plt.show()
print("CDS-vs-intergenic burden effect on founder log-fitness, per SD of burden:")
print(pd.DataFrame(KRD).round(4).to_string(index=False))
print("\nthe same effects per VARIANT, with what a founder actually carries:")
print(S_.reset_index(drop=True)[["label","carried_median","called_median","rate_mean","variants_per_SD",
      "beta_per_SD","pct_per_SD","beta_per_variant","pct_per_variant"]].round(6).to_string(index=False))
iv=S_.beta_per_variant
print(f"\nper-variant ratios:  insertion / deletion = {iv['ins']/iv['del']:.2f}x   "
      f"insertion / short indel = {iv['ins']/iv['indel']:.1f}x   insertion / SNP = {iv['ins']/iv['snp']:.0f}x")
print("(a local slope of a model linear in the burden RATE -- do not multiply up to a total load)")
print("\nconfounding checks on the burden rates across the 231 founders:")
print(S_.reset_index(drop=True)[["label","r_callrate_sv","r_own_intergenic"]].round(3).to_string(index=False))
RR={KL[k_]:BD[f"{k_}_shared_CDS_all"].to_numpy(float)/BD[f"called_{k_}_shared_CDS_all"].to_numpy(float) for k_ in KO}
print("\ncorrelation between the CDS burden rates of the four kinds:")
print(pd.DataFrame(RR).corr().round(3).to_string())
JT=pd.read_csv(f"{G}/sv_burden_lmm_joint_event.csv")
print("\njoint models on mean fitness -- does the insertion effect survive conditioning on SNP / indel CDS load?")
print(JT.round(4).to_string(index=False))
"""


md_del_read = """## What deletions say

**CDS-deletion burden also predicts lower founder fitness, but about 2.6× more weakly than
insertions.** Pooled β₁ = **−0.078 ± 0.017** (p = 4.4 × 10⁻⁶), negative in 28 of 31 gardens, against
−0.201 ± 0.024 for insertions. Both classes of coding SV look deleterious at the founder level, with
insertions the stronger of the two — the ordering the founder-panel spectra in Figure 1 imply.

**No climate gradient here either, and if anything it leans the other way.** I² = 0%, *Q* = 22.4 on
30 df (p = 0.84) — the per-garden deletion effects are even more uniform than the insertion ones.
β₁ vs garden bio1: Spearman ρ = **+0.25** (p = 0.17), i.e. the point estimate is *less* negative in
warmer gardens, the opposite sign to insertions. Neither is significant. Two uniform effects with
opposite-leaning non-significant slopes is what no climate gradient looks like.

**The deletion result contradicts the variant-level panel, and that is worth confronting.** At the
variant level, CDS deletions *rise* relative to matched SNPs (+0.7% [+0.1, +1.2], and +1.0% at
MAF ≥ 5% with null p = 0.028). At the founder level, founders carrying more coding deletions have
*lower* fitness. Both cannot be straightforward selection on coding deletions — if deletion-carrying
founders lose, coding deletion alleles should fall, not rise. The two contrasts are not identical
(matched-SNP versus CDS-vs-intergenic-within-deletions), but the sign flip is real, and it goes the same
way as everything else in this notebook: the level that controls relatedness says "deleterious", the
level that cannot says the opposite.

**Scale check — per SD of burden versus per variant, which order differently.** Per SD the effect is
largest for **SNPs** (−0.759 ± 0.068, negative in 31/31), then short indels (−0.519), insertions
(−0.201), deletions (−0.078). But an SD of burden is a very different number of variants in each class:
the median founder carries ≈79,000 of 376,000 called CDS SNPs and an SD of that rate is ≈10,000 SNPs,
whereas it carries ≈46 of 347 called CDS insertions and an SD is ≈9.5 insertions. Rescaling to one
additional variant reverses the order completely:

| CDS variant | per SD of burden | 1 SD in variants | per extra variant |
|---|---|---|---|
| SNP | −0.759 | ≈10,050 | −0.000076 |
| short indel | −0.519 | ≈279 | −0.00186 |
| **deletion** | −0.078 | ≈8.5 | **−0.00913** |
| **insertion** | −0.201 | ≈9.5 | **−0.0211** |

So per variant a coding insertion costs ≈2.3× a coding deletion, ≈11× a short indel and ≈280× a SNP —
a monotone ordering in how much coding sequence the mutation disrupts, which is what a functional
interpretation predicts. In log-fitness terms, one extra coding insertion is ≈2% lower lineage growth
over the experiment and a founder at +1 SD of coding-insertion burden is ≈18% lower.

**And the insertion effect is not a proxy for coding SNP or indel load.** The burden rates are
correlated across founders (insertion vs SNP *r* = 0.74, vs indel 0.79, vs deletion 0.47), so this has
to be tested rather than assumed. The joint models: conditioning on SNP CDS burden leaves
β = −0.187 (p = 0.0046), and on both SNP and indel CDS burden leaves β = −0.175 (p = 0.0079). No
equivalent joint model has been run for deletions, or for insertions conditional on deletions.

**One asymmetry to keep in view.** The CDS insertion rate correlates with the founder's SV call rate
at *r* = +0.45, while the deletion rate does not (−0.01). Call rate is a covariate in the model, so it
is adjusted for linearly, but the stronger of the two effects is also the one more entangled with how
well that founder is genotyped — a residual-confounding risk that does not apply to deletions.

*Caveats.* Gardens share one founder pool, so "28/31" is descriptive rather than 28 independent tests;
insertion and deletion burdens are correlated across founders, so their β₁ are not independent of each
other; and the SNP burden here is a random 10% of records, not the full SNP set."""

# ======================================================================================================
# Founder-level notebook: ONE full model (SNPs + insertions + deletions together). This is what
# `sv_climate_maf.ipynb` is built from. The variant-level sections above (md_fig ... md_par_read) and the
# earlier per-kind burden sections (md_burden ... md_del_read) are kept for reference but are not built.
# All numbers come from saved results -- the slow permutations are in `_sv_burden_full_model.py`.
# ======================================================================================================

md_title_fm = """# Founder burden and fitness in GrENE-Net — SNPs, insertions and deletions in one model

**Question.** Do founders that carry more variants in coding sequence do worse in the experiment — and
is that true for structural variants, not only for SNPs? One mixed model answers it for SNPs, insertions
and deletions **together**, so each class's effect is estimated with the other two held fixed, and with
relatedness between founders in the model.

**Contents.** *Figure 1* — the effect of each class, per carried variant. *Figure 2* — the same effect
garden by garden, and against garden temperature. *Supplementary Figure S1* — how the model was
simplified. *Supplementary Figure S2* — the full model on the scale it is fitted (per SD). *Table S1* —
joint models from one class up to all three. *Key numbers.* *Supplementary Figures S3–S4* — the same
model for promoter burden; *S5* — its intergenic terms garden by garden. *One-stage model* — all 7,161
founder × garden observations in one mixed model (Figures O1–O3). The final models are in
`sv_burden_final_models.ipynb`.

**Unit and filters.** SV **event** (arch3 alleles of a pangenome bubble merged with `truvari collapse`
defaults). Every variant called in ≥ 90% of the 231 founders and carried by ≥ 2 of them; **shared**
variants only (carried by both cold- and warm-origin founders); all SNP records. Context,
first match wins: CDS > UTR > intron > reference TE > promoter (≤ 2 kb upstream) > intergenic."""

code_setup_fm = r"""
import os
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({"figure.dpi":120,"font.size":8,"axes.linewidth":0.6,"savefig.bbox":"tight","legend.frameon":False})
G=os.path.abspath(os.path.join(os.getcwd(),"..","r1_sv_negative_selection","results","sv_adaptive")); os.makedirs(f"{G}/plots",exist_ok=True)
KINDS=("snp","ins","del"); KC={"snp":"0.15","ins":"#D95F02","del":"#1B9E77"}; YL={"snp":"SNPs","ins":"insertions","del":"deletions"}
def style(a,axis="x"):
    a.set_axisbelow(True); a.grid(color="0.9",lw=0.6,axis=axis)
    for s_ in a.spines.values(): s_.set_visible(False)
def panel(a,letter,x=-0.08,y=1.04): a.annotate(letter,xy=(x,y),xycoords="axes fraction",fontsize=11,fontweight="bold",ha="right")
M="−"                                                            # typographic minus for labels
def pct(v,d=3): return f"{v:+.{d}g}%".replace("-",M)
def pfmt(p,n=None):
    if n is not None and p==0: return f"p < {1/n:.0e}".replace("e-0","e-")
    return "p < 0.001" if p<0.001 else f"p = {p:.2g}"
def save(fig,name): fig.savefig(f"{G}/plots/{name}.png",dpi=200,bbox_inches="tight"); fig.savefig(f"{G}/plots/{name}.pdf",bbox_inches="tight")
FM=pd.read_csv(f"{G}/sv_burden_full_model_event.csv").set_index("term")          # mean fitness, all terms (_sv_burden_full_model.py)
FG=pd.read_csv(f"{G}/sv_burden_full_model_by_garden_event.csv")                   # the same model per garden
SIMP=pd.read_csv(f"{G}/sv_burden_model_simplify_event.csv")                        # model ladder (_sv_burden_model_simplify.py)
JOINT=pd.read_csv(f"{G}/sv_burden_joint_snp_ins_del_event.csv")                    # one class up to all three
PERM=np.load(f"{G}/sv_burden_full_model_perm_event.npz")                          # null draws per coding term (per garden, per SD)
def ew_slope(t):
    # equal-weights (OLS) slope of the per-garden coefficients on garden bio1, per SD; permutation p from the saved null
    # per-garden coefficients (same founder-label shuffle in every garden, so the dependence between gardens is kept)
    g0=FG[FG.term==t]; x=g0.bio1.to_numpy(); xc=x-x.mean(); b=g0.beta_sd.to_numpy()
    obs=(xc@(b-b.mean()))/(xc@xc); NGt=PERM[f"{t}_null_garden"].astype(float); nl=((NGt-NGt.mean(1,keepdims=True))@xc)/(xc@xc)
    return obs,(np.abs(nl)>=abs(obs)).mean()
NF=len(pd.read_csv(f"{G}/founder_burden_event.csv",usecols=["founder"])); NG=FG.site.nunique(); NPERM=int(FM.n_perm.dropna().iloc[0])
print(f"{NF} founders, {NG} gardens, {NPERM:,} founder-label permutations per tested coding term; kinship share of variance {FM.h2.iloc[0]:.2f}")
"""

md_fm = """## Figure 1 — the effect of each variant class, per carried variant

**Model.** For founder *f* and garden *s*, *y*<sub>f</sub> = log[(*H*<sub>0,f</sub> + Δ*H*<sub>s,f</sub> + 10⁻⁵) / (*H*<sub>0,f</sub> + 10⁻⁵)]
is the founder's log change in frequency from the seed mix to the last generation, averaged over the
31 gardens:

> *y* = β₀ + Σ<sub>k ∈ {SNP, ins, del}</sub> [ β<sub>k,coding</sub>·*r*<sub>k,coding</sub> + β<sub>k,intergenic</sub>·*r*<sub>k,intergenic</sub> ]
> + γ·[assembly route] + *u* + ε,   *u* ~ N(0, σ²<sub>g</sub>**K**)

*r* = the variants of that class and context the founder carries ÷ those it is called at; **K** = SNP
kinship. Each class enters twice, coding and intergenic, so its coding coefficient asks: at the same
overall load of that class, does carrying more of it **in coding sequence** change fitness? Assembly
route (cactus vs short-read) is the single technical term; the SV and SNP call rates were redundant with
it and with each other, and dropping them leaves the coding-insertion coefficient unchanged
(Supplementary Figure S1).

**Per carried variant.** The model is fitted on standardised rates; each coefficient is then divided by
the number of variants one SD of that rate represents, giving the % change in founder fitness for **one
more carried variant**, all else fixed. Only this scale compares classes: one SD of coding burden is
≈ 10,000 SNPs but ≈ 10 insertions. It is the *average* over the class's variants — most coding SNPs are
synonymous or harmless — and a local slope, not something to multiply up to a total load. The per-SD
version is Supplementary Figure S2.

**(a)** Each class on its own x-range, so every estimate and interval is readable. **(b)** The same six
values on one shared axis, so the classes compare directly. Circles = coding, squares = intergenic;
filled = p < 0.05, where coding terms use the founder-label permutation p (3,000 draws: the term's rate is
residualised on all other terms and the residual shuffled across founders) and intergenic terms the Wald p."""

code_fm = r"""
fig=plt.figure(figsize=(15.5,4.4)); gsf=fig.add_gridspec(1,2,width_ratios=[1.35,1],wspace=0.42)
def pval(t): return FM.loc[t,"p_perm"] if t.endswith("CDS") else FM.loc[t,"p_wald"]
def plab(t): return pfmt(pval(t),NPERM) if t.endswith("CDS") else pfmt(pval(t))
# (a) one small panel per class, own symmetric x-range
gsa=gsf[0,0].subgridspec(1,3,wspace=0.15)
for ci,k in enumerate(KINDS):
    a=fig.add_subplot(gsa[0,ci]); style(a); a.axvline(0,color="0.4",lw=0.9,ls="--")
    ext=max(abs(100*FM.loc[f"{k}_{x}","beta_per_variant"])+2.3*100*FM.loc[f"{k}_{x}","se_per_variant"] for x in ("CDS","intergenic"))
    for yi,(x,mk) in enumerate((("CDS","o"),("intergenic","s"))):
        t=f"{k}_{x}"; b=100*FM.loc[t,"beta_per_variant"]; s_=100*FM.loc[t,"se_per_variant"]; sig=pval(t)<0.05
        a.errorbar(b,yi,xerr=1.96*s_,fmt=mk,ms=7 if x=="CDS" else 6,color=KC[k],mfc=KC[k] if sig else "white",mew=1.4,capsize=3,elinewidth=1.5)
        a.annotate(f"{pct(b)}\n{plab(t)}",xy=(b,yi-0.17),ha="center",va="bottom",fontsize=6.8,color="0.3")
    a.set_xlim(-ext*1.1,ext*1.1); a.set_ylim(1.55,-0.75); a.set_yticks([0,1])
    a.set_yticklabels(["coding","intergenic"] if ci==0 else ["",""],fontsize=8.5)
    if ci>0: a.tick_params(axis="y",left=False)
    a.tick_params(axis="x",labelsize=7); a.annotate(YL[k],xy=(0.5,1.02),xycoords="axes fraction",ha="center",fontsize=9,color=KC[k])
    if ci==0: panel(a,"a",x=-0.32)
    if ci==1: a.set_xlabel("% change in founder fitness per carried variant (± 95% CI); each class on its own x-range",fontsize=8)
# (b) all six on one shared axis
a=fig.add_subplot(gsf[0,1]); style(a); panel(a,"b",x=-0.3); a.axvline(0,color="0.4",lw=0.9,ls="--")
rows_=[(f"{k}_{x}",k,x) for k in KINDS for x in ("CDS","intergenic")]
for yi,(t,k,x) in enumerate(rows_):
    b=100*FM.loc[t,"beta_per_variant"]; s_=100*FM.loc[t,"se_per_variant"]; sig=pval(t)<0.05
    a.errorbar(b,yi,xerr=1.96*s_,fmt="o" if x=="CDS" else "s",ms=7 if x=="CDS" else 6,color=KC[k],mfc=KC[k] if sig else "white",mew=1.4,capsize=3,elinewidth=1.5)
    a.annotate(f"{pct(b)}   {plab(t)}",xy=(1.02,yi),xycoords=("axes fraction","data"),va="center",fontsize=7,color="0.3",annotation_clip=False)
    if yi in (1,3): a.axhline(yi+0.5,color="0.88",lw=0.8)
a.set_yticks(range(len(rows_))); a.set_yticklabels([f"{YL[k]}, {'coding' if x=='CDS' else 'intergenic'}" for _,k,x in rows_],fontsize=8)
for tl,(_,k,_) in zip(a.get_yticklabels(),rows_): tl.set_color(KC[k])
a.set_ylim(len(rows_)-0.5,-0.6); a.set_xlabel("% change in founder fitness per carried variant (± 95% CI); one shared axis",fontsize=8)
a.errorbar([],[],fmt="o",color="0.3",ms=6,label="p < 0.05"); a.errorbar([],[],fmt="o",color="0.3",mfc="white",ms=6,label="n.s.")
a.legend(fontsize=7,loc="lower left",handletextpad=0.2)
save(fig,"burden_full_model_effects"); plt.show()
T=pd.DataFrame({"term":[f"{YL[t.split('_')[0]]}, {'coding' if t.endswith('CDS') else 'intergenic'}" for t in FM.index[:6]],
                "% per variant":[pct(100*FM.loc[t,'beta_per_variant']) for t in FM.index[:6]],
                "95% CI":[f"[{pct(100*(FM.loc[t,'beta_per_variant']-1.96*FM.loc[t,'se_per_variant']))}, {pct(100*(FM.loc[t,'beta_per_variant']+1.96*FM.loc[t,'se_per_variant']))}]" for t in FM.index[:6]],
                "per SD":[f"{FM.loc[t,'beta_sd']:+.3f}".replace("-",M) for t in FM.index[:6]],
                "variants per SD":[f"{FM.loc[t,'variants_per_sd']:,.1f}" for t in FM.index[:6]],
                "Wald p":[f"{FM.loc[t,'p_wald']:.2g}" for t in FM.index[:6]],
                "permutation p":[pfmt(FM.loc[t,'p_perm'],NPERM) if t.endswith('CDS') else "" for t in FM.index[:6]]})
print("Figure 1 — full model on mean fitness (all terms fitted together):"); print(T.to_string(index=False))
print(f"assembly route (cactus − short-read): {FM.loc['cactus','beta_sd']:+.3f} ± {FM.loc['cactus','se_sd']:.3f} (Wald p = {FM.loc['cactus','p_wald']:.2g})")
"""

md_gd = """## Figure 2 — garden by garden

The same full model fitted in each of the 31 gardens; one row per coding term (SNPs, insertions,
deletions). All panels share **one per-variant scale**, so the classes compare directly — the SNP
estimates are real but tiny on it (values in the labels and the table below).

**Left:** each garden's coefficient ± 95% CI, gardens sorted by effect (most negative at the bottom) and
coloured by garden temperature. Line and band = pooled coefficient ± 1.96 SE, with the SE from the
**founder-label permutation**: the 31 fits share the same founders, so the textbook inverse-variance SE
(printed for comparison) is 2–2.5× too small. **Right:** the same coefficients against garden
temperature (grey bars = 95% CI), with an ordinary least-squares line, every garden weighted equally; its
p-value is from the same founder-label permutation (the saved null coefficients of all 31 gardens), which
keeps the dependence between gardens. A leave-one-garden-out check of the slope is printed."""

code_gd = r"""
B1=FG.drop_duplicates("site").bio1; NORM=mpl.colors.Normalize(B1.min(),B1.max()); CMAP=plt.get_cmap("coolwarm")
def mreg(b,v,x):
    w=1/v; X=np.column_stack([np.ones(b.size),x]); XtW=X.T*w; bh=np.linalg.solve(XtW@X,XtW@b); return bh,np.sqrt(np.diag(np.linalg.inv(XtW@X)))
fig=plt.figure(figsize=(13.5,12.5)); gsg=fig.add_gridspec(3,2,width_ratios=[1.3,1],wspace=0.25,hspace=0.45); OUTG=[]; AXL=[]; AXR=[]
for r_,k in enumerate(KINDS):
    t=f"{k}_CDS"; g=FG[FG.term==t].sort_values("beta_sd").reset_index(drop=True); r=FM.loc[t]; sc=r.beta_per_variant/r.beta_sd
    b=100*g.beta_per_variant.to_numpy(); se=100*g.se_per_variant.to_numpy(); x=g.bio1.to_numpy(); v=se**2
    po=100*r.pooled_sd*sc; pse=100*r.pooled_se_perm*sc; pse_ind=100*r.pooled_se_independent*sc
    a=fig.add_subplot(gsg[r_,0]); AXL.append(a); style(a); panel(a,"abc"[r_],x=-0.07); a.axvline(0,color="0.4",lw=0.9,ls="--")
    a.axvspan(po-1.96*pse,po+1.96*pse,color=KC[k],alpha=0.14,lw=0); a.axvline(po,color=KC[k],lw=1.4)
    for yi in range(len(g)):
        c_=CMAP(NORM(g.bio1[yi])); a.errorbar(b[yi],yi,xerr=1.96*se[yi],fmt="o",ms=4.2,color=c_,mec="0.25",mew=0.5,ecolor=c_,elinewidth=1.0,capsize=1.6)
    a.set_yticks(range(len(g))); a.set_yticklabels(g.site.astype(int),fontsize=5.5); a.set_ylim(-1,len(g)); a.set_ylabel("garden (site), sorted by effect",fontsize=7.5)
    a.annotate(f"{YL[k]}, coding", xy=(0.0,1.075),xycoords="axes fraction",fontsize=9,color=KC[k],va="bottom",annotation_clip=False)
    a.annotate(f"pooled {pct(po)} ± {pse:.2g} (permutation SE), {pfmt(r.p_perm_pooled,NPERM)}; negative in {int(r.negative_in)}/{NG} gardens",
               xy=(0.0,1.01),xycoords="axes fraction",ha="left",va="bottom",fontsize=7,color="0.3",annotation_clip=False)
    if r_==2: a.set_xlabel("% change in founder fitness per carried coding variant (± 95% CI)",fontsize=8)
    a=fig.add_subplot(gsg[r_,1]); AXR.append(a); style(a,axis="both"); a.axhline(0,color="0.4",lw=0.9,ls="--"); a.axhline(po,color=KC[k],lw=1.1,alpha=0.7)
    a.errorbar(x,b,yerr=1.96*se,fmt="none",ecolor="0.82",elinewidth=0.8,zorder=1)
    a.scatter(x,b,s=34,c=x,cmap=CMAP,norm=NORM,linewidth=0.4,edgecolor="0.3",zorder=3)
    su,_=mreg(b,np.ones_like(v),x); ews,ewp=ew_slope(t); xx=np.linspace(x.min(),x.max(),50)
    a.plot(xx,su[0]+su[1]*xx,color="0.2",lw=1.5,label=f"linear fit: {pct(su[1])}/°C, permutation {pfmt(ewp,NPERM)}")
    a.set_ylabel("% per variant, per garden",fontsize=7.5); a.legend(fontsize=6.8,loc="lower left")
    if r_==2: a.set_xlabel("garden temperature  bio1 (°C)",fontsize=8)
    LOO=[]
    for i in range(len(b)):
        m_=np.ones(len(b),bool); m_[i]=False; lr=stats.linregress(x[m_],b[m_]); LOO.append((int(g.site[i]),lr.slope,lr.pvalue))
    L_=pd.DataFrame(LOO,columns=["site","slope","p"]).sort_values("p",ascending=False).iloc[0]
    OUTG.append({"coding term":YL[k],"pooled % per variant":pct(po),"SE (permutation)":f"{pse:.2g}","SE (textbook)":f"{pse_ind:.2g}",
                 "pooled p (perm.)":pfmt(r.p_perm_pooled,NPERM),"gardens < 0":f"{int(r.negative_in)}/{NG}",
                 "slope %/°C (equal weights)":pct(su[1]),"slope p (perm.)":pfmt(ewp,NPERM),"worst leave-one-out p (OLS)":f"{L_.p:.2g} (site {int(L_.site)})"})
xl=[a.get_xlim() for a in AXL]; yl=[a.get_ylim() for a in AXR]                        # one shared per-variant scale
for a in AXL: a.set_xlim(min(v_[0] for v_ in xl),max(v_[1] for v_ in xl))
for a in AXR: a.set_ylim(min(v_[0] for v_ in yl),max(v_[1] for v_ in yl))
for a in AXL[:2]+AXR[:2]: a.tick_params(labelbottom=True)
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORM,cmap=CMAP),ax=AXL+AXR,pad=0.015,fraction=0.016,aspect=40); cb.set_label("garden temperature  bio1 (°C)",fontsize=8)
save(fig,"burden_full_model_by_garden"); plt.show()
print("Figure 2 — per-garden summary, coding terms:"); print(pd.DataFrame(OUTG).T.to_string(header=False))
"""

md_sup = """## Supplement

**Supplementary Figure S1 — how the model was simplified.** Insertions alone (before the joint model).
**(a)** Every coefficient, model by model, from the original five covariates down to the coding rate
alone; kinship in all; AIC from maximum-likelihood refits. The coding-insertion coefficient stays at
−0.21 to −0.23 per SD in every model that keeps the intergenic rate, and collapses without it — the
contrast with the founder's overall insertion load is what carries the effect. The two call rates and
the assembly route are interchangeable: whichever is left in picks up the same technical difference, so
one is enough. **(b)** The coding-insertion coefficient with and without the kinship term: without it
the coefficient inflates and the fit is much worse (AIC ≈ 410 vs ≈ 355).

**Supplementary Figure S2 — the full model per SD.** The scale the model is fitted on, including the
assembly route. Here SNPs look largest because one SD of coding-SNP burden is ≈ 10,000 SNPs against
≈ 10 coding insertions; Figure 1 divides that out.

**Table S1 — joint models.** From one class alone up to all three: the coding-insertion coefficient
is stable, deletions are null throughout, and the insertion coding-minus-intergenic contrast weakens
once SNPs are in, because the intergenic insertion and SNP rates are highly correlated (r = 0.93).

*Caveats.* The model form (coding and intergenic rates, one technical term) was settled after looking
at these fits, which is why the whole ladder is shown. The coding-insertion effect depends on the
call-rate floor: at a 70% floor it weakens to −0.11 per SD (p = 0.11; `sv_main_figures_cr70.ipynb`), and
the coding-insertion rate correlates with founder SV call rate (r = 0.45)."""

code_sup = r"""
MODL={"full":"full: + SV call rate + SNP call rate + route","- callrate_snp":"− SNP call rate","- cactus":"− assembly route",
      "- both call rates":"− both call rates","intergenic + callrate_sv":"intergenic + SV call rate","intergenic only":"intergenic rate only",
      "CDS only":"coding rate only"}
TERML={"CDS":("coding rate",KC["ins"]),"intergenic":("intergenic rate","#7570B3"),"callrate_sv":("SV call rate","0.45"),
       "callrate_snp":("SNP call rate","0.72"),"cactus":("assembly route","#E7298A")}
MOD=list(MODL); I=SIMP[(SIMP.kind=="ins")&SIMP.fit.str.startswith("LMM")]
fig,ax=plt.subplots(1,2,figsize=(13.5,4.6),gridspec_kw=dict(width_ratios=[1.25,1],wspace=0.75))
a=ax[0]; style(a); panel(a,"a",x=-0.55); a.axvline(0,color="0.4",lw=0.9,ls="--")
for i,m in enumerate(MOD):
    d=I[I.model==m].set_index("term")
    for j,t in enumerate(TERML):
        if t in d.index:
            r=d.loc[t]; a.errorbar(r.beta,i+(j-2)*0.14,xerr=1.96*r.se,fmt="o",ms=4.5,color=TERML[t][1],mfc=TERML[t][1] if r.p<0.05 else "white",mew=1.1,capsize=1.5)
    a.annotate(f"{d.aic.iloc[0]:.1f}",xy=(1.03,i),xycoords=("axes fraction","data"),va="center",fontsize=7,color="0.35",annotation_clip=False)
    if i<len(MOD)-1: a.axhline(i+0.5,color="0.92",lw=0.8)
a.annotate("AIC",xy=(1.03,-0.75),xycoords=("axes fraction","data"),va="center",fontsize=7,color="0.35",annotation_clip=False)
a.set_yticks(range(len(MOD))); a.set_yticklabels([MODL[m] for m in MOD],fontsize=7.5); a.set_ylim(len(MOD)-0.5,-0.9)
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([],[],marker="o",ls="",color=c_,ms=5,label=l_) for l_,c_ in TERML.values()],fontsize=7,loc="upper center",bbox_to_anchor=(0.5,-0.16),ncol=5,handletextpad=0.2,columnspacing=0.9)
a.set_xlabel("coefficient per SD (± 95% CI); insertions alone, kinship in all; filled = p < 0.05",fontsize=8)
a=ax[1]; style(a); panel(a,"b",x=-0.55); a.axvline(0,color="0.4",lw=0.9,ls="--")
for fit,mk,al,dy in (("LMM (kinship)","o",1,-0.12),("OLS (no kinship)","x",0.55,0.12)):
    d=SIMP[(SIMP.kind=="ins")&(SIMP.term=="CDS")&(SIMP.fit==fit)].set_index("model").reindex(MOD)
    a.errorbar(d.beta,np.arange(len(MOD))+dy,xerr=1.96*d.se,fmt=mk,ms=5,color=KC["ins"],alpha=al,capsize=1.5,
               label=f"{'with' if mk=='o' else 'without'} kinship")
    for i,m in enumerate(MOD): a.annotate(f"{d.aic[m]:.0f}",xy=(1.03,i+dy),xycoords=("axes fraction","data"),va="center",fontsize=6.5,color="0.35",alpha=al,annotation_clip=False)
a.annotate("AIC",xy=(1.03,-0.75),xycoords=("axes fraction","data"),va="center",fontsize=7,color="0.35",annotation_clip=False)
a.set_yticks(range(len(MOD))); a.set_yticklabels([MODL[m] for m in MOD],fontsize=7.5); a.set_ylim(len(MOD)-0.5,-0.9)
a.legend(fontsize=7,loc="upper center",bbox_to_anchor=(0.5,-0.16),ncol=2)
a.set_xlabel("coding-insertion coefficient per SD (± 95% CI)",fontsize=8)
save(fig,"burden_model_ladder"); plt.show()

# Supplementary Figure S2 -- the full model per SD, with the assembly route
fig,a=plt.subplots(figsize=(7.5,3.9)); style(a); a.axvline(0,color="0.4",lw=0.9,ls="--")
rows_=[(f"{k}_{x}",f"{YL[k]}, {'coding' if x=='CDS' else 'intergenic'}",KC[k],"o" if x=="CDS" else "s") for k in KINDS for x in ("CDS","intergenic")]+[("cactus","assembly route (cactus − short-read)","#E7298A","D")]
for yi,(t,lab,col,mk) in enumerate(rows_):
    r=FM.loc[t]; p=r.p_perm if t.endswith("CDS") else r.p_wald
    a.errorbar(r.beta_sd,yi,xerr=1.96*r.se_sd,fmt=mk,ms=6,color=col,mfc=col if p<0.05 else "white",mew=1.3,capsize=2.5,elinewidth=1.4)
    if yi in (1,3,5): a.axhline(yi+0.5,color="0.88",lw=0.8)
a.set_yticks(range(len(rows_))); a.set_yticklabels([l for _,l,_,_ in rows_],fontsize=8); a.set_ylim(len(rows_)-0.5,-0.6)
for tl,(_,_,col,_) in zip(a.get_yticklabels(),rows_): tl.set_color(col)
a.set_xlabel("effect per SD of the rate (assembly route: per group) ± 95% CI\nfull model, the scale it is fitted on; filled = p < 0.05",fontsize=8)
save(fig,"burden_full_model_per_sd"); plt.show()

# Table S1 -- joint models
JM={"SNP alone":"SNPs","insertions alone":"ins","deletions alone":"del","SNP + insertions":"SNPs + ins","insertions + deletions":"ins + del","SNP + insertions + deletions":"all three"}
TL={"snp_CDS":"SNPs, coding","snp_intergenic":"SNPs, intergenic","ins_CDS":"insertions, coding","ins_intergenic":"insertions, intergenic",
    "del_CDS":"deletions, coding","del_intergenic":"deletions, intergenic","cactus":"assembly route",
    "snp: CDS - intergenic":"SNPs: coding − intergenic","ins: CDS - intergenic":"insertions: coding − intergenic","del: CDS - intergenic":"deletions: coding − intergenic"}
J=JOINT.copy(); J["cell"]=J.beta.map(lambda v:f"{v:+.3f}".replace("-",M))+" ("+J.p.map(lambda p:f"{p:.2g}")+")"
TS=J.pivot(index="term",columns="model",values="cell").reindex(list(TL))[list(JM)].fillna("").rename(index=TL,columns=JM)
TS.loc["AIC"]=[f"{J[J.model==m].aic.iloc[0]:.1f}" for m in JM]
print("Table S1 — joint models, coefficient per SD (p), kinship and assembly route in all:"); print(TS.to_string())
"""

md_key = """## Key numbers

Printed from the saved results, so they always match the figures above."""

code_key = r"""
def line(t):
    r=FM.loc[t]; s=r.beta_per_variant/r.beta_sd
    return (f"{pct(100*r.beta_per_variant)} per carried variant [{pct(100*(r.beta_per_variant-1.96*r.se_per_variant))}, {pct(100*(r.beta_per_variant+1.96*r.se_per_variant))}], "
            f"permutation {pfmt(r.p_perm,NPERM)};  per garden: pooled {pct(100*r.pooled_sd*s)}, {pfmt(r.p_perm_pooled,NPERM)}, negative in {int(r.negative_in)}/{NG}, "
            f"slope on garden temperature {pct(100*ew_slope(t)[0]*s)}/°C (equal weights), {pfmt(ew_slope(t)[1],NPERM)}")
for k in KINDS: print(f"coding {YL[k]:<11} {line(f'{k}_CDS')}")
print(f"\nper carried coding variant, insertion / SNP = {FM.loc['ins_CDS','beta_per_variant']/FM.loc['snp_CDS','beta_per_variant']:.0f}x (point estimates; the SNP average includes the many harmless coding SNPs)")
print(f"kinship share of variance {FM.h2.iloc[0]:.2f}; assembly route {FM.loc['cactus','beta_sd']:+.3f} (Wald p = {FM.loc['cactus','p_wald']:.2g})")
"""

md_prom = """## Supplementary Figures S3–S4 — the same model for promoter burden

The full model refitted with **promoter** burden (≤ 2 kb upstream of a gene, not overlapping a reference TE,
UTR or intron) in place of coding burden, for SNPs, insertions and deletions, each with its intergenic
rate; same covariates, kinship, 3,000 founder-label permutations (`_sv_burden_full_model.py`,
`CONTEXT=promoter`, outputs tagged `_promoter`). Same layout as Figures 1 and 2. Promoter and intergenic
rates are more correlated than coding and intergenic (insertions r = 0.94, deletions 0.83, SNPs 0.98), so
the average effects are less well separated from the intergenic baseline than in the coding model; the
garden-temperature slopes use the permutation p. A burden test averages over all promoter SVs a founder
carries (≈ 3,000 insertions, ≈ 1,000 deletions called per founder); it says nothing about which loci drive
an effect."""

code_prom_load = r"""
FM=pd.read_csv(f"{G}/sv_burden_full_model_event_promoter.csv").set_index("term")
FG=pd.read_csv(f"{G}/sv_burden_full_model_by_garden_event_promoter.csv")
PERM=np.load(f"{G}/sv_burden_full_model_perm_event_promoter.npz")                 # ew_slope() reads FG / PERM
print(f"promoter model: kinship share of variance {FM.h2.iloc[0]:.2f}; permutation draws {int(FM.n_perm.dropna().iloc[0]):,}")
"""


def _promoter(code):
    for a, c in (('save(fig,"burden_full_model_effects")', 'save(fig,"burden_promoter_model_effects")'),
                 ('save(fig,"burden_full_model_by_garden")', 'save(fig,"burden_promoter_model_by_garden")'),
                 ('print("Figure 1 —', 'print("Supplementary Figure S3 —'), ('print("Figure 2 —', 'print("Supplementary Figure S4 —'),
                 ("CDS", "promoter"), ("coding", "promoter")):
        code = code.replace(a, c)
    return code


md_prom_ig = """### Supplementary Figure S5 — the intergenic terms of the same promoter model, garden by garden

The intergenic insertion and deletion coefficients from the promoter model (the within-class baseline of
Supplementary Figures S3–S4), in the Figure 2 layout, on the same per-variant scale. Only the promoter
terms were permuted, so here the pooled estimate carries the textbook inverse-variance SE (too small by
the same ≈ 2–2.5× as in Figure 2) and the slope an ordinary least-squares p — both anti-conservative, so
a flat result is not at risk from them. Compare with Supplementary Figure S4c (promoter deletions)."""

code_prom_ig = r"""
fig=plt.figure(figsize=(13.5,8.6)); gsg=fig.add_gridspec(2,2,width_ratios=[1.3,1],wspace=0.25,hspace=0.45); AXL=[]; AXR=[]; OUT=[]
for r_,k in enumerate(("ins","del")):
    t=f"{k}_intergenic"; g=FG[FG.term==t].sort_values("beta_sd").reset_index(drop=True)
    b=100*g.beta_per_variant.to_numpy(); se=100*g.se_per_variant.to_numpy(); x=g.bio1.to_numpy(); w=1/se**2
    po=(w*b).sum()/w.sum(); pse=np.sqrt(1/w.sum()); pp_=2*stats.norm.sf(abs(po/pse))
    a=fig.add_subplot(gsg[r_,0]); AXL.append(a); style(a); panel(a,"ab"[r_],x=-0.07); a.axvline(0,color="0.4",lw=0.9,ls="--")
    a.axvspan(po-1.96*pse,po+1.96*pse,color=KC[k],alpha=0.14,lw=0); a.axvline(po,color=KC[k],lw=1.4)
    for yi in range(len(g)):
        c_=CMAP(NORM(g.bio1[yi])); a.errorbar(b[yi],yi,xerr=1.96*se[yi],fmt="o",ms=4.2,color=c_,mec="0.25",mew=0.5,ecolor=c_,elinewidth=1.0,capsize=1.6)
    a.set_yticks(range(len(g))); a.set_yticklabels(g.site.astype(int),fontsize=5.5); a.set_ylim(-1,len(g)); a.set_ylabel("garden (site), sorted by effect",fontsize=7.5)
    a.annotate(f"{YL[k]}, intergenic (promoter model)",xy=(0,1.075),xycoords="axes fraction",fontsize=9,color=KC[k],va="bottom",annotation_clip=False)
    a.annotate(f"pooled {pct(po)} ± {pse:.2g} (textbook SE), p = {pp_:.2g}; negative in {(b<0).sum()}/{NG} gardens",
               xy=(0,1.01),xycoords="axes fraction",fontsize=7,color="0.3",va="bottom",annotation_clip=False)
    if r_==1: a.set_xlabel("% change in founder fitness per carried intergenic variant (± 95% CI)",fontsize=8)
    a=fig.add_subplot(gsg[r_,1]); AXR.append(a); style(a,axis="both"); a.axhline(0,color="0.4",lw=0.9,ls="--"); a.axhline(po,color=KC[k],lw=1.1,alpha=0.7)
    a.errorbar(x,b,yerr=1.96*se,fmt="none",ecolor="0.82",elinewidth=0.8,zorder=1); a.scatter(x,b,s=34,c=x,cmap=CMAP,norm=NORM,linewidth=0.4,edgecolor="0.3",zorder=3)
    lr=stats.linregress(x,b); xx=np.linspace(x.min(),x.max(),50)
    a.plot(xx,lr.intercept+lr.slope*xx,color="0.2",lw=1.5,label=f"linear fit: {pct(lr.slope)}/°C, OLS p = {lr.pvalue:.2g}")
    a.set_ylabel("% per variant, per garden",fontsize=7.5); a.legend(fontsize=6.8,loc="lower left")
    if r_==1: a.set_xlabel("garden temperature  bio1 (°C)",fontsize=8)
    OUT.append({"intergenic term":YL[k],"pooled % per variant":pct(po),"SE (textbook)":f"{pse:.2g}","p":f"{pp_:.2g}",
                "gardens < 0":f"{(b<0).sum()}/{NG}","slope %/°C":pct(lr.slope),"slope p (OLS)":f"{lr.pvalue:.2g}"})
xl=[a.get_xlim() for a in AXL]; yl=[a.get_ylim() for a in AXR]
for a in AXL: a.set_xlim(min(v_[0] for v_ in xl),max(v_[1] for v_ in xl))
for a in AXR: a.set_ylim(min(v_[0] for v_ in yl),max(v_[1] for v_ in yl))
for a in AXL[:1]+AXR[:1]: a.tick_params(labelbottom=True)
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORM,cmap=CMAP),ax=AXL+AXR,pad=0.015,fraction=0.016,aspect=40); cb.set_label("garden temperature  bio1 (°C)",fontsize=8)
save(fig,"burden_promoter_model_by_garden_intergenic"); plt.show()
print("Supplementary Figure S5 — intergenic terms of the promoter model, per garden:"); print(pd.DataFrame(OUT).T.to_string(header=False))
"""

code_prom_fm = _promoter(code_fm)
code_prom_gd = _promoter(code_gd)
md_os = """## One-stage model — all founder × garden observations in one mixed model

Everything above is fitted on 231 founder values (the mean over gardens) or garden by garden. Here the
same question is asked of all **231 founders × 31 gardens = 7,161 observations** at once, so the overall
effect, the garden-to-garden variation and the dependence on garden temperature come out of one model.

**Data.** One row per founder *f* and garden *s*: *y*<sub>fs</sub> = log[(*H*<sub>0,f</sub> + Δ*H*<sub>s,f</sub> + 10⁻⁵)/(*H*<sub>0,f</sub> + 10⁻⁵)]
(the same fitness as above, not averaged); the founder's nine standardised burden rates (SNPs, insertions,
deletions × coding, promoter, intergenic) and assembly route; garden bio1, centred on the mean of the 31
gardens (so interaction terms are per °C and main effects are at the average garden). No missing cells.

**Models** (R formulas exactly as fitted; `kin` carries the kinship matrix, see below):

- **A — overall effects:** `y ~ <9 burden rates> + cactus + (1 | garden) + (1 | kin) + (1 | fid)`
- **B — the climate test:** A + `<tested> : bio1c` for the six coding and promoter terms + `(0 + <tested> | garden)`,
  one independent random garden slope per tested term
- **C — B without the random garden slopes** (to show what they do)

`garden` = random garden intercept; `kin` = founder genetic effect with covariance σ²<sub>kin</sub>**K**
(**K** = the 231 × 231 SNP kinship used throughout); `fid` = an independent founder effect, because each
founder is measured in 31 gardens. The random garden slope lets each garden have its own burden effect
around the overall one; without it, the 31 gardens would count as independent evidence for every burden
effect even though they share one founder set.

**How kinship enters.** lme4 cannot take a known covariance matrix, so the `(1 | kin)` design is replaced
by *Z*·*L*, with **K** = *LL*ᵀ from the eigendecomposition of **K** — lme4's documented modular interface
(`lFormula` → `mkLmerDevfun` → `optimizeLmer`; `_sv_burden_onestage.R`). This is what the lme4qtl package
does internally; the same models refitted with `lme4qtl::relmatLmer` (`_sv_burden_onestage_lme4qtl_check.R`)
give identical estimates (checked below). REML fits; Wald z-tests; for each climate interaction also a
likelihood-ratio test (ML refits of B with and without that one interaction).

*Caveats.* Frequencies within a garden compete (a founder can only rise if others fall), which none of these
models represent, so these p-values are model-based; the founder-label permutations of the two-stage
analysis remain the reference test. The independent founder effect `fid` is estimated at zero (kinship
absorbs all founder-level variation), which triggers lme4's "singular fit" note; it is harmless."""

code_os_load = r"""
OS=f"{G}/onestage"
oA=pd.read_csv(f"{OS}/coef_A.csv").set_index("term"); oB=pd.read_csv(f"{OS}/coef_B.csv").set_index("term"); oC=pd.read_csv(f"{OS}/coef_C.csv").set_index("term")
qA=pd.read_csv(f"{OS}/lme4qtl_coef_A.csv").set_index("term"); qB=pd.read_csv(f"{OS}/lme4qtl_coef_B.csv").set_index("term")
VC=pd.read_csv(f"{OS}/varcomp.csv"); LRT=pd.read_csv(f"{OS}/lrt_B.csv").set_index("term") if os.path.exists(f"{OS}/lrt_B.csv") else None
GB_=pd.read_csv(f"{OS}/garden_B.csv"); VB=pd.read_csv(f"{OS}/vcov_B.csv",index_col=0)
TS=pd.read_csv(f"{G}/sv_burden_full_model_event_CDS+promoter.csv").set_index("term")          # the two-stage fit of the same terms
LONG=pd.read_csv(f"{OS}/long.csv",usecols=["founder","garden","bio1"])
TERMS9=[f"{k}_{x}" for k in KINDS for x in ("CDS","promoter","intergenic")]; TESTED6=[f"{k}_{x}" for k in KINDS for x in ("CDS","promoter")]
print(f"data: {len(LONG):,} rows = {LONG.founder.nunique()} founders x {LONG.garden.nunique()} gardens; garden bio1 centred on {LONG.drop_duplicates('garden').bio1.mean():.2f} °C")
fx=" + ".join(TERMS9+["cactus"]); inter=" + ".join(f"{t}:bio1c" for t in TESTED6); sl=" + ".join(f"(0 + {t} | garden)" for t in TESTED6)
print(f"\nmodel A:  y ~ {fx} + (1 | garden) + (1 | kin) + (1 | fid)")
print(f"model B:  y ~ {fx} + {inter} + (1 | garden) + (1 | kin) + (1 | fid) + {sl}")
print("model C:  model B without the (0 + term | garden) slopes")
print(f"\nlme4qtl check: max |difference| in estimates  A {(oA.estimate-qA.estimate).abs().max():.1e},  B {(oB.estimate-qB.estimate).abs().max():.1e};"
      f"  in SEs  A {(oA.se-qA.se).abs().max():.1e},  B {(oB.se-qB.se).abs().max():.1e}")
print("\nvariance components, model B:"); print(VC.assign(vcov=VC.vcov.round(4)).rename(columns={"grp":"group","var1":"term","vcov":"variance"}).to_string(index=False))
"""

md_os1 = """### Figure O1 — overall effects: one model on 7,161 observations vs the founder-mean model

Model A against the two-stage full model with the same nine burden terms (`CONTEXT=CDS,promoter`), per SD
of each rate. Because every founder is measured in every garden, a random garden intercept cannot change a
founder-level effect, so the two should agree."""

code_os1 = r"""
LAB=lambda t: f"{YL[t.split('_')[0]]}, {'coding' if t.endswith('CDS') else t.split('_')[1]}" if t!="cactus" else "assembly route"
rows_=TERMS9+["cactus"]
fig,a=plt.subplots(figsize=(8.5,5.2)); style(a); a.axvline(0,color="0.4",lw=0.9,ls="--")
for yi,t in enumerate(rows_):
    col=KC[t.split("_")[0]] if t!="cactus" else "#E7298A"
    a.errorbar(TS.loc[t,"beta_sd"],yi-0.15,xerr=1.96*TS.loc[t,"se_sd"],fmt="s",ms=5,color=col,mfc="white",mew=1.2,capsize=2,label="founder means (231 values)" if yi==0 else None)
    a.errorbar(oA.loc[t,"estimate"],yi+0.15,xerr=1.96*oA.loc[t,"se"],fmt="o",ms=5.5,color=col,capsize=2,label="one-stage, model A (7,161 values)" if yi==0 else None)
    if yi in (2,5,8): a.axhline(yi+0.5,color="0.88",lw=0.8)
a.set_yticks(range(len(rows_))); a.set_yticklabels([LAB(t) for t in rows_],fontsize=8); a.set_ylim(len(rows_)-0.5,-0.7)
for tl,t in zip(a.get_yticklabels(),rows_): tl.set_color(KC[t.split("_")[0]] if t!="cactus" else "#E7298A")
a.set_xlabel("effect per SD of the rate (± 95% CI; Wald)",fontsize=8); a.legend(fontsize=7.5,loc="lower right")
save(fig,"onestage_main_effects"); plt.show()
print(pd.DataFrame({"term":[LAB(t) for t in rows_],"one-stage":[f"{oA.loc[t,'estimate']:+.3f} ± {oA.loc[t,'se']:.3f} (p {oA.loc[t,'p']:.2g})" for t in rows_],
                    "founder means":[f"{TS.loc[t,'beta_sd']:+.3f} ± {TS.loc[t,'se_sd']:.3f} (p {TS.loc[t,'p_wald']:.2g})" for t in rows_]}).to_string(index=False))
"""

md_os2 = """### Figure O2 — the climate test: burden × garden temperature

The six interaction coefficients of model B (per SD of burden per °C of garden bio1): the change in a
term's effect per degree of garden temperature. Filled = with the random garden slopes (model B, the
honest version); open = without them (model C), which shows how much the garden-to-garden variation
matters. The likelihood-ratio p (ML, model B with vs without that one interaction) is printed at the
right when available."""

code_os2 = r"""
fig,a=plt.subplots(figsize=(8.5,4.2)); style(a); a.axvline(0,color="0.4",lw=0.9,ls="--")
for yi,t in enumerate(TESTED6):
    col=KC[t.split("_")[0]]; i=f"{t}:bio1c"
    a.errorbar(oB.loc[i,"estimate"],yi-0.14,xerr=1.96*oB.loc[i,"se"],fmt="o",ms=6,color=col,capsize=2.5,elinewidth=1.4,label="with random garden slopes (B)" if yi==0 else None)
    a.errorbar(oC.loc[i,"estimate"],yi+0.14,xerr=1.96*oC.loc[i,"se"],fmt="o",ms=6,color=col,mfc="white",mew=1.2,capsize=2.5,elinewidth=1.0,label="without (C)" if yi==0 else None)
    txt=f"B: p = {oB.loc[i,'p']:.2g}"+(f";  LRT p = {LRT.loc[i,'p']:.2g}" if LRT is not None and i in LRT.index else "")+f"   |   C: p = {oC.loc[i,'p']:.2g}"
    a.annotate(txt,xy=(1.02,yi),xycoords=("axes fraction","data"),va="center",fontsize=7,color="0.3",annotation_clip=False)
    if yi in (1,3): a.axhline(yi+0.5,color="0.88",lw=0.8)
a.set_yticks(range(6)); a.set_yticklabels([LAB(t) for t in TESTED6],fontsize=8); a.set_ylim(5.6,-0.7)
for tl,t in zip(a.get_yticklabels(),TESTED6): tl.set_color(KC[t.split("_")[0]])
a.set_xlabel("burden × garden temperature: change in effect per SD of burden per °C (± 95% CI)",fontsize=8); a.legend(fontsize=7.5,loc="lower left")
save(fig,"onestage_climate_interactions"); plt.show()
if LRT is None: print("likelihood-ratio tests not yet available (lrt_B.csv) — rebuild when `_sv_burden_onestage.R` has finished")
"""

md_os3 = """### Figure O3 — each garden's effect, from the one model

For the four SV terms, each garden's total effect from model B — the overall effect plus the temperature
term plus that garden's own random slope (β + γ·bio1 + *b*<sub>s</sub>) — against garden temperature,
with the model's fitted line (β + γ·bio1) and its 95% confidence band. Bars = ± 1.96 conditional SD of the
garden's random slope. Unlike Figure 2, the gardens are not fitted separately: noisy gardens are pulled
towards the line (shrinkage). Converted to % change in fitness per carried variant with the same
per-variant scale as Figure 1 (one SD of the rate = the number of variants given in Figure 1)."""

code_os3 = r"""
SV4=[t for t in TESTED6 if not t.startswith("snp")]
fig,axs=plt.subplots(2,2,figsize=(11.5,7.6),sharex=True); B1s=LONG.drop_duplicates("garden").bio1; NORM_=mpl.colors.Normalize(B1s.min(),B1s.max()); CM=plt.get_cmap("coolwarm")
bc0=B1s.mean(); xx=np.linspace(B1s.min(),B1s.max(),60); xc=xx-bc0
for a,t,letter in zip(axs.T.ravel(),["ins_CDS","ins_promoter","del_CDS","del_promoter"],"abcd"):
    style(a,axis="both"); a.axhline(0,color="0.4",lw=0.9,ls="--"); panel(a,letter,x=-0.1)
    sc=100/TS.loc[t,"variants_per_sd"]                                          # per SD -> % per carried variant
    g=GB_[GB_.term==t]; a.errorbar(g.bio1,sc*g.total,yerr=1.96*sc*g.cond_sd,fmt="none",ecolor="0.8",elinewidth=0.9,zorder=1)
    a.scatter(g.bio1,sc*g.total,c=g.bio1,cmap=CM,norm=NORM_,s=36,linewidth=0.4,edgecolor="0.3",zorder=3)
    b0,b1=oB.loc[t,"estimate"],oB.loc[f"{t}:bio1c","estimate"]; i0,i1=t,f"{t}:bio1c"
    v=VB.loc[i0,i0]+xc**2*VB.loc[i1,i1]+2*xc*VB.loc[i0,i1]; fit=b0+b1*xc
    a.fill_between(xx,sc*(fit-1.96*np.sqrt(v)),sc*(fit+1.96*np.sqrt(v)),color=KC[t.split("_")[0]],alpha=0.15,lw=0)
    a.plot(xx,sc*fit,color=KC[t.split("_")[0]],lw=1.8)
    lr=f";  LRT p = {LRT.loc[i1,'p']:.2g}" if LRT is not None and i1 in LRT.index else ""
    a.annotate(f"{LAB(t)}\nslope {pct(sc*b1)}/°C, p = {oB.loc[i1,'p']:.2g}{lr}",xy=(0.02,0.97),xycoords="axes fraction",va="top",fontsize=7.5,color=KC[t.split("_")[0]])
    a.set_ylabel("% per carried variant, per garden",fontsize=7.5)
    if letter in "bd": a.set_xlabel("garden temperature  bio1 (°C)",fontsize=8)
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORM_,cmap=CM),ax=axs,pad=0.015,fraction=0.02); cb.set_label("garden temperature  bio1 (°C)",fontsize=8)
save(fig,"onestage_garden_effects"); plt.show()
T_=pd.DataFrame([{"term":LAB(t),"overall (per SD)":f"{oB.loc[t,'estimate']:+.3f} ± {oB.loc[t,'se']:.3f}","× bio1 (per SD per °C)":f"{oB.loc[t+':bio1c','estimate']:+.4f} ± {oB.loc[t+':bio1c','se']:.4f}",
                  "Wald p":f"{oB.loc[t+':bio1c','p']:.2g}","LRT p":(f"{LRT.loc[t+':bio1c','p']:.2g}" if LRT is not None and t+':bio1c' in LRT.index else "pending"),
                  "random-slope SD":f"{np.sqrt(VC[(VC.var1==t)].vcov.iloc[0]):.3f}"} for t in TESTED6])
print("model B, tested terms:"); print(T_.to_string(index=False))
"""
md_final = """## Final models — fitted with lme4qtl

Two mixed models on all **7,161 founder × garden observations** (231 founders × 31 gardens), fitted with
`lme4qtl::relmatLmer` (Ziyatdinov et al. 2018, *BMC Bioinformatics* 19:68; built on lme4, Bates et al.
2015, *J. Stat. Softw.* 67), kinship passed as `relmat` (`_sv_burden_final_models.R`, cluster array
`run_final_models.sbatch`).

**Model 1 — is each class good or bad for fitness overall?**

```
y ~ snp_CDS + ins_CDS + ins_promoter + ins_intergenic + del_CDS + del_promoter + del_intergenic + cactus
    + (1 | garden) + (1 | kin) + (1 | fid)
    + (0 + ins_CDS | garden) + (0 + ins_promoter | garden) + (0 + del_CDS | garden) + (0 + del_promoter | garden)
```

**Model 2 — does the effect depend on garden temperature?** Model 1 + `ins_CDS:bio1c + ins_promoter:bio1c + del_CDS:bio1c + del_promoter:bio1c`.

*y* = founder log-fitness in the garden, log[(*H*<sub>0</sub> + Δ*H*)/(*H*<sub>0</sub>)] (+10⁻⁵ in both);
burden terms = standardised shared-variant rates, reported here as % change in fitness **per carried
variant** (one SD = the number of variants in Figure 1); `snp_CDS` = coding SNP load (control for
genomic divergence); intergenic terms = within-class baselines; `cactus` = assembly route (technical);
`bio1c` = garden bio1 centred on the mean of the 31 gardens (main effects in Model 2 are at the average
garden). Random effects: `garden` intercept; `kin` = founder genetic effect with covariance σ²**K** (SNP
kinship); `fid` = independent founder effect; `(0 + term | garden)` = each garden's own effect of the four
SV terms of interest, so the overall effects generalise across gardens and the temperature interactions
are tested against real between-garden variation. REML fits, Wald z-tests; for the interactions also a
likelihood-ratio test (ML refit of Model 2 with vs without that one interaction), Bonferroni-corrected for
the four tests. The controls (`snp_CDS`, intergenic terms, `cactus`) are shown for completeness, not
interpreted."""

code_final_load = r"""
FD=f"{G}/final_models"
F1=pd.read_csv(f"{FD}/coef_M1_reml.csv").set_index("term"); F2=pd.read_csv(f"{FD}/coef_M2_reml.csv").set_index("term")
VF1=pd.read_csv(f"{FD}/varcomp_M1_reml.csv"); VF2=pd.read_csv(f"{FD}/varcomp_M2_reml.csv")
GF2=pd.read_csv(f"{FD}/garden_M2_reml.csv"); VV2=pd.read_csv(f"{FD}/vcov_M2_reml.csv",index_col=0)
LL=pd.concat([pd.read_csv(f) for f in sorted(__import__("glob").glob(f"{FD}/loglik_*.csv"))]).set_index("model")
SVT=["ins_CDS","ins_promoter","del_CDS","del_promoter"]; ROWS7=["snp_CDS","ins_CDS","ins_promoter","ins_intergenic","del_CDS","del_promoter","del_intergenic"]
VPS=pd.read_csv(f"{G}/sv_burden_full_model_event_CDS+promoter.csv").set_index("term").variants_per_sd     # one SD of each rate, in variants
PV=lambda t: 100/VPS[t.split(":")[0]]                                                               # per SD -> % per carried variant
LRTF=pd.DataFrame([{"term":f"{t}:bio1c","chisq":2*(LL.loc["M2_ml","logLik"]-LL.loc[f"M2_ml_drop_{t}","logLik"])} for t in SVT]).set_index("term")
LRTF["p"]=stats.chi2.sf(LRTF.chisq,1); LRTF["p_bonferroni"]=np.minimum(LRTF.p*4,1)
LABF=lambda t: {"snp":"SNPs","ins":"insertions","del":"deletions"}[t.split("_")[0]]+", "+{"CDS":"coding","promoter":"promoter","intergenic":"intergenic"}[t.split("_")[1].split(":")[0]]
MKF={"CDS":"o","promoter":"^","intergenic":"s"}
print(f"{len(LL)} fits; singular (founder term at zero) in: {', '.join(LL.index[LL.singular])}")
print("\nvariance components, Model 2:"); print(VF2.assign(vcov=VF2.vcov.round(4)).rename(columns={'grp':'group','var1':'term','vcov':'variance'}).to_string(index=False))
"""

md_final1 = """### Figure M1 — Model 1: overall effect of every term

Every fixed effect of Model 1 on one axis, % change in founder fitness per carried variant (± 95% CI).
Circles = coding, triangles = promoter, squares = intergenic; filled = p < 0.05 (Wald). Assembly route
is per group, not per variant, and is printed below the figure."""

code_final1 = r"""
def coef_forest(a,C,rows,interaction=False):
    style(a); a.axvline(0,color="0.4",lw=0.9,ls="--")
    for yi,t in enumerate(rows):
        k=t.split("_")[0]; ctx=t.split("_")[1].split(":")[0]; sc=PV(t); b=sc*C.loc[t,"estimate"]; se=sc*C.loc[t,"se"]; p=C.loc[t,"p"]
        a.errorbar(b,yi,xerr=1.96*se,fmt=MKF[ctx],ms=7,color=KC[k],mfc=KC[k] if p<0.05 else "white",mew=1.4,capsize=3,elinewidth=1.5)
        lab=f"{pct(b)}{'/°C' if interaction else ''}   p = {p:.2g}"
        if interaction: lab+=f"   (LRT {LRTF.loc[t,'p']:.2g}; Bonferroni {LRTF.loc[t,'p_bonferroni']:.2g})"
        a.annotate(lab,xy=(1.02,yi),xycoords=("axes fraction","data"),va="center",fontsize=7,color="0.3",annotation_clip=False)
        if yi<len(rows)-1 and rows[yi+1].split("_")[0]!=k: a.axhline(yi+0.5,color="0.88",lw=0.8)
    a.set_yticks(range(len(rows))); a.set_yticklabels([LABF(t)+(" × temperature" if interaction else "") for t in rows],fontsize=8)
    for tl,t in zip(a.get_yticklabels(),rows): tl.set_color(KC[t.split("_")[0]])
    a.set_ylim(len(rows)-0.5,-0.6)
fig,a=plt.subplots(figsize=(8,4.6)); coef_forest(a,F1,ROWS7)
a.set_xlabel("% change in founder fitness per carried variant (± 95% CI) — Model 1",fontsize=8)
a.errorbar([],[],fmt="o",color="0.3",ms=6,label="p < 0.05"); a.errorbar([],[],fmt="o",color="0.3",mfc="white",ms=6,label="n.s."); a.legend(fontsize=7,loc="lower left")
save(fig,"final_model1_coefficients"); plt.show()
print(f"assembly route (cactus − short-read): {F1.loc['cactus','estimate']:+.3f} ± {F1.loc['cactus','se']:.3f} log-fitness, p = {F1.loc['cactus','p']:.2g}")
print(pd.DataFrame({"term":[LABF(t) for t in ROWS7],"per SD":[f"{F1.loc[t,'estimate']:+.3f} ± {F1.loc[t,'se']:.3f}" for t in ROWS7],
                    "% per variant":[f"{pct(PV(t)*F1.loc[t,'estimate'])} [{pct(PV(t)*(F1.loc[t,'estimate']-1.96*F1.loc[t,'se']))}, {pct(PV(t)*(F1.loc[t,'estimate']+1.96*F1.loc[t,'se']))}]" for t in ROWS7],
                    "p":[f"{F1.loc[t,'p']:.2g}" for t in ROWS7]}).to_string(index=False))
"""

md_final2 = """### Figure M2 — Model 2: overall effects at the average garden, and the temperature interactions

**(a)** The same terms as Figure M1, now from Model 2 (effects at the average garden temperature).
**(b)** The four burden × garden temperature interactions: the change in % fitness per carried variant
per °C of garden bio1, with the Wald p and, at the right, the likelihood-ratio p and its Bonferroni
correction for the four tests."""

code_final2 = r"""
fig,ax=plt.subplots(2,1,figsize=(8,7.2),gridspec_kw=dict(height_ratios=[7,4],hspace=0.45))
coef_forest(ax[0],F2,ROWS7); panel(ax[0],"a",x=-0.3); ax[0].set_xlabel("% change in founder fitness per carried variant (± 95% CI) — Model 2, at the average garden",fontsize=8)
coef_forest(ax[1],F2,[f"{t}:bio1c" for t in SVT],interaction=True); panel(ax[1],"b",x=-0.3)
ax[1].set_xlabel("change in % fitness per carried variant, per °C of garden temperature (± 95% CI)",fontsize=8)
save(fig,"final_model2_coefficients"); plt.show()
print("temperature interactions, Model 2:")
print(pd.DataFrame([{"term":LABF(t)+" × temperature","per SD per °C":f"{F2.loc[t+':bio1c','estimate']:+.4f} ± {F2.loc[t+':bio1c','se']:.4f}",
                     "% per variant per °C":pct(PV(t)*F2.loc[t+':bio1c','estimate']),"Wald p":f"{F2.loc[t+':bio1c','p']:.2g}",
                     "LRT chi²":f"{LRTF.loc[t+':bio1c','chisq']:.2f}","LRT p":f"{LRTF.loc[t+':bio1c','p']:.2g}","Bonferroni (×4)":f"{LRTF.loc[t+':bio1c','p_bonferroni']:.2g}"} for t in SVT]).to_string(index=False))
"""

md_final3 = """### Figure M3 — Model 2, garden by garden, for each SV term

For each of the four SV terms, each garden's effect from Model 2 — the overall effect, plus the
temperature term, plus that garden's own random slope — in % per carried variant. **Left:** gardens
sorted by effect, coloured by garden temperature, bars = ± 1.96 conditional SD of the garden's slope;
line and band = the term's effect at the average garden ± 95% CI. **Right:** the same against garden
temperature, with Model 2's fitted line and its 95% confidence band. All panels on one shared scale.
Because gardens are estimated within one model, noisy gardens are pulled towards the fitted line."""

code_final3 = r"""
B1s=GF2.drop_duplicates("garden").bio1; NORMF=mpl.colors.Normalize(B1s.min(),B1s.max()); CMF=plt.get_cmap("coolwarm"); bc0=B1s.mean()
fig=plt.figure(figsize=(13.5,16)); gsf=fig.add_gridspec(4,2,width_ratios=[1.3,1],wspace=0.25,hspace=0.45); AXL=[]; AXR=[]
for r_,t in enumerate(SVT):
    k=t.split("_")[0]; sc=PV(t); g=GF2[GF2.term==t].sort_values("total").reset_index(drop=True)
    b=sc*g.total.to_numpy(); e=1.96*sc*g.cond_sd.to_numpy(); x=g.bio1.to_numpy()
    b0=sc*F2.loc[t,"estimate"]; s0=sc*F2.loc[t,"se"]; i0,i1=t,f"{t}:bio1c"
    a=fig.add_subplot(gsf[r_,0]); AXL.append(a); style(a); panel(a,"abcd"[r_],x=-0.07); a.axvline(0,color="0.4",lw=0.9,ls="--")
    a.axvspan(b0-1.96*s0,b0+1.96*s0,color=KC[k],alpha=0.14,lw=0); a.axvline(b0,color=KC[k],lw=1.4)
    for yi in range(len(g)):
        c_=CMF(NORMF(x[yi])); a.errorbar(b[yi],yi,xerr=e[yi],fmt="o",ms=4.2,color=c_,mec="0.25",mew=0.5,ecolor=c_,elinewidth=1.0,capsize=1.6)
    a.set_yticks(range(len(g))); a.set_yticklabels(g.garden.astype(int),fontsize=5.5); a.set_ylim(-1,len(g)); a.set_ylabel("garden (site), sorted by effect",fontsize=7.5)
    a.annotate(LABF(t),xy=(0,1.075),xycoords="axes fraction",fontsize=9,color=KC[k],va="bottom",annotation_clip=False)
    a.annotate(f"at the average garden {pct(b0)} ± {1.96*s0:.2g} (95% CI), p = {F2.loc[t,'p']:.2g}; negative in {(b<0).sum()}/{len(g)} gardens",
               xy=(0,1.01),xycoords="axes fraction",fontsize=7,color="0.3",va="bottom",annotation_clip=False)
    if r_==3: a.set_xlabel("% change in founder fitness per carried variant, per garden (Model 2)",fontsize=8)
    a=fig.add_subplot(gsf[r_,1]); AXR.append(a); style(a,axis="both"); a.axhline(0,color="0.4",lw=0.9,ls="--")
    a.errorbar(x,b,yerr=e,fmt="none",ecolor="0.82",elinewidth=0.8,zorder=1); a.scatter(x,b,s=34,c=x,cmap=CMF,norm=NORMF,linewidth=0.4,edgecolor="0.3",zorder=3)
    xx=np.linspace(x.min(),x.max(),60); xc=xx-bc0; fit=F2.loc[i0,"estimate"]+F2.loc[i1,"estimate"]*xc
    v=VV2.loc[i0,i0]+xc**2*VV2.loc[i1,i1]+2*xc*VV2.loc[i0,i1]
    a.fill_between(xx,sc*(fit-1.96*np.sqrt(v)),sc*(fit+1.96*np.sqrt(v)),color=KC[k],alpha=0.15,lw=0); a.plot(xx,sc*fit,color=KC[k],lw=1.8)
    a.annotate(f"slope {pct(sc*F2.loc[i1,'estimate'])}/°C;  LRT p = {LRTF.loc[i1,'p']:.2g} (Bonferroni {LRTF.loc[i1,'p_bonferroni']:.2g})",
               xy=(0.02,0.03),xycoords="axes fraction",fontsize=7.2,color="0.25")
    a.set_ylabel("% per variant, per garden",fontsize=7.5)
    if r_==3: a.set_xlabel("garden temperature  bio1 (°C)",fontsize=8)
xl=[a.get_xlim() for a in AXL]; yl=[a.get_ylim() for a in AXR]
for a in AXL: a.set_xlim(min(v_[0] for v_ in xl),max(v_[1] for v_ in xl))
for a in AXR: a.set_ylim(min(v_[0] for v_ in yl),max(v_[1] for v_ in yl))
for a in AXL[:3]+AXR[:3]: a.tick_params(labelbottom=True)
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORMF,cmap=CMF),ax=AXL+AXR,pad=0.015,fraction=0.016,aspect=50); cb.set_label("garden temperature  bio1 (°C)",fontsize=8)
save(fig,"final_model2_by_garden"); plt.show()
"""

md_methods_fm = """## Methods, as run

- **Fitness.** Founder frequencies *H* per pool from kMate (`founder_site_dh.npz`): *H*<sub>0</sub> = mean of
  the 8 seed-mix estimates; *H*<sub>0</sub> + Δ*H*<sub>s</sub> = last-generation frequency in garden *s*
  (plots equally weighted); *y* as in Figure 1, ε = 10⁻⁵.
- **Burden rates** (`_sv_burden_build.py` → `founder_burden_event.csv`): per founder and class, variants
  carried ÷ variants called; shared variants, ≥ 2 carriers, called in ≥ 90% of founders; all SNP records;
  SVs as Truvari events.
- **Model** (`_sv_burden_full_model.py`): linear mixed model, REML by eigendecomposition of the kinship matrix
  (standardised GRM from a 2% draw of SNPs with ≥ 5 carriers; `LMM` class of `_sv_burden_lmm.py`); rates
  standardised for the fit; Wald z-tests. Per-variant coefficient = coefficient per SD ÷ (SD of the rate ×
  median number of that class's variants called per founder).
- **Permutation** (3,000 draws per coding term): the term's rate regressed (OLS) on all other terms, the
  residual shuffled across founders, the same shuffled vector used in the mean-fitness fit and in all 31
  per-garden fits; two-sided p = share of draws at least as extreme. Pooled per-garden coefficient =
  inverse-variance mean, with SE and p from the permutation spread; slope on garden temperature = OLS of
  the 31 per-garden coefficients on garden bio1 (equal weights), p from the same permutation.
- **Supplement** (`_sv_burden_model_simplify.py`, joint models): AIC from maximum-likelihood refits (REML
  likelihoods are not comparable across fixed effects)."""


if __name__ == "__main__":
    build("sv_climate_maf.ipynb", [
        new_markdown_cell(md_title_fm), new_code_cell(code_setup_fm),
        new_markdown_cell(md_fm), new_code_cell(code_fm),
        new_markdown_cell(md_gd), new_code_cell(code_gd),
        new_markdown_cell(md_sup), new_code_cell(code_sup),
        new_markdown_cell(md_key), new_code_cell(code_key),
        new_markdown_cell(md_prom), new_code_cell(code_prom_load), new_code_cell(code_prom_fm), new_code_cell(code_prom_gd),
        new_markdown_cell(md_prom_ig), new_code_cell(code_prom_ig),
        new_markdown_cell(md_os), new_code_cell(code_os_load), new_markdown_cell(md_os1), new_code_cell(code_os1),
        new_markdown_cell(md_os2), new_code_cell(code_os2), new_markdown_cell(md_os3), new_code_cell(code_os3),
        new_markdown_cell(md_methods_fm),
    ])
