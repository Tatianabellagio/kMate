#!/usr/bin/env python
"""Build+execute `sv_burden_model.ipynb` -- THE founder-burden notebook.

ONE model, reported on its own. No specification comparisons: every alternative arm ever fitted is in
`results/sv_adaptive/model_specs.csv` / `model_registry.csv`, which is where they belong. Earlier notebooks
are in `archive/burden_superseded_2026-09-29/`.

  y ~ snp_CDS + sv_CDS + sv_TE + sv_promoter + sv_intergenic + cactus
      + sv_CDS:bio1c + sv_TE:bio1c + sv_promoter:bio1c
      + obio1c + obio1c:bio1c
      + (1 | garden) + (1 | kin) + (0 + obio1c | garden)
      + (0 + sv_CDS | garden) + (0 + sv_TE | garden) + (0 + sv_promoter | garden)

SVs pooled over insertions and deletions (that split is defined against Col-0, not biology), no origin
filter, reference TE kept as its own context, founder origin climate and its garden-climate interaction
always in. The independent founder intercept `fid` is pruned (RE_DROP=fid): it estimated 4.9e-07, so
removing it costs one parameter, leaves no variance component at the boundary, and changes the
log-likelihood and every coefficient and standard error by nothing to five decimals.

Diagnostic sections reuse the tested code from `_build_sv_burden_influence_nb.py`, extracted by
`_extract_sections.py` into `_sv_burden_kept_sections.json` so this builder stays self-contained.
`basic` env.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_fig_common import build, new_markdown_cell, new_code_cell, code_setup_burden as code_setup

_H = os.path.dirname(os.path.abspath(__file__))
_K = json.load(open(os.path.join(_H, "_sv_burden_kept_sections.json")))
_ns = {}
exec(open(os.path.join(_H, "_sv_burden_strata_sections.py")).read(), _ns)

md_title = """# The founder SV-burden model

**One model, reported on its own.** Alternative specifications — the `shared` origin filter, dropping the
reference-TE term, omitting the local-adaptation control, the split insertion/deletion parameterisation,
single-technology panels — were all fitted, and are indexed in `results/sv_adaptive/model_specs.csv` and
`model_registry.csv`. None is shown here: carrying retired arms through the figures made the result harder
to read, not easier.

## The model, exactly as run

```r
y ~ snp_CDS + sv_CDS + sv_TE + sv_promoter + sv_intergenic + cactus
    + sv_CDS:bio1c + sv_TE:bio1c + sv_promoter:bio1c
    + obio1c + obio1c:bio1c
    + (1 | garden) + (1 | kin) + (0 + obio1c | garden)
    + (0 + sv_CDS | garden) + (0 + sv_TE | garden) + (0 + sv_promoter | garden)
```

Fitted with `lme4qtl::relmatLmer`, kinship passed as `relmat`, on **all 231 founders × 31 gardens =
7,161 rows** — 151 short-read plus 80 cactus-assembly founders, with assembly route entering as the
covariate `cactus` rather than as a subset.

**Burdens.** Truvari-collapsed SV **events**, insertions and deletions pooled, one genomic context each by
first match (CDS > UTR > intron > reference TE > promoter ≤ 2 kb > intergenic). No origin filter — every
variant passing call rate ≥ 90 % and ≥ 2 carriers enters:

```
sv_<ctx>  =  z( (ins carried + del carried) / (ins called + del called) )
```

| symbol | meaning |
|---|---|
| `bio1c` | the **garden's** mean annual temperature, centred on the mean of the 31 gardens |
| `obio1c` | the **founder's origin** bio1, centred over the 231 founders |
| `obio1c:bio1c` | **local adaptation**: origin climate × garden climate, containing no SV information |
| `snp_CDS` | matching SNP coding burden · `cactus` = assembly route (0/1) |
| `kin` | founder random effect with covariance σ²K, K = SNP kinship |
| `(0 + sv_ctx \\| garden)` | independent random garden slope per modelled burden |
| `(0 + obio1c \\| garden)` | local adaptation gets its own strength in each garden |

Insertions and deletions are **pooled** because which label a variant receives depends on the Col-0
reference and is not a biological distinction — the same reason the site-frequency spectra in Figure 1 are
folded.

The independent founder intercept `fid` is **pruned**. In the unpruned fit its variance was 4.9 × 10⁻⁷ —
it carried nothing, all founder-level variation being absorbed by kinship — so removing it costs one
parameter and leaves no component at the boundary. The two fits agree in log-likelihood (−11501.0028) and
in every coefficient and standard error to five decimal places.

## Reading the coefficients

`bio1c` is centred, so each burden splits into two numbers that are read together:

- **`sv_ctx`** — the burden's effect **at a garden of average temperature** (13.5 °C): *are SVs in this
  context deleterious?*
- **`sv_ctx:bio1c`** — how that effect **changes per °C** of garden temperature.

Both come out of this one fit; no separate main-effects model is needed or used.

⚠ `obio1c:bio1c` is itself an interaction, so it is near-orthogonal to the burden main effects: it tests
the `sv_ctx:bio1c` terms hard and the `sv_ctx` terms only mildly. §13 returns to this."""

code_load = r"""
PA=f"{G}/final_models_pooled_all"
SUF="_x_obio1cobio1cbio1c_re0obio1cgarden_refid"           # THE model: + local adaptation per garden, fid pruned
SCS="_x_obio1cobio1cbio1c"                                 # RE_DROP does not change the burdens, so the
                                                           # per-variant scale file carries no _refid tag
SVT=["sv_CDS","sv_TE","sv_promoter"]                       # burdens carrying a temperature interaction
BURD=["snp_CDS","sv_CDS","sv_TE","sv_promoter","sv_intergenic"]
PAL2={"snp_CDS":"0.25","sv_CDS":"#4A3B6B","sv_TE":"#6B5A8E","sv_promoter":"#9084AE",
      "sv_intergenic":"#C6BCD8","cactus":"0.6","obio1c":"#8C6D1F"}
PC2=lambda t: PAL2[("obio1c" if "obio1c" in t else t.split(":")[0])]
LAB2=lambda t: {"snp_CDS":"SNPs, coding","sv_CDS":"SVs, coding","sv_TE":"SVs, reference TE",
                "sv_promoter":"SVs, promoter","sv_intergenic":"SVs, intergenic","cactus":"assembly route",
                "obio1c":"founder origin climate"}[("obio1c" if "obio1c" in t else t.split(":")[0])]
CO=pd.read_csv(f"{PA}/coef_M2p_reml{SUF}.csv").set_index("term")
VC=pd.read_csv(f"{PA}/varcomp_M2p_reml{SUF}.csv")
GD=pd.read_csv(f"{PA}/garden_M2p_reml{SUF}.csv")
VV=pd.read_csv(f"{PA}/vcov_M2p_reml{SUF}.csv",index_col=0)
LLR=pd.read_csv(f"{PA}/loglik_M2p_reml{SUF}.csv").iloc[0]
# LRTs from the NAMED loglik files -- this directory also holds SUBSET=, RE_DROP= and other EXTRA_FIXED
# runs, and a glob would silently mix likelihoods from models fitted on different data
ll=lambda n: pd.read_csv(f"{PA}/loglik_{n}{SUF}.csv").logLik[0]
LRT={t: stats.chi2.sf(2*(ll("M2p_ml")-ll(f"M2p_ml_drop_{t}")),1) for t in SVT}
OBX=next(t for t in CO.index if set(t.split(":"))=={"obio1c","bio1c"})   # lme4 spells it bio1c:obio1c

SCL=pd.read_csv(f"{PA}/pooled_rate_scale{SCS}.csv").set_index("context")
VPS={f"sv_{c}": SCL.loc[c,"sd_rate"]*SCL.loc[c,"med_called"] for c in SCL.index}
VPS["snp_CDS"]=pd.read_csv(f"{G}/sv_burden_full_model_event_CDS+promoter.csv").set_index("term").variants_per_sd["snp_CDS"]
PV=lambda t: 100/VPS[t.split(":")[0]]                      # per SD -> % per carried variant (display only)

FB=pd.read_csv(f"{G}/founder_burden_event.csv"); sd_ob=FB.bio1.std(ddof=1)
gbn=pd.read_csv(f"{G}/onestage_all/long.csv").drop_duplicates("garden").bio1
cac=(FB.is_cactus.astype(str)=="True")
print(f"founders {len(FB)} = {int((~cac).sum())} short-read + {int(cac.sum())} cactus   "
      f"gardens {len(gbn)}   rows {len(FB)*len(gbn):,}")
print(f"origin bio1 {FB.bio1.min():.1f} to {FB.bio1.max():.1f} degC (SD {sd_ob:.2f});  "
      f"garden bio1 {gbn.min():.1f} to {gbn.max():.1f} degC (mean {gbn.mean():.2f} = where bio1c is zero)")
print(f"\nlogLik {LLR.logLik:.4f}   df {int(LLR.df)}   singular {LLR.singular}")
print("\nvariance components (none at the boundary):"); print(VC.to_string(index=False))
print("\nvariants per SD of each burden (the per-carried-variant divisor):")
for c in SCL.index: print(f"  sv_{c:<11} {VPS['sv_'+c]:8.1f} variants")
"""

md_counts = """## 1 · What enters the model

Every truvari-collapsed SV event passing call rate ≥ 90 % and ≥ 2 carriers, split by class so the quantity
the pooling merges is visible: the burden for a context is the green plus the orange segment.

`UTR` and `intron` carry no term in this model; they are shown so the context accounting is complete."""

code_counts = r"""
CAT=pd.read_csv(f"{G}/sv_event_catalog.csv"); AP=CAT[CAT.passes]
CTX6=["CDS","UTR","intron","TE","promoter","intergenic"]; MODC={"CDS","TE","promoter","intergenic"}
CI,CD="#2F4A1F","#B04A00"
fig,a=plt.subplots(figsize=(7.6,3.4)); style(a)
for i,c in enumerate(CTX6):
    ni=int(((AP.context==c)&(AP.kind=="ins")).sum()); nd=int(((AP.context==c)&(AP.kind=="del")).sum())
    a.barh(i,ni,height=0.55,color=CI,lw=0); a.barh(i,nd,height=0.55,left=ni,color=CD,lw=0)
    a.annotate(f"{ni+nd:,}",xy=(ni+nd,i),xytext=(5,0),textcoords="offset points",va="center",fontsize=7.5)
    if i<len(CTX6)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(CTX6)))
a.set_yticklabels([c+("" if c in MODC else "  (not a model term)") for c in CTX6],fontsize=8.5)
a.invert_yaxis(); a.set_xlabel("SV events",fontsize=8.5); a.set_xlim(0,8200)
a.legend([plt.Rectangle((0,0),1,1,color=CI),plt.Rectangle((0,0),1,1,color=CD)],
         ["insertions","deletions"],fontsize=8,loc="lower right",handlelength=1.2)
save(fig,"model_inputs",nb=True); plt.show()
T=pd.DataFrame({"insertions":[int(((AP.context==c)&(AP.kind=="ins")).sum()) for c in CTX6],
                "deletions":[int(((AP.context==c)&(AP.kind=="del")).sum()) for c in CTX6]},index=CTX6)
T["pooled total"]=T.insertions+T.deletions; T["in the model"]=[c in MODC for c in CTX6]
print("SV events entering the pooled burdens, by context:"); print(T.to_string())
print(f"\nTOTAL {int(T['pooled total'].sum()):,} events")
"""

md_model = """## 2 · The burden effects

The effect of each burden **at a garden of average temperature**, in % fitness per carried variant
(± 95 % CI). `bio1c` is centred, so this is what the `sv_ctx` coefficients mean; how each effect changes
with temperature is §3.

Filled markers are p < 0.05."""

code_model = r"""
fig,a=plt.subplots(figsize=(7.4,3.2)); style(a); a.axvline(0,color="0.35",lw=1,ls="--")
for i,t in enumerate(BURD):
    e,se=CO.loc[t,"estimate"]*PV(t), CO.loc[t,"se"]*PV(t)
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=PC2(t),mfc=PC2(t) if CO.loc[t,"p"]<0.05 else "white",
               mew=1.2,ms=7,capsize=2.5,zorder=3)
    if i<len(BURD)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(BURD))); a.set_yticklabels([LAB2(t) for t in BURD],fontsize=8.5); a.invert_yaxis()
x0,x1=a.get_xlim(); a.set_xlim(x0,x1+0.32*(x1-x0))
for i,t in enumerate(BURD):
    a.annotate(f"p = {CO.loc[t,'p']:.3g}",xy=(0.985,i),xycoords=("axes fraction","data"),
               ha="right",va="center",fontsize=7,color="0.35")
a.set_xlabel("% fitness per carried variant, at an average-temperature garden  (± 95% CI)",fontsize=8.5)
save(fig,"model_burden_effects",nb=True); plt.show()
print("BURDEN EFFECTS at an average-temperature garden, % fitness per carried variant:\n")
print(pd.DataFrame([{"term":LAB2(t),"%/variant":f"{CO.loc[t,'estimate']*PV(t):+.4f}",
    "95% CI":f"[{(CO.loc[t,'estimate']-1.96*CO.loc[t,'se'])*PV(t):+.3f}, "
             f"{(CO.loc[t,'estimate']+1.96*CO.loc[t,'se'])*PV(t):+.3f}]",
    "per SD":f"{CO.loc[t,'estimate']:+.4f} ± {CO.loc[t,'se']:.4f}",
    "Wald p":f"{CO.loc[t,'p']:.4g}"} for t in BURD]).to_string(index=False))
print(f"\nnuisance terms: intercept {CO.loc['(Intercept)','estimate']:+.4f}, "
      f"cactus {CO.loc['cactus','estimate']:+.4f} (p = {CO.loc['cactus','p']:.3g})")
"""

md_inter = """## 3 · The garden-temperature interactions

How each burden's effect **changes per °C** of garden temperature, in % per carried variant per °C
(± 95 % CI). Significance is the likelihood-ratio test — ML refits with and without that single
interaction, 1 df. Three interactions are fitted; p-values are reported as they come, uncorrected.

Filled markers are p < 0.05."""

code_inter = r'''
fig,a=plt.subplots(figsize=(7.4,2.9)); style(a); a.axvline(0,color="0.35",lw=1,ls="--")
for i,t in enumerate(SVT):
    n=f"{t}:bio1c"; e,se=CO.loc[n,"estimate"]*PV(t), CO.loc[n,"se"]*PV(t)
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=PC2(t),mfc=PC2(t) if LRT[t]<0.05 else "white",
               mew=1.2,ms=7,capsize=2.5,zorder=3)
    if i<len(SVT)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(SVT))); a.set_yticklabels([LAB2(t) for t in SVT],fontsize=8.5); a.invert_yaxis()
a.set_ylim(len(SVT)-0.5,-0.5)
x0,x1=a.get_xlim(); a.set_xlim(x0,x1+0.34*(x1-x0))
for i,t in enumerate(SVT):
    a.annotate(f"LRT p = {LRT[t]:.3g}",xy=(0.985,i),xycoords=("axes fraction","data"),
               ha="right",va="center",fontsize=7,color="0.35")
a.set_xlabel("change in % fitness per carried variant,\nper °C of garden temperature  (± 95% CI)",fontsize=8.5)
save(fig,"model_interactions",nb=True); plt.show()
print("BURDEN x GARDEN TEMPERATURE, % per carried variant per degC:\n")
print(pd.DataFrame([{"term":LAB2(t),"%/variant/degC":f"{CO.loc[t+':bio1c','estimate']*PV(t):+.5f}",
    "95% CI":f"[{(CO.loc[t+':bio1c','estimate']-1.96*CO.loc[t+':bio1c','se'])*PV(t):+.4f}, "
             f"{(CO.loc[t+':bio1c','estimate']+1.96*CO.loc[t+':bio1c','se'])*PV(t):+.4f}]",
    "Wald p":f"{CO.loc[t+':bio1c','p']:.3g}","LRT p":f"{LRT[t]:.3g}",
    "p < 0.05":LRT[t]<0.05} for t in SVT]).to_string(index=False))
print(f"\n  {sum(LRT[t]<0.05 for t in SVT)} of {len(SVT)} interactions reach p < 0.05.")
'''

md_scale = (_K["md_scale"].replace("## 4 ·", "## 4 ·")
            .replace("All values are from the FULL fit.", "All values are from the fitted model."))

code_scale = r"""
TT=[(OBX,sd_ob)]+[(f"{t}:bio1c",1.0) for t in SVT]
fig,a=plt.subplots(figsize=(6.8,2.7)); style(a); a.axvline(0,color="0.35",lw=1,ls="--")
for i,(n,k) in enumerate(TT):
    e,se=CO.loc[n,"estimate"]*k, CO.loc[n,"se"]*k
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=PC2(n),ms=6.5,mew=1.2,capsize=2,zorder=3)
    # annotate the same statistic section 3 reports -- the likelihood-ratio p -- so the same term cannot
    # appear to differ between sections. The local-adaptation term has no LRT refit, so it carries its Wald p.
    lab=(f"LRT p = {LRT[n.split(':')[0]]:.2g}" if n.split(":")[0] in LRT
         else f"Wald p = {CO.loc[n,'p']:.2g}")
    a.annotate(f"{e:+.4f}   {lab}",xy=(1.03,i),xycoords=("axes fraction","data"),
               va="center",fontsize=7.2,color="0.25",annotation_clip=False)
    if i<len(TT)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(TT))); a.set_yticklabels([LAB2(n)+" × garden temp." for n,_ in TT],fontsize=8)
a.invert_yaxis()
a.set_xlabel("change in log-fitness per SD of the predictor,\nper °C of garden temperature  (± 95% CI)",fontsize=7.5)
save(fig,"model_local_adaptation_scale",nb=True); plt.show()
ob=CO.loc[OBX,"estimate"]*sd_ob
print("all x temperature terms on one scale (per SD of predictor, per degC of garden):")
for n,k in TT:
    e=CO.loc[n,"estimate"]*k
    print(f"  {LAB2(n)+' x garden temp.':<42} {e:+.5f}   ({abs(e/ob):.2f}x the local-adaptation term)")
"""

md_garden = """## 5 · The fitted effect in each garden

Each garden's own total effect of the coding-SV burden, **from the fitted model** — the fixed part, plus
the temperature term, plus that garden's random slope — in % per carried variant. Every covariate is in:
kinship, the SNP coding burden, assembly route, origin climate and the other burdens. These are not
separate per-garden regressions. **Left:** gardens sorted by effect, coloured by
garden temperature; bars are ± 1.96 conditional SD of the slope, and the line and band are the effect at
the average garden ± 95 % CI. **Right:** the same against garden temperature, with the model's fitted line.

⚠ The points are **conditional modes from the same model that fitted the line**, so noisy gardens are
shrunk towards it. Their agreement with the line is partly by construction; the test is the fixed-effect
slope and its LRT, not the scatter."""

code_garden = r"""
T="sv_CDS"; sc=PV(T)
g=GD[GD.term==T]; b=sc*g.total.to_numpy(); e=1.96*sc*g.cond_sd.to_numpy(); xg=g.bio1.to_numpy()
NORM=mpl.colors.Normalize(xg.min(),xg.max()); CM=plt.get_cmap("coolwarm"); bc=xg.mean()
fig=plt.figure(figsize=(11.8,3.9)); gs=fig.add_gridspec(1,2,width_ratios=[1.3,1],wspace=0.26)
a=fig.add_subplot(gs[0,0]); style(a); panel(a,"a",x=-0.07); a.axvline(0,color="0.4",lw=0.9,ls="--")
o=np.argsort(b); b0_=sc*CO.loc[T,"estimate"]; s0_=sc*CO.loc[T,"se"]
a.axvspan(b0_-1.96*s0_,b0_+1.96*s0_,color=PC2(T),alpha=0.14,lw=0); a.axvline(b0_,color=PC2(T),lw=1.4)
for yi,j in enumerate(o):
    a.errorbar(b[j],yi,xerr=e[j],fmt="o",ms=4.2,color=CM(NORM(xg[j])),mec="0.25",mew=0.5,
               ecolor=CM(NORM(xg[j])),elinewidth=1.0,capsize=1.6)
a.set_yticks(range(len(g))); a.set_yticklabels(g.garden.to_numpy()[o].astype(int),fontsize=5.5)
a.set_ylim(-1,len(g)); a.set_ylabel("site",fontsize=7.5)
a.annotate(f"{pct(b0_)} ± {1.96*s0_:.2g} (95% CI), p = {CO.loc[T,'p']:.2g};  "
           f"negative in {(b<0).sum()}/{len(g)} gardens",xy=(0,1.02),xycoords="axes fraction",
           fontsize=7,color="0.3",va="bottom")
a.set_xlabel("% fitness per carried coding SV, per garden",fontsize=8)
a=fig.add_subplot(gs[0,1]); style(a,axis="both"); panel(a,"b",x=-0.10); a.axhline(0,color="0.4",lw=0.9,ls="--")
a.errorbar(xg,b,yerr=e,fmt="none",ecolor="0.85",elinewidth=0.8,zorder=1)
a.scatter(xg,b,s=34,c=xg,cmap=CM,norm=NORM,linewidth=0.4,edgecolor="0.3",zorder=3)
xx=np.linspace(xg.min(),xg.max(),60); xc=xx-bc; i0,i1=T,f"{T}:bio1c"
fit=CO.loc[i0,"estimate"]+CO.loc[i1,"estimate"]*xc
v=VV.loc[i0,i0]+xc**2*VV.loc[i1,i1]+2*xc*VV.loc[i0,i1]
a.fill_between(xx,sc*(fit-1.96*np.sqrt(v)),sc*(fit+1.96*np.sqrt(v)),color=PC2(T),alpha=0.15,lw=0)
a.plot(xx,sc*fit,color=PC2(T),lw=1.8)
a.annotate(f"slope {pct(sc*CO.loc[i1,'estimate'])}/°C,  LRT p = {LRT[T]:.2g}",
           xy=(0.02,0.04),xycoords="axes fraction",fontsize=7.2,color="0.25")
a.set_xlabel("garden temperature  bio1 (°C)",fontsize=8)
a.set_ylabel("% per carried variant, per garden",fontsize=7.5)
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORM,cmap=CM),ax=a,pad=0.02,fraction=0.03,aspect=26)
cb.set_label("garden bio1 (°C)",fontsize=7.5)
save(fig,"model_by_garden",nb=True); plt.show()
print(f"coding-SV effect per garden: {b.min():+.2f}% .. {b.max():+.2f}%, "
      f"negative in {(b<0).sum()}/{len(g)} gardens")
"""

md_garden2 = """### Every burden, garden by garden

The same quantity for each burden that carries a temperature interaction, so the coding result can be read
against the others. Each point is one garden's **total** effect from the full model — fixed part + its
temperature term + that garden's random slope — with every covariate held in: kinship, the SNP coding
burden, assembly route, origin climate, local adaptation fitted per garden, and the other burdens
simultaneously. **Red** marks gardens where the effect is positive. The line is the model's fitted
temperature trend with its 95 % CI."""

code_garden2 = r'''
fig,axs=plt.subplots(1,len(SVT),figsize=(4.2*len(SVT),3.5))
for a,T in zip(np.atleast_1d(axs),SVT):
    style(a,axis="both"); a.axhline(0,color="0.4",lw=0.9,ls="--")
    sc=PV(T); g=GD[GD.term==T]; b=sc*g.total.to_numpy(); xg=g.bio1.to_numpy()
    er=1.96*sc*g.cond_sd.to_numpy(); pos=b>0
    a.errorbar(xg,b,yerr=er,fmt="none",ecolor="0.85",elinewidth=0.8,zorder=1)
    a.scatter(xg[~pos],b[~pos],s=34,color=PC2(T),edgecolor="0.25",linewidth=0.4,zorder=3)
    a.scatter(xg[pos],b[pos],s=34,color="#C0392B",edgecolor="0.25",linewidth=0.4,zorder=3)
    xx=np.linspace(xg.min(),xg.max(),60); xc=xx-xg.mean(); i1=f"{T}:bio1c"
    fit=CO.loc[T,"estimate"]+CO.loc[i1,"estimate"]*xc
    v=VV.loc[T,T]+xc**2*VV.loc[i1,i1]+2*xc*VV.loc[T,i1]
    a.fill_between(xx,sc*(fit-1.96*np.sqrt(v)),sc*(fit+1.96*np.sqrt(v)),color=PC2(T),alpha=0.15,lw=0)
    a.plot(xx,sc*fit,color=PC2(T),lw=1.8,zorder=2)
    a.annotate(f"{LAB2(T)}\n{(b>0).sum()}/{len(g)} gardens positive\n"
               f"slope {pct(sc*CO.loc[i1,'estimate'])}/degC, LRT p = {LRT[T]:.2g}",
               xy=(0.03,0.97),xycoords="axes fraction",va="top",ha="left",fontsize=7.2,color="0.25")
    y0,y1=a.get_ylim(); a.set_ylim(y0,y1+0.30*(y1-y0))
    a.set_xlabel("garden temperature  bio1 (degC)",fontsize=8)
    a.set_ylabel("% per carried variant, per garden",fontsize=7.5)
fig.subplots_adjust(wspace=0.32)
save(fig,"model_by_garden_all",nb=True); plt.show()
for T in SVT:
    g=GD[GD.term==T]; b=PV(T)*g.total.to_numpy()
    print(f"{LAB2(T):20s} per garden {b.min():+.3f}% .. {b.max():+.3f}%   "
          f"positive in {(b>0).sum():2d}/{len(g)} gardens")
'''

md_pergarden = """### Garden by garden, as the model fits it

The panel you asked for, one per garden, ordered by garden temperature — but the y axis is **not** raw
ecotype fitness and the line is **not** a per-garden regression of fitness on burden.

Each point is a **partial residual** from the full model: that founder's fitness in that garden with every
other covariate already removed — kinship, the SNP coding burden, assembly route, origin climate, local
adaptation fitted per garden, and the other burden. The coloured line is that garden's **conditional**
effect from the model, fixed part + temperature term + its random slope; by construction it is the slope
through those points. The grey dashed line is the raw within-garden regression of fitness on burden, with
nothing controlled, so the two can be read against each other.

The number in each panel is the model's slope for that garden, in % per carried variant, printed **red**
where it is positive."""

code_pergarden = r'''
RS=pd.read_csv(f"{PA}/resid_M2p_reml{SUF}.csv",dtype={"founder":str,"garden":str})
GDX=GD; TT=SVT; LBP={t:LAB2(t) for t in SVT}
gb=GDX.drop_duplicates("garden")[["garden","bio1"]].copy(); gb["garden"]=gb.garden.astype(str)
order=gb.sort_values("bio1").garden.to_list(); bio=dict(zip(gb.garden,gb.bio1))
NORM=mpl.colors.Normalize(min(bio.values()),max(bio.values())); CM=plt.get_cmap("coolwarm")
NEGC,POSC="#1F4E79","#C0392B"                      # negative slope / positive slope
fig,axs=plt.subplots(len(TT),len(order),figsize=(0.80*len(order),2.0*len(TT)),
                     sharex=True,sharey="row",gridspec_kw=dict(hspace=0.34,wspace=0.05))
axs=np.atleast_2d(axs)
for ri,T in enumerate(TT):
    sc=PV(T); tot={str(r.garden):r.total for r in GDX[GDX.term==T].itertuples()}
    pr=RS.resid+RS.garden.map(tot).to_numpy()*RS[T].to_numpy()      # component + residual
    ylo,yhi=np.percentile(pr,[1,99]); xlo,xhi=RS[T].min(),RS[T].max()
    for ci,gname in enumerate(order):
        a=axs[ri,ci]; m=(RS.garden==gname).to_numpy()
        x=RS[T].to_numpy()[m]; yv=pr.to_numpy()[m]
        lc=NEGC if tot[gname]<0 else POSC                # the line is coloured by the SIGN of the slope
        a.scatter(x,yv,s=2.2,color="0.45",alpha=0.13,linewidth=0,rasterized=True)
        xx=np.array([x.min(),x.max()])
        a.plot(xx,tot[gname]*(xx-x.mean())+yv.mean(),color=lc,lw=1.8,zorder=3)
        rw=np.polyfit(x,RS.y.to_numpy()[m],1)[0]                    # raw, nothing controlled
        a.plot(xx,rw*(xx-x.mean())+yv.mean(),color="0.45",lw=0.9,ls="--",zorder=2)
        b=sc*tot[gname]
        a.annotate(f"{b:+.2f}",xy=(0.5,-0.02),xycoords="axes fraction",ha="center",va="top",
                   fontsize=5.8,color=lc)
        a.set_ylim(ylo,yhi); a.set_xlim(xlo-0.15,xhi+0.15)
        a.set_xticks([]); a.set_yticks([]); [sp.set_linewidth(0.4) for sp in a.spines.values()]
        if ri==0: a.annotate(gname,xy=(0.5,1.04),xycoords="axes fraction",ha="center",fontsize=5.8,
                             color=tuple(0.80*np.array(CM(NORM(bio[gname]))[:3])),weight="bold")
        if ci==0: a.set_ylabel(f"partial residual\n{LBP[T]}",fontsize=7)
fig.supxlabel("burden (z-scored rate) -- one panel per garden, ordered by garden temperature",
              fontsize=8.5,y=0.02)
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORM,cmap=CM),ax=axs,pad=0.006,fraction=0.012,aspect=34)
cb.set_label("garden bio1 (degC) -- colour of the garden number",fontsize=7.5)
cb.ax.tick_params(labelsize=6)
save(fig,"model_per_garden_panels",nb=True); plt.show()
print("per-garden CONDITIONAL slope (full model) vs RAW within-garden slope, % per carried variant\n")
for T in TT:
    sc=PV(T); tot={str(r.garden):r.total for r in GDX[GDX.term==T].itertuples()}
    rw=np.array([sc*np.polyfit(RS[T].to_numpy()[(RS.garden==g).to_numpy()],
                               RS.y.to_numpy()[(RS.garden==g).to_numpy()],1)[0] for g in order])
    cd=np.array([sc*tot[g] for g in order])
    print(f"{LBP[T]:18s} conditional {cd.min():+.2f}..{cd.max():+.2f}, positive in {(cd>0).sum():2d}/{len(order)}"
          f"   |   raw {rw.min():+.2f}..{rw.max():+.2f}, positive in {(rw>0).sum():2d}/{len(order)}"
          f"   |   r(raw, conditional) = {np.corrcoef(rw,cd)[0,1]:+.2f}")
'''
md_freq = _ns["md_freq"].replace("## 7 ·", "## 6 ·")
md_size = _ns["md_size"].replace("## 8 ·", "## 7 ·")
code_freq, code_size = _ns["code_freq"], _ns["code_size"]

md_avp  = _K["md_avp"].replace("## 1 ·", "## 8 ·")
md_infl = _K["md_infl"].replace("## 2 ·", "## 9 ·")
md_deld = _K["md_deld"].replace("## 3 ·", "## 10 ·")
code_load_diag, code_avp, code_infl, code_deld = (_K["code_load"], _K["code_avp"],
                                                  _K["code_infl"], _K["code_deld"])

md_lare = """## 11 · Every term of the model, in one figure

The whole fit in one place: every fixed effect with its 95 % confidence interval. Filled markers are
p < 0.05. Three panels because the scales differ by more than an order of magnitude — on a shared axis the
interaction intervals collapse onto the zero line and you cannot see whether they cross it.

**(a)** burden effects at an average-temperature garden. **(b)** how each changes per °C of garden
temperature, tested by likelihood ratio. **(c)** the two local-adaptation terms, on the model's own scale —
these have no per-carried-variant conversion.

Local adaptation is fitted **per garden** here, `(0 + obio1c | garden)`: the single `obio1c:bio1c`
coefficient forces it to scale linearly with garden temperature, and the per-garden slopes run from −0.15
to +0.19, so it does not. Everything above is fitted this way too."""

code_lare = r'''
LRS=SUF
CO2=pd.read_csv(f"{PA}/coef_M2p_reml{LRS}.csv").set_index("term")
VC2=pd.read_csv(f"{PA}/varcomp_M2p_reml{LRS}.csv")
LL2=pd.read_csv(f"{PA}/loglik_M2p_reml{LRS}.csv").iloc[0]
ll2=lambda n: pd.read_csv(f"{PA}/loglik_{n}{LRS}.csv").logLik[0]
LRT2={t: stats.chi2.sf(2*(ll2("M2p_ml")-ll2(f"M2p_ml_drop_{t}")),1) for t in ["sv_CDS","sv_TE","sv_promoter"]}
OB2=next(t for t in CO2.index if set(t.split(":"))=={"obio1c","bio1c"})
PVX={t:PV(t) for t in VPS}; CLX=PAL2
LBL={**{t:LAB2(t) for t in VPS},"obio1c":"origin climate","bio1c:obio1c":"origin x garden climate"}
BT=["snp_CDS","sv_CDS","sv_TE","sv_promoter","sv_intergenic"]; IT=["sv_CDS","sv_TE","sv_promoter"]
BLK=[[(t,t,CO2.loc[t,"p"],PVX[t],CLX[t],LBL[t]) for t in BT],
     [(f"{t}:bio1c",t,LRT2[t],PVX[t],CLX[t],LBL[t]+" x temp.") for t in IT],
     [("obio1c","obio1c",CO2.loc["obio1c","p"],1.0,"#8C6D1F",LBL["obio1c"]),
      (OB2,OB2,CO2.loc[OB2,"p"],1.0,"#8C6D1F",LBL["bio1c:obio1c"])]]
XL=["% fitness per carried variant\nat an average-temperature garden",
    "change in % fitness per carried variant\nper degC of garden temperature",
    "log-fitness per degC"]
NR=max(len(b) for b in BLK)
fig,axs=plt.subplots(1,3,figsize=(14.2,0.42*NR+1.9),
                     gridspec_kw=dict(width_ratios=[1.1,1,0.9],wspace=1.45))
for a,rows,xl,pl in zip(axs,BLK,XL,"abc"):
    style(a); panel(a,pl,x=-0.42 if pl!="a" else -0.50); a.axvline(0,color="0.35",lw=1,ls="--")
    for i,(term,base,pv,k,c,lb) in enumerate(rows):
        e,se=CO2.loc[term,"estimate"]*k,CO2.loc[term,"se"]*k
        a.errorbar(e,i,xerr=1.96*se,fmt="o",color=c,mfc=c if pv<0.05 else "white",
                   mew=1.3,ms=7.5,capsize=2.5,zorder=3)
        if i<len(rows)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
    a.set_yticks(range(len(rows))); a.set_yticklabels([r[5] for r in rows],fontsize=8.5)
    pad=(NR-len(rows))/2.0; a.set_ylim(len(rows)-0.5+pad,-0.5-pad)
    for i,(term,base,pv,k,c,lb) in enumerate(rows):
        tag=("LRT p = " if term.endswith(":bio1c") else "p = ")+f"{pv:.3g}"
        a.annotate(tag,xy=(1.03,i),xycoords=("axes fraction","data"),ha="left",va="center",
                   fontsize=7.0,color="0.35",annotation_clip=False)
    a.set_xlabel(xl,fontsize=8)
save(fig,"model_fullmodel_allterms",nb=True); plt.show()
print("THE FULL MODEL, local adaptation fitted per garden -- all fixed effects\n")
out=[]
for t in CO2.index:
    base=t.split(":")[0]; k=PVX.get(base)
    out.append({"term":t,"estimate":f"{CO2.loc[t,'estimate']:+.6f}","se":f"{CO2.loc[t,'se']:.6f}",
                "95% CI":f"[{CO2.loc[t,'estimate']-1.96*CO2.loc[t,'se']:+.4f}, "
                         f"{CO2.loc[t,'estimate']+1.96*CO2.loc[t,'se']:+.4f}]",
                "Wald p":f"{CO2.loc[t,'p']:.4g}",
                "LRT p":(f"{LRT2[base]:.4g}" if t.endswith(":bio1c") and base in LRT2 else ""),
                "%/carried variant":(f"{CO2.loc[t,'estimate']*k:+.4f}" if k else "")})
print(pd.DataFrame(out).to_string(index=False))
print("\nINTERACTIONS -- are they significant?")
print(pd.DataFrame([{"interaction":f"{t}:bio1c","%/variant/degC":f"{CO2.loc[t+':bio1c','estimate']*PVX[t]:+.5f}",
                     "LRT p":f"{LRT2[t]:.4g}","p < 0.05":LRT2[t]<0.05} for t in IT]).to_string(index=False))
print("\nRANDOM EFFECTS"); print(VC2.to_string(index=False))
print(f"\nlogLik {LL2.logLik:.4f}   df {int(LL2.df)}   singular {LL2.singular}")
'''

md_tab = """## 12 · The complete coefficient table

Every fixed effect on the model's own scale — per SD of the standardised rate, which is what R reports —
rather than the per-variant display scale used in the figures. Nothing is omitted."""

code_tab = r"""
print("="*84); print("FIXED EFFECTS")
print(pd.DataFrame([{"term":t,"estimate":f"{CO.loc[t,'estimate']:+.5f}","se":f"{CO.loc[t,'se']:.5f}",
                     "z":f"{CO.loc[t,'z']:+.3f}","Wald p":f"{CO.loc[t,'p']:.3g}"} for t in CO.index]
                  ).to_string(index=False))
print(f"\nLIKELIHOOD-RATIO TESTS for the burden x temperature terms (1 df each, {len(SVT)} fitted):")
print(pd.DataFrame([{"term":t,"LRT p":f"{LRT[t]:.5g}","p < 0.05":LRT[t]<0.05} for t in SVT]).to_string(index=False))
print("\nRANDOM EFFECTS"); print(VC.to_string(index=False))
print("\nPER-VARIANT SCALING (display only; enters no fit, changes no p-value):")
print(pd.DataFrame([{"context":c,"sd_rate":f"{SCL.loc[c,'sd_rate']:.6f}",
                     "median called":f"{SCL.loc[c,'med_called']:,.0f}",
                     "variants per SD":f"{VPS['sv_'+c]:.1f}"} for c in SCL.index]).to_string(index=False))
"""

md_read = """## 13 · Reading it

### The result

**Carrying a structural variant in coding sequence costs about 1.5 % of founder fitness per variant** —
−1.49 % at a garden of average temperature, 95 % CI [−2.34, −0.65], p = 5.6 × 10⁻⁴, with founder origin
climate and its garden-climate interaction in the model.

It is **broad, not concentrated**:

| | |
|---|---|
| frequency (§6) | rare −2.49 %, common −0.97 % per variant; contrast p = 0.08 |
| size (§7) | −1.02 → −1.25 → −2.24 per variant; omnibus across bins p = 0.56 |
| ecotypes (§9) | 231 leave-one-out refits, all negative, never within 2.9 SE of zero |
| dropping 20 % at random (§10) | 600 refits, 100 % negative, 99 % still p < 0.05 |
| gardens (§5) | negative in 27 of 31 gardens |

Neither stratification isolates a driver. Estimates order as purifying selection predicts — rarer and
larger variants worse per variant — but in neither case can the strata be told apart, and the effect
survives every partition precisely because it does not rest on a narrow subset.

### What does not hold

**The temperature interactions are weak** (§3). Coding SVs give the only one near significance,
LRT p = 0.043; reference TE and promoter are at 0.093 and 0.095. Fitted *without* the local-adaptation
terms all three are clearly significant; adding origin climate × garden climate attenuates the estimates by
29–50 % while leaving the standard errors essentially unchanged — the signature of confounding rather than
of lost precision. A ~100× shift in the coding interaction (0.0004 → 0.043) is the finding here, wherever
one puts a threshold. Those uncontrolled fits are in the registry.

### What this does not establish

⚠ **Robustness is not identification.** §§9–10 show the estimate does not depend on particular ecotypes.
Whether something *correlated with* SV burden — rather than the SVs — causes the fitness difference is a
separate question. The local-adaptation control is the relevant test and leaves `sv_CDS` almost untouched
(0.22 SE), but `obio1c` is a single summary of provenance, so the control is incomplete.

⚠ **The control tests the interactions much harder than the main effects.** `obio1c:bio1c` is an
interaction and therefore near-orthogonal to the burden main effects; the main-effect counterpart is
`obio1c` itself (p ≈ 0.04). The coding-SV effect is robust to provenance *as far as this design can tell*,
not established to the standard by which the interactions are refuted.

⚠ **Private variants are an assembly artifact** (§6). They are reported only within the cactus panel, and
even there the non-monotone result is a caution about the bin rather than a finding.

⚠ The per-garden panels (§5, §12) come from the fitted model, so they carry every covariate. A temperature
trend in them is the §3 interaction, which is marginal at best.

### Suggested reporting

> Carrying a structural variant in coding sequence costs about 1.5 % of founder fitness per variant. The
> effect holds with founder origin climate and its garden-climate interaction in the model, appears in both
> sequencing panels independently, is not concentrated in rare or in large variants, and is not driven by
> any subset of ecotypes. No climate-dependent SV burden effect survives the local-adaptation control."""

md_methods = """## As run

**The fit.** `_sv_burden_pooled_models.R` with
`POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all`,
`EXTRA_FIXED="obio1c + obio1c:bio1c"`, `RE_DROP=fid`, tasks 0–5 (`_submit_main_refid.sh`). Outputs carry
the suffix `_x_obio1cobio1cbio1c_refid`. Nothing is refitted in this notebook; every number is read from
the R output.

- **Burdens.** `founder_burden_event.csv` — truvari-collapsed SV **events**, insertions and deletions summed
  within a context, origin class `all`, z-scored across the 231 founders. The stratified sections read
  `founder_burden_event_{freq,size}.csv` from `_sv_burden_strata.py`, which lowers the carrier floor to 1
  for the frequency bins and holds it at 2 for the size bins so the singleton class cannot leak in.
- **`RE_DROP=fid`.** The independent founder intercept estimated 4.9 × 10⁻⁷ in the unpruned fit. Removing
  it leaves the log-likelihood and every coefficient and standard error unchanged to five decimals, costs
  one parameter, and leaves no variance component at the boundary.
- **`obio1c` vs `bio1c`.** `obio1c` is the **founder's origin** bio1, centred over the founders and *not*
  standardised (hence the rescaling in §3); `bio1c` is the **garden's** bio1, centred over the 31 gardens.
  lme4 names the interaction `bio1c:obio1c` by model-frame order; the load cell resolves the spelling.
- **Per-variant scaling.** `sd_rate × med_called` from the arm's `pooled_rate_scale*.csv`. A display scaling
  only. `RE_DROP` does not change the burdens, so the reported fit shares the scale file of its unpruned
  twin. Each *stratified* arm has its own file: before 2026-09-29 the name was tagged only with `SUBSET`
  and stratified arms overwrote one another; `_regen_rate_scales.py` rebuilt the lost ones, validated by
  reproducing the survivors to machine precision.
- **LRTs.** Read from the named `loglik_*.csv` files, never by globbing: this directory also holds
  `SUBSET=`, `RE_DROP=` and other `EXTRA_FIXED` runs.
- **Diagnostics (§§7–10).** Computed here, not in R. The response is averaged over the 31 gardens, reducing
  the mixed model to a 231-observation GLS whose main effects reproduce the LMM's to machine precision (§7
  verifies it) because `bio1c` is centred and the design is balanced. Leave-one-out and delete-d refits
  *condition* on the retained founders — the covariance is re-inverted on the retained submatrix of `V`,
  not row-deleted from `V⁻¹`, which would be a different model.
- **Every other specification** ever fitted — `shared` origin filter, TE dropped, no local-adaptation
  control, split insertions/deletions, single-technology panels — is indexed in
  `results/sv_adaptive/model_specs.csv` and `model_registry.csv`, each fit carrying a `spec_*.csv` with its
  literal formula, founder counts, singularity flag and git commit."""

if __name__ == "__main__":
    build("sv_burden_model.ipynb", [
        new_markdown_cell(md_title), new_code_cell(code_setup), new_code_cell(code_load),
        new_markdown_cell(md_counts), new_code_cell(code_counts),
        new_markdown_cell(md_model),  new_code_cell(code_model),
        new_markdown_cell(md_inter),  new_code_cell(code_inter),
        new_markdown_cell(md_scale),  new_code_cell(code_scale),
        new_markdown_cell(md_garden), new_code_cell(code_garden),
        new_markdown_cell(md_garden2), new_code_cell(code_garden2),
        new_markdown_cell(md_pergarden), new_code_cell(code_pergarden),
        new_markdown_cell(md_freq),   new_code_cell(code_freq),
        new_markdown_cell(md_size),   new_code_cell(code_size),
        new_markdown_cell(md_avp),    new_code_cell(code_load_diag), new_code_cell(code_avp),
        new_markdown_cell(md_infl),   new_code_cell(code_infl),
        new_markdown_cell(md_deld),   new_code_cell(code_deld),
        new_markdown_cell(md_lare),   new_code_cell(code_lare),
        new_markdown_cell(md_tab),    new_code_cell(code_tab),
        new_markdown_cell(md_read),
        new_markdown_cell(md_methods),
    ])
