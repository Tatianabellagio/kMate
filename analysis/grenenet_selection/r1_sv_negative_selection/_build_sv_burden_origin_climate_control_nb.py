#!/usr/bin/env python
"""Build+execute `sv_burden_origin_climate_control.ipynb`: THE FULL POOLED, NO-ORIGIN-FILTER burden model
-- SV burdens + burden x garden-temperature interactions + a founder-origin-climate (local-adaptation)
term -- reported in full, and compared against the same model with the local-adaptation terms removed.

Only two fits are named anywhere in the notebook:
  REDUCED  = burdens + burden:bio1c                               (coef_M2p_reml.csv)
  FULL     = burdens + burden:bio1c + obio1c + obio1c:bio1c       (coef_M2p_reml_x_obio1c...csv)
The main-effects-only fit (M1p) is deliberately NOT read: burden main effects are reported from the FULL
model, where they are the burden effect at an average-temperature garden (bio1c is centred).

Primary arm: POOL_SET=sv6 BURDEN_ORIGIN=all -> final_models_pooled_all (reference TE kept).
Sensitivity arm: POOL_SET=sv4 -> final_models_pooled_all_noTE. Control submitted by
`_submit_pooled_sv6_origctrl.sh`. Setup cell reused from `_main_fig_common.code_setup_burden`. `basic` env.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_fig_common import build, new_markdown_cell, new_code_cell, code_setup_burden as code_setup

md_title = """# The full pooled SV burden model — interactions and local adaptation together

## The two fits, and nothing else

Only **two** fits appear in this notebook. They are named **REDUCED** and **FULL** throughout; no other
model is referred to anywhere.

| name | what it contains |
|---|---|
| **REDUCED** | SV burdens + burden × garden-temperature interactions |
| **FULL** | the same, **plus** founder origin climate and its interaction with garden climate |

FULL is the model being reported. REDUCED exists only to show what the local-adaptation terms change.

```r
# REDUCED
y ~ snp_CDS + sv_CDS + sv_TE + sv_promoter + sv_intergenic + cactus
    + sv_CDS:bio1c + sv_TE:bio1c + sv_promoter:bio1c
    + (1 | garden) + (1 | kin) + (1 | fid)
    + (0 + sv_CDS | garden) + (0 + sv_TE | garden) + (0 + sv_promoter | garden)

# FULL   ( = REDUCED + the two boxed terms )
y ~ snp_CDS + sv_CDS + sv_TE + sv_promoter + sv_intergenic + cactus
    + sv_CDS:bio1c + sv_TE:bio1c + sv_promoter:bio1c
    + obio1c + obio1c:bio1c                          # <-- local adaptation
    + (1 | garden) + (1 | kin) + (1 | fid)
    + (0 + sv_CDS | garden) + (0 + sv_TE | garden) + (0 + sv_promoter | garden)
```

Both transcribed verbatim from the run logs. Fitted with `lme4qtl::relmatLmer`, kinship passed as
`relmat`, on **all 231 founders** × 31 gardens = **7,161 rows**.

**Both fits use the whole panel — 151 short-read plus 80 cactus-assembly founders.** Assembly route enters
as the fixed effect `cactus`, it is not a subset: `_sv_burden_pooled_models.R` *removes* that term whenever
it is fitting one technology panel, so `cactus` appearing in the formulas above is itself the proof that
neither fit is restricted. Single-panel refits do exist in this section (files suffixed `_cactus` and
`_shortread`), and **none of them is read by this notebook.**

## The data

Insertions and deletions are **pooled** into one SV burden per genomic context, and **no origin filter** is
applied — every variant passing the technical filters enters, not only those carried by both cold- and
warm-origin founders. Unit is truvari-collapsed SV **events**, one context each by first match
(CDS > UTR > intron > reference TE > promoter ≤ 2 kb > intergenic):

```
sv_<ctx>  =  z( (ins carried + del carried) / (ins called + del called) )
```

| symbol | meaning |
|---|---|
| `bio1c` | the **garden's** mean annual temperature, centred on the mean of the 31 gardens |
| `obio1c` | the **founder's origin** bio1 — temperature where the accession comes from — centred over the 231 founders |
| `obio1c:bio1c` | **local adaptation**: origin climate × garden climate. Contains no SV information at all |
| `snp_CDS` | matching SNP coding burden · `cactus` = assembly route (0/1) |
| `kin` / `fid` | founder random effects, with covariance σ²K (SNP kinship) and independent |
| `(0 + sv_ctx \\| garden)` | independent random garden slope per modelled burden |

## Reading the FULL model's coefficients

Because `bio1c` is **centred**, each burden term splits into two numbers that must be read together:

- **`sv_ctx`** — the burden's effect **at a garden of average temperature** (13.5 °C). This is the
  "are SVs in this context deleterious?" number.
- **`sv_ctx:bio1c`** — how that effect **changes per °C** of garden temperature.

Both come out of the one FULL fit. A separate main-effects-only model is not needed and is not used here.

## Why the local-adaptation terms belong in the model

A `burden × garden-temperature` term asks whether carrying more SVs changes fitness as gardens warm. But
founders differ in **where they came from**, and origin × garden climate is plain local adaptation — an
effect that exists whether or not SVs do anything. `sv_founder_confound.ipynb` showed founder SV content
*is* related to origin climate (inserted kb vs origin bio1, ρ = −0.542), so a burden × temperature term is
exposed to picking up provenance instead.

Dropping the `shared` origin filter does not address this: it only stopped a single *variant* from being a
lineage marker, never the *aggregate burden* from tracking provenance, and it put nothing about origin into
the model. `obio1c + obio1c:bio1c` does.

⚠ `obio1c:bio1c` is itself an interaction, so it is near-orthogonal to the burden main effects: it is a
**strong test of the `sv_ctx:bio1c` terms and a weak test of the `sv_ctx` terms**. §8 returns to this."""

code_load = r"""
PA=f"{G}/final_models_pooled_all"                           # primary: pooled, no origin filter, TE kept (sv6)
PN=f"{G}/final_models_pooled_all_noTE"                      # sensitivity: same, TE term dropped (sv4)
SUF="_x_obio1cobio1cbio1c"                                  # the FULL fit's filename suffix
SVT6=["sv_CDS","sv_TE","sv_promoter"]                       # burdens carrying a temperature interaction
SVT4=["sv_CDS","sv_promoter"]
BURD=["snp_CDS","sv_CDS","sv_TE","sv_promoter","sv_intergenic"]    # burdens with a per-variant scale
PAL2={"snp_CDS":"0.25","sv_CDS":"#4A3B6B","sv_TE":"#6B5A8E","sv_promoter":"#9084AE",
      "sv_intergenic":"#C6BCD8","cactus":"0.6","obio1c":"#8C6D1F"}
PC2=lambda t: PAL2[("obio1c" if "obio1c" in t else t.split(":")[0])]
LAB2=lambda t: {"snp_CDS":"SNPs, coding","sv_CDS":"SVs, coding","sv_TE":"SVs, reference TE",
                "sv_promoter":"SVs, promoter","sv_intergenic":"SVs, intergenic","cactus":"assembly route",
                "obio1c":"founder origin climate"}[("obio1c" if "obio1c" in t else t.split(":")[0])]

# Only the M2p fits are read. M1p (the main-effects-only model) is deliberately not loaded: burden main
# effects are reported from the FULL fit, where bio1c is centred so they are the effect at an
# average-temperature garden.
def load(d,terms,suf):
    o={"co":pd.read_csv(f"{d}/coef_M2p_reml{suf}.csv").set_index("term"),
       "VC":pd.read_csv(f"{d}/varcomp_M2p_reml{suf}.csv"),
       "G":pd.read_csv(f"{d}/garden_M2p_reml{suf}.csv"),
       "V":pd.read_csv(f"{d}/vcov_M2p_reml{suf}.csv",index_col=0)}
    # named loglik files, never a glob: these directories also hold SUBSET= and other EXTRA_FIXED runs,
    # and a glob would silently mix likelihoods from models fitted on different data
    ll=lambda n: pd.read_csv(f"{d}/loglik_{n}{suf}.csv").logLik[0]
    o["LRT"]={t: stats.chi2.sf(2*(ll("M2p_ml")-ll(f"M2p_ml_drop_{t}")),1) for t in terms}
    o["BONF"]={t: min(o["LRT"][t]*len(terms),1.0) for t in terms}
    return o
RED =load(PA,SVT6,"")   ; FULL =load(PA,SVT6,SUF)           # primary arm
REDN=load(PN,SVT4,"")   ; FULLN=load(PN,SVT4,SUF)           # sensitivity arm (TE dropped)
OBX=next(t for t in FULL["co"].index if set(t.split(":"))=={"obio1c","bio1c"})   # lme4 spells it bio1c:obio1c

# per-carried-variant conversion: coefficients are per SD of the standardised rate, so
# estimate / (SD of raw rate x median events called) x 100 = % fitness per carried variant.
# Display scaling only -- enters no fit, changes no p-value. The two arms' scale files are byte-identical.
SCL=pd.read_csv(f"{PA}/pooled_rate_scale.csv").set_index("context")
VPS={f"sv_{c}": SCL.loc[c,"sd_rate"]*SCL.loc[c,"med_called"] for c in SCL.index}
VPS["snp_CDS"]=pd.read_csv(f"{G}/sv_burden_full_model_event_CDS+promoter.csv").set_index("term").variants_per_sd["snp_CDS"]
PV=lambda t: 100/VPS[t.split(":")[0]]

FB=pd.read_csv(f"{G}/founder_burden_event.csv"); sd_ob=FB.bio1.std(ddof=1)
gb=pd.read_csv(f"{G}/onestage_all/long.csv").drop_duplicates("garden").bio1
# BOTH fits use the whole panel. `cactus` being present as a fixed effect proves it: _sv_burden_pooled_models.R
# drops that term whenever SUBSET is set, because it is constant within a single technology panel.
nc=int((FB.is_cactus==True).sum()|0) if FB.is_cactus.dtype==bool else int((FB.is_cactus.astype(str)=="True").sum())
assert "cactus" in RED["co"].index and "cactus" in FULL["co"].index, "a SUBSET= fit got loaded by mistake"
print(f"FOUNDERS IN BOTH FITS: all {len(FB)}  =  {len(FB)-nc} short-read + {nc} cactus assembly")
print(f"  assembly route is a COVARIATE, not a subset: cactus {FULL['co'].loc['cactus','estimate']:+.4f} "
      f"(p = {FULL['co'].loc['cactus','p']:.3g}) in the FULL fit")
print(f"  (single-technology refits exist elsewhere as _cactus / _shortread files; none is read here)\n")
print(f"founders {len(FB)}   origin bio1 {FB.bio1.min():.1f} to {FB.bio1.max():.1f} degC (SD {sd_ob:.2f})")
print(f"gardens  {len(gb)}   garden bio1 {gb.min():.1f} to {gb.max():.1f} degC "
      f"(mean {gb.mean():.2f} = where bio1c is zero)")
print(f"rows     {len(FB)*len(gb):,}\n")
print("variants per SD of each burden rate (the per-carried-variant divisor):")
for c in SCL.index:
    print(f"  sv_{c:<11} {VPS['sv_'+c]:8.1f} variants")
print(f"\nthe local-adaptation terms in the FULL fit (lme4 spells the interaction {OBX!r}):")
for t in ("obio1c",OBX):
    r=FULL["co"].loc[t]; print(f"  {t:14s} {r.estimate:+.6f} +- {r.se:.6f}   p = {r.p:.3g}")
"""

md_counts = """## 1 · What enters the model

Every truvari-collapsed SV event passing call rate ≥ 90 % and ≥ 2 carriers — **no origin filter**. Bars are
split by class so the quantity the pooling merges is visible: the pooled burden for a context is the green
plus the orange segment. The `shared`-filtered subset earlier models used is drawn open, for scale only —
it is not used anywhere in this notebook.

`UTR` and `intron` carry no term in this model; they are shown so the context accounting is complete."""

code_counts = r"""
CAT=pd.read_csv(f"{G}/sv_event_catalog.csv"); AP=CAT[CAT.passes]; SP=AP[AP.origin=="shared"]
CTX6=["CDS","UTR","intron","TE","promoter","intergenic"]; MODC={"CDS","TE","promoter","intergenic"}
CI,CD="#2F4A1F","#B04A00"
fig,a=plt.subplots(figsize=(7.8,3.6)); style(a); hh=0.36
for i,c in enumerate(CTX6):
    for off,src,fill in ((-hh/2,AP,True),(hh/2,SP,False)):
        ni=int(((src.context==c)&(src.kind=="ins")).sum()); nd=int(((src.context==c)&(src.kind=="del")).sum())
        a.barh(i+off,ni,height=hh,color=CI if fill else "white",edgecolor=CI,lw=1.0)
        a.barh(i+off,nd,height=hh,left=ni,color=CD if fill else "white",edgecolor=CD,lw=1.0)
        if fill: a.annotate(f"{ni+nd:,}",xy=(ni+nd,i+off),xytext=(4,0),textcoords="offset points",
                            va="center",fontsize=7)
    if i<len(CTX6)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(CTX6)))
a.set_yticklabels([c+("" if c in MODC else "  (not a model term)") for c in CTX6],fontsize=8)
a.invert_yaxis(); a.set_xlabel("SV events",fontsize=8.5); a.set_xlim(0,8400)
h_=[plt.Rectangle((0,0),1,1,color=CI),plt.Rectangle((0,0),1,1,color=CD),
    plt.Rectangle((0,0),1,1,facecolor="white",edgecolor="0.45",lw=1.0)]
a.legend(h_,["insertions","deletions","shared-filtered subset (not used here)"],
         fontsize=7.5,loc="lower right",handlelength=1.2)
save(fig,"fullmodel_inputs",nb=True); plt.show()
T=pd.DataFrame({"insertions":[int(((AP.context==c)&(AP.kind=="ins")).sum()) for c in CTX6],
                "deletions":[int(((AP.context==c)&(AP.kind=="del")).sum()) for c in CTX6]},index=CTX6)
T["pooled total"]=T.insertions+T.deletions; T["in the model"]=[c in MODC for c in CTX6]
print("SV events entering the pooled burdens, by context (no origin filter):")
print(T.to_string())
print(f"\nTOTAL {int(T['pooled total'].sum()):,} events  |  the `shared` subset would have been "
      f"{len(SP):,} ({len(SP)/len(AP):.0%}) -- this model uses all of them")
"""

md_full = """## 2 · The FULL model

Every fixed effect of the FULL fit, in one figure. The three blocks have **different units**, so they get
their own panels — read each against its own axis.

**(a) Burden effect at an average-temperature garden** (`sv_ctx`), in % fitness per carried variant.
**(b) How that effect changes with temperature** (`sv_ctx:bio1c`), in % per carried variant per °C.
Filled markers clear Bonferroni over the three interactions tested.
**(c) The local-adaptation terms**, which are not burdens and have no per-variant scale — `obio1c` is per
°C of origin climate, `obio1c:bio1c` per °C of origin × °C of garden.

This is the whole model. Nothing else is fitted."""

code_full = r"""
fig=plt.figure(figsize=(13.2,3.5)); gs=fig.add_gridspec(1,3,width_ratios=[1.25,1,0.95],wspace=0.50)
a=fig.add_subplot(gs[0,0]); style(a); panel(a,"a",x=-0.34); a.axvline(0,color="0.35",lw=1,ls="--")
for i,t in enumerate(BURD):
    e,se=FULL["co"].loc[t,"estimate"]*PV(t), FULL["co"].loc[t,"se"]*PV(t)
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=PC2(t),mfc=PC2(t) if FULL["co"].loc[t,"p"]<0.05 else "white",
               mew=1.2,ms=6.5,capsize=2,zorder=3)
    if i<len(BURD)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(BURD))); a.set_yticklabels([LAB2(t) for t in BURD]); a.invert_yaxis()
x0,x1=a.get_xlim(); a.set_xlim(x0,x1+0.42*(x1-x0))          # room for the p-values inside the panel
for i,t in enumerate(BURD):
    a.annotate(f"p = {FULL['co'].loc[t,'p']:.2g}",xy=(0.985,i),xycoords=("axes fraction","data"),
               ha="right",va="center",fontsize=6.6,color="0.35")
a.set_xlabel("% fitness per carried variant,\nat an average-temperature garden  (± 95% CI)",fontsize=7.5)
a=fig.add_subplot(gs[0,1]); style(a); panel(a,"b",x=-0.42); a.axvline(0,color="0.35",lw=1,ls="--")
for i,t in enumerate(SVT6):
    n=f"{t}:bio1c"; e,se=FULL["co"].loc[n,"estimate"]*PV(t), FULL["co"].loc[n,"se"]*PV(t)
    sig=FULL["BONF"][t]<0.05
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=PC2(t),mfc=PC2(t) if sig else "white",
               mew=1.2,ms=6.5,capsize=2,zorder=3)
    if i<len(SVT6)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(SVT6))); a.set_yticklabels([LAB2(t) for t in SVT6]); a.invert_yaxis()
a.set_ylim(len(SVT6)-0.5,-0.5)
x0,x1=a.get_xlim(); a.set_xlim(x0,x1+0.50*(x1-x0))
for i,t in enumerate(SVT6):
    a.annotate(f"Bonf {FULL['BONF'][t]:.2g}",xy=(0.985,i),xycoords=("axes fraction","data"),
               ha="right",va="center",fontsize=6.6,color="0.35")
a.set_xlabel("change in % fitness per carried variant,\nper °C of garden temperature  (± 95% CI)",fontsize=7.5)
a=fig.add_subplot(gs[0,2]); style(a); panel(a,"c",x=-0.44); a.axvline(0,color="0.35",lw=1,ls="--")
LA=[("obio1c","origin climate"),(OBX,"origin × garden climate")]
for i,(n,lb) in enumerate(LA):
    e,se=FULL["co"].loc[n,"estimate"], FULL["co"].loc[n,"se"]
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=PAL2["obio1c"],mew=1.2,ms=6.5,capsize=2,zorder=3)
    a.annotate(f"p = {FULL['co'].loc[n,'p']:.2g}",xy=(0.5,i+0.30),xycoords=("axes fraction","data"),
               ha="center",fontsize=6.8,color="0.3")
    if i<len(LA)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(LA))); a.set_yticklabels([lb for _,lb in LA],fontsize=8); a.invert_yaxis()
a.set_ylim(len(LA)-0.5,-0.5); a.set_xlabel("log-fitness per °C  (± 95% CI)",fontsize=7.5)
save(fig,"fullmodel_coefficients",nb=True); plt.show()
print("THE FULL MODEL — burdens + burden×temperature + local adaptation\n")
print("(a) burden effect at an average-temperature garden, % fitness per carried variant:")
print(pd.DataFrame([{"term":LAB2(t),"%/variant":f"{FULL['co'].loc[t,'estimate']*PV(t):+.4f}",
    "95% CI":f"[{(FULL['co'].loc[t,'estimate']-1.96*FULL['co'].loc[t,'se'])*PV(t):+.3f}, "
             f"{(FULL['co'].loc[t,'estimate']+1.96*FULL['co'].loc[t,'se'])*PV(t):+.3f}]",
    "Wald p":f"{FULL['co'].loc[t,'p']:.4g}"} for t in BURD]).to_string(index=False))
print(f"\n(b) burden × garden temperature, % per carried variant per °C (LRT, Bonferroni ×{len(SVT6)}):")
print(pd.DataFrame([{"term":LAB2(t),"%/variant/°C":f"{FULL['co'].loc[t+':bio1c','estimate']*PV(t):+.4f}",
    "Wald p":f"{FULL['co'].loc[t+':bio1c','p']:.3g}","LRT p":f"{FULL['LRT'][t]:.3g}",
    "Bonferroni":f"{FULL['BONF'][t]:.3g}","survives":FULL["BONF"][t]<0.05} for t in SVT6]).to_string(index=False))
print("\n(c) local adaptation (not burdens, no per-variant scale):")
print(pd.DataFrame([{"term":n,"estimate":f"{FULL['co'].loc[n,'estimate']:+.6f}",
    "se":f"{FULL['co'].loc[n,'se']:.6f}","Wald p":f"{FULL['co'].loc[n,'p']:.3g}"} for n,_ in LA]).to_string(index=False))
print(f"\nnon-burden nuisance terms: intercept {FULL['co'].loc['(Intercept)','estimate']:+.4f}, "
      f"cactus {FULL['co'].loc['cactus','estimate']:+.4f} (p = {FULL['co'].loc['cactus','p']:.3g})")
"""

md_cmp = """## 3 · What the local-adaptation terms change — FULL vs REDUCED

The same coefficients in the two fits. **Open** markers are REDUCED (no local-adaptation terms), **filled**
markers FULL. Panel **(c)** separates the two ways a coefficient can move when a term is added:

- **the estimate shifting** (plotted in units of its own REDUCED standard error) — that is *confounding*,
  the new term taking back signal the burden term had been carrying;
- **the standard error changing** — that would be *variance competition*, two predictors sharing
  information and each widening the other.

Which of the two dominates is the whole question, and the answer differs sharply between the burden main
effects and the burden × temperature interactions."""

code_cmp = r"""
fig=plt.figure(figsize=(13.4,3.6)); gs=fig.add_gridspec(1,3,width_ratios=[1.1,0.9,1.5],wspace=0.55); OFF=0.17
a=fig.add_subplot(gs[0,0]); style(a); panel(a,"a",x=-0.38); a.axvline(0,color="0.35",lw=1,ls="--")
for i,t in enumerate(BURD):
    for o,(src,hollow) in zip((OFF,-OFF),((RED,True),(FULL,False))):
        e,se=src["co"].loc[t,"estimate"]*PV(t), src["co"].loc[t,"se"]*PV(t)
        a.errorbar(e,i+o,xerr=1.96*se,fmt="o",color=PC2(t),mfc="white" if hollow else PC2(t),
                   mew=1.2,ms=6,capsize=2,zorder=3)
    if i<len(BURD)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(BURD))); a.set_yticklabels([LAB2(t) for t in BURD]); a.invert_yaxis()
a.set_xlabel("% fitness per carried variant,\nat an average-temperature garden",fontsize=7.5)
a=fig.add_subplot(gs[0,1]); style(a); panel(a,"b",x=-0.46); a.axvline(0,color="0.35",lw=1,ls="--")
for i,t in enumerate(SVT6):
    n=f"{t}:bio1c"
    for o,(src,hollow) in zip((OFF,-OFF),((RED,True),(FULL,False))):
        e,se=src["co"].loc[n,"estimate"]*PV(t), src["co"].loc[n,"se"]*PV(t)
        a.errorbar(e,i+o,xerr=1.96*se,fmt="o",color=PC2(t),mfc="white" if hollow else PC2(t),
                   mew=1.2,ms=6,capsize=2,zorder=3)
    if i<len(SVT6)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(SVT6))); a.set_yticklabels([LAB2(t) for t in SVT6]); a.invert_yaxis()
a.set_ylim(len(SVT6)-0.5,-0.5)
a.set_xlabel("change per °C of\ngarden temperature",fontsize=7.5)
TERMS=[(t,"burden") for t in BURD]+[(f"{t}:bio1c","interaction") for t in SVT6]
sh=[(FULL["co"].loc[t,"estimate"]-RED["co"].loc[t,"estimate"])/RED["co"].loc[t,"se"] for t,_ in TERMS]
ds=[100*(FULL["co"].loc[t,"se"]/RED["co"].loc[t,"se"]-1) for t,_ in TERMS]
a=fig.add_subplot(gs[0,2]); style(a,axis="y"); panel(a,"c",x=-0.07)
xs=np.arange(len(TERMS)); w=0.38
a.bar(xs-w/2,sh,width=w,color=[PC2(t) for t,_ in TERMS],lw=0,zorder=3)
a.bar(xs+w/2,[d/10 for d in ds],width=w,color="white",edgecolor=[PC2(t) for t,_ in TERMS],lw=1.1,zorder=3)
for i,(q,d) in enumerate(zip(sh,ds)):
    a.annotate(f"{q:+.2f}",xy=(i-w/2,q),xytext=(0,3 if q>=0 else -9),textcoords="offset points",
               ha="center",fontsize=6.0,color="0.25")
    a.annotate(f"{d:+.0f}%",xy=(i+w/2,d/10),xytext=(0,3 if d>=0 else -9),textcoords="offset points",
               ha="center",fontsize=6.0,color="0.5")
a.axhline(0,color="0.35",lw=1); a.axvline(len(BURD)-0.5,color="0.55",lw=1.0,ls="--")
m_=max(max(np.abs(sh)),max(np.abs(ds))/10)*1.30
a.set_ylim(-m_,m_)
a.set_xticks(xs)
a.set_xticklabels([LAB2(t).replace(", ","\n")+("\n× temp." if k=="interaction" else "") for t,k in TERMS],fontsize=5.6)
a.set_ylabel("shift in estimate (base SEs)",fontsize=8)
a2=a.twinx(); a2.set_ylim(-m_*10,m_*10); a2.set_ylabel("change in SE (%)",fontsize=8,color="0.45")
a2.tick_params(labelsize=6.5,colors="0.45")
for sp in a2.spines.values(): sp.set_visible(False)
a.annotate("burden effects",xy=((len(BURD)-1)/2/len(TERMS),1.03),xycoords="axes fraction",
           ha="center",va="bottom",fontsize=7.5,color="0.35",annotation_clip=False)
a.annotate("burden × temperature",xy=((len(BURD)+len(TERMS)-1)/2/len(TERMS),1.03),xycoords="axes fraction",
           ha="center",va="bottom",fontsize=7.5,color="0.35",annotation_clip=False)
h_=[a.plot([],[],"o",color="0.35",mfc="white",mew=1.2,ms=6)[0],a.plot([],[],"o",color="0.35",ms=6)[0],
    plt.Rectangle((0,0),1,1,color="0.45"),plt.Rectangle((0,0),1,1,facecolor="white",edgecolor="0.45",lw=1.1)]
fig.legend(h_,["REDUCED","FULL","shift in estimate","change in SE"],loc="lower center",
           bbox_to_anchor=(0.5,-0.17),ncol=4,fontsize=7.5,columnspacing=1.5,handlelength=1.3)
save(fig,"fullmodel_vs_reduced",nb=True); plt.show()
print("FULL vs REDUCED — every burden coefficient:\n")
print(pd.DataFrame([{"term":LAB2(t)+(" × temp." if k=="interaction" else ""),"block":k,
    "REDUCED":f"{RED['co'].loc[t,'estimate']*PV(t):+.4f}","FULL":f"{FULL['co'].loc[t,'estimate']*PV(t):+.4f}",
    "shift (base SEs)":f"{s:+.2f}","Δ SE":f"{d:+.1f}%"}
    for (t,k),s,d in zip(TERMS,sh,ds)]).to_string(index=False))
nb_=len(BURD)
print(f"\nmean |shift|:  burden effects {np.mean(np.abs(sh[:nb_])):.2f} SE      "
      f"burden × temperature {np.mean(np.abs(sh[nb_:])):.2f} SE")
print(f"mean |Δ SE| :  burden effects {np.mean(np.abs(ds[:nb_])):.1f}%          "
      f"burden × temperature {np.mean(np.abs(ds[nb_:])):.1f}%")
print("\nestimates move, standard errors do not -> confounding, not variance competition.")
print(f"\ninteraction evidence, REDUCED -> FULL (LRT, Bonferroni ×{len(SVT6)}):")
print(pd.DataFrame([{"term":LAB2(t),"REDUCED Bonf":f"{RED['BONF'][t]:.4g}","FULL Bonf":f"{FULL['BONF'][t]:.4g}",
    "survives REDUCED":RED["BONF"][t]<0.05,"survives FULL":FULL["BONF"][t]<0.05} for t in SVT6]).to_string(index=False))
"""

md_scale = """## 4 · How big is local adaptation, next to the burden interactions?

Burden coefficients are per SD of the standardised rate; `obio1c` is entered **centred but not
standardised**, so its coefficient is per °C of origin climate. Multiplying it by the SD of founder origin
bio1 puts both on one scale: change in log-fitness **per SD of the predictor, per °C of garden
temperature**. All values are from the FULL fit."""

code_scale = r"""
TT=[(OBX,sd_ob)]+[(f"{t}:bio1c",1.0) for t in SVT6]
fig,a=plt.subplots(figsize=(6.8,2.7)); style(a); a.axvline(0,color="0.35",lw=1,ls="--")
for i,(n,k) in enumerate(TT):
    e,se=FULL["co"].loc[n,"estimate"]*k, FULL["co"].loc[n,"se"]*k
    a.errorbar(e,i,xerr=1.96*se,fmt="o",color=PC2(n),ms=6.5,mew=1.2,capsize=2,zorder=3)
    a.annotate(f"{e:+.4f}   p = {FULL['co'].loc[n,'p']:.2g}",xy=(1.03,i),xycoords=("axes fraction","data"),
               va="center",fontsize=7.2,color="0.25",annotation_clip=False)
    if i<len(TT)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(TT))); a.set_yticklabels([LAB2(n)+" × garden temp." for n,_ in TT],fontsize=8)
a.invert_yaxis()
a.set_xlabel("change in log-fitness per SD of the predictor,\nper °C of garden temperature  (± 95% CI)",fontsize=7.5)
save(fig,"fullmodel_local_adaptation_scale",nb=True); plt.show()
ob=FULL["co"].loc[OBX,"estimate"]*sd_ob
print("all × temperature terms on one scale (per SD of predictor, per degC of garden), FULL fit:")
for n,k in TT:
    e=FULL["co"].loc[n,"estimate"]*k
    print(f"  {LAB2(n)+' x garden temp.':<42} {e:+.5f}   ({abs(e/ob):.2f}x the local-adaptation term)")
print(f"\nvariance components (FULL fit):"); print(FULL["VC"].to_string(index=False))
"""

md_garden = """## 5 · Garden by garden

Each garden's own total effect of the coding-SV burden — the fixed part, plus the temperature term, plus
that garden's random slope — in % per carried variant. **Left** REDUCED, **right** FULL, on shared axes.
Points are coloured by garden temperature; the line and band are the model's fit ± 95 % CI.

The coding burden is negative in **every** garden in both panels — that is the `sv_CDS` effect showing
through, and it is what the local-adaptation terms leave alone. What they flatten is the *tilt* of the
cloud against temperature, which is `sv_CDS:bio1c`.

⚠ The points are **conditional modes from the same model that fitted the line**, so noisy gardens are shrunk
towards it. Their agreement with the line is partly by construction; the test is the fixed-effect slope and
its LRT, not the scatter."""

code_garden = r"""
T="sv_CDS"; sc=PV(T)
b1a=RED["G"].drop_duplicates("garden").bio1
NORM=mpl.colors.Normalize(b1a.min(),b1a.max()); CM=plt.get_cmap("coolwarm"); bc=b1a.mean()
fig,axs=plt.subplots(1,2,figsize=(11.6,3.8),sharey=True,gridspec_kw=dict(wspace=0.08))
for j,(lab,src) in enumerate((("REDUCED",RED),("FULL",FULL))):
    a=axs[j]; style(a,axis="both"); panel(a,"ab"[j],x=-0.06)
    g=src["G"][src["G"].term==T]; b=sc*g.total.to_numpy(); e=1.96*sc*g.cond_sd.to_numpy(); xg=g.bio1.to_numpy()
    a.axhline(0,color="0.4",lw=0.9,ls="--")
    a.errorbar(xg,b,yerr=e,fmt="none",ecolor="0.85",elinewidth=0.8,zorder=1)
    a.scatter(xg,b,s=34,c=xg,cmap=CM,norm=NORM,linewidth=0.4,edgecolor="0.3",zorder=3)
    xx=np.linspace(xg.min(),xg.max(),60); xc=xx-bc; V=src["V"]; i0,i1=T,f"{T}:bio1c"
    fit=src["co"].loc[i0,"estimate"]+src["co"].loc[i1,"estimate"]*xc
    v=V.loc[i0,i0]+xc**2*V.loc[i1,i1]+2*xc*V.loc[i0,i1]
    a.fill_between(xx,sc*(fit-1.96*np.sqrt(v)),sc*(fit+1.96*np.sqrt(v)),color=PC2(T),alpha=0.15,lw=0)
    a.plot(xx,sc*fit,color=PC2(T),lw=1.8)
    a.annotate(lab,xy=(0,1.03),xycoords="axes fraction",fontsize=9.5,color=PC2(T),va="bottom",fontweight="bold")
    a.annotate(f"slope {pct(sc*src['co'].loc[i1,'estimate'])}/°C,  Bonferroni {src['BONF'][T]:.2g}\n"
               f"negative in {(b<0).sum()}/{len(g)} gardens",
               xy=(0.02,0.04),xycoords="axes fraction",fontsize=7.2,color="0.25")
    a.set_xlabel("garden temperature  bio1 (°C)",fontsize=8)
    if j==0: a.set_ylabel("% fitness per carried coding SV, per garden",fontsize=8)
cb=fig.colorbar(mpl.cm.ScalarMappable(norm=NORM,cmap=CM),ax=list(axs),pad=0.015,fraction=0.02,aspect=32)
cb.set_label("garden temperature  bio1 (°C)",fontsize=8)
save(fig,"fullmodel_by_garden",nb=True); plt.show()
for lab,src in (("REDUCED",RED),("FULL",FULL)):
    g=src["G"][src["G"].term==T]; b=sc*g.total.to_numpy()
    print(f"  {lab:<9} coding-SV effect per garden: {b.min():+.2f}% .. {b.max():+.2f}%, "
          f"negative in {(b<0).sum()}/{len(g)} gardens")
"""

md_te = """## 6 · Sensitivity — dropping the reference-TE term

The same two fits with `sv_TE` removed entirely, so TE-context variants enter no burden at all (Bonferroni
over two interactions instead of three). Shown because the promoter interaction moves sharply between the
specifications: `sv_TE:bio1c` is itself non-zero and same-signed, so removing `sv_TE` hands its signal to
the correlated promoter term — a bias, not a gain in precision.

Four fits: TE kept / TE dropped × REDUCED / FULL. The thing to read is that the **coding-SV burden effect
is stable across all four**, while the promoter interaction depends on both choices."""

code_te = r"""
ARMS=[("TE kept, REDUCED",RED,"o",True),("TE kept, FULL",FULL,"o",False),
      ("TE dropped, REDUCED",REDN,"s",True),("TE dropped, FULL",FULLN,"s",False)]
RT=["sv_CDS","sv_promoter"]
fig,axs=plt.subplots(1,2,figsize=(11.6,3.4),gridspec_kw=dict(wspace=0.34))
a=axs[0]; style(a); panel(a,"a",x=-0.22); a.axvline(0,color="0.35",lw=1,ls="--")
for i,t in enumerate(RT):
    for j,(lab,src,mk,hollow) in enumerate(ARMS):
        e,se=src["co"].loc[t,"estimate"]*PV(t), src["co"].loc[t,"se"]*PV(t)
        a.errorbar(e,i+0.27-0.18*j,xerr=1.96*se,fmt=mk,color=PC2(t),
                   mfc="white" if hollow else PC2(t),mew=1.1,ms=5.5,capsize=2,zorder=3)
    if i<len(RT)-1: a.axhline(i+0.5,color="0.93",lw=0.8)
a.set_yticks(range(len(RT))); a.set_yticklabels([LAB2(t) for t in RT]); a.invert_yaxis()
a.set_xlabel("burden effect at an average-temperature garden\n(% per carried variant, ± 95% CI)",fontsize=7.5)
a=axs[1]; style(a,axis="y"); panel(a,"b",x=-0.16)
xs=np.arange(len(RT)); w=0.2
for j,(lab,src,mk,hollow) in enumerate(ARMS):
    a.bar(xs+(j-1.5)*w,[src["BONF"][t] for t in RT],width=w,
          color=["white" if hollow else PC2(t) for t in RT],edgecolor=[PC2(t) for t in RT],lw=1.1,
          hatch=None if mk=="o" else "////",zorder=3)
    for i,t in enumerate(RT):
        a.annotate(f"{src['BONF'][t]:.2g}",xy=(i+(j-1.5)*w,src["BONF"][t]),xytext=(2,4),
                   textcoords="offset points",ha="center",fontsize=5.6,color="0.25",rotation=90)
a.set_yscale("log"); a.axhline(0.05,color="0.35",lw=1,ls="--"); a.set_ylim(top=30)
a.annotate("α = 0.05",xy=(0.005,0.05),xycoords=("axes fraction","data"),xytext=(0,3),
           textcoords="offset points",ha="left",fontsize=7,color="0.35")
a.set_xticks(xs); a.set_xticklabels([LAB2(t).replace(", ","\n")+"\n× temp." for t in RT],fontsize=6.5)
a.set_ylabel("burden × temperature:\nLRT p, Bonferroni",fontsize=8.5)
h_=[plt.Line2D([],[],marker=mk,color="0.4",mfc="white" if hollow else "0.4",ls="none",ms=5.5)
    for _,_,mk,hollow in ARMS]
fig.legend(h_,[l for l,_,_,_ in ARMS],loc="lower center",bbox_to_anchor=(0.5,-0.17),ncol=4,fontsize=7.5,
           handlelength=1.2,columnspacing=1.4)
save(fig,"fullmodel_te_sensitivity",nb=True); plt.show()
print("all four fits (pooled, no origin filter):\n")
print(pd.DataFrame([dict([("fit",lab)]+
    [(f"{LAB2(t)}",f"{src['co'].loc[t,'estimate']*PV(t):+.3f} (p={src['co'].loc[t,'p']:.3g})") for t in RT]+
    [(f"{LAB2(t)} × temp. Bonf",f"{src['BONF'][t]:.4g}") for t in RT])
    for lab,src,_,_ in ARMS]).to_string(index=False))
print(f"\nwhat sv_TE carries in the TE-kept FULL fit:")
print(f"  sv_TE          {FULL['co'].loc['sv_TE','estimate']*PV('sv_TE'):+.3f} %/variant "
      f"(p = {FULL['co'].loc['sv_TE','p']:.3g})")
print(f"  sv_TE:bio1c    {FULL['co'].loc['sv_TE:bio1c','estimate']*PV('sv_TE'):+.4f} %/variant/degC "
      f"(Bonferroni {FULL['BONF']['sv_TE']:.3g})")
"""

md_tab = """## 7 · The complete coefficient tables

Every fixed effect of both fits, on the model's own scale — per SD of the standardised rate, which is what
R reports — rather than the per-variant display scale used in the figures. Nothing is omitted: intercept,
SNP burden, all four SV contexts, assembly route, every interaction, the local-adaptation terms, the LRTs,
and the same for the TE-dropped arm."""

code_tab = r"""
def dump(red,full,terms,label,nb):
    print("="*98); print(label)
    ix=list(dict.fromkeys(list(red["co"].index)+list(full["co"].index)))
    g=lambda src,t,f: (f"{src['co'].loc[t,f]:+.5f}" if f!="p" else f"{src['co'].loc[t,'p']:.3g}") \
                      if t in src["co"].index else "—"
    print(pd.DataFrame([{"term":t,
        "REDUCED est":g(red,t,"estimate"),"REDUCED se":g(red,t,"se"),"REDUCED p":g(red,t,"p"),
        "FULL est":g(full,t,"estimate"),"FULL se":g(full,t,"se"),"FULL p":g(full,t,"p")}
        for t in ix]).to_string(index=False))
    print(f"\nLikelihood-ratio tests for the burden × temperature terms (1 df each, Bonferroni ×{nb}):")
    print(pd.DataFrame([{"term":t,"REDUCED LRT p":f"{red['LRT'][t]:.5g}","REDUCED Bonf":f"{red['BONF'][t]:.5g}",
                         "FULL LRT p":f"{full['LRT'][t]:.5g}","FULL Bonf":f"{full['BONF'][t]:.5g}"}
                        for t in terms]).to_string(index=False))
    print()
dump(RED,FULL,SVT6,"PRIMARY ARM — pooled, no origin filter, reference TE kept",len(SVT6))
dump(REDN,FULLN,SVT4,"SENSITIVITY ARM — pooled, no origin filter, reference TE dropped",len(SVT4))
"""

md_read = """## 8 · Reading it

In the FULL model — burdens, burden × temperature and local adaptation fitted together — the two kinds of
burden term behave completely differently.

**The coding-SV burden effect stands.** −1.49 % of founder fitness per carried coding SV at an
average-temperature garden, 95 % CI [−2.34, −0.65], p = 5.5 × 10⁻⁴. Adding the local-adaptation terms moved
it by 0.22 of a standard error. It is negative in **31 of 31 gardens** (§5) and stable across all four
specifications in §6. This is the result of the model.

**No burden × temperature interaction stands.** All three fall below significance once local adaptation is
in the model: coding Bonferroni 0.00039 → 0.13, reference TE 0.055 → 0.28, promoter 0.050 → 0.28, with
estimates dropping 29–50 %.

**And they fall because the estimates move, not because the errors grow.** §3c separates the two: mean
absolute shift is **0.13 SE for the burden effects and 1.16 SE for the burden × temperature terms** — an
order of magnitude apart — while every standard error changes by about 1 %. If local adaptation were
merely competing for variance, the standard errors would widen and the estimates would hold. The opposite
happens, which is the signature of confounding: the local-adaptation terms take back signal the
interactions had been carrying.

**Local adaptation is the largest effect in the model.** `obio1c:bio1c` has p ≈ 3 × 10⁻³², and on a common
per-SD scale (§4) it is larger than any burden interaction.

**The promoter interaction survives only when `sv_TE` is dropped** (Bonferroni 0.012 there, 0.28 here).
Since `sv_TE:bio1c` is non-zero and same-signed, dropping `sv_TE` reallocates its signal to the correlated
promoter term. The TE-kept specification is the better-behaved one, and the promoter interaction does not
stand up in it.

### What this does and does not establish

⚠ **The local-adaptation terms test the interactions much harder than they test the burden effects.**
`obio1c:bio1c` is an interaction, so it is near-orthogonal to the burden main effects and can only absorb
founder × garden-climate structure. The main-effect counterpart is `obio1c` itself, which *is* in the FULL
fit (p = 0.036) and which moves `sv_CDS` by only 0.22 SE. That is a real test, but a milder one: the coding
effect is robust to founder provenance *as far as this design can tell*, not established to the same
standard as the interactions are refuted.

**This design cannot separate confounding from mechanism.** "The burden × temperature terms are confounded
by provenance" and "the burden is one of the routes by which provenance produces local adaptation" make
identical predictions here — in both, adding `obio1c:bio1c` absorbs the burden signal. Nothing fitted at
the founder level can distinguish them, because both live entirely at the founder level. Breaking that
needs evidence at a different unit — per-variant, or within-founder.

**`obio1c` is one summary of origin.** Mean annual temperature is a single axis; provenance is
multidimensional, so the control is incomplete and a residual gradient can remain in either direction.

### Suggested reporting

From the FULL model: **carrying an SV in coding sequence costs about 1.5 % of founder fitness per variant,
and this holds with founder origin climate and its garden-climate interaction in the model.** Do not report
a climate-dependent SV burden effect from this model — no burden × temperature term survives once local
adaptation is fitted alongside it."""

md_marg = """## 9 · Appendix — why `obio1c` is kept alongside `obio1c:bio1c`

If the point of the control is the *interaction*, does the model also need origin climate as a main effect?
Refitting with `EXTRA_FIXED="obio1c:bio1c"` alone answers it directly.

**Nothing that this notebook concludes depends on the choice.** The design is balanced — all 231 founders in
all 31 gardens — so `obio1c`, which is constant within a founder, is **orthogonal to every within-founder
contrast**. Dropping it therefore cannot move a single interaction: the local-adaptation term and all three
burden × temperature terms come back bit-identical. Only the burden *main effects* shift, and by at most a
quarter of a standard error.

**It is kept for three reasons.** (i) Marginality: fitting an interaction without its main effect constrains
that main effect to exactly zero, which here is a substantive claim — that origin climate has no effect on
fitness at an average-temperature garden — and the data mildly reject it (LRT p = 0.036). (ii) It costs one
degree of freedom. (iii) It is the **conservative** choice for this notebook's headline: `obio1c` absorbs a
little between-founder signal, so `sv_CDS` is *smaller* with it in (−1.49 %) than without (−1.59 %).
Reporting the number with `obio1c` in the model is the harder test of the coding-SV claim."""

code_marg = r'''
SUF2="_x_obio1cbio1c"                                       # EXTRA_FIXED="obio1c:bio1c" -- interaction only
ONLY=load(PA,SVT6,SUF2)
llml=lambda suf: pd.read_csv(f"{PA}/loglik_M2p_ml{suf}.csv").logLik[0]
chi=2*(llml(SUF)-llml(SUF2)); pm=stats.chi2.sf(chi,1)
print("FULL (obio1c + obio1c:bio1c)   vs   INTERACTION-ONLY (obio1c:bio1c, no main effect)\n")
print("burden effect at an average-temperature garden, % per carried variant:")
print(pd.DataFrame([{"term":LAB2(t),
    "FULL":f"{FULL['co'].loc[t,'estimate']*PV(t):+.4f}","FULL p":f"{FULL['co'].loc[t,'p']:.3g}",
    "interaction-only":f"{ONLY['co'].loc[t,'estimate']*PV(t):+.4f}","p":f"{ONLY['co'].loc[t,'p']:.3g}",
    "shift (FULL SEs)":f"{(ONLY['co'].loc[t,'estimate']-FULL['co'].loc[t,'estimate'])/FULL['co'].loc[t,'se']:+.3f}"}
    for t in BURD]).to_string(index=False))
print(f"\nburden × temperature (Bonferroni ×{len(SVT6)}):")
print(pd.DataFrame([{"term":LAB2(t),
    "FULL %/°C":f"{FULL['co'].loc[t+':bio1c','estimate']*PV(t):+.5f}","FULL Bonf":f"{FULL['BONF'][t]:.4g}",
    "interaction-only %/°C":f"{ONLY['co'].loc[t+':bio1c','estimate']*PV(t):+.5f}","Bonf":f"{ONLY['BONF'][t]:.4g}"}
    for t in SVT6]).to_string(index=False))
print(f"\nthe local-adaptation interaction itself:")
for lab,src in (("FULL",FULL),("interaction-only",ONLY)):
    print(f"  {lab:<18} {OBX} {src['co'].loc[OBX,'estimate']:+.6f} +- {src['co'].loc[OBX,'se']:.6f}  "
          f"p = {src['co'].loc[OBX,'p']:.3g}")
print(f"  obio1c main effect (FULL only): {FULL['co'].loc['obio1c','estimate']:+.6f} "
      f"+- {FULL['co'].loc['obio1c','se']:.6f}  p = {FULL['co'].loc['obio1c','p']:.3g}")
mx=max(abs(ONLY["co"].loc[t,"estimate"]-FULL["co"].loc[t,"estimate"])/FULL["co"].loc[t,"se"]
       for t in BURD+[f"{u}:bio1c" for u in SVT6])
print(f"\nlargest change in ANY burden coefficient when obio1c is dropped: {mx:.3f} standard errors")
print(f"largest change among the INTERACTIONS only: "
      f"{max(abs(ONLY['co'].loc[u+':bio1c','estimate']-FULL['co'].loc[u+':bio1c','estimate'])/FULL['co'].loc[u+':bio1c','se'] for u in SVT6):.2e} SEs")
print("  -> exactly zero to numerical precision: a between-founder covariate cannot move a within-founder")
print("     contrast in a balanced 231 x 31 design.")
print(f"\nis dropping obio1c justified? LRT, FULL vs interaction-only (ML, 1 df): "
      f"chi2 = {chi:.3f}, p = {pm:.4g}")
print("  -> keeping it is mildly favoured, costs 1 df, and is the conservative choice for sv_CDS.")
'''

md_methods = """## As run

- **Burdens.** `results/sv_adaptive/founder_burden_event.csv` — truvari-collapsed SV **events**, insertions
  and deletions summed within a context, origin class `all` (no `shared` filter), frequency class `all`.
  Carried = `ins_all_<ctx>_all + del_all_<ctx>_all`; called = the matching `called_*` columns. Rates are
  z-scored across the 231 founders inside `_sv_burden_pooled_models.R`.
- **REDUCED.** `_sv_burden_pooled_models.R` with
  `POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all`, tasks 1–5
  (M2p REML, M2p ML, and the three single-interaction drops for the LRTs).
- **FULL.** The same plus `EXTRA_FIXED="obio1c + obio1c:bio1c"`, submitted by
  `_submit_pooled_sv6_origctrl.sh`. Its outputs carry the `_x_obio1cobio1cbio1c` suffix, so REDUCED is not
  overwritten. `EXTRA_FIXED` is empty by default, so every other invocation of the script is unchanged.
- **Sensitivity arm.** Both fits again with `POOL_SET=sv4 FINAL_OUT=final_models_pooled_all_noTE` (two
  interactions, Bonferroni ×2).
- **The main-effects-only fit (`M1p`) is not read.** Burden effects are reported from the FULL model, where
  `bio1c` is centred on the mean of the 31 gardens, so `sv_ctx` is the burden's effect at an
  average-temperature garden (13.5 °C). Reporting it from a separate no-interaction fit would be a
  different model and is not done here.
- **`obio1c` vs `bio1c`.** `obio1c` is the **founder's origin** bio1 from `founder_burden_event.csv`,
  centred over the 231 founders and **not** standardised — hence the rescaling by the founder SD in §4.
  `bio1c` is the **garden's** bio1, centred over the 31 gardens. lme4 names the interaction `bio1c:obio1c`
  by model-frame order; the load cell resolves the spelling rather than hard-coding it.
- **LRTs.** Read from the named `loglik_*.csv` files, never by globbing: these directories also hold
  `SUBSET=` and other `EXTRA_FIXED` runs, and a glob would mix likelihoods across different models.
- **Per-variant scaling.** `sd_rate × med_called` from `pooled_rate_scale.csv`, ×100. A display scaling
  only — it enters no fit and changes no p-value. The two arms' scale files are byte-identical, so one
  scale serves both and the arms are directly comparable.
- **Nothing is fitted in this notebook** — every number is read from the R output."""

if __name__ == "__main__":
    build("sv_burden_origin_climate_control.ipynb", [
        new_markdown_cell(md_title), new_code_cell(code_setup), new_code_cell(code_load),
        new_markdown_cell(md_counts), new_code_cell(code_counts),
        new_markdown_cell(md_full), new_code_cell(code_full),
        new_markdown_cell(md_cmp), new_code_cell(code_cmp),
        new_markdown_cell(md_scale), new_code_cell(code_scale),
        new_markdown_cell(md_garden), new_code_cell(code_garden),
        new_markdown_cell(md_te), new_code_cell(code_te),
        new_markdown_cell(md_tab), new_code_cell(code_tab),
        new_markdown_cell(md_read),
        new_markdown_cell(md_marg), new_code_cell(code_marg),
        new_markdown_cell(md_methods),
    ])
