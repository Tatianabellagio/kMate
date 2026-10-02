#!/usr/bin/env python
"""Build+execute `sv_burden_diagnostics.ipynb`: what the model's terms actually look like.

Not a results notebook. It asks the questions that sit underneath the fitted coefficients:

  1. How are the six burdens distributed across founders? z-scoring centres and scales them, it does
     not make them normal. Skew, tails and outliers are shown as they are.
  2. How do the burdens relate to each other pairwise, and is the relationship linear?
  3. How do the FITTED per-garden effects of the six contexts relate to each other, and to garden
     temperature?

The third asks a trap question, and it is handled explicitly: each garden's TOTAL effect is
`fixed_part + slope_blup`, and `fixed_part` is a deterministic straight line in bio1c. Correlating two
contexts' totals therefore mostly correlates two straight lines in the same variable. The pairwise panels
use the BLUP deviations, which is the part that is actually free to vary per garden; the vs-temperature
panels use the totals, because there the line is the question.

Reads the reported all6 fit. `basic` env.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_fig_common import build, new_markdown_cell, new_code_cell, code_setup_burden as code_setup

md_title = """# What the model's terms actually look like

The fitted coefficients assume things about their inputs. This notebook looks at the inputs.

Three questions:

1. **Are the burdens normally distributed?** They are z-scored, which centres and scales them — it does
   **not** make them normal. If a burden has a long tail, a handful of founders carry a lot of the
   leverage behind its coefficient.
2. **How do the burdens relate pairwise?** We know they correlate 0.83–0.97. Whether that relationship is
   *linear* is a separate question, and a curved one would mean the linear model is mis-specified.
3. **How do the fitted per-garden effects behave** — against each other, and against garden temperature?

⚠ **One trap, handled explicitly in §3.** Each garden's total effect is `fixed_part + slope_blup`, and
`fixed_part` is by construction a straight line in `bio1c`. Correlating two contexts' *totals* therefore
mostly correlates two straight lines in the same variable, and will look impressive whatever is true. The
pairwise panels below use the **BLUP deviations** — the part genuinely free to vary per garden. The
vs-temperature panels use the totals, because there the line is exactly what is being asked about."""

code_load = r"""
PA=f"{G}/final_models_pooled_all"
SUF="_x_obio1cobio1cbio1c_re0obio1cgarden_all6_refid_intall"
CTX=["CDS","UTR","intron","TE","promoter","intergenic"]
CC={"CDS":"#2F4A1F","UTR":"#4E6B3A","intron":"#7E9B6B","TE":"#C46A10",
    "promoter":"#5B7A99","intergenic":"#B9B9B9"}
CO=pd.read_csv(f"{PA}/coef_M2p_reml{SUF}.csv").set_index("term")
GD=pd.read_csv(f"{PA}/garden_M2p_reml{SUF}.csv")
SC=pd.read_csv(f"{PA}/pooled_rate_scale_x_obio1cobio1cbio1c_all6.csv").set_index("context")
PVv={f"sv_{c}":100/(SC.loc[c,"sd_rate"]*SC.loc[c,"med_called"]) for c in SC.index}
B=pd.read_csv(f"{G}/founder_burden_event.csv",dtype={"founder":str})
L=pd.read_csv(f"{G}/onestage_all/long.csv",dtype={"founder":str})
TECH=L.drop_duplicates("founder").set_index("founder").cactus
car={c:(B[f"ins_all_{c}_all"]+B[f"del_all_{c}_all"]).astype(float) for c in CTX}
cal={c:(B[f"called_ins_all_{c}_all"]+B[f"called_del_all_{c}_all"]).astype(float) for c in CTX}
RAW={c:(car[c]/cal[c]).to_numpy() for c in CTX}                    # the carried/called rate
zz=lambda v:(v-v.mean())/v.std(ddof=1)
Z={c:zz(RAW[c]) for c in CTX}                                      # exactly what the model fits
CARR={c:car[c].to_numpy() for c in CTX}
print(f"{len(B):,} founders, {len(CTX)} burdens")
print("the model fits the z-scored rate; z-scoring centres and scales, it does not normalise")
"""

md_dist = """## 1 · How is each burden distributed?

Each panel is one burden as the model sees it — the z-scored `carried / called` rate across the 231
founders. The dashed curve is the normal with the same mean and SD, so departures are visible directly.
The annotation gives skew, excess kurtosis, and a Shapiro–Wilk test of normality.

The second row is the same burdens as **raw carried counts**, which is what a biologist would picture, and
is on a very different scale per context."""

code_dist = r"""
fig,axs=plt.subplots(2,len(CTX),figsize=(2.42*len(CTX),4.9),
                     gridspec_kw=dict(hspace=0.62,wspace=0.30))
for j,c in enumerate(CTX):
    a=axs[0,j]; style(a,axis="both")
    v=Z[c]
    a.hist(v,bins=26,density=True,color=CC[c],alpha=0.75,lw=0)
    xs=np.linspace(v.min()-0.5,v.max()+0.5,200)
    a.plot(xs,stats.norm.pdf(xs,v.mean(),v.std(ddof=1)),color="0.25",lw=1.2,ls="--")
    sk=stats.skew(v); ku=stats.kurtosis(v); W,pw=stats.shapiro(v)
    a.annotate(f"skew {sk:+.2f}\nkurt {ku:+.2f}\nSW p {pw:.2g}",xy=(0.96,0.95),
               xycoords="axes fraction",ha="right",va="top",fontsize=6.6,
               color="#B04A00" if pw<0.05 else "0.35")
    a.annotate(c,xy=(0.5,1.10),xycoords="axes fraction",ha="center",fontsize=8.5,
               color=CC[c],fontweight="bold")
    a.set_yticks([]); a.tick_params(labelsize=6.5)
    if j==0: a.set_ylabel("z-scored rate\n(what the model fits)",fontsize=7.5)
    a=axs[1,j]; style(a,axis="both")
    a.hist(CARR[c],bins=26,color=CC[c],alpha=0.75,lw=0)
    a.annotate(f"median {np.median(CARR[c]):,.0f}",xy=(0.96,0.95),xycoords="axes fraction",
               ha="right",va="top",fontsize=6.6,color="0.35")
    a.set_yticks([]); a.tick_params(labelsize=6.5)
    a.set_xlabel("events carried",fontsize=7)
    if j==0: a.set_ylabel("raw carried count",fontsize=7.5)
save(fig,"diag_burden_distributions",nb=True); plt.show()
print("DISTRIBUTION OF EACH BURDEN ACROSS THE 231 FOUNDERS\n")
rows=[]
for c in CTX:
    v=Z[c]; W,pw=stats.shapiro(v)
    rows.append({"context":c,"skew":f"{stats.skew(v):+.3f}","excess kurtosis":f"{stats.kurtosis(v):+.3f}",
                 "Shapiro-Wilk p":f"{pw:.3g}","normal at 0.05":"no" if pw<0.05 else "yes",
                 "min z":f"{v.min():+.2f}","max z":f"{v.max():+.2f}",
                 "|z|>3":int((np.abs(v)>3).sum()),
                 "median carried":f"{np.median(CARR[c]):,.0f}"})
print(pd.DataFrame(rows).to_string(index=False))
print("\n  Skew > 0 means a right tail: a few founders carry far more than the rest. Those founders")
print("  have the most leverage on that burden's coefficient -- which is what the jackknife and")
print("  delete-d sections of sv_burden_model.ipynb test directly.")
"""

md_pairs = """## 2 · How do the burdens relate to each other?

Lower triangle: every pair of burdens, z-scored, with the fitted line. Upper triangle: the correlation.
Diagonal: the burden against itself as a reference.

The question is not *whether* they correlate — we know they do — but whether the relationship is
**linear**. A curved cloud would mean the model's additive, linear form is wrong, and no amount of
re-parameterising would fix it. Points are coloured by sequencing technology, because that is the
covariate most likely to produce structure in these clouds."""

code_pairs = r"""
n=len(CTX); tc=TECH.loc[B.founder].to_numpy(float)
fig,axs=plt.subplots(n,n,figsize=(1.62*n,1.62*n),gridspec_kw=dict(hspace=0.12,wspace=0.12))
for i in range(n):
    for j in range(n):
        a=axs[i,j]; a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values(): sp.set_linewidth(0.4); sp.set_color("0.8")
        ci,cj=CTX[i],CTX[j]
        if i==j:
            a.hist(Z[ci],bins=18,color=CC[ci],alpha=0.8,lw=0)
            a.annotate(ci,xy=(0.5,0.5),xycoords="axes fraction",ha="center",va="center",
                       fontsize=8,color=CC[ci],fontweight="bold")
        elif i>j:
            a.scatter(Z[cj][tc==0],Z[ci][tc==0],s=5,color="0.55",alpha=0.55,linewidth=0)
            a.scatter(Z[cj][tc==1],Z[ci][tc==1],s=5,color="#C46A10",alpha=0.65,linewidth=0)
            m,b0=np.polyfit(Z[cj],Z[ci],1); xs=np.array([Z[cj].min(),Z[cj].max()])
            a.plot(xs,m*xs+b0,color="0.2",lw=0.9)
        else:
            r=np.corrcoef(Z[ci],Z[cj])[0,1]
            a.annotate(f"{r:.2f}",xy=(0.5,0.5),xycoords="axes fraction",ha="center",va="center",
                       fontsize=9+6*(r-0.8),color="0.2")
            a.set_facecolor(plt.get_cmap("Oranges")((r-0.78)/0.24))
        if j==0: a.set_ylabel(ci,fontsize=7,color=CC[ci])
        if i==n-1: a.set_xlabel(cj,fontsize=7,color=CC[cj])
fig.legend([plt.Line2D([],[],marker="o",ls="none",color="0.55",ms=5),
            plt.Line2D([],[],marker="o",ls="none",color="#C46A10",ms=5)],
           ["short-read","cactus assembly"],frameon=False,fontsize=7.5,
           loc="upper center",bbox_to_anchor=(0.5,0.055),ncol=2)
save(fig,"diag_burden_pairs",nb=True); plt.show()
print("PAIRWISE BURDEN CORRELATIONS (z-scored rate, 231 founders)\n")
M=pd.DataFrame({c:Z[c] for c in CTX}).corr()
print(M.round(3).to_string())
print("\n\nIS ANY PAIR NONLINEAR?  (Spearman vs Pearson -- a gap means curvature or outlier leverage)\n")
rows=[]
for i in range(n):
    for j in range(i+1,n):
        pe=np.corrcoef(Z[CTX[i]],Z[CTX[j]])[0,1]
        sp=stats.spearmanr(Z[CTX[i]],Z[CTX[j]]).correlation
        rows.append({"pair":f"{CTX[i]} vs {CTX[j]}","Pearson":f"{pe:+.3f}",
                     "Spearman":f"{sp:+.3f}","difference":f"{sp-pe:+.3f}"})
D_=pd.DataFrame(rows); D_["absdiff"]=D_.difference.astype(float).abs()
print(D_.sort_values("absdiff",ascending=False).drop(columns="absdiff").head(8).to_string(index=False))
print(f"\n  largest |Spearman - Pearson| across all {len(rows)} pairs: "
      f"{D_.absdiff.max():.3f}  -- small gaps mean the linear form is not obviously wrong.")
"""

md_garden = """## 3 · The fitted per-garden effects

Each garden gets its own effect for each burden: `total = fixed_part + slope_blup`, where `fixed_part` is
the fixed effect plus its temperature term and `slope_blup` is that garden's random deviation.

**(a) Against garden temperature**, using the **totals**. This is the interaction, drawn: a sloped line
means the burden's effect depends on temperature.

**(b) Against each other**, using the **BLUP deviations only**. This is the part that matters, and it is
why the panel is not drawn from the totals: two totals both contain a straight line in `bio1c`, so
correlating them would mostly measure that shared line rather than anything about the gardens. The
deviations have the fixed part removed, so a correlation between them means two contexts genuinely move
together across gardens beyond what temperature explains."""

code_garden = r"""
P=GD[GD.term.str.startswith("sv_")].copy()
P["ctx"]=P.term.str.replace("sv_","",regex=False)
TOT=P.pivot(index="garden",columns="ctx",values="total")
BLU=P.pivot(index="garden",columns="ctx",values="slope_blup")
BIO=P.drop_duplicates("garden").set_index("garden").bio1.reindex(TOT.index)
for c in CTX: TOT[c]=TOT[c]*PVv[f"sv_{c}"]; BLU[c]=BLU[c]*PVv[f"sv_{c}"]
fig,axs=plt.subplots(2,3,figsize=(11.8,6.2),gridspec_kw=dict(hspace=0.42,wspace=0.34))
for k,c in enumerate(CTX):
    a=axs[k//3,k%3]; style(a,axis="both")
    a.scatter(BIO,TOT[c],s=26,color=CC[c],edgecolor="0.3",linewidth=0.4,zorder=3)
    m,b0=np.polyfit(BIO,TOT[c],1); xs=np.array([BIO.min(),BIO.max()])
    a.plot(xs,m*xs+b0,color=CC[c],lw=1.6)
    a.axhline(0,color="0.5",lw=0.9,ls="--")
    r,pv=stats.pearsonr(BIO,TOT[c])
    a.annotate(f"{c}\nslope {m:+.4f} %/degC\nr = {r:+.2f}",xy=(0.03,0.97),xycoords="axes fraction",
               va="top",fontsize=7.2,color="0.25")
    a.set_xlabel("garden bio1 (degC)",fontsize=7.5)
    a.set_ylabel("% per carried variant",fontsize=7.5)
save(fig,"diag_garden_vs_temp",nb=True); plt.show()
n=len(CTX)
fig,axs=plt.subplots(n,n,figsize=(1.62*n,1.62*n),gridspec_kw=dict(hspace=0.12,wspace=0.12))
for i in range(n):
    for j in range(n):
        a=axs[i,j]; a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values(): sp.set_linewidth(0.4); sp.set_color("0.8")
        ci,cj=CTX[i],CTX[j]
        if i==j:
            a.annotate(ci,xy=(0.5,0.5),xycoords="axes fraction",ha="center",va="center",
                       fontsize=8,color=CC[ci],fontweight="bold")
        elif i>j:
            a.scatter(BLU[cj],BLU[ci],s=13,color="0.4",alpha=0.75,linewidth=0)
            if BLU[cj].std()>1e-12 and BLU[ci].std()>1e-12:
                m,b0=np.polyfit(BLU[cj],BLU[ci],1); xs=np.array([BLU[cj].min(),BLU[cj].max()])
                a.plot(xs,m*xs+b0,color="#B04A00",lw=0.9)
        else:
            if BLU[ci].std()<1e-12 or BLU[cj].std()<1e-12:
                a.annotate("n/a",xy=(0.5,0.5),xycoords="axes fraction",ha="center",va="center",
                           fontsize=7,color="0.6")
            else:
                r=np.corrcoef(BLU[ci],BLU[cj])[0,1]
                a.annotate(f"{r:+.2f}",xy=(0.5,0.5),xycoords="axes fraction",ha="center",va="center",
                           fontsize=9,color="0.2")
                a.set_facecolor(plt.get_cmap("RdBu_r")(0.5+r/2))
        if j==0: a.set_ylabel(ci,fontsize=7,color=CC[ci])
        if i==n-1: a.set_xlabel(cj,fontsize=7,color=CC[cj])
save(fig,"diag_garden_blup_pairs",nb=True); plt.show()
print("PER-GARDEN EFFECT vs GARDEN TEMPERATURE  (totals -- the interaction, drawn)\n")
print(pd.DataFrame([{"context":c,"slope %/degC":f"{np.polyfit(BIO,TOT[c],1)[0]:+.5f}",
                     "r":f"{stats.pearsonr(BIO,TOT[c])[0]:+.3f}",
                     "p":f"{stats.pearsonr(BIO,TOT[c])[1]:.3g}",
                     "range across gardens":f"{TOT[c].min():+.3f} .. {TOT[c].max():+.3f}"}
                    for c in CTX]).to_string(index=False))
print("\n  These r values are near 1 by construction -- `total` contains the fixed temperature term,")
print("  which IS a straight line in bio1c. The slope is the content; the correlation is not.\n")
print("\nPER-GARDEN BLUP DEVIATIONS, pairwise  (the fixed line removed)\n")
sd=BLU.std()
dead=[c for c in CTX if sd[c]<1e-9]
if dead: print(f"  at the variance boundary, no deviation to correlate: {', '.join(dead)}\n")
live=[c for c in CTX if c not in dead]
print(BLU[live].corr().round(3).to_string())
print("\n  A positive entry means two contexts are jointly stronger in the same gardens, beyond what")
print("  garden temperature already explains.")
"""

md_read = """## 4 · Reading it

**On the distributions (§1).** The burdens are not normal, and the model does not require them to be —
linear mixed models assume normality of the *residuals*, not of the covariates. A skewed covariate is
fine; what it costs is leverage concentration, since the founders in the tail carry disproportionate
weight in that coefficient. That is exactly what §10–§11 of `sv_burden_model.ipynb` test by refitting
without them, and the estimate does not move.

**On the pairwise structure (§2).** The correlations are high but the relationships are close to linear —
Spearman and Pearson agree closely across all 15 pairs. If they had diverged, the additive linear model
would have been the wrong functional form and no re-parameterisation would have rescued it. They don't,
so the issue with these burdens is collinearity, not mis-specification. Those are different problems with
different fixes.

**On the per-garden effects (§3).** The vs-temperature panels are the interactions drawn out, and their
correlations are near 1 **by construction** — the total contains the fixed temperature term, which is a
straight line in `bio1c`. Only the slope carries information there; the correlation is an artifact of what
is being plotted, and would look equally impressive under a null.

The pairwise BLUP panel is the one that can surprise. Contexts whose garden deviations correlate are
jointly stronger in the same gardens for a reason temperature does not capture — which is a hint about an
unmodelled garden-level variable, not about the contexts themselves.

⚠ Burdens whose garden random slope sat at the variance boundary have no deviation to correlate, and are
marked `n/a` rather than being shown as a spurious zero."""

md_methods = """## As run

Everything is read from the reported fit — `coef_M2p_reml`, `garden_M2p_reml` with suffix
`_x_obio1cobio1cbio1c_re0obio1cgarden_all6_refid_intall` — and from `founder_burden_event.csv`. Nothing
is refitted.

- **Burdens** are `carried / called` per context, z-scored across the 231 founders, exactly as the model
  receives them. The raw carried counts are shown alongside because the z-score hides the scale.
- **Per-garden effects** are `total = fixed_part + slope_blup` from the R output, converted to % per
  carried variant with the arm's `pooled_rate_scale`.
- **Why the pairwise panel uses BLUPs and the temperature panel uses totals** is in §3: `fixed_part` is a
  deterministic line in `bio1c`, so a total-vs-total correlation is largely that line appearing twice.
- Shapiro–Wilk at n = 231 is sensitive — it will reject mild departures. The skew and kurtosis values are
  the more useful summary; the test is reported for completeness, not as a gate."""

if __name__ == "__main__":
    build("sv_burden_diagnostics.ipynb", [
        new_markdown_cell(md_title), new_code_cell(code_setup), new_code_cell(code_load),
        new_markdown_cell(md_dist),   new_code_cell(code_dist),
        new_markdown_cell(md_pairs),  new_code_cell(code_pairs),
        new_markdown_cell(md_garden), new_code_cell(code_garden),
        new_markdown_cell(md_read),   new_markdown_cell(md_methods),
    ])
