#!/usr/bin/env python
"""Build + execute the 'predict ecotype performance: kinship vs climate-match' notebook.
Standalone; load-only (precomputed json/npz/csv from predict_ecotype_performance.py) -> runs in
the `basic` env. Model math + variance partition + leave-one-site-out CV + permutation null."""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/notebooks/predict_ecotype_performance.ipynb"

md_title = """# Predicting ecotype performance: "who you are" (kinship) vs climate-match

Can we predict how a founder ecotype performs in a common garden from **who it is** (genome-wide
relatedness / generalist genetic merit) alone, or does **climate match** (the founder's home
climate vs the garden's climate) add predictive power — the signature of climate-driven local
adaptation *beyond* clade?"""

md_model = r"""## The model

A **GxE reaction-norm variance partition**. The trait is per-founder selection
$S_{f,\text{site}}$ — the flower-weighted, plot-averaged slope of the founder's genome-wide
frequency $h_f$ over generations $0\!-\!3$ at a site (large $\Rightarrow$ that ecotype *won*).

$$S_{f,\text{site}}=\underbrace{\text{site}_\text{site}}_{\text{garden mean}}
+\underbrace{u_f}_{\text{kinship: "who you are"}}
+\underbrace{\gamma\cdot\text{match}(o_f,c_\text{site})}_{\text{local adaptation}}
+\varepsilon,\qquad u_f\sim\mathcal N(0,\sigma_g^2\mathbf K).$$

- **$\mathbf K$** = founder genomic-relationship matrix, built from **common haplotype clusters
  only (MAF $\ge$ 5%)** — without the floor the $1/\sqrt{p(1-p)}$ standardisation lets a *private*
  haplotype (1 founder) dominate $\mathbf K$ and inflate relatedness.
- **$u_f$** is shared across all sites, so it encodes a single **generalist** value — the same
  ranking at every garden. It can only capture "who wins everywhere".
- **match** $=-(z(o_f)-z(c_\text{site}))^2$ (home-field advantage; $\gamma>0$ = founders win in
  gardens that resemble home). Computed for **bio1** (annual temp) and a **multivariate** version
  over all 19 bioclim variables. It is the *only* term that can re-rank founders for a new climate.
- Origin climate $o_f$ = founder home WorldClim bio1 (all 231 founders matched by ecotype id);
  garden climate $c_\text{site}$ = site ERA5 bio1. Low-rank REML (Woodbury; founder effect shared
  across sites).

**The decisive test — leave-one-site-out CV.** Hold out an entire garden; predict its founders'
performance from the other 16. *Model A* (kinship only) predicts the generalist BLUP $\hat g_f$ —
the same ranking everywhere. *Model B* adds $\hat\gamma\cdot\text{match}(o_f,c_\text{held-out})$ —
re-ranking founders for the held-out climate. If **B > A** at a *never-seen* garden, climate
adaptation is real beyond clade. A permutation null (shuffle founder origin-climate labels:
preserves kinship & marginal performance, breaks only the match) calibrates the climate gain."""

code_setup = '''import json, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
H="/global/scratch/users/tbellg/kmate/results/grenenet_gea/hapfreq"
d=json.load(open(f"{H}/predict_ecotype_performance.json"))
arr=np.load(f"{H}/predict_ecotype_performance_arrays.npz")
tab=pd.read_csv(f"{H}/predict_ecotype_performance.csv")
print(f"{d['n_sites']} gardens x {d['n_founders']} founders")'''

code_partition = '''vp=d["variance_partition"]
rows=[]
for tag in ["bio1","multivar"]:
    v=vp[tag]
    rows.append([tag, v["h2_kinship_nomatch"], v["h2_kinship"], v["gamma"], v["gamma_z"], v["gamma_p"]])
part=pd.DataFrame(rows,columns=["climate","h2_kinship(no match)","h2_kinship(+match)",
                                "gamma(match)","gamma z","gamma p (in-sample)"])
print("IN-SAMPLE variance partition (site fixed effects + kinship random + climate-match):")
display(part.round(3))
print("\\nKinship explains ~half the founder-performance variance; the home-field match term is")
print("strongly significant IN-SAMPLE (origin-matched founders over-perform within a garden).")'''

code_cv = '''cv=d["cross_validation"]; pm=d["permutation"]
print("LEAVE-ONE-SITE-OUT CV  (predict a NEVER-SEEN garden):")
for tag in ["bio1","multivar"]:
    c=cv[tag]; p=pm[tag]
    print(f"  {tag:9s}: kinship r={c['r_kinship']:.3f}  ->  +climate-match r={c['r_kinship_match']:.3f}"
          f"   (Δr={c['delta']:+.3f})")
    print(f"             across-fold paired t p={c['paired_p']:.3f} (n={c['n_sites']} folds, noisy)  |  "
          f"origin-permutation p={p['perm_p']:.3f}  [REAL vs shuffled origin]")
print("\\n=> Kinship alone predicts a held-out garden at r~0.69 (it IS mostly 'who you are').")
print("   Climate-match adds only a SMALL increment (Δr +0.006 bio1 / +0.013 multivar) -- but the")
print("   permutation test shows that small gain is STATISTICALLY REAL (real origin climate beats")
print("   shuffled labels). The across-fold paired t is underpowered at 17 folds, hence the two")
print("   p-values diverge; the permutation null is the appropriate, better-powered test here.")'''

code_fig = '''fig,ax=plt.subplots(1,3,figsize=(15,4.6))
# (1) generalist value vs origin climate
gv=tab["generalist_value"].values; o1=tab["origin_bio1"].values
ok=np.isfinite(gv)&np.isfinite(o1); rr=stats.pearsonr(o1[ok],gv[ok])[0]
ax[0].scatter(o1,gv,s=14,c="#34495e",alpha=.6,edgecolors="none")
ax[0].set_xlabel("founder origin bio1 (home °C)"); ax[0].set_ylabel("generalist value (mean selection)")
ax[0].set_title(f"who wins overall vs home climate  r={rr:.2f}",fontsize=10,loc="left")
# (2) CV bars
tags=["bio1","multivar"]; xb=np.arange(2)
ax[1].bar(xb-0.2,[cv[t]["r_kinship"] for t in tags],0.4,label="kinship only",color="#95a5a6")
ax[1].bar(xb+0.2,[cv[t]["r_kinship_match"] for t in tags],0.4,label="kinship + climate-match",color="#c0392b")
ax[1].set_xticks(xb); ax[1].set_xticklabels(tags); ax[1].axhline(0,color="k",lw=.5)
ax[1].set_ylabel("LOSO CV predictive r (held-out garden)")
ax[1].set_title("does climate add power for a NEW garden?",fontsize=10,loc="left")
ax[1].legend(fontsize=8,frameon=False)
for t,x in zip(tags,xb):
    ax[1].text(x,max(cv[t]["r_kinship"],cv[t]["r_kinship_match"])+.02,f"Δ={cv[t]['delta']:+.3f}\\np={cv[t]['paired_p']:.2f}",
               ha="center",fontsize=7)
# (3) permutation null (multivariate climate = stronger signal)
null=arr["null_multivar"]; obs=d["permutation"]["multivar"]["obs_delta"]; pp=d["permutation"]["multivar"]["perm_p"]
ax[2].hist(null,bins=30,color="#bdc3c7",edgecolor="none")
ax[2].axvline(obs,color="#c0392b",lw=2,label=f"observed Δ={obs:+.3f}\\nperm p={pp:.3f}")
ax[2].set_xlabel("climate gain Δr (kinship+match − kinship)"); ax[2].set_ylabel("permutations")
ax[2].set_title("origin-label permutation null (multivar)",fontsize=10,loc="left"); ax[2].legend(fontsize=8,frameon=False)
for a in ax: a.spines[["top","right"]].set_visible(False)
fig.suptitle(f"Predicting ecotype performance: 'who you are' (kinship) vs climate-match  "
             f"({d['n_sites']} gardens, {d['n_founders']} founders)",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.96]); plt.show()'''

md_interp = """## Interpretation

**It is mostly "who you are," with a small but real climate component on top.**

**Kinship dominates.** Genome-wide relatedness (generalist genetic merit) predicts a founder's
selection response at a **never-seen** garden with $r\\approx0.69$ — the dominant,
out-of-sample-validated signal. A handful of broadly-fit ecotypes rise in almost every climate.

**Climate adaptation is statistically real but small.** The home-field match term is strongly
significant *in-sample* (origin-matched founders out-perform their kinship-predicted value within a
garden). Out-of-sample it adds only a small increment ($\\Delta r\\approx+0.006$ for bio1,
$+0.013$ for multivariate climate), but the **origin-label permutation test shows that gain is not
chance** — real origin climate re-ranks a held-out garden better than shuffled origin labels
($\\text{perm }p\\approx0.005$). The across-fold paired $t$ is non-significant only because 17
gardens give too few folds to see so small an effect through fold-to-fold noise; the permutation
null, which holds the CV structure fixed, is the appropriate and better-powered test.

**Why the effect is small.** Founder origin climate is confounded with kinship/clade
(climatically-similar ecotypes are genetically similar), so most of what climate "knows" about who
wins is already encoded in kinship; only a thin, independent slice is separable — and with ~17
gardens there isn't much independent climate contrast to resolve it. This is the same
clade-vs-climate collinearity wall seen in the haploblock GEA, here quantified predictively:
**performance is predictable, that prediction is mostly generalist genetic merit, and climate
matching contributes a small but genuine extra.**"""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_markdown_cell(md_model),
    new_markdown_cell("## Load the precomputed solution"), new_code_cell(code_setup),
    new_markdown_cell("## In-sample variance partition"), new_code_cell(code_partition),
    new_markdown_cell("## Leave-one-site-out cross-validation (the decisive test)"), new_code_cell(code_cv),
    new_markdown_cell("## Figure"), new_code_cell(code_fig),
    new_markdown_cell(md_interp),
])
os.makedirs(os.path.dirname(OUT), exist_ok=True)
ExecutePreprocessor(timeout=300, kernel_name="python3", startup_timeout=180).preprocess(nb, {"metadata": {"path": ROOT}})
nbf.write(nb, OUT)
print(f"[done] executed -> {OUT}")
