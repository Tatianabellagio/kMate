#!/usr/bin/env python
"""Build the RDA variance-partition notebook: explicit methods + validated results."""
import os, nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

C = []
C.append(new_markdown_cell(
    "# Do SVs / indels add climate-adaptive information beyond SNPs?\n"
    "## Partial-RDA variance partitioning by variant class (SNP / indel / SV / ALL)\n\n"
    "Instead of hunting individual loci (null under isolation-by-environment), we ask "
    "**how much of each variant class's among-site allele-frequency variance is explained "
    "by climate, net of population structure**, and whether adding indels/SVs to SNPs "
    "explains any climate variance the SNPs can't."))

C.append(new_markdown_cell(
    "## Exactly what we ran\n\n"
    "**Method:** partial Redundancy Analysis (pRDA) variance partitioning "
    "([Capblancq & Forester 2021, *Meth. Ecol. Evol.*](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722)), "
    "implemented in numpy (`rdavp_varpart.py`; `vegan` not used — hangs on Savio NFS).\n\n"
    "| component | choice |\n"
    "|---|---|\n"
    "| **Response Y** | site allele frequencies, `[31 sites × L]` per class (plot AF aggregated to site by **flower-weighted mean**), each variant **column-centered** — *no* Hellinger, *no* scaling, *no* depth-weighting (kMate AF is calibrated / equal-precision) |\n"
    "| **Predictor X (climate)** | PCA of the 19 standardized bioclim vars across the 31 sites; **3 PCs** kept by **broken-stick** (cum. var 0.82) |\n"
    "| **Covariate Z (structure)** | top **3 SNP-PCA axes** (PCA of the SNP site-AF matrix) — the *same* covariate conditions **every** class, so the confound cancels in the class contrast |\n"
    "| **Classes** | SNP · smallindel · SV · ALL (SNP+indel+SV stacked) |\n"
    "| **MAF filter** | site-level MAF ≥ **0.05** (primary) + ≥ **0.01** (sensitivity), identical for all classes; invariant sites dropped |\n"
    "| **Statistics** | **pure-climate R²** = constrained variance of Y (residualized on Z) explained by X (residualized on Z); **adjusted R²** (Ezekiel, n=31, conditioned df); **confounded** = full-climate R² − pure-climate R²; **permutation p** = permute the 31 climate rows ×2000, recompute constrained variance |\n"
    "| **n** | **31 sites** is the true sample size (not the 355 plots — pseudoreplication) |\n\n"
    "**Constrained variance** uses the identity `SS_constrained = trace(Hx · G)` with "
    "`G = Yz·Yzᵀ` (31×31) precomputed and `Hx` the climate hat-matrix — so 2000 permutations "
    "are near-instant.\n\n"
    "**Code audit (`rdavp_audit.py`) — passed:** Gram-trick identity exact; **positive "
    "control** (response built from climate) → R²=0.99, p=0.0005 (detects real signal); "
    "**negative control** (random response) → R²=0.100 = null expectation p/(n−1), uniform "
    "permutation p; site-AF aggregation matches manual; confounded ≥ 0."))

C.append(new_markdown_cell("## Results"))
C.append(new_code_cell(
    "import pandas as pd, numpy as np, matplotlib.pyplot as plt, matplotlib as mpl\n"
    "mpl.rcParams.update({'figure.dpi':110,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})\n"
    "OUT='/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/results/rda_varpart'\n"
    "df=pd.read_csv(f'{OUT}/varpart_by_class.csv')\n"
    "df"))

C.append(new_markdown_cell(
    "### Variance-partition figure\n"
    "Each bar = climate-explained among-site variance, split into the **confounded** "
    "(climate∩structure) fraction and the **pure-climate** fraction (net of structure). "
    "Adjusted R² and permutation p annotated. MAF≥0.05."))
C.append(new_code_cell(
    "d=df[df.maf==0.05].set_index('cls').loc[['snp','smallindel','sv','all']]\n"
    "lab={'snp':'SNP','smallindel':'indel','sv':'SV','all':'ALL'}\n"
    "kc=int(d.kclim.iloc[0]); kz=int(d.kstruct.iloc[0]); n=31\n"
    "null_exp=kc/(n-kz-1)   # E[pure-climate R²] under NO signal (overfitting expectation)\n"
    "x=np.arange(len(d)); fig,ax=plt.subplots(figsize=(8,4.8))\n"
    "conf=d.confounded.values; pure=d.pure_climate_R2.values\n"
    "ax.bar(x,conf,color='#c9ccd1',label='confounded (climate ∩ structure)')\n"
    "ax.bar(x,pure,bottom=conf,color='#3b7dd8',alpha=0.9,label='pure-climate raw R² (= null expectation, n.s.)')\n"
    "# dashed line = the pure-climate R² EXPECTED under no signal; bars sit on it -> not real\n"
    "for i in range(len(d)):\n"
    "    ax.plot([i-0.42,i+0.42],[conf[i]+null_exp]*2,ls='--',lw=1.3,color='#b8860b',\n"
    "            label=('pure-climate null expectation (no signal)' if i==0 else None))\n"
    "for i,(c,r) in enumerate(d.iterrows()):\n"
    "    ax.text(i, conf[i]+pure[i]+0.005, f\"adjR²={r.pure_climate_R2adj:+.3f}\\np={r.perm_p:.2f}\",\n"
    "            ha='center',va='bottom',fontsize=8)\n"
    "ax.set_xticks(x); ax.set_xticklabels([f'{lab[c]}\\n(L={int(d.loc[c,\"L\"]):,})' for c in d.index])\n"
    "ax.set_ylabel('climate-explained among-site variance (R²)')\n"
    "ax.set_ylim(0, d.full_climate_R2.max()*1.4)\n"
    "ax.set_title('pRDA climate variance by variant class — MAF≥0.05, 31 sites, 3 climate PCs | 3 structure PCs\\n'\n"
    "             'pure-climate bars sit ON the null-expectation line → adjR²≈0, not significant',fontsize=10)\n"
    "ax.legend(loc='upper right',fontsize=7.5,frameon=False)\n"
    "fig.tight_layout(); fig.savefig(f'{OUT}/varpart_by_class.png',dpi=150,bbox_inches='tight')\n"
    "print(f'null expectation for pure-climate R² (no signal) = kc/(n-kz-1) = {null_exp:.3f}')\n"
    "print('wrote',f'{OUT}/varpart_by_class.png'); plt.show()"))

C.append(new_code_cell(
    "# headline: does adding indels+SVs to SNPs add pure-climate information?\n"
    "for thr in [0.05,0.01]:\n"
    "    t=df[df.maf==thr].set_index('cls')\n"
    "    add=t.loc['all','pure_climate_R2adj']-t.loc['snp','pure_climate_R2adj']\n"
    "    print(f'MAF>={thr}: pure-climate adjR²  SNP={t.loc[\"snp\",\"pure_climate_R2adj\"]:+.4f}  '\n"
    "          f'indel={t.loc[\"smallindel\",\"pure_climate_R2adj\"]:+.4f}  '\n"
    "          f'SV={t.loc[\"sv\",\"pure_climate_R2adj\"]:+.4f}  ALL={t.loc[\"all\",\"pure_climate_R2adj\"]:+.4f}  '\n"
    "          f'| ALL−SNP={add:+.4f}')"))

C.append(new_markdown_cell(
    "## Read-out\n"
    "1. **All classes are interchangeable** — SNP, indel, SV and ALL explain the *same* "
    "climate variance (full ≈0.176, confounded ≈0.064, pure adjusted ≈0) at both MAF cuts.\n"
    "2. **Adding indels+SVs to SNPs adds nothing** — `ALL − SNP ≈ 0` (pure-climate adjusted R²). "
    "In the variance-explained framing, **SVs/indels carry no climate information beyond SNPs** "
    "(fully redundant / SNP-tagged — consistent with the SV-heritability literature and prior "
    "'SVs are passengers' results).\n"
    "3. **Net of structure there is no significant climate signal for any class** — pure-climate "
    "adjusted R² ≈ 0, permutation **p ≈ 0.44–0.54**. Of the ~17.6% of among-site variance climate "
    "'explains', ~6.4% is confounded with structure (isolation-by-environment) and the rest is "
    "indistinguishable from chance at n=31.\n\n"
    "**Answer to the project question:** in a structure-robust, polygenic-native variance framing, "
    "structural variants and indels **do not add climate-adaptive information beyond SNPs** — a "
    "clean, quantified negative (SVs add −0.08% adjusted climate variance, p≈0.5)."))

nb=new_notebook(); nb['cells']=C
nb.metadata['kernelspec']={'name':'python3','display_name':'Python 3','language':'python'}
os.makedirs('notebooks',exist_ok=True)
with open('notebooks/rda_varpart_by_class.ipynb','w') as f: nbformat.write(nb,f)
print('wrote notebooks/rda_varpart_by_class.ipynb')
