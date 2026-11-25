#!/usr/bin/env python
"""Build+execute `sv_main_fig3.ipynb` (which SVs are purified). Shared setup/data code in `_main_fig_common.py`. `basic` env."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_fig_common import code_setup, code_grenenet, build, new_markdown_cell, new_code_cell

md_title = """# Selection on structural variants — Figure 3: which SVs are purified?

Condensed from `sv_selection_key_results.ipynb` (full per-context panels, tables and robustness there); Figures 1–3 are built separately (`_build_main_fig{1,2,3}_nb.py`), shared code in `_main_fig_common.py`.

**Unit and filters** (as in the key-results notebook): the SV **event** (arch3 alleles of a pangenome
bubble merged with `truvari collapse` defaults); called in ≥ 90% of the 231 founders throughout;
≥ 2 carrier founders for the GrENE-Net analyses and the composition panels, singletons included in the
spectra (Figure 1a–b, 3a). Context, first match wins: CDS > UTR > intron > reference TE > promoter
(≤ 2 kb upstream) > intergenic. Colours: insertions orange, deletions green, SNPs black, short indels grey."""

md_fig3 = """## Figure 3 — Which SVs are purified?

**(a)** Nature against experiment, one point per class of CDS insertion. *x*: Tajima's *D* in the
founders (95% bootstrap CI; vertical lines = synonymous and nonsynonymous SNP references). *y*: GrENE-Net load
frequency change of the class's events minus that of their matched CDS SNPs (Figure 2c statistic, % of start), mean ±
95% CI over the 31 gardens. Inserted-sequence classes from BLAST of the inserted sequence against
TAIR10 (`insertion_tair10_class.csv`): TE-derived, duplicated gene sequence, other (dispersed / local)
duplication, novel (< 50% of the sequence aligns). **(b)** The experimental side of (a) garden by
garden: frequency change of each class relative to matched CDS SNPs, one small dot per garden, large marker
= mean ± 95% CI. **(c)** What insertions are made of, by context: composition of insertion events by
inserted-sequence class."""

code_fig3 = r"""
I=E[E.kind=="ins"].copy(); C=I[I.context=="CDS"].copy(); iC=np.where((E.kind=="ins")&(E.context=="CDS"))[0]
snpCDS=(KD=="snp")&(CX==0)
SUB=[("TE-derived",["insCDS_te_derived"],(C.seq_class=="te_derived").to_numpy(),"#B2182B"),
     ("other duplication",["insCDS_dispersed_dup","insCDS_local_dup"],C.seq_class.isin(["dispersed_dup","local_dup"]).to_numpy(),"#EF8A62"),
     ("gene duplicate",["insCDS_gene_dup"],(C.seq_class=="gene_dup").to_numpy(),"#67A9CF"),
     ("novel",["insCDS_novel"],(C.seq_class=="novel").to_numpy(),"#2166AC"),
     ("all CDS insertions",["ins_CDS"],np.ones(len(C),bool),"0.2")]
fig,ax=plt.subplots(1,3,figsize=(16,4.3),gridspec_kw=dict(width_ratios=[1,0.8,0.9]))
a=ax[0]; style(a); panel(a,"a"); a.axhline(0,color="0.35",lw=1,ls="--"); rows=[]; DD={}; PN={}
for name,lab,ls in (("snp_syn","synonymous SNPs","--"),("snp_nonsyn","nonsynonymous SNPs",":")):
    d0=tajD(name)[0]; a.axvline(d0,color="0.35",lw=1,ls=ls); a.annotate(lab,xy=(d0,0.075),fontsize=7,ha="center",va="bottom",color="0.35")
for lab,names,m,col in SUB:
    d,lo,hi=tajD_ci(names)
    idx=iC[m]; ms=np.zeros(KD.size,bool); ms[CATZ[np.where(ok)[0][idx]]]=True; dd=load_diff(ms,snpCDS); mu,se=dd.mean(),dd.std(ddof=1)/np.sqrt(dd.size); DD[lab]=(dd,col)
    a.errorbar(d,mu,xerr=[[d-lo],[hi-d]],yerr=1.96*se,fmt="o",ms=8 if lab!="all CDS insertions" else 6,color=col,capsize=2,mfc=col if lab!="all CDS insertions" else "white",mew=1.5)
    OFF={"TE-derived":((-9,-18),"right"),"other duplication":((9,-16),"left"),"gene duplicate":((9,9),"left"),"novel":((9,9),"left"),"all CDS insertions":((-9,8),"right")}[lab]
    a.annotate(f"{lab}\n(n = {int(m.sum())})",xy=(d,mu),xytext=OFF[0],textcoords="offset points",ha=OFF[1],fontsize=7.2,color=col if col!="0.2" else "0.3")
    nl=matched_null(ms,snpCDS); PN[lab]=null_p(mu,nl)
    rows.append(dict(cls=lab,n=int(m.sum()),TajD=d,lo=lo,hi=hi,load=mu,load_se=se,gardens_below=int((dd<0).sum()),p_matched_null=PN[lab],null_sd=nl.std()))
a.set_xlabel("Tajima's D in the founders  (← rarer in nature)"); a.set_ylabel("frequency change minus matched CDS SNPs, % of start\n(← lost during the experiment)")
a.set_ylim(-0.16,0.09); a.set_xlim(-3.1,-0.4); a.yaxis.set_major_formatter(PercentFormatter(1.0,decimals=0))
# (b) per garden by sub-class
a=ax[1]; style(a,axis="x"); panel(a,"b"); a.axvline(0,color="0.35",lw=1,ls="--")
for j,(lab,(dd,col)) in enumerate(DD.items()):
    a.scatter(dd,j+rng.uniform(-0.18,0.18,dd.size),color=col,s=9,alpha=0.45,linewidth=0,zorder=2)
    a.errorbar(dd.mean(),j,xerr=1.96*dd.std(ddof=1)/np.sqrt(dd.size),fmt="o",ms=8,color=col,mfc=col if lab!="all CDS insertions" else "white",mew=1.5,capsize=3,zorder=3)
    a.annotate(f"{(dd<0).sum()}/31   {pfmt(PN[lab])}",xy=(0.99,j),xycoords=("axes fraction","data"),ha="right",va="center",fontsize=7,color="0.35")
a.set_yticks(range(len(DD))); a.set_yticklabels(list(DD)); a.set_ylim(len(DD)-0.5,-0.5); a.set_xlabel("frequency change minus matched CDS SNPs, per garden (% of start)"); a.xaxis.set_major_formatter(PercentFormatter(1.0,decimals=0))
a.annotate("gardens < 0   matched-SNP null",xy=(0.99,-0.45),xycoords=("axes fraction","data"),ha="right",va="center",fontsize=7,color="0.35")
# (c) composition
a=ax[2]; style(a,axis="y"); panel(a,"c")
SC=[("te_derived","TE-derived","#B2182B"),("dispersed_dup","dispersed duplication","#EF8A62"),("local_dup","local duplication","#FDDBC7"),("gene_dup","gene duplicate","#67A9CF"),("novel","novel","#2166AC")]
comp=pd.crosstab(I.context,I.seq_class,normalize="index").reindex(CTX); bottom=np.zeros(6)
for key,lab,col in SC:
    v=comp[key].to_numpy(); a.bar(range(6),v,bottom=bottom,color=col,width=0.7,label=lab,edgecolor="white",lw=0.5); bottom+=v
for i,c in enumerate(CTX): a.annotate(f"n = {int((I.context==c).sum()):,}",xy=(i,1.01),ha="center",va="bottom",fontsize=7,color="0.35")
a.set_xticks(range(6)); a.set_xticklabels(CTX); a.set_ylim(0,1.08); a.set_ylabel("share of insertion events"); a.legend(frameon=False,fontsize=7,loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=3)
fig.tight_layout(); fig.savefig(f"{G}/plots/main_fig3.png",dpi=200); fig.savefig(f"{G}/plots/main_fig3.pdf"); plt.show()
print("Figure 3a:"); print(pd.DataFrame(rows).round(3).to_string(index=False))
print("\nFigure 3c  composition by context:"); print(comp.rename(columns=dict((k,l) for k,l,_ in SC)).round(3).to_string())
"""

md_methods = """## Models and statistics, as run

- **Spectra (Fig. 1a, d; 3a).** Per variant, carriers *c* among *k* called founders (k ≥ 208 of 231,
  1 ≤ c ≤ k − 1; singletons included); folded MAF = min(c/k, 1 − c/k). Two-sample Kolmogorov–Smirnov on
  MAF against synonymous SNPs (effect size *D* reported; p-values are driven by *n*).
  Tajima's *D* = (π − S/a₁) / √(e₁S + e₂S(S−1)) with the Tajima (1989) constants for n = median called
  founders of the class, π = Σ<sub>v</sub> 2c<sub>v</sub>(k<sub>v</sub> − c<sub>v</sub>)/(k<sub>v</sub>(k<sub>v</sub> − 1)),
  S = number of variants; CI from bootstrap resampling of variants (300 draws; 40 for the SNP classes,
  whose CI is negligible). No down-projection: variants are called in a median of 216–231 of 231
  founders, and projecting to a common sample of 200 changes *D* by < 0.03 in every class. Degeneracy: standard code, reference codon, `.1` gene models; sites with conflicting
  models dropped.
- **Allele frequencies in GrENE-Net (Fig. 2c–d, 3a–b).** kMate's stored per-pool allele frequencies
  (`common/results/pool_matrices/pool_gen9_*`: every plot at its last sampled generation, 352 pools) and
  the seed-mix founding frequency (mean of 8 runs, `af_store/p0_*`). Garden frequency = mean over the
  garden's plots. SV event = sum of its alleles' frequencies, capped at 1. SNPs a random 10% of records,
  short indels a random 10%; filters and labels as in the founder panel (`_sv_garden_af.py`).
- **Load (Fig. 2c, 3a–b).** *L* = Σ<sub>v</sub> *p*<sub>v</sub> over all events of the class (any
  climate of origin); matched SNP set = CDS SNPs re-weighted by (SV count / SNP count) per carrier-count ×
  carrier-temperature cell. **Test:** the 31 gardens share one founder pool, so garden counts ("22/31")
  are descriptive only; the p-value is empirical — 1,000 random sets of CDS SNPs, one SNP per event
  drawn within the event's cell, same statistic, computed in this notebook. The supplementary cell
  splits the classes by climate of origin against same-origin SNPs.
- **Robustness (`_sv_main_figures_checks.py`).** Burden model with founder origin bio1 (and bio1²) as
  extra covariates: β₁ unchanged (−0.22); the burden–fitness relation holds within origin tertiles
  (cold ρ = −0.51, warm ρ = −0.27). Spectra without singletons (minor-allele count ≥ 2): same ordering
  (values printed under Figure 1)."""

build("sv_main_fig3.ipynb", [new_markdown_cell(md_title), new_code_cell(code_setup), new_code_cell(code_grenenet), new_markdown_cell(md_fig3), new_code_cell(code_fig3), new_markdown_cell(md_methods)])
