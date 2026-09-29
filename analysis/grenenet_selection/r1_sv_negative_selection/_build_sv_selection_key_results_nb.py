#!/usr/bin/env python
"""Build+execute `sv_selection_key_results.ipynb`: the key results on selection on SVs, on Truvari events.
Load-only: sv_event_catalog.csv(+_p.npy), sv_sfs_context_origin_event.npz, context_genome_share.csv, grenenet_sfs_event.npz,
sv_burden_lmm{,_joint,_by_garden}_event.csv, sv_burden_lmm_perm_event.npz, cds_shift_null_event.npz.
`basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "notebooks", "sv_selection_key_results.ipynb")

md_title = """# Selection on structural variants — key results

Three questions, in order:

1. **Evolutionary history** — in the 231 founders, do SVs carry the signature of purifying selection?
2. **GrENE-Net** — is that selection still acting during the experiment?
3. **Which SVs** — what are the variants being purified: where do they land, how big are they, what do
   they contain?

**Unit: the SV event.** arch3 stores one record per assembly path, so one insertion appears as several
near-identical "alleles". Alleles of a pangenome bubble are merged with `truvari collapse` defaults
(size and sequence similarity ≥ 0.95, ≤ 500 bp apart; `_build_event_map.py`): 226k SV records → 89.9k
events (50.2k insertions, 39.7k deletions). A founder carries an event if it carries any of its
alleles. **Filters throughout:** ≥ 2 carrier founders (private variants excluded), called in ≥ 90% of
the 231 founders (the spectra of section 1 include singletons). SNPs and short indels (1–50 bp) per
record, same filters.
Context, first match wins: CDS > UTR > intron > reference TE > promoter (≤ 2 kb upstream) > intergenic.

Full analyses: `sv_purging_fixed.ipynb` (per-variant trajectories, climate of origin), `sv_burden_lmm.ipynb`
(founder burden, population load, robustness); allele-level twins `*_allele.ipynb`."""

code_setup = r"""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'figure.dpi':110,'font.size':8,'axes.linewidth':0.6})
G=os.path.abspath(os.path.join(os.getcwd(),"..","r1_sv_negative_selection","results","sv_adaptive"))
os.makedirs(f"{G}/plots",exist_ok=True)
CTX=["CDS","UTR","intron","TE","promoter","intergenic"]
KC={"snp":"0.2","indel":"0.6","ins":"#D95F02","del":"#1B9E77"}; KL={"snp":"SNPs","indel":"short indels","ins":"insertions","del":"deletions"}
def style(ax,axis="both"):
    ax.set_axisbelow(True); ax.grid(color="0.88",lw=0.7,axis=axis)
    for sp in ax.spines.values(): sp.set_visible(False)
rng=np.random.default_rng(0)
D=pd.read_csv(f"{G}/sv_event_catalog.csv"); PALL=np.load(f"{G}/sv_event_catalog_p.npy")
ok=D.passes.to_numpy(); E=D[ok].reset_index(drop=True); PE=PALL[ok]
S=np.load(f"{G}/sv_sfs_context_origin_event.npz")
# one frame per class: carriers, context, origin counts (>= 2 carriers)
def cls_frame(k):
    if k in ("ins","del"):
        d=E[E.kind==k]; return pd.DataFrame(dict(n=d.n_carriers.to_numpy(),ctx=d.context.map({c:i for i,c in enumerate(CTX)}).to_numpy(),
                                                  nc=d.n_cold.to_numpy(),nw=d.n_warm.to_numpy()))
    m=S[f"{k}_n"]>=2; return pd.DataFrame(dict(n=S[f"{k}_n"][m],ctx=S[f"{k}_ctx"][m],nc=S[f"{k}_nc"][m],nw=S[f"{k}_nw"][m]))
FR={k:cls_frame(k) for k in ("snp","indel","ins","del")}
print({KL[k]:f"{len(v):,}" for k,v in FR.items()}, "(SNPs and indels: random 5% of records)")
"""

md_1 = """## 1. Evolutionary history — the founder panel

The standard test: the **site-frequency spectrum** of a class of variants against a putatively neutral
reference class. Purifying selection keeps deleterious variants rare, so a selected class shows an
excess of low-frequency variants and a more negative Tajima's *D*.

- **Neutral reference:** SNPs at **4-fold degenerate (synonymous) sites**. **Known-selected
  comparator:** SNPs at **0-fold degenerate (nonsynonymous) sites**. Degeneracy from the TAIR10
  sequence and the representative gene models (`_sv_classic_sfs.py`).
- **Spectra:** folded (minor-allele frequency; insertions and deletions are defined against the
  reference, not polarised), 231 inbred founders as a haploid sample, variants called in ≥ 90% of
  founders, **singletons included**. SVs per event.
- **Tests vs 4-fold SNPs:** two-sample Kolmogorov–Smirnov and one-sided Mann–Whitney *U* on
  minor-allele frequency. **Tajima's *D*** per class from the spectrum projected to a common sample of
  200 (hypergeometric down-projection for missing calls), bootstrap 95% CI. Because SNP spectra
  themselves differ by genomic context, the same two tests are repeated against **SNPs of the same
  context** (`_ctx` columns; SNPs per context also plotted in 1b).
- **1c–d, supporting:** where each class is found, as its share in a context divided by that
  context's share of genomic length (`_context_genome_share.py`, same first-match rule; 1 = uniform
  placement), SNPs alongside; and the share of insertions and deletions whose length is a multiple of
  3 (1/3 expected outside coding sequence)."""

code_1 = r"""
CS=np.load(f"{G}/classic_sfs.npz"); M=int(CS["M"])
def maf(name): c=CS[f"{name}_c"].astype(float); k=CS[f"{name}_k"].astype(float); f=c/k; return np.minimum(f,1-f)
def tajD(sfs):
    n=M; j=np.arange(1,n); S_=sfs.sum(); pi=(sfs*2*j*(n-j)/(n*(n-1))).sum()
    a1=(1/j).sum(); a2=(1/j**2).sum(); b1=(n+1)/(3*(n-1)); b2=2*(n*n+n+3)/(9*n*(n-1))
    c1=b1-1/a1; c2=b2-(n+2)/(a1*n)+a2/a1**2; e1=c1/a1; e2=c2/(a1**2+a2)
    return (pi-S_/a1)/np.sqrt(e1*S_+e2*S_*(S_-1))
def tajD_ci(sfs,B=300):
    S_=int(round(sfs.sum())); p=sfs/sfs.sum(); bs=[tajD(rng.multinomial(S_,p).astype(float)) for _ in range(B)]
    return tajD(sfs),np.percentile(bs,2.5),np.percentile(bs,97.5)
REF=maf("snp_4fold"); bins=np.linspace(0,0.5,11); rows=[]
def row(name,label,ctx=None):
    v=maf(name); d,lo,hi=tajD_ci(CS[f"{name}_sfs"]); ks=stats.ks_2samp(v,REF); mw=stats.mannwhitneyu(v,REF,alternative="less")
    o=dict(cls=label,n=v.size,singletons=(CS[f"{name}_c"]==1).mean(),maf_lt5=(v<0.05).mean(),mean_maf=v.mean(),
           KS_D=ks.statistic,KS_p=ks.pvalue,MWU_p=mw.pvalue,TajD=d,lo=lo,hi=hi,KS_D_ctx=np.nan,MWU_p_ctx=np.nan)
    if ctx:                                  # same tests against SNPs of the same genomic context
        r=maf(f"snp_{ctx}"); o["KS_D_ctx"]=stats.ks_2samp(v,r).statistic; o["MWU_p_ctx"]=stats.mannwhitneyu(v,r,alternative="less").pvalue
    return o
rows.append(row("snp_4fold","SNPs 4-fold (neutral ref.)")); rows.append(row("snp_0fold","SNPs 0-fold (nonsyn.)"))
for k in ("snp","ins","del","indel"):
    for x in CTX: rows.append(row(f"{k}_{x}",f"{KL[k]} {x}",None if k=="snp" else x))
T1=pd.DataFrame(rows)
# 1a folded SFS per context
fig,ax=plt.subplots(2,3,figsize=(16,7),sharex=True,sharey=True); w=0.2; xb=np.arange(10)
for a,x in zip(ax.ravel(),CTX):
    style(a,axis="y")
    for j,(name,lab,col) in enumerate((("snp_4fold","SNPs 4-fold (neutral)","0.75"),("snp_0fold","SNPs 0-fold (nonsyn.)","0.35"),
                                      (f"del_{x}","deletions","#1B9E77"),(f"ins_{x}","insertions","#D95F02"))):
        h=np.histogram(maf(name),bins)[0]; a.bar(xb+(j-1.5)*w,h/h.sum(),width=w*0.95,color=col,label=lab)
    ki=stats.ks_2samp(maf(f"ins_{x}"),REF); kd=stats.ks_2samp(maf(f"del_{x}"),REF)
    a.annotate(f"{x}\ninsertions n = {maf(f'ins_{x}').size:,}  KS D = {ki.statistic:.2f}, p = {ki.pvalue:.1g}\n"
               f"deletions n = {maf(f'del_{x}').size:,}  KS D = {kd.statistic:.2f}, p = {kd.pvalue:.1g}",
               xy=(0.97,0.95),xycoords="axes fraction",ha="right",va="top",fontsize=7.8,color="0.2")
for a in ax[1]: a.set_xticks(xb); a.set_xticklabels([f"{b:.2f}" for b in bins[1:]],fontsize=7); a.set_xlabel("minor-allele frequency (bin upper edge)")
for a in ax[:,0]: a.set_ylabel("proportion of variants")
h_,l_=ax[0,0].get_legend_handles_labels(); fig.legend(h_,l_,loc="lower center",bbox_to_anchor=(0.5,0.99),ncol=4,frameon=False,fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_key_1a_sfs.png",dpi=130,bbox_inches="tight"); plt.show()
# 1b Tajima's D, 1c density, 1d frame
fig,ax=plt.subplots(1,3,figsize=(16,4),gridspec_kw=dict(width_ratios=[1.5,1,0.8]))
a=ax[0]; style(a,axis="y")
for lab,ls in (("SNPs 4-fold (neutral ref.)","--"),("SNPs 0-fold (nonsyn.)",":")):
    a.axhline(T1.set_index("cls").loc[lab,"TajD"],color="0.35",lw=1,ls=ls); a.annotate(lab,xy=(-0.55,T1.set_index("cls").loc[lab,"TajD"]),fontsize=7,va="bottom",ha="left",color="0.3")
for j,k in enumerate(("snp","indel","del","ins")):
    d=T1[T1.cls.isin([f"{KL[k]} {x}" for x in CTX])]; a.errorbar(np.arange(6)+(j-1.5)*0.2,d.TajD,yerr=[d.TajD-d.lo,d.hi-d.TajD],fmt="o",ms=5,color=KC[k],capsize=2,label=KL[k]+(" (same context)" if k=="snp" else ""))
a.set_xticks(range(6)); a.set_xticklabels(CTX); a.set_ylabel("Tajima's D (sample of 200)"); a.legend(frameon=False,fontsize=7.5,loc="lower right")
a.annotate("1b  Tajima's D by class and context",xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5)
a=ax[1]; style(a,axis="y"); a.axhline(1,color="0.35",lw=1,ls="--")
GSH=pd.read_csv(f"{G}/context_genome_share.csv").set_index("context").share.reindex(CTX).to_numpy(); dens={}
for j,k in enumerate(("snp","indel","del","ins")):
    n=np.bincount(FR[k].ctx,minlength=6); sh=n/n.sum(); r=sh/GSH; dens[k]=r
    a.errorbar(np.arange(6)+(j-1.5)*0.2,r,yerr=1.96*np.sqrt(sh*(1-sh)/n.sum())/GSH,fmt="o",ms=5,color=KC[k],capsize=2,label=KL[k])
a.set_xticks(range(6)); a.set_xticklabels(CTX,fontsize=7.5); a.set_ylabel("share in context / share of genomic length"); a.legend(frameon=False,fontsize=7.5,loc="upper left")
a.annotate("1c  where they are found (dashed = uniform placement)",xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5)
a=ax[2]; style(a,axis="y"); a.axhline(1/3,color="0.35",lw=1,ls="--"); fr=[]
for j,k in enumerate(("del","ins")):
    Ik=E[E.kind==k]
    for i,c in enumerate(CTX):
        d=Ik[Ik.context==c]; x_=int((d.size_bp%3==0).sum()); n=len(d); ci=stats.binomtest(x_,n,1/3).proportion_ci()
        a.errorbar([i+(j-0.5)*0.24],[x_/n],yerr=[[x_/n-ci.low],[ci.high-x_/n]],fmt="o",ms=6,color=KC[k],capsize=2,label=KL[k] if i==0 else None)
        fr.append(dict(kind=KL[k],context=c,n=n,in_frame=x_/n,p_vs_third=stats.binomtest(x_,n,1/3).pvalue))
a.set_xticks(range(6)); a.set_xticklabels(CTX,fontsize=7.5); a.set_ylabel("share with length multiple of 3"); a.legend(frameon=False,fontsize=7.5,loc="upper right")
a.annotate("1d  reading frame (dashed = 1/3)",xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_key_1_founders.png",dpi=130,bbox_inches="tight"); plt.show()
pd.set_option("display.width",220)
print("1a-b  folded SFS vs 4-fold SNPs: KS and one-sided Mann-Whitney (class MAF lower), Tajima's D [95% CI]; _ctx = same tests vs SNPs of the same context")
o=T1.copy(); o["Tajima D"]=o.TajD.map("{:+.2f}".format)+" ["+o.lo.map("{:+.2f}".format)+", "+o.hi.map("{:+.2f}".format)+"]"
print(o.drop(columns=["TajD","lo","hi"]).to_string(index=False,formatters={"singletons":"{:.3f}".format,"maf_lt5":"{:.3f}".format,"mean_maf":"{:.4f}".format,
      "KS_D":"{:.3f}".format,"KS_p":"{:.2g}".format,"MWU_p":"{:.2g}".format,"KS_D_ctx":"{:.3f}".format,"MWU_p_ctx":"{:.2g}".format,"n":"{:,}".format}))
print("\n1c  share in context relative to the context's share of genomic length (1 = uniform)")
print(pd.DataFrame(dens,index=CTX).rename(columns=KL).T.round(2).to_string()); print("    genomic share:", dict(zip(CTX,GSH.round(3))))
print("    insertions relative to SNPs:", dict(zip(CTX,(dens["ins"]/dens["snp"]).round(2))))
print("\n1d  in-frame share by context"); print(pd.DataFrame(fr).round(4).to_string(index=False))
Ic=E[(E.kind=="ins")&(E.context=="CDS")]
print("\n    CDS insertions, in-frame share by inserted sequence:")
print(Ic.groupby("seq_class").apply(lambda d: pd.Series(dict(n=len(d),in_frame=(d.size_bp%3==0).mean(),p_vs_third=stats.binomtest(int((d.size_bp%3==0).sum()),len(d),1/3).pvalue))).sort_values("in_frame",ascending=False).round(4).to_string())
"""

md_1_read = """**Reading it.** The classic footprint of purifying selection, in the classic order.

- **The reference behaves.** 4-fold SNPs have the least skewed spectrum of any class (Tajima's
  *D* = −0.58, 23% singletons). 0-fold SNPs are shifted towards rare (*D* = −1.40, 33% singletons,
  KS *D* = 0.14): the known-selected class separates from the neutral one as it should.
- **Insertions in coding sequence are the most skewed class in the genome:** 43% singletons, 84% below
  MAF 0.05, KS *D* = 0.28 against 4-fold SNPs (p = 3e-87), Tajima's *D* = −2.09 [−2.20, −1.99] —
  beyond nonsynonymous SNPs (−1.40) and beyond short indels in CDS (−1.82). UTR insertions follow
  (−1.87), then promoter (−1.80).
- **Read like for like.** SNPs are not uniform across the genome: SNPs in reference TEs (−1.88) and
  intergenic sequence (−1.51) are *more* skewed than nonsynonymous SNPs, so part of the skew of TE and
  intergenic insertions is the region, not the insertion. Against SNPs of the same context the excess
  of rare insertions is graded exactly by function — the gap in *D* is 1.10 in CDS (KS *D* = 0.21,
  p = 5e-49), 0.69 in UTRs, 0.47 in promoters, 0.34 in introns, 0.25 in intergenic sequence and 0.05
  in TEs (KS p = 0.07): insertions into TE sequence are indistinguishable from SNPs there.
- **Deletions** carry the signal only in CDS (*D* = −1.72 against −0.99 for CDS SNPs; KS *D* = 0.15,
  p = 6e-23). Everywhere else they sit at or *above* the SNP background (one-sided Mann–Whitney
  p ≈ 1), most clearly in reference TEs (−1.03 against −1.88): a "deletion" there is mostly the
  absence of a TE copy the reference happens to carry.
- **Missing from genes (1c).** Coding sequence is 28% of the genome. Insertions reach a fifth of
  that (0.20× uniform placement; SNPs 0.79×, short indels 0.13×, deletions 0.46×) and are also
  under-represented in introns (0.74×; SNPs 1.31×). In UTRs they sit at expectation (0.97×), in TE
  sequence at the SNP level (0.79× vs 0.76×), and they accumulate in promoters (2.29×) and
  intergenic sequence (1.56×). Relative to SNPs the CDS depletion is 4×. Variants called in < 90% of
  founders are excluded here as everywhere; their exclusion is conservative for this panel.
- **Frame (1d).** 44% of CDS insertions and 44% of CDS deletions have lengths that are multiples of
  3, against 33–36% in every other context (both p < 1e-8). Among insertions the excess is carried by
  duplicated sequence: gene duplicates 59%, local duplications 53%, novel 40%, dispersed duplications
  38%, **TE-derived 33%** — the TE insertions that survive in coding sequence are not frame-selected,
  while duplicated coding sequence is.
- **Caveats.** Every class has negative *D*: the 231 founders are a structured, selfing, range-wide
  sample and assembly-based calling adds singletons, so absolute values are not interpretable — only
  the ranking of classes and the within-context contrasts are. Alleles are defined against the
  reference, not polarised, hence folded spectra. So over evolutionary time SVs **are** under
  purifying selection, graded by how much function they disrupt: CDS > UTR > promoter > intron >
  intergenic > TE, and insertions > deletions."""

md_2 = """## 2. GrENE-Net — is it still happening?

Pool frequencies in GrENE-Net are founder frequencies projected onto founder genotypes, so the
question is asked at the level where it is identifiable:

- **2a. Founder burden.** Do founders carrying more events of a class lose frequency, beyond what
  relatedness predicts? Mixed model: founder mean log-fitness over the 31 gardens ~ burden in the
  context (rate: carried / called) + burden in intergenic sequence of the same kind + call rates +
  sequencing route + SNP-kinship random effect. Coefficient per SD of burden ± 95% CI; **shared**
  variants. SNPs and short indels alongside.
- **2b. Population load.** Expected number of a class's variants per genome, start vs last generation,
  per garden: log change of shared insertions **minus** that of SNPs matched on carrier count and
  carriers' origin climate. Below 0 = depleted beyond matched SNPs. One dot per garden."""

code_2 = r"""
R=pd.read_csv(f"{G}/sv_burden_lmm_event.csv"); J=pd.read_csv(f"{G}/sv_burden_lmm_joint_event.csv")
GB=pd.read_csv(f"{G}/sv_burden_lmm_by_garden_event.csv"); PM=np.load(f"{G}/sv_burden_lmm_perm_event.npz")
GS=np.load(f"{G}/grenenet_sfs_event.npz",allow_pickle=True); PV=GS["p"]; KD=GS["kind"]; CX=GS["ctx"]; ZN=GS["zone"]
CB=np.array([2,3,4,5,6,8,11,16,26,51,1000]); BQ=np.quantile(GS["b1"][KD=="snp"],[0.2,0.4,0.6,0.8])
CELL=np.digitize(GS["n"],CB)*10+np.digitize(GS["b1"],BQ); SH=(GS["nc"]>0)&(GS["nw"]>0)
def load_diff(p_sv,cell_sv,m_snp):
    cs=pd.Series(cell_sv).value_counts(); cn=pd.Series(CELL[m_snp]).value_counts(); kp=np.isin(cell_sv,cn.index)
    Ls=p_sv[kp].sum(0); w=(cs/cn).reindex(CELL[m_snp]).fillna(0).to_numpy(); Ln=(PV[m_snp]*w[:,None]).sum(0)
    return np.log(Ls[1:]/Ls[0])-np.log(Ln[1:]/Ln[0])
fig,ax=plt.subplots(1,2,figsize=(15,4.2),gridspec_kw=dict(width_ratios=[1,1.25]))
a=ax[0]; style(a,axis="x"); a.axvline(0,color="0.35",lw=1,ls="--"); C5=["CDS","UTR","intron","TE","promoter"]
for j,k in enumerate(("ins","del","indel","snp")):
    d=R[(R.freq=="all")&(R.kind==k)&(R.origin=="shared")].set_index("context").reindex(C5); y=np.arange(5)+(j-1.5)*0.18
    a.errorbar(d.beta,y,xerr=1.96*d.se,fmt="none",ecolor=KC[k],elinewidth=1.2)
    sig=(d.p<0.05).to_numpy(); a.scatter(d.beta[sig],y[sig],color=KC[k],s=28,zorder=3,label=KL[k]); a.scatter(d.beta[~sig],y[~sig],facecolor="white",edgecolor=KC[k],s=28,zorder=3)
a.set_yticks(range(5)); a.set_yticklabels(C5); a.set_ylim(4.6,-0.6); a.set_xlabel("Δ founder mean log-fitness per SD of burden (filled: p < 0.05)")
a.legend(frameon=False,fontsize=7.5,loc="lower left"); a.annotate("2a  founder burden, shared variants",xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5)
a=ax[1]; style(a,axis="y"); a.axhline(0,color="0.35",lw=1,ls="--"); rows=[]
for i,c in enumerate(CTX):
    for j,k in enumerate(("ins","del")):
        ms=SH&(KD==k)&(CX==i); d=load_diff(PV[ms],CELL[ms],SH&(KD=="snp")&(CX==i)); x=i+(j-0.5)*0.34
        zc=np.where(ZN=="cold","#2166AC",np.where(ZN=="hot","#B2182B","0.6"))
        a.scatter(x+rng.uniform(-0.08,0.08,d.size),d,c=zc,s=9,alpha=0.65 if k=="ins" else 0.3,linewidth=0,zorder=2)
        a.boxplot(d,positions=[x],widths=0.28,patch_artist=True,showfliers=False,manage_ticks=False,boxprops=dict(facecolor="none",edgecolor=KC[k],linewidth=1.5 if k=="ins" else 0.9),
                  medianprops=dict(color=KC[k],lw=2),whiskerprops=dict(color=KC[k]),capprops=dict(color=KC[k]),zorder=3)
        rows.append(dict(kind=KL[k],context=c,n_events=int(ms.sum()),mean_diff=d.mean(),gardens_below=int((d<0).sum()),p=stats.wilcoxon(d).pvalue))
a.set_xticks(range(6)); a.set_xticklabels(CTX); a.set_ylim(-0.11,0.09)
a.set_ylabel("log load change: SV class minus matched SNPs")
a.annotate("2b  population load, shared insertions (bold) and deletions; dots = gardens (blue cold, grey mid, red hot)",xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_key_2_grenenet.png",dpi=130,bbox_inches="tight"); plt.show()
h=R[(R.kind=="ins")&(R.origin=="shared")&(R.context=="CDS")].set_index("freq")
print(f"2a  shared CDS insertion burden: {h.loc['all','beta']:+.3f} per SD (p = {h.loc['all','p']:.2g}; permutation p "
      f"{'< 0.001' if (np.abs(PM['null'])>=abs(PM['obs'])).mean()==0 else (np.abs(PM['null'])>=abs(PM['obs'])).mean()}); "
      f"negative in {(GB.beta<0).sum()}/31 gardens; rare events (2-10 carriers) alone {h.loc['rare','beta']:+.3f} (p = {h.loc['rare','p']:.2g})")
j=J[J.term=="ins_shared_CDS_all"].set_index("model"); print("    with SNP CDS load in the model: %+.3f (p = %.2g); with SNP + indel CDS load: %+.3f (p = %.2g)"%(
      j.loc["ins + SNP","beta"],j.loc["ins + SNP","p"],j.loc["ins + SNP + indel","beta"],j.loc["ins + SNP + indel","p"]))
t=R[(R.freq=="all")&(R.origin=="shared")].copy(); t["c"]=t.beta.map("{:+.2f}".format)+" ("+t.p.map("{:.2g}".format)+")"
print(t.pivot(index="kind",columns="context",values="c")[C5].rename(index=KL).to_string())
print("\n2b  load change vs matched SNPs"); print(pd.DataFrame(rows).round(4).to_string(index=False))
SHN=np.load(f"{G}/cds_shift_null_event.npz")
print(f"\nrobustness: same coefficient with stricter merging (-0.22) and looser (-0.24); without merging -0.14."
      f"\nrandom-CDS control (200 shifted annotations): burden p = {float(SHN['p_burden_matched']):.3f} size-matched "
      f"({float(SHN['p_burden']):.2f} unmatched); load p = {float(SHN['p_load_matched']):.3f} size-matched ({float(SHN['p_load']):.3f} unmatched)")
"""

md_2_read = """**Reading it.** Yes — weakly, and in the same class that nature purifies.

- **Founders carrying more shared CDS insertions lose frequency**: −0.22 per SD of burden
  (p = 6e-4, permutation p < 0.001), negative in 29 of 31 gardens, present for rare events alone, and
  still −0.18 (p = 0.008) with the same founders' SNP and indel CDS load in the model. UTR (−0.16)
  and intron (−0.20) follow; TE-sequence load predicts nothing — for any variant class.
- **In the population, shared CDS insertion load falls 1.5% beyond matched SNPs**, in 27 of 31 gardens,
  cold and hot alike (p = 2e-4). No other context, and deletions not at all.
- **Limits.** The effect is small. Genic load of *every* class (SNPs, indels) predicts lower founder
  fitness, so part of 2a is a genome-wide property of founders, only partly specific to SVs. Against
  randomly shifted annotations the CDS result is borderline (size-matched p = 0.02 for burden, 0.055
  for load). 2a and 2b use the same founder frequencies and are not independent of each other;
  section 1 is.
- No class shows a consistent **rise**: there is no sign of positive directional selection on SVs."""


md_2c = """### 2c. Climate of origin × garden climate

Each variant is labelled from its own founder carriers: **cold-specific** (carried by ≥ 2 cold-origin
ecotypes and no warm-origin ones), **warm-specific** (the mirror), **shared** (both); origin tertiles
of ecotype bio1. Its per-garden statistic is the residual *s*: the variant's temporal selection
coefficient in that garden (slope of logit pool frequency over generations) minus the median *s* of
SNPs matched on starting frequency **and** carrier-origin climate (`sv_purging_fixed.ipynb`).

- **Left, per garden:** median residual of the group's variants in each garden against garden bio1,
  ordinary least-squares line with 95% band; slope per °C, Spearman ρ. Insertions in colour, **short
  indels in grey as the null** — same estimator, same matching, same labelling.
- **Right, per variant:** each variant's own slope of residual on garden bio1 (> 0 = it does better
  in warmer gardens). Cold- vs warm-specific by Wilcoxon rank-sum, per class.
- **Test of SV-specificity:** does the slope of cold-specific insertions differ from that of
  cold-specific indels? Garden-level interaction (origin or class × bio1 on the 31 garden medians) and
  per-variant rank-sum. Unit: the SV event; the allele unit is reported for comparison."""

code_2c = r"""
import statsmodels.api as sm, statsmodels.formula.api as smf
ORD=["cold-specific","shared","warm-specific"]; OC={"cold-specific":"#2166AC","shared":"0.45","warm-specific":"#B2182B"}
KEY={"insertion":"ins","deletion":"del","indel":"indel"}; CL={"insertion":"insertions","deletion":"deletions","indel":"indels (null)"}
def load_unit(tag):
    V=pd.read_csv(f"{G}/sv_purging_annotated{tag}.csv"); Ds=pd.read_csv(f"{G}/sv_temporal_climate_matched_sites{tag}.csv").sort_values("bio1").reset_index(drop=True)
    return V,Ds,np.load(f"{G}/sv_temporal_climate_matched_persite{tag}.npz")
def gmed(V,Ds,Ps,cls,o):
    ri=V[(V.cls==cls)&(V.origin==o)].ridx.to_numpy()
    return np.array([np.nanmedian(Ps[f"{s}_res_{KEY[cls]}_matched"][ri]) for s in Ds.site.astype(int)])
VP,DS,PS=load_unit("_event2"); x=DS.bio1.to_numpy(); xx=np.linspace(x.min()-0.5,x.max()+0.5,60)
def ols(y):
    m=sm.OLS(y,sm.add_constant(x)).fit(); pr=m.get_prediction(sm.add_constant(xx)); return m,pr.predicted_mean,pr.conf_int(alpha=0.05)
fig,ax=plt.subplots(1,3,figsize=(15,4.3),sharey=True); rows=[]
for a,o in zip(ax,ORD):
    style(a); a.axhline(0,color="0.3",lw=1)
    for cls,col,mk,z in (("indel","0.6","s",2),("insertion",OC[o],"o",3)):
        y=gmed(VP,DS,PS,cls,o); m,yy,ci=ols(y); rho=stats.spearmanr(x,y); neg=y<0
        a.fill_between(xx,ci[:,0],ci[:,1],color=col,alpha=0.12 if cls=="insertion" else 0.08,lw=0)
        a.plot(xx,yy,color=col,lw=1.6 if cls=="insertion" else 1.2,ls="-" if cls=="insertion" else "--")
        a.scatter(x[neg],y[neg],s=30 if cls=="insertion" else 18,marker=mk,color=col,edgecolor="0.2",linewidth=0.4,zorder=z)
        a.scatter(x[~neg],y[~neg],s=30 if cls=="insertion" else 18,marker=mk,facecolor="white",edgecolor=col,linewidth=1,zorder=z)
        rows.append(dict(origin=o,cls=CL[cls],n=int(((VP.cls==cls)&(VP.origin==o)).sum()),gardens_below0=int(neg.sum()),
                         slope_per_C=m.params[1],se=m.bse[1],p=m.pvalues[1],rho=rho.statistic,p_rho=rho.pvalue,wilcoxon_p=stats.wilcoxon(y).pvalue))
    r=pd.DataFrame(rows).set_index(["origin","cls"])
    a.annotate(f"{o}\ninsertions (n = {r.loc[(o,'insertions'),'n']:,}): slope {r.loc[(o,'insertions'),'slope_per_C']:+.4f} ± {r.loc[(o,'insertions'),'se']:.4f} per °C, p = {r.loc[(o,'insertions'),'p']:.2g}, ρ = {r.loc[(o,'insertions'),'rho']:+.2f}\n"
               f"indels (n = {r.loc[(o,'indels (null)'),'n']:,}): slope {r.loc[(o,'indels (null)'),'slope_per_C']:+.4f} ± {r.loc[(o,'indels (null)'),'se']:.4f} per °C, p = {r.loc[(o,'indels (null)'),'p']:.2g}, ρ = {r.loc[(o,'indels (null)'),'rho']:+.2f}",
               xy=(0.03,0.04),xycoords="axes fraction",va="bottom",fontsize=7.3,color="0.25",bbox=dict(facecolor="white",alpha=0.85,edgecolor="none",pad=1.5))
    a.set_xlabel("garden temperature  bio1 (°C)")
ax[0].set_ylabel("median residual vs matched SNPs\n(< 0 = declines)"); ax[0].set_ylim(-0.36,0.3)
ax[0].annotate("2c  per garden: insertions (colour, solid) and short indels (grey, dashed); filled = below 0",xy=(0.0,1.02),xycoords="axes fraction",fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_key_2c_origin_by_garden.png",dpi=130,bbox_inches="tight"); plt.show()
T2c=pd.DataFrame(rows)
# per-variant slopes
fig,a=plt.subplots(figsize=(11,4)); style(a,axis="y"); a.axhline(0,color="0.35",lw=1,ls="--"); TL=[]; SL={}
for i,o in enumerate(ORD):
    for j,cls in enumerate(("insertion","deletion","indel")):
        v=VP[(VP.cls==cls)&(VP.origin==o)].slope_bio1.dropna().to_numpy(); SL[(o,cls)]=v; xp=i+(j-1)*0.27
        take=v if v.size<=300 else rng.choice(v,300,replace=False); col=OC[o] if cls!="indel" else "0.6"
        a.scatter(xp+rng.uniform(-0.09,0.09,take.size),take,s=6,color=col,alpha=0.35 if cls!="deletion" else 0.6,linewidth=0,zorder=2)
        a.boxplot(v,positions=[xp],widths=0.22,showfliers=False,manage_ticks=False,patch_artist=True,zorder=3,
                  boxprops=dict(facecolor="none",edgecolor="0.1",linewidth=1.1,linestyle="-" if cls=="insertion" else ":" if cls=="deletion" else "--"),
                  medianprops=dict(color="k",lw=1.8),whiskerprops=dict(color="0.2"),capprops=dict(color="0.2"))
        a.annotate(f"{np.median(v):+.4f}",(xp,-0.062),ha="center",fontsize=6.8,color="0.3")
    TL.append(f"{o}\n(ins n={SL[(o,'insertion')].size:,}, del n={SL[(o,'deletion')].size:,}, indel n={SL[(o,'indel')].size:,})")
a.set_xticks(range(3)); a.set_xticklabels(TL,fontsize=7.6); a.set_ylim(-0.07,0.07); a.set_ylabel("per-variant slope of residual vs garden bio1 (per °C)\n(> 0 = does better in warmer gardens)")
txt="\n".join(f"{CL[cls]}: cold- vs warm-specific rank-sum p = {stats.ranksums(SL[('cold-specific',cls)],SL[('warm-specific',cls)]).pvalue:.1g}" for cls in ("insertion","deletion","indel"))
a.annotate(txt,xy=(0.99,0.03),xycoords="axes fraction",ha="right",va="bottom",fontsize=7.5,color="0.2",bbox=dict(facecolor="white",alpha=0.9,edgecolor="none",pad=2))
a.annotate("boxes per origin: insertions (solid), deletions (dotted), indels (dashed, grey); numbers = medians",xy=(0.0,1.02),xycoords="axes fraction",fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_key_2c_origin_slopes.png",dpi=130,bbox_inches="tight"); plt.show()
print("2c  per garden (31 gardens): OLS slope of median residual on garden bio1")
print(T2c.round(5).to_string(index=False))
print("\n    per-variant slope, insertions vs indels within origin (rank-sum):")
for o in ORD: print(f"      {o:<14} medians {np.median(SL[(o,'insertion')]):+.4f} vs {np.median(SL[(o,'indel')]):+.4f}   p = {stats.ranksums(SL[(o,'insertion')],SL[(o,'indel')]).pvalue:.2g}")
def inter(frames,term):
    m=smf.ols(f"med ~ bio1*C({term})",pd.concat(frames)).fit(); k=[t for t in m.params.index if ":" in t][0]; return m.params[k],m.pvalues[k]
print("\n    garden-level interaction tests (OLS on the 31 garden medians per group):")
for cls in ("insertion","deletion","indel"):
    b,pv=inter([pd.DataFrame(dict(bio1=x,med=gmed(VP,DS,PS,cls,o),origin=o)) for o in ("cold-specific","warm-specific")],"origin")
    print(f"      {CL[cls]:<14} cold- vs warm-specific slope difference {b:+.4f} per °C, p = {pv:.2g}")
for o in ORD:
    b,pv=inter([pd.DataFrame(dict(bio1=x,med=gmed(VP,DS,PS,c,o),cls=c)) for c in ("insertion","indel")],"cls")
    print(f"      {o:<14} insertions vs indels slope difference {b:+.4f} per °C, p = {pv:.2g}")
VA,DA,PA=load_unit("_rec2"); xa=DA.bio1.to_numpy()
print("\n    allele unit (sv_purging_fixed_allele.ipynb), for comparison:")
for o in ORD:
    si=VA[(VA.cls=="insertion")&(VA.origin==o)].slope_bio1.dropna(); sn=VA[(VA.cls=="indel")&(VA.origin==o)].slope_bio1.dropna()
    gi=stats.linregress(xa,gmed(VA,DA,PA,"insertion",o)); gn=stats.linregress(xa,gmed(VA,DA,PA,"indel",o))
    print(f"      {o:<14} garden slope ins {gi.slope:+.4f} (p = {gi.pvalue:.2g}) vs indel {gn.slope:+.4f} (p = {gn.pvalue:.2g});  per-variant medians {si.median():+.4f} vs {sn.median():+.4f}, ins vs indel p = {stats.ranksums(si,sn).pvalue:.1g}")
"""

md_2c_read = """**Reading it.** The pattern is real, and it is not about SVs.

- **Variants private to cold-origin ecotypes lose ground in warm gardens** and hold in cold ones:
  cold-specific insertions fall by 0.0074 per °C of garden temperature (p = 0.02, ρ = −0.35), and
  their per-variant slopes are negative (median −0.007) where warm-specific insertions' are positive
  (+0.006; rank-sum p = 1e-19). Warm-specific variants show the mirror, more weakly at the garden
  level. That is the "cold SVs selected against in warm climates, warm SVs in cold" picture.
- **Indels and deletions do exactly the same.** Cold-specific indels: −0.0079 per °C (p = 0.01);
  cold- vs warm-specific rank-sum p = 2e-49; deletions p = 4e-9. Within the cold-specific group,
  insertions and indels have indistinguishable slopes (per-variant p = 0.74; garden-level
  interaction p = 0.9), and the same holds for shared and warm-specific variants.
- So this is the signature of **lineage turnover**: warm-origin ecotypes win in warm gardens,
  cold-origin ones in cold gardens, and every variant private to those lineages moves with them,
  SNP-matched or not. Matching SNPs on carrier-origin climate removes most of it, but a variant
  private to one origin group shares its founder background more completely than a SNP with the same
  mean carrier climate, so a residue remains. It is evidence of **local adaptation at the level of
  the founder ecotype**, read through their private variants — not of selection on the SVs themselves.
  The per-variant tests overstate the evidence (variants on the same lineages are not independent);
  the garden-level regressions on 31 points are the ones to quote.
- **Allele unit.** In `sv_purging_fixed_allele.ipynb` insertions look steeper than indels
  (cold-specific −0.012 vs −0.007 per °C), but the excess is the same in shared and warm-specific
  insertions too, i.e. uniform rather than origin-specific: it is the arch3 allele split (one
  insertion counted as several alleles with divided carriers), and it vanishes at the event unit."""

md_3 = """## 3. Which SVs are being purified?

The class under selection is insertions in genes and their regulatory flanks. What are they?

- **3a.** Per context: how many events and genes, how big, what the inserted sequence is
  (`insertion_tair10_class.csv`: TE-derived, duplicated gene sequence, other duplication, novel),
  how often frame-preserving, how often rare.
- **3b.** Inside CDS, which sub-classes carry the signal — in nature (share with ≤ 3 carriers) and in
  GrENE-Net (shared events: load change vs matched shared CDS SNPs, mean over gardens ± 95% CI).
- **3c.** The same sub-classes the classic way: folded SFS against 4-fold SNPs, KS / Mann–Whitney,
  Tajima's *D*."""

code_3 = r"""
I=E[E.kind=="ins"].copy(); I["inframe"]=I.size_bp%3==0
tab=I.groupby("context").apply(lambda d: pd.Series({"events":len(d),"genes":d.gene[d.gene.notna()&(d.gene!="")].nunique(),
     "median bp":d.size_bp.median(),"≥ 1 kb":(d.size_bp>=1000).mean(),"TE-derived":(d.seq_class=="te_derived").mean(),
     "gene dup.":(d.seq_class=="gene_dup").mean(),"novel":(d.seq_class=="novel").mean(),"in frame":d.inframe.mean(),
     "≤ 3 carriers":(d.n_carriers<=3).mean(),"shared":(d.origin=="shared").mean()})).reindex(CTX)
tab.loc[["TE","intergenic"],"genes"]=np.nan        # no gene assigned outside genes / promoters
print("3a  insertion events (>= 2 carriers) by context"); print(tab.round(3).to_string())
C=I[I.context=="CDS"].copy(); iC=np.where((E.kind=="ins")&(E.context=="CDS"))[0]
C["sub_seq"]=C.seq_class.map({"te_derived":"TE-derived","gene_dup":"gene duplicate","novel":"novel"}).fillna("other dup.")
C["sub_size"]=pd.cut(C.size_bp,[50,100,300,1000,1e7],labels=["51–100 bp","101–300 bp","301 bp–1 kb","> 1 kb"]).astype(str)
C["sub_frame"]=np.where(C.inframe,"in frame","frameshift")
snpCDS=SH&(KD=="snp")&(CX==0); cellE=np.digitize(E.n_carriers,CB)*10+np.digitize(E.carrier_bio1,BQ)
groups=[("all CDS insertions",np.ones(len(C),bool))]
for col in ("sub_seq","sub_size","sub_frame"):
    for v in (["TE-derived","gene duplicate","other dup.","novel"] if col=="sub_seq" else ["51–100 bp","101–300 bp","301 bp–1 kb","> 1 kb"] if col=="sub_size" else ["in frame","frameshift"]):
        groups.append((v,(C[col]==v).to_numpy()))
snp_r=(FR["snp"].n[FR["snp"].ctx==0]<=3).mean(); rows=[]
for lab,m in groups:
    d=C[m]; v=(d.n_carriers<=3).to_numpy(); bs=[rng.choice(v,v.size).mean() for _ in range(600)]
    sh=(d.origin=="shared").to_numpy(); idx=iC[m][sh]
    dd=load_diff(PE[idx],cellE[idx],snpCDS) if idx.size>=25 else np.full(31,np.nan)
    rows.append(dict(group=lab,n=len(d),rare_fold=v.mean()/snp_r,lo=np.percentile(bs,2.5)/snp_r,hi=np.percentile(bs,97.5)/snp_r,
                     n_shared=int(idx.size),load=np.nanmean(dd),lse=np.nanstd(dd,ddof=1)/np.sqrt(31),below=int(np.nansum(dd<0)),
                     p_load=stats.wilcoxon(dd).pvalue if idx.size>=25 else np.nan))
T3=pd.DataFrame(rows)
fig,ax=plt.subplots(1,2,figsize=(13,4.6),sharey=True); y=np.arange(len(T3)); sep=[0.5,4.5,8.5]
for a,(xc,lo,hi,ref,xl,ttl) in zip(ax,(("rare_fold","lo","hi",1,"share with ≤ 3 carriers, relative to CDS SNPs","3b  in nature"),
                                       ("load",None,None,0,"log load change vs matched shared CDS SNPs","3b  in GrENE-Net (shared events)"))):
    style(a,axis="x"); a.axvline(ref,color="0.35",lw=1,ls="--")
    err=[T3[xc]-T3[lo],T3[hi]-T3[xc]] if lo else 1.96*T3.lse
    a.errorbar(T3[xc],y,xerr=err,fmt="o",ms=6,color=KC["ins"],capsize=2)
    for s in sep: a.axhline(s,color="0.8",lw=0.8)
    a.set_xlabel(xl); a.annotate(ttl,xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5)
ax[0].set_yticks(y); ax[0].set_yticklabels([f"{g}  (n={n:,})" for g,n in zip(T3.group,T3.n)]); ax[0].set_ylim(len(T3)-0.5,-0.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_key_3_which.png",dpi=130,bbox_inches="tight"); plt.show()
print("\n3b  CDS insertion sub-classes"); print(T3.drop(columns=["lo","hi","lse"]).round(3).to_string(index=False))
# 3c  classic SFS of the CDS sub-classes against 4-fold SNPs
def maf2(names):
    c=np.concatenate([CS[f"{n}_c"] for n in names]).astype(float); k=np.concatenate([CS[f"{n}_k"] for n in names]).astype(float); f=c/k
    return np.minimum(f,1-f), np.sum([CS[f"{n}_sfs"] for n in names],0)
SUB=[("TE-derived",["insCDS_te_derived"],"#B2182B"),("other dup.",["insCDS_dispersed_dup","insCDS_local_dup"],"#EF8A62"),
     ("gene duplicate",["insCDS_gene_dup"],"#67A9CF"),("novel",["insCDS_novel"],"#2166AC")]
fig,a=plt.subplots(figsize=(9,3.8)); style(a,axis="y"); w=0.16; h=np.histogram(REF,bins)[0]
a.bar(xb-2*w,h/h.sum(),width=w*0.95,color="0.75",label="SNPs 4-fold (neutral)"); rows=[]
for j,(lab,names,col) in enumerate(SUB):
    v,sf=maf2(names); h=np.histogram(v,bins)[0]; a.bar(xb+(j-1)*w,h/h.sum(),width=w*0.95,color=col,label=f"{lab} (n={v.size:,})")
    ks=stats.ks_2samp(v,REF); d,lo,hi=tajD_ci(sf)
    rows.append(dict(cls=lab,n=v.size,singletons=(np.concatenate([CS[f"{n}_c"] for n in names])==1).mean(),maf_lt5=(v<0.05).mean(),
                     KS_D=ks.statistic,KS_p=ks.pvalue,MWU_p=stats.mannwhitneyu(v,REF,alternative="less").pvalue,TajD=f"{d:+.2f} [{lo:+.2f}, {hi:+.2f}]"))
a.set_xticks(xb); a.set_xticklabels([f"{b:.2f}" for b in bins[1:]],fontsize=7); a.set_xlabel("minor-allele frequency (bin upper edge)"); a.set_ylabel("proportion of variants")
a.legend(frameon=False,fontsize=7.5); a.annotate("3c  folded SFS, insertions in CDS by inserted sequence",xy=(0.02,1.02),xycoords="axes fraction",fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_key_3c_sfs_cds_subclasses.png",dpi=130,bbox_inches="tight"); plt.show()
print("\n3c  CDS insertion sub-classes vs 4-fold SNPs"); print(pd.DataFrame(rows).to_string(index=False,formatters={"singletons":"{:.3f}".format,"maf_lt5":"{:.3f}".format,"KS_D":"{:.3f}".format,"KS_p":"{:.2g}".format,"MWU_p":"{:.2g}".format}))
g=C.gene.value_counts(); print(f"\nCDS insertion events: {len(C)} in {g.size} genes; genes with >= 3 events: {(g>=3).sum()}")
print("most hit:", ", ".join(f"{k} ({v})" for k,v in g.head(10).items()))
"""

md_3_read = """**Reading it.**

- **Where.** Of 13.0k insertion events, only 715 sit in CDS (591 genes) and 790 in UTRs; 43% are in
  promoters. The purified class is small because most of it has already been removed (section 1a).
- **What they are.** CDS insertions are unlike the rest: a third are **duplicated gene sequence**
  (35%, vs 5–7% in promoter / intergenic / TE) and fewer are TE-derived (37% vs 46–63%); they are
  bimodal in size (many short, 38% ≥ 1 kb) and enriched for frame-preserving lengths (44%, carried
  by the duplicated-sequence classes).
- **Which ones carry the signal — the same sub-classes in nature and in GrENE-Net.** Inside CDS,
  **TE-derived** insertions are the rarest (2.04× vs CDS SNPs) and are the ones depleted in GrENE-Net
  (−8.5% beyond matched SNPs, 26/31 gardens); other (dispersed / local) duplications behave the same
  (1.83×; −9.2%, 27/31). **Duplicated gene sequence** is less rare (1.55×) and not depleted (−0.4%,
  14/31), and **novel sequence** looks neutral on both counts (1.05×; −0.1%).
- **The classic view agrees (3c).** Against 4-fold SNPs, TE-derived CDS insertions are the most skewed
  sub-class (49% singletons, Tajima's *D* = −2.46 [−2.55, −2.37], KS *D* = 0.36), other duplications
  next (−2.29), duplicated gene sequence intermediate (−1.88), novel sequence least (−1.18
  [−1.55, −0.76], KS *D* = 0.17) — the same ranking as the rarity ratio and as the GrENE-Net load
  change. Size: rarity rises with
  size in nature; in GrENE-Net the 301 bp–1 kb class falls most. In-frame and frameshifting insertions
  that survived are equally rare and equally depleted — the frame filter decided which insertions
  persist at all, not how they fare afterwards.
- **Deletions are different:** shared CDS deletions do not fall — they *gain* slightly on matched
  SNPs (+1.1%, 22/31 gardens above; p = 0.008).
- **Regulatory sequence.** UTR and promoter insertions are rarer than SNPs in nature (shared: 1.70×,
  1.41×) and UTR / intron burden predicts lower founder fitness, but their population load is not
  detectably depleted within three years.

## Bottom line

1. **Evolutionary history:** SVs carry the classic site-frequency footprint of purifying selection,
   graded by the function they disrupt. Insertions in coding sequence are the most skewed class in the
   genome (Tajima's *D* −2.09, against −1.40 for nonsynonymous and −0.58 for synonymous SNPs), sit in
   CDS at a fifth of the genomic expectation (SNPs: four fifths), and the duplicated-sequence
   insertions that survive there are frame-preserving; UTR and promoter insertions show a weaker
   version; insertions in TE and intergenic sequence are no more skewed than SNPs in the same
   sequence. Deletions: only in CDS.
2. **GrENE-Net:** the same class — insertions in coding sequence — keeps losing ground during the
   experiment: their carriers decline and their load falls ~1.5% beyond matched SNPs in 27/31 gardens.
   Small, not specific to a climate, and only partly separable from the general cost of genic
   variation. No positive selection on any SV class. Garden climate does sort the **founder
   lineages** — variants private to cold-origin ecotypes fall in warm gardens and vice versa — but
   SNP-matched indels and deletions do so equally, so that is local adaptation of the ecotypes seen
   through their private variants, not selection on SVs.
3. **Which SVs:** insertions into exons, and among them the **TE-derived and duplicated-repeat
   insertions** — the most skewed spectra in nature (Tajima's *D* −2.5 and −2.3, against −1.2 for
   novel sequence) and the only sub-classes measurably lost in GrENE-Net (−8 to −9%); duplicated gene
   fragments and novel sequence in exons are tolerated. UTR and promoter insertions
   carry a weaker, nature-only signal."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_setup),
    new_markdown_cell(md_1), new_code_cell(code_1), new_markdown_cell(md_1_read),
    new_markdown_cell(md_2), new_code_cell(code_2), new_markdown_cell(md_2_read),
    new_markdown_cell(md_2c), new_code_cell(code_2c), new_markdown_cell(md_2c_read),
    new_markdown_cell(md_3), new_code_cell(code_3), new_markdown_cell(md_3_read),
])
ep = ExecutePreprocessor(timeout=1800, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
