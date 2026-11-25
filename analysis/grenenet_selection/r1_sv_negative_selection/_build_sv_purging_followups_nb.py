#!/usr/bin/env python
"""Build+execute the follow-ups to `sv_purging_vs_hitchhiking.ipynb`:
  1. is the decline size-dependent?   2. candidate cold-specific purged insertions.
  3. coding vs non-coding.
Load-only from `_sv_purging_classify.py` (per-variant table). `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "notebooks", "sv_purging_followups.ipynb")

md_title = """# Follow-ups: which insertions are being purged, and why?

`sv_purging_vs_hitchhiking.ipynb` found that lineage-restricted insertions — cold-specific ones above
all — decline in GrENE-Net in almost every garden, including cold gardens where their carrier
lineages are winning, while shared insertions and equally-restricted indels do not. Three questions
follow, each a test of whether the declining insertions look deleterious:

1. **Size** — are bigger insertions purged more?
2. **Coding vs non-coding** — is the decline stronger for insertions in genes?
3. **Inserted sequence** — does it matter whether the insertion is TE-derived?
4. **Candidates** — which cold-specific insertions show the global decline, and what do they hit?
5. **Frequency spectra by context** — in the natural founder panel, are insertions in coding sequence
   shifted to rare frequencies more than elsewhere?

Everything uses the same per-variant statistics: each variant's residual `s` against SNPs matched on
starting frequency **and** carrier-origin climate, summarised over the 31 gardens as the **fraction
of gardens where it declines** (0.5 = no tendency)."""

code_setup = r"""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf
plt.rcParams.update({'figure.dpi':110,'font.size':8,'axes.linewidth':0.6})
G=os.path.abspath(os.path.join(os.getcwd(),"..","r1_sv_negative_selection","results","sv_adaptive"))
os.makedirs(f"{G}/plots",exist_ok=True)
def style(ax,axis="both"):
    ax.set_axisbelow(True); ax.grid(color="0.88",lw=0.7,axis=axis)
    for sp in ax.spines.values(): sp.set_visible(False)
def dotbox(a,x,v,color,rng,w=0.55,maxdots=400):
    take=v if v.size<=maxdots else rng.choice(v,maxdots,replace=False)
    a.scatter(x+rng.uniform(-w*0.36,w*0.36,take.size),take,s=7,color=color,alpha=0.55,linewidth=0,zorder=2)
    a.boxplot(v,positions=[x],widths=w,patch_artist=True,showfliers=False,
              boxprops=dict(facecolor="none",edgecolor="0.1",linewidth=1.1),medianprops=dict(color="k",lw=1.8),
              whiskerprops=dict(color="0.2",lw=0.9),capprops=dict(color="0.2",lw=0.9),manage_ticks=False,zorder=3)
V=pd.read_csv(f"{G}/sv_purging_classify_variants.csv")
V["lsize"]=np.log10(V.size_bp); V["lcar"]=np.log10(V.n_carriers)
ORD=["cold-specific","shared","warm-specific"]
DOT={"cold-specific":"#2166AC","shared":"0.55","warm-specific":"#B2182B"}
I=V[(V.cls=="insertion")&(V.origin.isin(ORD))].copy()
print(f"{len(I):,} tested insertions:", I.origin.value_counts().to_dict())
"""

md_size = """## 1. Is the decline size-dependent?

If the declining insertions are deleterious because of what they disrupt, bigger ones might be worse.
**Size is confounded with origin**, though: restricted insertions are about twice as long as shared
ones (medians 221 bp cold-specific, 156 bp warm-specific, 106 bp shared), and size also tracks
frequency (bigger insertions are rarer). So size has to be tested **within** each origin group, and
jointly with origin, carrier count and call rate."""

code_size = r"""
print("median insertion size (bp):", I.groupby("origin").size_bp.median().reindex(ORD).to_dict())
bins=[50,100,200,500,2000,1e7]; lab=["51-100","101-200","201-500","501-2k",">2k"]
I["sbin"]=pd.cut(I.size_bp,bins,labels=lab)
rng=np.random.default_rng(0)
fig,ax=plt.subplots(1,3,figsize=(13,4.3),sharey=True)
for a,o in zip(ax,ORD):
    style(a,axis="y"); a.axhline(0.5,color="0.35",lw=1.0,ls="--",zorder=1)
    d=I[I.origin==o]; TL=[]
    for i,b in enumerate(lab):
        v=d[d.sbin==b].frac_neg.dropna().to_numpy()
        if v.size>=5: dotbox(a,i,v,DOT[o],rng)
        TL.append(f"{b}\nn={v.size}")
    r=stats.spearmanr(d.lsize,d.frac_neg)
    a.annotate(f"{o}   within-group Spearman vs size: rho = {r.statistic:+.3f}, p = {r.pvalue:.2g}",
               xy=(0.02,1.01),xycoords="axes fraction",fontsize=8,color=DOT[o] if o!="shared" else "0.3")
    a.set_xticks(range(len(lab))); a.set_xticklabels(TL,fontsize=7.5); a.set_xlabel("insertion size (bp)")
    if o=="cold-specific": a.set_ylabel("fraction of gardens where residual < 0\n(0.5 = no tendency)")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_size.png",dpi=130,bbox_inches="tight"); plt.show()

print("confounds (Spearman with log size, all insertions):")
for c in ("p0","n_carriers","call_rate"):
    r=stats.spearmanr(I.lsize,I[c]); print(f"   {c:<11} rho={r.statistic:+.3f}  p={r.pvalue:.2g}")
m=smf.ols("frac_neg ~ lsize + C(origin, Treatment('shared')) + call_rate + lcar",data=I).fit()
print("\njoint model: frac_neg ~ log10 size + origin (vs shared) + call rate + log10 carriers")
print(m.summary().tables[1])
"""

md_size_read = """**Reading it.** No — bigger insertions are not purged more. Within cold-specific insertions the
decline is flat across sizes from 51 bp to over 2 kb (ρ ≈ +0.07, n.s.), shared insertions are flat
at ~0.48, and warm-specific ones go slightly the *other* way (bigger → declines less consistently).
In the joint model the size coefficient is small and negative while **cold-specific stays strongly
positive after size, carrier count and call rate are in the model**, so the cold-specific decline is
not a size effect even though those insertions are longer.

**Carrier count matters, and is worth keeping in view.** Insertions carried by fewer founders decline
more consistently, even after frequency matching. That is what purifying selection predicts —
deleterious variants are held rare — but it overlaps with "restricted". Cold-specific survives it;
warm-specific does not (its origin coefficient is n.s. once carriers are in the model)."""

md_ctx = """## 2. Coding vs non-coding

If the declining insertions are deleterious because they disrupt genes, the decline should be
strongest in coding sequence. Each tested insertion is given one context from the TAIR10 annotation
(`_sv_purging_annotate.py`), first match wins: **CDS > UTR > intron > reference TE > promoter
(≤ 2 kb upstream, strand-aware) > intergenic**. Boxes are drawn where a context holds ≥ 8
insertions; smaller groups are shown as dots only. Indels split the same way are in the table as the
null, not in the figure."""

code_ctx = r"""

A=pd.read_csv(f"{G}/sv_purging_annotated.csv")
CTX=["CDS","UTR","intron","TE","promoter","intergenic"]
AI=A[(A.cls=="insertion")&(A.origin.isin(ORD))].copy()
rng=np.random.default_rng(0)
fig,ax=plt.subplots(1,3,figsize=(13,4.3),sharey=True)
for a,o in zip(ax,ORD):
    style(a,axis="y"); a.axhline(0.5,color="0.35",lw=1.0,ls="--",zorder=1)
    d=AI[AI.origin==o]; TL=[]
    for i,c in enumerate(CTX):
        v=d[d.context==c].frac_neg.dropna().to_numpy()
        if v.size>=8: dotbox(a,i,v,DOT[o],rng)
        elif v.size: a.scatter(i+rng.uniform(-0.15,0.15,v.size),v,s=10,color=DOT[o],alpha=0.7,linewidth=0)
        TL.append(f"{c}\nn={v.size}")
    a.annotate(o,xy=(0.02,1.01),xycoords="axes fraction",fontsize=8.5,color=DOT[o] if o!="shared" else "0.3")
    a.set_xticks(range(len(CTX))); a.set_xticklabels(TL,fontsize=7.5)
    if o=="cold-specific": a.set_ylabel("fraction of gardens where residual < 0\n(0.5 = no tendency)")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_context.png",dpi=130,bbox_inches="tight"); plt.show()

for cls in ("insertion","indel"):
    d=A[(A.cls==cls)&A.origin.isin(ORD)]
    t=d.groupby(["context","origin"]).frac_neg.median().unstack().reindex(CTX)[ORD]
    n=d.groupby(["context","origin"]).size().unstack().reindex(CTX)[ORD]
    print(f"\n{cls}s: median fraction of gardens declining (n)")
    print(t.round(3).astype(str).add(" (").add(n.astype(str)).add(")").to_string())
print("\nshare of variants in CDS:")
for cls in ("insertion","deletion","indel"):
    d=A[(A.cls==cls)&A.origin.isin(ORD)]
    print(f"  {cls:<10}"+"  ".join(f"{o} {100*(d[d.origin==o].context=='CDS').mean():.1f}%" for o in ORD))
AI["ctx3"]=AI.context.map({"CDS":"coding","UTR":"coding","intron":"intron","TE":"outside","promoter":"outside","intergenic":"outside"})
AI["lcar"]=np.log10(AI.n_carriers)
m=smf.ols("frac_neg ~ C(ctx3, Treatment('outside')) * C(origin, Treatment('shared')) + lcar",data=AI).fit()
print("\njoint model, insertions: frac_neg ~ context (coding / intron / outside genes) x origin + log10 carriers")
print(m.summary().tables[1])
"""

md_ctx_read = """**Reading it.** The decline of cold-specific insertions is **not concentrated in genes**. It is
present in every context — UTR, intron, reference TE, promoter, intergenic all sit around 0.63–0.73 —
and CDS is, if anything, the lowest (0.55, on only 15 insertions). Indels split the same way stay
at 0.45–0.55 everywhere. The joint model finds no context effect and no context × origin
interaction; only cold-specific origin matters.

So "deleterious because they disrupt genes" is not supported. What the data point to is narrower:
**a cost that is specific to insertions but independent of where they land.** Most insertion cargo
is transposable-element sequence, so a position-independent TE cost is a plausible candidate — but
that is a hypothesis, and the direct test is to split by what was *inserted* (TE vs non-TE sequence)
rather than by where.

Cold-specific insertions are, separately, somewhat enriched in CDS (8.6% vs 5.0% shared, 2.3%
warm-specific), but that is 15 variants."""

md_seq = """## 3. What was inserted — TE-derived vs other sequence

Section 2 found no dependence on *where* an insertion lands. The next candidate was *what* it
carries: most insertion cargo is transposable-element sequence, and TE copies can carry costs
independent of position. Each tested insertion is classed by its inserted sequence
(`insertion_tair10_class.csv`, joined on the arch3 panel record index): **te_derived**, **novel**
(no TAIR10 hit — partly a reflection of alignment power, not guaranteed-novel sequence), and three
duplication classes."""

code_seq = r"""

S=pd.read_csv(f"{G}/sv_purging_insertions_seqclass.csv"); S=S[S.origin.isin(ORD)]
SQ=["te_derived","novel","dispersed_dup","gene_dup","local_dup"]
rng=np.random.default_rng(0)
fig,ax=plt.subplots(1,3,figsize=(13,4.3),sharey=True)
for a,o in zip(ax,ORD):
    style(a,axis="y"); a.axhline(0.5,color="0.35",lw=1.0,ls="--",zorder=1)
    d=S[S.origin==o]; TL=[]
    for i,c in enumerate(SQ):
        v=d[d.seq_class==c].frac_neg.dropna().to_numpy()
        if v.size>=8: dotbox(a,i,v,DOT[o],rng)
        elif v.size: a.scatter(i+rng.uniform(-0.15,0.15,v.size),v,s=10,color=DOT[o],alpha=0.7,linewidth=0)
        TL.append(f"{c.replace('_',chr(10))}\nn={v.size}")
    te=d[d.seq_class=="te_derived"].frac_neg; nt=d[d.seq_class!="te_derived"].frac_neg
    a.annotate(f"{o}   TE-derived vs rest: rank-sum p = {stats.ranksums(te,nt).pvalue:.2g}",
               xy=(0.02,1.01),xycoords="axes fraction",fontsize=8,color=DOT[o] if o!="shared" else "0.3")
    a.set_xticks(range(len(SQ))); a.set_xticklabels(TL,fontsize=7.3)
    if o=="cold-specific": a.set_ylabel("fraction of gardens where residual < 0\n(0.5 = no tendency)")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_purging_seqclass.png",dpi=130,bbox_inches="tight"); plt.show()
t=S.groupby(["seq_class","origin"]).frac_neg.agg(["median","size"])
print(t["median"].unstack().reindex(SQ)[ORD].round(3).astype(str).add(" (").add(t["size"].unstack().reindex(SQ)[ORD].astype(str)).add(")").to_string())
"""

md_seq_read = """**Reading it.** No difference either: TE-derived and non-TE insertions decline alike in every
origin group (cold-specific 0.63 vs 0.68, p = 0.75). So the TE-cost idea from section 2 is **not
supported**. Neither size, nor genomic context, nor inserted-sequence type predicts which insertions
decline — the only predictors so far are being lineage-restricted (cold-specific above all) and being
carried by few founders. A property of the insertion itself that explains the cost has not been
found."""

md_cand = """## 4. Candidates — cold-specific insertions with a global decline

Cold-specific insertions called **global** (sign test p < 0.05 over the 31 gardens), ranked by how
consistently they decline, then by mean residual. Gene = the overlapping gene, or the gene whose
2 kb promoter window holds the insertion; function from TAIR GO + UniProt
(`annotate_genes_tair_uniprot.py`).

**Read this as a candidate list, not a confirmed set.** Equally-restricted indels are called global
~32% of the time against ~49% for these insertions, so roughly a third of the 86 are expected by
chance. The top of the ranking is the most trustworthy part."""

code_cand = r"""

C=pd.read_csv(f"{G}/sv_purging_candidates.csv")
print(f"{len(C)} candidates")
print("  by genomic context:", C.context.value_counts().to_dict())
print("  by inserted sequence:", C.seq_class.value_counts().to_dict())
pd.set_option("display.width",260); pd.set_option("display.max_colwidth",70)
cols=["chrom","pos","size_bp","seq_class","te_family","context","gene","protein_name","n_cold","frac_neg","mean_res"]
show=C[[c for c in cols if c in C.columns]].copy()
show["frac_neg"]=show.frac_neg.round(3); show["mean_res"]=show.mean_res.round(3)
print("\nthe 5 candidates inside a coding sequence:")
print(show[show.context=="CDS"].to_string(index=False))
print("\ntop 25 by consistency:")
print(show.head(25).to_string(index=False))
print(f"\nfull list: {G}/sv_purging_candidates.csv")
"""

md_sfsc = """## 5. Frequency spectra by genomic context

The classic population-genetic test for deleteriousness: a class of variants under purifying
selection has a spectrum shifted towards **rare** variants relative to a neutral-ish class **in the
same context**. If insertions in coding sequence are deleterious, their excess of rare variants over
SNPs should be largest in CDS.

Built on the full 231-founder panel by `_sv_sfs_context.py`: every insertion and deletion record plus
a random 5% of SNPs and short indels, each given one context (same rules as section 2). Frequency =
carriers / called founders; kept if called in ≥ 90% of all 231 founders and carried by ≥ 1.
Singletons stay in. Insertions and deletions are shown as **bubbles** (alleles of one event merged) —
the unit to read, since allele-splitting inflates the rare end and does so unevenly across contexts —
with insertion alleles dashed for comparison.

Unlike sections 1–4, this is about **natural populations** (the founder panel), not GrENE-Net: it asks
whether insertions look purged over the species' history, which is independent evidence."""

code_sfsc = r"""

Z=np.load(f"{G}/sv_sfs_context.npz")
CTX=["CDS","UTR","intron","TE","promoter","intergenic"]
LN=[("snp","SNPs","0.35","-",1.6),("indel","short indels","0.7","-",1.4),
    ("ins_bubble","insertions (bubbles)","#D55E00","-",1.9),("ins_allele","insertions (alleles)","#D55E00","--",1.2),
    ("del_bubble","deletions (bubbles)","#009E73","-",1.7)]
fig,ax=plt.subplots(2,3,figsize=(13,7),sharex=True,sharey=True); ax=ax.ravel()
for a,x in zip(ax,CTX):
    style(a)
    for key,lab,col,ls,lw in LN:
        f=np.sort(Z[f"{key}_{x}"])
        if f.size<20: continue
        a.plot(f,np.arange(1,f.size+1)/f.size,color=col,ls=ls,lw=lw,label=lab)
    a.set_xlim(0,0.5); a.set_ylim(0,1.02)
    a.annotate(f"{x}\nins bubbles n={Z[f'ins_bubble_{x}'].size:,}  SNPs n={Z[f'snp_{x}'].size:,}",
               xy=(0.97,0.05),xycoords="axes fraction",ha="right",va="bottom",fontsize=8,color="0.25")
for a in ax[3:]: a.set_xlabel("ALT frequency among called founders")
for a in (ax[0],ax[3]): a.set_ylabel("cumulative distribution")
ax[0].legend(frameon=False,fontsize=7.5,loc="center right")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_sfs_context_ecdf.png",dpi=130,bbox_inches="tight"); plt.show()

# summary: share of rare variants (<= 5%) relative to SNPs in the same context, bootstrap 95% CI
def fold(f,s,B=1000,seed=0):
    r=np.random.default_rng(seed); est=(f<=0.05).mean()/(s<=0.05).mean()
    bs=[(r.choice(f,f.size,True)<=0.05).mean()/(r.choice(s,s.size,True)<=0.05).mean() for _ in range(B)]
    return est,np.percentile(bs,2.5),np.percentile(bs,97.5)
SUMC=[("indel","short indels","0.6"),("ins_bubble","insertions (bubbles)","#D55E00"),
      ("ins_allele","insertions (alleles)","#F4A582"),("del_bubble","deletions (bubbles)","#009E73")]
fig,a=plt.subplots(figsize=(10,4.3)); style(a,axis="y"); a.axhline(1,color="0.35",lw=1.0,ls="--")
rows=[]
for k,(key,lab,col) in enumerate(SUMC):
    for i,x in enumerate(CTX):
        f=Z[f"{key}_{x}"]; s=Z[f"snp_{x}"]
        if f.size<20: continue
        e,lo,hi=fold(f,s)
        xx=i+(k-1.5)*0.18
        a.errorbar([xx],[e],yerr=[[e-lo],[hi-e]],fmt="o",color=col,ms=6,capsize=2,lw=1.2,
                   label=lab if i==0 else None)
        rows.append(dict(cls=lab,context=x,n=f.size,rare=(f<=0.05).mean(),snp_rare=(s<=0.05).mean(),fold=e,lo=lo,hi=hi))
a.set_xticks(range(len(CTX))); a.set_xticklabels(CTX)
a.set_ylabel("rare share (<= 5%) relative to\nSNPs in the same context")
a.legend(frameon=False,fontsize=8,loc="upper right")
fig.tight_layout(); fig.savefig(f"{G}/plots/sv_sfs_context_fold.png",dpi=130,bbox_inches="tight"); plt.show()
R=pd.DataFrame(rows)
print(R.assign(rare=lambda d:(100*d.rare).round(1),snp_rare=lambda d:(100*d.snp_rare).round(1),
               fold=lambda d:d.fold.round(2),lo=lambda d:d.lo.round(2),hi=lambda d:d.hi.round(2)).to_string(index=False))
"""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_setup),
    new_markdown_cell(md_size), new_code_cell(code_size), new_markdown_cell(md_size_read),
    new_markdown_cell(md_ctx), new_code_cell(code_ctx), new_markdown_cell(md_ctx_read),
    new_markdown_cell(md_seq), new_code_cell(code_seq), new_markdown_cell(md_seq_read),
    new_markdown_cell(md_cand), new_code_cell(code_cand),
    new_markdown_cell(md_sfsc), new_code_cell(code_sfsc),
])
ep = ExecutePreprocessor(timeout=1800, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
