#!/usr/bin/env python
"""Characterize the non-SNP new-peak BLOCKS (not genes) — LFMM, raw, clq0.9 tiling.

The set: every clq0.9 tiling block that, on >=1 of the 20 climate axes, has a lead
(min-p) non-SNP (SV+small-indel) record clearing per-class Bonferroni while its lead
SNP record does NOT (SNP-absent blocks count as SNP-not-sig — the "SNPs missed it"
case). Input: production `results/multiaxis/tile_cr50` (MAF >= 0.05 + founder call rate
>= 0.50; small indels per record, SVs per Truvari event; unmerged tiling blocks).
Question: are these a special kind of genomic region —
TE/repeat-rich, non-genic, pericentromeric — i.e. places SNP calling is blind, or are
they scattered like the SNP peaks?

Three pillars, every enrichment against a MATCHED BACKGROUND = all tiling blocks that
carry a testable non-SNP record (the foreground is a subset; both counts printed in-notebook):
  A  feature/TE composition  (TAIR10 genes_transposons.gff: gene / TE / TE-gene /
     pseudogene / ncRNA coverage; TE-superfamily mix from the Alias field)
  B  structural context      (distance to centromere, block size, SNP-sparsity, SV size)
  C  signal-quality overlay  (lead MAF; per-axis new-peak count next to that axis' lambda)

All blocks are kept but split by recurrence: 1-axis (noisy tail) vs >=2-axis (robust core).
Raw uncalibrated p throughout, never GIF-corrected (project rule; pillar C states lambda as a
caveat instead of a "survives correction" test, which earlier versions computed). No chart titles
(project convention): markdown headers + axis labels + in-panel annotations carry context.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
import os

HERE = os.path.dirname(os.path.abspath(__file__))
# rendered notebooks live in the tree-level notebooks/, not beside their builder
OUT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "notebooks", "nonsnp_block_characterization.ipynb"))

C = []
md = lambda s: C.append(new_markdown_cell(s))
co = lambda s: C.append(new_code_cell(s))

md(r"""# What are the non-SNP new-peak blocks? — feature / structural / signal characterization

**The set.** clq0.9 **tiling** blocks where, on >=1 of the 20 climate axes (bio1-19 + pc1),
the lead (min-p) **non-SNP** (SV + small-indel) LFMM record clears per-class Bonferroni
(0.05/n) while the lead **SNP** record in that same block does **not** — SNP-absent blocks
count as SNP-not-significant (the "SNPs missed it" case, which is the whole point). gen9,
LFMM K=16, raw p.

**Input.** Production `results/multiaxis/tile_cr50` (2026-09-30 re-run): MAF >= 0.05 **and**
founder call rate >= 0.50. Small indels are tested per record, SVs per **Truvari event** (pooled
non-SNP = small-indel records + SV events; an event's pos / size come from its representative
allele). Blocks are the clq r2=0.9 tiling **without** `merge_small_blocks`, so a block id names
the same LD interval for SNPs and non-SNPs — the SNP-vs-non-SNP contrast within a block is
like-for-like (the old per-class merge mis-keyed ~27% of SV records).

**Question.** Are these a *special kind of genomic region* — TE/repeat-rich, non-genic,
pericentromeric (places short-read SNP calling goes blind) — or are they scattered like the
SNP peaks? Answered against a **matched background** = every tiling block carrying a testable
non-SNP record (the foreground is a subset; counts printed below). Whole-genome would be the
wrong denominator.

**Three pillars**
- **A — feature / TE composition:** per-block coverage by protein-coding gene, TE
  (`transposable_element` + `transposon_fragment`), TE-gene, pseudogene, ncRNA; TE-superfamily
  mix (Alias field); a genic / TE / intergenic classification.
- **B — structural context:** distance to centromere, block physical size, how many blocks are
  genuinely SNP-sparse/absent (makes the caveat quantitative), lead-SV size profile.
- **C — signal-quality overlay:** lead non-SNP MAF; per-axis new-peak count next to that axis'
  genomic-inflation lambda — ties the characterization to the raw-p inflation caveat.

All blocks kept, **split by recurrence**: `1-axis` (noisy tail) vs `>=2-axis` (robust core); the
split is printed below. Raw uncalibrated p throughout — never GIF-corrected (project rule).""")

# ---------------------------------------------------------------- setup ------
co(r"""import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
# NB: `__file__` is NOT defined in a Jupyter kernel; walk up from cwd to the tree root.
def _tree_root(start=None):
    d = os.path.abspath(start or os.getcwd())
    while d != os.path.dirname(d):
        if os.path.basename(d) == "grenenet_selection" and os.path.exists(os.path.join(d, "lib.py")):
            return d
        d = os.path.dirname(d)
    raise RuntimeError("could not locate the grenenet_selection tree root from " + os.getcwd())
sys.path.insert(0, _tree_root())
import lib
GEA = lib.GEA
# 2026-09-30 production re-run (MAF 0.05 + founder call rate 0.50, SV/non-SNP on truvari
# EVENTS, unmerged tiling blocks). Old input: _superseded_pre_cr50/wza_in_clq09_tile.
WZAIN = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/tile_cr50"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
OUTDIR = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
plt.rcParams.update({"figure.dpi": 110, "axes.spines.top": False, "axes.spines.right": False})

CHROM_LEN = {"Chr1":30427671,"Chr2":19698289,"Chr3":23459830,"Chr4":18585056,"Chr5":26975502}
# TAIR10 centromere midpoints (approx CEN, bp) — standard cited positions.
CEN = {"Chr1":15_086_000,"Chr2":3_607_000,"Chr3":13_799_000,"Chr4":3_956_000,"Chr5":11_725_000}

# colours: recurrence classes and background
COL = {"1-axis":"#9DB4C0", ">=2-axis":"#1B5E20", "background":"#CED4DA"}

def block_spans_tiling(r2=0.9):
    tag=f"clq{r2}"; rows=[]
    for ci in range(1,6):
        ch=f"Chr{ci}"
        g=pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv",sep="\t").sort_values("start_pos").reset_index(drop=True)
        ends=g.end_pos.to_numpy(np.int64)
        for i in range(len(g)):
            lo=1 if i==0 else int(ends[i-1])+1
            hi=int(ends[i]) if i<len(g)-1 else CHROM_LEN[ch]     # cap last block at chrom length
            rows.append((f"{ch}_{i}",ch,lo,hi))
    return pd.DataFrame(rows,columns=["block","chrom","start","end"]).set_index("block")
SPANS = block_spans_tiling()
SPANS["len"] = SPANS.end - SPANS.start + 1
print("tiling blocks:", len(SPANS))""")

# ------------------------------------------------ feature GFF + coverage -----
co(r"""# --- load TAIR10 genes+transposons GFF into per-class interval tables ---
GFF = lib.TAIR10_GENES_TE
cols = ["chrom","src","feat","start","end","score","strand","frame","attr"]
gff = pd.read_csv(GFF, sep="\t", comment="#", header=None, names=cols,
                  usecols=["chrom","feat","start","end","attr"])
gff = gff[gff.chrom.isin(CHROM_LEN)].copy()

FEATSETS = {
    "gene":        ["gene"],                                   # protein-coding (this GFF)
    "te":          ["transposable_element","transposon_fragment"],
    "te_gene":     ["transposable_element_gene"],
    "pseudogene":  ["pseudogene"],
    "ncRNA":       ["ncRNA","miRNA","tRNA","snoRNA","snRNA","rRNA"],
}
FEAT = {k: gff[gff.feat.isin(v)][["chrom","start","end"]].copy() for k,v in FEATSETS.items()}
for k,v in FEAT.items(): print(f"  {k:11s} {len(v):7d} features")

# TE superfamily from the Alias= tag on transposable_element rows (coarse, name-prefix).
te = gff[gff.feat=="transposable_element"].copy()
te["alias"] = te.attr.str.extract(r"Alias=([^;]+)")[0].fillna("")
def te_superfamily(a):
    a=a.upper()
    if "COPIA" in a: return "LTR/Copia"
    if "GYPSY" in a or "HILA" in a or a.startswith("ATLANTYS"): return "LTR/Gypsy"
    if "LINE" in a: return "LINE"
    if "SINE" in a: return "SINE"
    if "ENSPM" in a or "CACTA" in a: return "DNA/EnSpm-CACTA"
    if a.startswith("VANDAL") or "MUDR" in a or a.startswith("ATMU"): return "DNA/MuDR"
    if a.startswith("ATREP") or "HELITRON" in a: return "DNA/Helitron"
    if "HAT" in a or "POGO" in a or "MARINER" in a or "TAG1" in a: return "DNA/hAT-Tc1"
    if "HARB" in a: return "DNA/Harbinger"
    return "other/unknown"
te["superfam"] = te.alias.map(te_superfamily)
print("\nTE superfamily counts (all genome TEs):")
print(te.superfam.value_counts().to_string())""")

co(r"""# --- fast interval-union coverage: covered bp of each block by a feature class ---
def _merge(s, e):
    o=np.argsort(s); s=s[o].astype(np.int64); e=e[o].astype(np.int64)
    ms=[s[0]]; me=[e[0]]
    for i in range(1,len(s)):
        if s[i] <= me[-1]+1: me[-1]=max(me[-1],e[i])
        else: ms.append(s[i]); me.append(e[i])
    return np.array(ms), np.array(me)

def _cov_upto(x, ms, me, prefix):
    x=np.asarray(x,np.int64); res=np.zeros(len(x),np.int64)
    pos=np.searchsorted(ms, x, side="right")            # intervals [0,pos) have start<=x
    m=pos>0; idx=pos[m]-1
    last=np.minimum(me[idx], x[m]) - ms[idx] + 1
    res[m]=prefix[idx] + np.maximum(last,0)
    return res

def coverage_bp(blocks, feat):
    # covered bp of each block (rows of `blocks`: chrom,start,end) by union of `feat` intervals
    out=np.zeros(len(blocks),np.int64)
    bidx=np.arange(len(blocks))
    for ch in CHROM_LEN:
        bm = blocks.chrom.values==ch
        if not bm.any(): continue
        fm = feat[feat.chrom==ch]
        if len(fm)==0: continue
        ms,me=_merge(fm.start.values, fm.end.values)
        prefix=np.concatenate([[0], np.cumsum(me-ms+1)])[:-1]   # prefix[k]=cov bp intervals 0..k-1
        bs=blocks.start.values[bm]; be=blocks.end.values[bm]
        cov=_cov_upto(be,ms,me,prefix) - _cov_upto(bs-1,ms,me,prefix)
        out[bidx[bm]]=cov
    return out

def cov_superfam(blocks):
    # per-block TE bp broken down by superfamily (returns dict superfam->bp array)
    return {sf: coverage_bp(blocks, te[te.superfam==sf][["chrom","start","end"]])
            for sf in te.superfam.unique()}
print("coverage helpers ready")""")

# ------------------------------------------------ foreground ----------------
co(r"""# --- recompute the non-SNP new-peak blocks over 20 axes (slow: loads 40 CSVs) ---
def load(cls, axis):
    f=f"{WZAIN}/lfmm_{cls}_gen9_{axis}.csv"
    if not os.path.exists(f): return None
    d=pd.read_csv(f); d=d[d.MAF>0.05].copy()
    d["nlp"]=-np.log10(d.pval.clip(lower=1e-300)); return d
def leads(d):
    d=d[d.block.notna()&(d.block!="")]; return d.loc[d.groupby("block").nlp.idxmax()]
def gif_lambda(p):
    p=np.clip(np.asarray(p,float),1e-300,1.0)
    return float(np.median(stats.chi2.isf(p,1))/stats.chi2.isf(0.5,1))

rows=[]; LAM={}
for a in AXES:
    s,c=load("snp",a),load("nonsnp",a)
    if s is None or c is None: print("missing",a); continue
    LAM[a]=gif_lambda(c.pval)                                  # per-axis genomic inflation (nonsnp)
    bs,bc=-np.log10(0.05/len(s)),-np.log10(0.05/len(c))
    ls_,lc_=leads(s),leads(c)
    m=lc_.merge(ls_[["block","nlp"]],on="block",suffixes=("","_snp"),how="left")
    m["nlp_snp"]=m["nlp_snp"].fillna(0.0)
    npk=m[(m.nlp>bc)&(m.nlp_snp<=bs)].copy()
    npk["axis"]=a; npk["bonf_c"]=bc
    rows.append(npk[["block","chrom","pos","ref_len","alt_len","MAF","nlp","nlp_snp","axis","bonf_c"]])
PK=pd.concat(rows,ignore_index=True)
PK["size"]=(PK.alt_len-PK.ref_len).abs()
PK["lam"]=PK["axis"].map(LAM)                                  # stated as a caveat only; p stays raw
print(f"(block,axis) new peaks: {len(PK)} | unique blocks: {PK.block.nunique()}")
print("per-axis genomic inflation lambda (nonsnp; stated, never applied): min %.2f  median %.2f  max %.2f"%(min(LAM.values()),np.median(list(LAM.values())),max(LAM.values())))""")

co(r"""# --- collapse to one row per block; attach recurrence, best signal, SNP presence/density ---
# variant-density + SNP presence from bio1 (block membership is axis-invariant)
snp_b=load("snp","bio1"); nsnp_b=load("nonsnp","bio1")
n_snp   = snp_b.groupby("block").size().rename("n_snp")
n_nonsnp= nsnp_b.groupby("block").size().rename("n_nonsnp")
snp_blocks=set(snp_b.block.unique())
BG_blocks = sorted(nsnp_b.block.dropna().unique())             # matched background (count printed below)

g=PK.sort_values("nlp",ascending=False)
best=g.drop_duplicates("block").set_index("block")
FG=pd.DataFrame(index=sorted(PK.block.unique()))
FG["n_axes"]=PK.groupby("block")["axis"].nunique()
FG["best_axis"]=best["axis"]; FG["best_nlp"]=best.nlp
FG["lead_size"]=best["size"]; FG["lead_MAF"]=best.MAF; FG["lead_pos"]=best.pos
FG["recur"]=np.where(FG.n_axes>=2,">=2-axis","1-axis")
FG=FG.join(n_snp).join(n_nonsnp)
FG["n_snp"]=FG.n_snp.fillna(0).astype(int); FG["n_nonsnp"]=FG.n_nonsnp.fillna(0).astype(int)
FG["snp_present"]=FG.index.isin(snp_blocks)
FG=FG.join(SPANS[["chrom","start","end","len"]])
FG["dist_cen"]=FG.apply(lambda r:abs((r.start+r.end)/2 - CEN[r.chrom]),axis=1)
print("foreground blocks:",len(FG),"| matched background blocks:",len(BG_blocks))
print(FG.recur.value_counts().to_string())
print("SNP-absent blocks:",int((~FG.snp_present).sum()),f"({100*(~FG.snp_present).mean():.0f}%)")
display(FG.head())""")

co(r"""# --- feature coverage for foreground AND background; classify each block ---
BG=SPANS.loc[[b for b in BG_blocks if b in SPANS.index]].reset_index()[["block","chrom","start","end","len"]]
def add_cov(df):
    df=df.copy()
    for k in FEATSETS:
        df[f"cov_{k}"]=coverage_bp(df, FEAT[k])
    L=df["len"].values
    for k in FEATSETS: df[f"f_{k}"]=df[f"cov_{k}"]/L
    # TE-any = UNION of TE + TE-gene bp (not sum of fractions, which double-counts overlaps)
    df["f_te_any"]=coverage_bp(df, pd.concat([FEAT["te"],FEAT["te_gene"]]))/L
    df["f_intergenic"]=1-np.clip(df["f_gene"]+df["f_te_any"]+df["f_pseudogene"]+df["f_ncRNA"],0,1)
    def klass(r):
        if r.f_te_any>=0.5: return "TE"
        if r.f_gene>=0.5:   return "genic"
        if r.f_te_any>=0.2: return "TE-mixed"
        if r.f_gene>=0.2:   return "genic-mixed"
        return "intergenic"
    df["klass"]=df.apply(klass,axis=1)
    return df
FGc=add_cov(FG.reset_index().rename(columns={"index":"block"}))
BGc=add_cov(BG)
FGc["recur"]=FGc.block.map(FG.recur)
FGc.to_csv(f"{OUTDIR}/nonsnp_block_characterization.csv",index=False)
print("wrote nonsnp_block_characterization.csv ;", FGc.shape)
print("\nclass composition (%):")
comp=pd.concat([
    FGc[FGc.recur=="1-axis"].klass.value_counts(normalize=True).rename("FG 1-axis"),
    FGc[FGc.recur==">=2-axis"].klass.value_counts(normalize=True).rename("FG >=2-axis"),
    BGc.klass.value_counts(normalize=True).rename("background"),
],axis=1).fillna(0)*100
print(comp.round(1).to_string())""")

# ------------------------------------------------ Pillar A -------------------
md(r"""## Pillar A — feature / TE composition

Per-block coverage fractions (foreground split by recurrence vs matched background). If the
non-SNP-only signal lives where SNP calling is blind, we expect **higher TE coverage / lower
gene coverage** in the foreground than background. **ECDFs** (not histograms) so the two
foreground classes and the background overlay cleanly; enrichment by Mann-Whitney (coverage
fraction) + Fisher (block overlaps feature yes/no). Right: TE-superfamily bp mix; far right:
genic/TE/intergenic class shares.""")

co(r"""def ecdf(ax, vals, label, color, lw=1.8):
    v=np.sort(np.asarray(vals,float)); y=np.arange(1,len(v)+1)/len(v)
    ax.plot(v,y,label=f"{label} (n={len(v)})",color=color,lw=lw)

def panel_frac(ax, key, xlabel):
    for r,cc in [("1-axis",COL["1-axis"]),(">=2-axis",COL[">=2-axis"])]:
        ecdf(ax, FGc[FGc.recur==r][key], f"FG {r}", cc)
    ecdf(ax, BGc[key], "background", COL["background"], lw=1.4)
    ax.set_xlabel(xlabel,fontsize=9); ax.set_ylabel("ECDF",fontsize=9)
    ax.legend(frameon=False,fontsize=7,loc="lower right")

def enr(key):
    fg=FGc[key].values; bg=BGc[key].values
    u,p=stats.mannwhitneyu(fg,bg,alternative="greater")
    # Fisher on "overlaps feature" (frac>0)
    a=int((fg>0).sum()); b=len(fg)-a; c=int((bg>0).sum()); d=len(bg)-c
    orr,pf=stats.fisher_exact([[a,b],[c,d]])
    return dict(key=key, fg_med=np.median(fg), bg_med=np.median(bg), mwu_p=p,
                fg_ovl=a/len(fg), bg_ovl=c/len(bg), fisher_or=orr, fisher_p=pf)

fig,axs=plt.subplots(1,3,figsize=(16,4.2))
panel_frac(axs[0],"f_te_any","TE coverage fraction")
panel_frac(axs[1],"f_gene","protein-coding gene coverage fraction")
panel_frac(axs[2],"f_intergenic","intergenic fraction")
fig.tight_layout(); plt.show()

print("enrichment (foreground vs background; MWU alt=greater on fraction, Fisher on overlap y/n):")
print(pd.DataFrame([enr("f_te_any"),enr("f_gene"),enr("f_pseudogene"),enr("f_ncRNA"),
                    enr("f_intergenic")]).set_index("key").round(4).to_string())""")

co(r"""# TE-superfamily bp mix (foreground vs background) + class-share bars
sf_fg=cov_superfam(FGc); sf_bg=cov_superfam(BGc)
def sf_frac(covd):
    tot=max(sum(v.sum() for v in covd.values()),1)
    return {sf: covd[sf].sum()/tot for sf in covd}
FGsf=sf_frac(sf_fg); BGsf=sf_frac(sf_bg)
order=[k for k,_ in sorted(BGsf.items(),key=lambda kv:-kv[1])]
fig,axs=plt.subplots(1,2,figsize=(14,4.2))
x=np.arange(len(order)); w=0.38
axs[0].bar(x-w/2,[FGsf[k] for k in order],w,label=f"foreground ({len(FGc)})",color=COL[">=2-axis"])
axs[0].bar(x+w/2,[BGsf[k] for k in order],w,label="background",color=COL["background"])
axs[0].set_xticks(x); axs[0].set_xticklabels(order,rotation=40,ha="right",fontsize=7)
axs[0].set_ylabel("share of TE bp",fontsize=9); axs[0].legend(frameon=False,fontsize=8)
# genic/TE/intergenic class shares
ks=["genic","genic-mixed","TE","TE-mixed","intergenic"]
share=pd.concat([
    FGc[FGc.recur=="1-axis"].klass.value_counts(normalize=True).rename("FG 1-axis"),
    FGc[FGc.recur==">=2-axis"].klass.value_counts(normalize=True).rename("FG >=2-axis"),
    BGc.klass.value_counts(normalize=True).rename("background")],axis=1).reindex(ks).fillna(0)
share.T.plot(kind="barh",stacked=True,ax=axs[1],
             color=["#2E7D32","#A5D6A7","#8E44AD","#CBA6D6","#CED4DA"],width=0.7)
axs[1].set_xlabel("share of blocks",fontsize=9)
axs[1].legend(frameon=False,fontsize=7,loc="center left",bbox_to_anchor=(1.0,0.5))
fig.tight_layout(); plt.show()
print("share of TE bp by superfamily (foreground vs background):")
print(pd.DataFrame({"foreground":FGsf,"background":BGsf}).reindex(order).round(3).to_string())""")

# ------------------------------------------------ Pillar B -------------------
md(r"""## Pillar B — structural context

Are the blocks pericentromeric, physically large, and SNP-sparse? **Left:** distance from
block midpoint to its chromosome's centromere (foreground vs background — a leftward shift =
pericentromeric enrichment). **Middle:** block length. **Right:** SNP density (# SNP records
per block) with the SNP-absent fraction called out — the quantitative form of the "SNPs missed
it" caveat. Below: lead-variant size (SV vs small-indel) split by recurrence.""")

co(r"""fig,axs=plt.subplots(1,3,figsize=(16,4.2))
# distance to centromere (Mb)
for r,cc in [("1-axis",COL["1-axis"]),(">=2-axis",COL[">=2-axis"])]:
    ecdf(axs[0], FGc[FGc.recur==r].block.map(FG.dist_cen)/1e6, f"FG {r}", cc)
ecdf(axs[0], BG.apply(lambda r:abs((r.start+r.end)/2-CEN[r.chrom]),axis=1)/1e6,"background",COL["background"],1.4)
axs[0].set_xlabel("distance to centromere (Mb)",fontsize=9); axs[0].set_ylabel("ECDF",fontsize=9)
axs[0].legend(frameon=False,fontsize=7,loc="lower right")
# block length (log)
for r,cc in [("1-axis",COL["1-axis"]),(">=2-axis",COL[">=2-axis"])]:
    ecdf(axs[1], FGc[FGc.recur==r].len, f"FG {r}", cc)
ecdf(axs[1], BGc.len,"background",COL["background"],1.4)
axs[1].set_xscale("log"); axs[1].set_xlabel("block length (bp)",fontsize=9); axs[1].set_ylabel("ECDF",fontsize=9)
axs[1].legend(frameon=False,fontsize=7,loc="lower right")
# SNP density
for r,cc in [("1-axis",COL["1-axis"]),(">=2-axis",COL[">=2-axis"])]:
    ecdf(axs[2], FG[FG.recur==r].n_snp, f"FG {r}", cc)
ecdf(axs[2], n_snp.reindex(BG_blocks).fillna(0),"background",COL["background"],1.4)
axs[2].set_xscale("symlog"); axs[2].set_xlabel("# SNP records in block",fontsize=9); axs[2].set_ylabel("ECDF",fontsize=9)
axs[2].legend(frameon=False,fontsize=7,loc="lower right")
axs[2].annotate(f"SNP-absent: FG {100*(~FG.snp_present).mean():.0f}%  vs  bg "
                f"{100*(n_snp.reindex(BG_blocks).isna()).mean():.0f}%",
                xy=(0.03,0.5),xycoords="axes fraction",fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2",fc="white",ec="0.7"))
fig.tight_layout(); plt.show()

# stats
d1=stats.mannwhitneyu(FG.dist_cen, BG.apply(lambda r:abs((r.start+r.end)/2-CEN[r.chrom]),axis=1),alternative="less")
print(f"dist-to-centromere FG<bg  MWU p={d1.pvalue:.2e}  (FG median {FG.dist_cen.median()/1e6:.2f} Mb vs bg {(BG.apply(lambda r:abs((r.start+r.end)/2-CEN[r.chrom]),axis=1)).median()/1e6:.2f} Mb)")
sA=int((~FG.snp_present).sum()); print(f"SNP-absent foreground blocks: {sA}/{len(FG)} ({100*sA/len(FG):.0f}%)")
_bg_nsnp=n_snp.reindex(BG_blocks).fillna(0)
print(f"block length median: FG {FG.len.median():.0f} bp vs bg {BG.len.median():.0f} bp  "
      f"(MWU FG>bg p={stats.mannwhitneyu(FG.len,BG.len,alternative='greater').pvalue:.1e})")
print(f"SNP records per block median: FG {FG.n_snp.median():.0f} vs bg {_bg_nsnp.median():.0f}  "
      f"(MWU FG>bg p={stats.mannwhitneyu(FG.n_snp,_bg_nsnp,alternative='greater').pvalue:.1e})")""")

co(r"""# lead-variant size split by recurrence (SV >50 bp vs small-indel <=50 bp)
fig,ax=plt.subplots(figsize=(7,4))
for r,cc in [("1-axis",COL["1-axis"]),(">=2-axis",COL[">=2-axis"])]:
    ecdf(ax, FG[FG.recur==r].lead_size.clip(lower=1), f"FG {r}", cc)
ax.axvline(50,ls=":",c="0.5",lw=1); ax.set_xscale("log")
ax.annotate("50 bp\nSV / indel",xy=(50,0.05),fontsize=7,color="0.4",ha="left")
ax.set_xlabel("lead non-SNP |alt-ref| size (bp)",fontsize=9); ax.set_ylabel("ECDF",fontsize=9)
ax.legend(frameon=False,fontsize=8,loc="lower right")
plt.show()
print("lead is SV (>50bp): FG 1-axis %.0f%%  vs  >=2-axis %.0f%%"%(
    100*(FG[FG.recur=="1-axis"].lead_size>50).mean(),
    100*(FG[FG.recur==">=2-axis"].lead_size>50).mean()))""")

# ------------------------------------------------ Pillar C -------------------
md(r"""## Pillar C — signal-quality overlay (the calibration caveat, quantified)

These are **raw, uncalibrated** LFMM p, and they stay raw: the project does not GIF-correct GEA
p-values (lambda here is partly real climate signal correlated with structure, so dividing it out
removes signal with the artifact). Lambda is **stated** as a caveat instead. **(left)** lead
non-SNP MAF — low-MAF leads are the fragile ones; **(right)** per-axis count of new-peak blocks
against that axis' non-SNP lambda — if the axes with the most peaks are the most inflated ones,
absolute counts (not the class contrast) carry that caveat. Candidates are strengthened with
independent evidence downstream, not by a correction.""")

co(r"""fig,axs=plt.subplots(1,2,figsize=(12,4.2))
bg_lead=nsnp_b.loc[nsnp_b.groupby("block")["pval"].idxmin(),["block","MAF"]].set_index("block")
for r,cc in [("1-axis",COL["1-axis"]),(">=2-axis",COL[">=2-axis"])]:
    ecdf(axs[0], FG[FG.recur==r].lead_MAF, f"FG {r}", cc)
ecdf(axs[0], bg_lead.MAF.reindex(BG_blocks).dropna(),"background lead MAF",COL["background"],1.4)
axs[0].set_xlabel("lead non-SNP MAF",fontsize=9); axs[0].set_ylabel("ECDF",fontsize=9)
axs[0].legend(frameon=False,fontsize=7,loc="lower right")

# per-axis new-peak block count vs that axis' non-SNP lambda (lambda stated, never applied)
cnt=PK.groupby("axis").block.nunique().reindex(AXES).fillna(0)
lam=pd.Series(LAM).reindex(AXES)
axs[1].scatter(lam,cnt,s=22,color=COL[">=2-axis"],zorder=3)
for a in AXES:
    axs[1].annotate(a,(lam[a],cnt[a]),xytext=(3,2),textcoords="offset points",fontsize=6.5,color="0.3")
axs[1].axvline(1,ls=":",c="0.5",lw=1)
axs[1].set_xlabel("genomic inflation lambda (non-SNP, per axis)",fontsize=9)
axs[1].set_ylabel("# new-peak blocks on that axis",fontsize=9)
fig.tight_layout(); plt.show()

rho=stats.spearmanr(lam,cnt)
print(f"lead MAF median: FG {FG.lead_MAF.median():.3f}  vs  background lead {bg_lead.MAF.reindex(BG_blocks).median():.3f}")
print(f"non-SNP lambda per axis: min {lam.min():.2f}  median {lam.median():.2f}  max {lam.max():.2f} ({lam.idxmax()})")
print(f"Spearman(axis lambda, axis new-peak count) = {rho.statistic:.2f}  (p={rho.pvalue:.1e}, n={len(AXES)} axes)")""")

md(r"""## Summary table + master CSV

`nonsnp_block_characterization.csv` — one row per block: chrom, span, len, recurrence class,
best axis/p, lead size/MAF, SNP presence & density, distance-to-centromere, per-class coverage
fractions, and genic/TE/intergenic class. Top rows below by recurrence then signal.""")

co(r"""cols=["block","chrom","start","end","len","recur","n_axes","best_axis","best_nlp",
      "lead_size","lead_MAF","snp_present","n_snp","n_nonsnp","dist_cen",
      "f_gene","f_te_any","f_pseudogene","f_intergenic","klass"]
# FGc already carries every column (block stats + coverage); just select + sort.
show=FGc[[c for c in cols if c in FGc.columns]].sort_values(
    ["recur","best_nlp"],ascending=[True,False])
display(show.head(40))
print("full table:", f"{OUTDIR}/nonsnp_block_characterization.csv")""")

md(r"""## Conclusions — what are these blocks?

A **structurally distinct region class**; the raw-p caveats below say how far to trust it as
a list of adaptive loci. Numbers are the printed outputs above (rebuilt on `tile_cr50`).

1. **The SNP-vs-non-SNP contrast is not SNP-absence.** No foreground block lacks SNPs — every
   one carries SNP records that are simply *not significant* while a non-SNP is. The blocks are
   SNP-**richer** and physically **larger** than background (Pillar B prints), so the "SNPs
   missed it" caveat from the gene-level view does not bite at block level. Caveat: a larger,
   denser block holds more tests, which by itself raises the chance that some lead clears
   Bonferroni — part of the size/density shift is expected from that alone.
2. **TE-enriched, pseudogene-enriched, gene-depleted** — TE overlap roughly 2x background with
   Fisher OR ~3, pseudogene OR >2, and a smaller genic share than background (Pillar A tables).
   The TE-bp mix is shifted toward LTR/Gypsy relative to background (superfamily shares printed
   under the Pillar A bar chart; not formally tested).
3. **Pericentromere-leaning** — median distance to the centromere ~1.3 Mb shorter than
   background (MWU, Pillar B).
4. **A small-indel story, not SV** — only a few percent of leads are >50 bp, in both recurrence
   classes (Pillar B).
5. **Signal-quality caveats (stated, not corrected).** Leads are lower-MAF than background leads,
   and the per-axis new-peak count tracks that axis' lambda closely (Spearman printed in Pillar
   C): the most inflated axes (bio15, bio10, bio1) contribute most blocks. That is what inflation
   alone would produce, and it is also what climate signal correlated with structure would
   produce — so it qualifies absolute counts, not the class contrast.

**Read.** The non-SNP-only peaks concentrate in the compartment where indels are dense *and*
LFMM inflation is highest (TE-rich, pericentromeric, low-MAF leads, high-lambda axes). Pillars
A-B say *where* such signals arise. Per the project rule the p-values stay raw — this is a
candidate net, not a calibrated test — so individual blocks need independent evidence
(garden-level aggregation, a matched empirical background, locus context, cross-scan
recurrence) before being read as adaptive loci.""")

nb = new_notebook(); nb["cells"] = C
ep = ExecutePreprocessor(timeout=5400, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
