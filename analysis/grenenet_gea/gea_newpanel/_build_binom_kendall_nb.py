#!/usr/bin/env python
"""Build the new-panel SNP-vs-nonSNP notebook for the OTHER two phase-1 tests:
raw **binomial-regression** and raw **Kendall-tau** climate-GEA (bio1, gen9),
per-record (RAW frequencies, no WZA). SNP (1.99M) vs non-SNP=indel+SV (693k).

kMate env builds it; execute in the `basic` env (matplotlib + nbconvert).
"""
import os
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

C = []
C.append(new_markdown_cell(
    "# New-panel climate-GEA — SNP vs non-SNP: **binomial regression** + **Kendall-tau** (raw, no WZA)\n\n"
    "The other two phase-1 GEA tests (companions to the LFMM notebook), run on the "
    "new arch3 panel and split **SNP** (1.99M records) vs **non-SNP = indel+SV pooled** "
    "(693k). Both are **raw per-record** scans across the 355 gen9 `site_gen_plot` pools "
    "vs **bio1**, MAF≥0.05, flower-weighted Δp. **No WZA** — raw frequencies only.\n\n"
    "- **Binomial** (`run_binomial.py`): per-record GLM `[alt,ref] ~ const + z(bio1)`, "
    "counts = AF×flowers×2; slope + Wald p.\n"
    "- **Kendall** (`run_kendall.py`): per-record tau-b between per-pool AF and bio1.\n\n"
    "Pools within a site share a climate value (pseudoreplication) — same as phase-1; "
    "these raw scans do **not** correct structure (that is LFMM/WZA's job). So the "
    "comparison here is SNP-vs-nonSNP *relative* signal under identical, uncorrected tests."))

C.append(new_code_cell(
    "import sys, os\n"
    "PROJ='/global/scratch/users/tbellg/kmate'\n"
    "sys.path.insert(0, f'{PROJ}/analysis/grenenet_gea')            # lib.py\n"
    "import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl\n"
    "import scipy.stats as st, lib\n"
    "mpl.rcParams.update({'figure.dpi':110,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})\n"
    "GNP=f'{PROJ}/results/grenenet_gea/gea_newpanel'\n"
    "OUT=f'{GNP}/binom_kendall'; os.makedirs(OUT, exist_ok=True)\n"
    "CHR=[f'Chr{i}' for i in range(1,6)]\n"
    "ACC={'snp':'#c1443c','nonsnp':'#2e7d5b'}; NAME={'snp':'SNP','nonsnp':'non-SNP (indel+SV)'}\n"
    "TESTS={'binomial':'binomial','kendall':'kendall'}"))

C.append(new_markdown_cell("## Load per-record raw scans (binomial + kendall, both classes), build genome coords"))
C.append(new_code_cell(
    "dat={}   # dat[(test,cls)] -> DataFrame\n"
    "for test in ['binomial','kendall']:\n"
    "    for cls in ['snp','nonsnp']:\n"
    "        d=pd.read_csv(f'{GNP}/{test}/{test}_{cls}_gen9_bio1.csv')\n"
    "        d=d[np.isfinite(d.pval)].copy(); d['nlp']=-np.log10(d.pval.clip(1e-300))\n"
    "        dat[(test,cls)]=d\n"
    "# genome offsets from max pos per chrom across everything\n"
    "chrlen={}\n"
    "for c in CHR:\n"
    "    chrlen[c]=int(max(D.loc[D.chrom==c,'pos'].max() for D in dat.values()))\n"
    "off={}; run=0\n"
    "for c in CHR: off[c]=run; run+=chrlen[c]+int(2e6)\n"
    "ticks=[off[c]+chrlen[c]/2 for c in CHR]\n"
    "for D in dat.values(): D['gx']=D.pos.astype(float)+D.chrom.map(off)\n"
    "for test in ['binomial','kendall']:\n"
    "    for cls in ['snp','nonsnp']:\n"
    "        d=dat[(test,cls)]\n"
    "        print(f'{test:9s} {cls:7s}: {len(d):>9,} records tested | min p={d.pval.min():.1e} '\n"
    "              f'| p<1e-5: {(d.pval<1e-5).sum():>4d} | p<0.05: {(d.pval<0.05).sum():>8,}')"))

C.append(new_markdown_cell(
    "## Raw per-record Manhattans\n"
    "Each panel = one class; top row SNP, bottom non-SNP. Dashed = per-panel Bonferroni "
    "(0.05/n_tested); dotted gold = BH FDR q<0.05 line if any pass. Points thinned by "
    "uniform subsample (real hits nlp>2 always kept). These are raw, structure-uncorrected "
    "scans — read them as relative SNP-vs-nonSNP signal, not calibrated significance."))
C.append(new_code_cell(
    "cols={'Chr1':'#3b5b92','Chr2':'#8bb0d0','Chr3':'#3b5b92','Chr4':'#8bb0d0','Chr5':'#3b5b92'}\n"
    "def _bh(p):\n"
    "    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)\n"
    "    q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q,0,1)\n"
    "def manhattan(test, fname, pcol='pval', tag='RAW '):\n"
    "    fig,axes=plt.subplots(2,1,figsize=(13,6.5),sharex=True)\n"
    "    for ax,cls in zip(axes,['snp','nonsnp']):\n"
    "        d=dat[(test,cls)].copy()\n"
    "        d['_nlp']=-np.log10(d[pcol].clip(1e-300))\n"
    "        pv=d[pcol].to_numpy(float); pv=pv[np.isfinite(pv)]; nt=len(pv)\n"
    "        bonf=-np.log10(0.05/nt); pas=pv[_bh(pv)<0.05]\n"
    "        fdr=float(-np.log10(pas.max())) if pas.size else None\n"
    "        base=d.sample(frac=min(1.0,250000/len(d)),random_state=0) if len(d)>250000 else d\n"
    "        t=pd.concat([base,d[d._nlp>2]]).drop_duplicates(subset=['gx','_nlp'])\n"
    "        ax.scatter(t.gx,t._nlp,s=3,c=t.chrom.map(cols),rasterized=True,linewidths=0)\n"
    "        ax.axhline(bonf,ls='--',lw=1.0,color='#444',label=f'Bonferroni ({bonf:.1f})')\n"
    "        if fdr is not None: ax.axhline(fdr,ls=':',lw=1.3,color='#b8860b',label=f'FDR q<.05 ({fdr:.1f})')\n"
    "        else: ax.plot([],[],' ',label='FDR q<.05: none pass')\n"
    "        top=d.nlargest(1,'_nlp').iloc[0]\n"
    "        ax.scatter([top.gx],[top._nlp],s=42,facecolors='none',edgecolors=ACC[cls],linewidths=1.6,zorder=5)\n"
    "        ax.annotate(f'{top.chrom}:{int(top.pos):,} {top.block}',(top.gx,top._nlp),\n"
    "                    xytext=(6,-1),textcoords='offset points',fontsize=8,color=ACC[cls])\n"
    "        ax.set_ylabel('-log10 p'+(' (GIF-corr)' if pcol!='pval' else ''))\n"
    "        ax.text(0.995,0.90,NAME[cls],transform=ax.transAxes,ha='right',fontweight='bold',color=ACC[cls])\n"
    "        ax.set_ylim(-0.4, max(float(d._nlp.max()),bonf)*1.10)\n"
    "        ax.legend(loc='upper left',fontsize=7.5,frameon=False,handlelength=1.6)\n"
    "    axes[1].set_xticks(ticks); axes[1].set_xticklabels(CHR); axes[1].set_xlabel('genome position')\n"
    "    fig.suptitle(f'{tag}{test} climate-GEA (bio1, gen9) — SNP vs non-SNP',y=0.98); fig.tight_layout()\n"
    "    fig.savefig(f'{OUT}/{fname}',dpi=150,bbox_inches='tight'); plt.show()\n"
    "    print('wrote',fname)"))

C.append(new_markdown_cell("### Binomial-regression Manhattan"))
C.append(new_code_cell("manhattan('binomial','manhattan_binomial_snp_vs_nonsnp.png')"))
C.append(new_markdown_cell("### Kendall-tau Manhattan"))
C.append(new_code_cell("manhattan('kendall','manhattan_kendall_snp_vs_nonsnp.png')"))

C.append(new_markdown_cell(
    "## Raw top hits per test/class — gene-annotated\n"
    "Top records straight from each test, with the overlapping TAIR10 gene (±2 kb)."))
C.append(new_code_cell(
    "genes=lib.load_genes()\n"
    "def annotate(t):\n"
    "    a=lib.annotate_svs(t.copy(), flank=2000, genes=genes)\n"
    "    t=t.copy(); t['gene_name']=a.gene_name.values; t['genes']=a.genes_all.values; return t\n"
    "def top(test,cls,n=20):\n"
    "    d=dat[(test,cls)]\n"
    "    stat='slope' if test=='binomial' else 'tau'\n"
    "    t=d.nsmallest(n,'pval')[['chrom','pos','ref_len','alt_len','MAF','block',stat,'pval']]\n"
    "    return annotate(t).reset_index(drop=True)\n"
    "for test in ['binomial','kendall']:\n"
    "    for cls in ['snp','nonsnp']:\n"
    "        top(test,cls).to_csv(f'{OUT}/top_{test}_{cls}_annotated.csv',index=False)\n"
    "pd.set_option('display.width',220,'display.max_colwidth',38)\n"
    "print('saved top-record tables to',OUT)"))

C.append(new_markdown_cell("### Binomial — top SNP records"))
C.append(new_code_cell("top('binomial','snp')"))
C.append(new_markdown_cell("### Binomial — top non-SNP records"))
C.append(new_code_cell("top('binomial','nonsnp')"))
C.append(new_markdown_cell("### Kendall — top SNP records"))
C.append(new_code_cell("top('kendall','snp')"))
C.append(new_markdown_cell("### Kendall — top non-SNP records"))
C.append(new_code_cell("top('kendall','nonsnp')"))

C.append(new_markdown_cell(
    "## Concordance — do the tests and classes agree?\n"
    "(a) **Within class, across tests**: Spearman of -log10 p between binomial and kendall "
    "(records shared by chrom:pos). (b) **Within test, across classes**: SNP-vs-nonSNP peak "
    "concordance at the phase-1 LD-block level (block peak = min per-record p), matching the "
    "LFMM notebook's block-Spearman. (c) Class-specific top-p records that the other class misses."))
C.append(new_code_cell(
    "# unique record key: chrom,pos,ref_len,alt_len (pos alone collides on multi-allelic sites)\n"
    "def key(d): return (d.chrom.astype(str)+':'+d.pos.astype(str)+':'\n"
    "                    +d.ref_len.astype(str)+':'+d.alt_len.astype(str))\n"
    "# (a) within-class across-test Spearman on -log10 p\n"
    "print('=== (a) within-class, binomial vs kendall -log10p Spearman ===')\n"
    "for cls in ['snp','nonsnp']:\n"
    "    b=dat[('binomial',cls)].assign(k=lambda x:key(x)).drop_duplicates('k')[['k','nlp']].rename(columns={'nlp':'nlp_b'})\n"
    "    kd=dat[('kendall',cls)].assign(k=lambda x:key(x)).drop_duplicates('k')[['k','nlp']].rename(columns={'nlp':'nlp_k'})\n"
    "    m=b.merge(kd,on='k')\n"
    "    rho=st.spearmanr(m.nlp_b,m.nlp_k).correlation\n"
    "    print(f'  {cls:7s}: n_shared={len(m):>9,}  Spearman rho={rho:.3f}')\n"
    "# (b) within-test across-class LD-block peak concordance\n"
    "print('\\n=== (b) within-test, SNP vs non-SNP LD-block peak (min-p) Spearman ===')\n"
    "def block_peak(test,cls):\n"
    "    d=dat[(test,cls)]; d=d[d.block!='']\n"
    "    return d.groupby('block').pval.min().rename('p')\n"
    "for test in ['binomial','kendall']:\n"
    "    ps=block_peak(test,'snp'); pn=block_peak(test,'nonsnp')\n"
    "    j=pd.concat([ps,pn],axis=1,join='inner',keys=['snp','non'])\n"
    "    rho=st.spearmanr(-np.log10(j['snp'].clip(1e-300)),-np.log10(j['non'].clip(1e-300))).correlation\n"
    "    print(f'  {test:9s}: shared LD-blocks={len(j):>6,}  peak Spearman rho={rho:.3f}')"))

C.append(new_markdown_cell(
    "### Class-specific top records (raw)\n"
    "Records in one class's top-p tail whose *nearest* record in the other class (same LD block, "
    "min p) is unremarkable — i.e. signal one class carries and the other misses."))
C.append(new_code_cell(
    "def class_specific(test, lead='nonsnp', other='snp', topn=40, other_thresh=1e-3):\n"
    "    d=dat[(test,lead)]\n"
    "    stat='slope' if test=='binomial' else 'tau'\n"
    "    cand=d.nsmallest(topn,'pval').copy()\n"
    "    other_blockmin=dat[(test,other)][dat[(test,other)].block!=''].groupby('block').pval.min()\n"
    "    cand['p_other_block']=cand.block.map(other_blockmin)\n"
    "    spec=cand[(cand.block!='')&((cand.p_other_block.isna())|(cand.p_other_block>other_thresh))]\n"
    "    spec=spec[['chrom','pos','ref_len','alt_len','MAF','block',stat,'pval','p_other_block']]\n"
    "    return annotate(spec).reset_index(drop=True)\n"
    "for test in ['binomial','kendall']:\n"
    "    cs=class_specific(test)\n"
    "    cs.to_csv(f'{OUT}/{test}_nonsnp_specific.csv',index=False)\n"
    "    print(f'{test}: {len(cs)} non-SNP-specific top records (block p_snp>1e-3 or absent)')"))
C.append(new_markdown_cell("#### Binomial — non-SNP-specific top records"))
C.append(new_code_cell("class_specific('binomial')"))
C.append(new_markdown_cell("#### Kendall — non-SNP-specific top records"))
C.append(new_code_cell("class_specific('kendall')"))

C.append(new_markdown_cell(
    "## Controlling the inflation — genomic control (GIF), the LFMM recalibration\n"
    "The raw binomial/Kendall p above are wildly inflated: the binomial GLM treats "
    "AF×flowers×2 (up to ~1350 genomes/pool) as that many independent draws (over-precision), "
    "and both tests soak up genome-wide structure and the pseudoreplication of 355 pools "
    "sharing only 31 climate values. **LFMM controls the analogous inflation with genomic "
    "control (`calibrate='gif'`)**: rescale the genome-wide χ² so the null is flat — "
    "λ = median(χ²_obs)/median(χ²_null,1df), then divide every statistic by λ and re-derive p. "
    "That recalibration is method-agnostic, so we apply the identical move here: recover a "
    "χ²₁ from each two-sided p (`chi2.isf(p,1)`), deflate by λ, back to `pval_gc`. "
    "**Caveat:** unlike LFMM (where the GIF only mops up residual inflation *after* K latent "
    "factors), here one global λ is absorbing structure + over-precision + pseudoreplication "
    "at once — it flattens the bulk but under-corrects a heavy tail, and it does NOT remove "
    "pseudoreplication at its source (that needs the 31-site honest unit, as in the LFMM notebook)."))
C.append(new_code_cell(
    "from scipy.stats import chi2\n"
    "NULLMED=chi2.ppf(0.5,1)   # 0.4549 = median of a 1-df chi-square (the GIF null)\n"
    "def gif(d):\n"
    "    p=d.pval.to_numpy(float).clip(1e-300,1.0)\n"
    "    x2=chi2.isf(p,1)                    # observed 1-df chi-square from the two-sided p\n"
    "    lam=np.median(x2)/NULLMED           # genomic inflation factor (LFMM's GIF)\n"
    "    return lam, chi2.sf(x2/lam,1)       # deflate by lambda -> recalibrated p\n"
    "print('=== genomic-control (GIF) recalibration — the LFMM calibrate=\"gif\" move ===')\n"
    "print(f'{\"scan\":16s} {\"GIF\":>6s} | {\"raw p<1e-5\":>10s} {\"GC p<1e-5\":>10s} | '\n"
    "      f'{\"raw FDR<.05\":>11s} {\"GC FDR<.05\":>11s}')\n"
    "for test in ['binomial','kendall']:\n"
    "    for cls in ['snp','nonsnp']:\n"
    "        d=dat[(test,cls)]; lam,p_gc=gif(d)\n"
    "        d['pval_gc']=p_gc\n"
    "        r5=int((d.pval<1e-5).sum()); g5=int((p_gc<1e-5).sum())\n"
    "        rf=int((_bh(d.pval.to_numpy(float))<0.05).sum()); gf=int((_bh(p_gc)<0.05).sum())\n"
    "        print(f'{test+\"/\"+cls:16s} {lam:6.1f} | {r5:>10,} {g5:>10,} | {rf:>11,} {gf:>11,}')"))

C.append(new_markdown_cell(
    "### QQ plots — raw vs GIF-corrected\n"
    "Grey = raw (line lifts far off the diagonal = inflation); coloured = after GIF deflation "
    "(bulk collapses onto the diagonal; residual tail excursions = candidate real signal). "
    "λ annotated per panel."))
C.append(new_code_cell(
    "fig,axes=plt.subplots(2,2,figsize=(11,9))\n"
    "combos=[(t,c) for t in ['binomial','kendall'] for c in ['snp','nonsnp']]\n"
    "for ax,(test,cls) in zip(axes.ravel(),combos):\n"
    "    d=dat[(test,cls)]; lam,_=gif(d)\n"
    "    for col,color,label in [('pval','#b0b0b0','raw'),('pval_gc',ACC[cls],'GIF-corrected')]:\n"
    "        p=np.sort(d[col].to_numpy(float)); n=len(p)\n"
    "        exp=-np.log10((np.arange(1,n+1)-0.5)/n)\n"
    "        obs=-np.log10(np.clip(p,1e-300,1))\n"
    "        idx=np.unique(np.r_[np.linspace(0,n-1,2000).astype(int),np.arange(max(0,n-500),n)])\n"
    "        ax.plot(exp[idx],obs[idx],'.',ms=3,color=color,label=label)\n"
    "    mx=float(max(exp.max(),obs.max()))\n"
    "    ax.plot([0,mx],[0,mx],'k--',lw=0.8)\n"
    "    ax.set_title(f'{test} — {NAME[cls]}  (GIF={lam:.1f})',fontsize=9)\n"
    "    ax.set_xlabel('expected -log10 p'); ax.set_ylabel('observed -log10 p')\n"
    "    ax.legend(fontsize=7.5,frameon=False)\n"
    "fig.suptitle('QQ: raw vs genomic-control (GIF) recalibration',y=0.997); fig.tight_layout()\n"
    "fig.savefig(f'{OUT}/qq_gif_binom_kendall.png',dpi=140,bbox_inches='tight'); plt.show()\n"
    "print('wrote qq_gif_binom_kendall.png')"))

C.append(new_markdown_cell("### GIF-corrected Manhattans\nSame scans after genomic control; the flat baseline is now honest, hits stand above it."))
C.append(new_code_cell(
    "manhattan('binomial','manhattan_binomial_gc_snp_vs_nonsnp.png',pcol='pval_gc',tag='GIF-corrected ')\n"
    "manhattan('kendall','manhattan_kendall_gc_snp_vs_nonsnp.png',pcol='pval_gc',tag='GIF-corrected ')"))

C.append(new_markdown_cell(
    "## Read-out\n"
    "- Two raw scans (binomial GLM + Kendall tau-b), each SNP vs non-SNP, on identical gen9 "
    "pools vs bio1 — no WZA. Raw p are massively inflated (binomial min p→0, 85% p<0.05).\n"
    "- **(a)** binomial↔kendall Spearman shows the two raw tests agree within each class (ρ≈0.71); "
    "**(b)** the LD-block peak Spearman shows SNP↔nonSNP concordance per test (ρ≈0.66–0.71), "
    "comparable to the LFMM notebook's clq0.9 block ρ≈0.55–0.66.\n"
    "- **Inflation control:** genomic control (GIF) — the same recalibration LFMM applies with "
    "`calibrate='gif'` — flattens the QQ bulk onto the diagonal and collapses the hit counts; "
    "the printed GIF/λ per scan quantifies the inflation, and the GC Manhattans/QQ show what "
    "survives. This is the honest across-method comparison; the residual pseudoreplication "
    "(31 climate values) still calls for the site-level honest unit used in the LFMM notebook.\n"
    "- The non-SNP-specific tables name the indel/SV-led records (and genes) the SNP scan misses "
    "— candidate 'signal SNPs miss', to cross-check against the LFMM/site-level PC1 hits "
    "(e.g. `Chr1_9323`) before any selection claim."))

nb = new_notebook(cells=C)
nb.metadata['kernelspec'] = {'name':'python3','display_name':'Python 3','language':'python'}
outp = '/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/notebooks/binom_kendall_snp_vs_nonsnp.ipynb'
os.makedirs(os.path.dirname(outp), exist_ok=True)
with open(outp,'w') as f:
    nbformat.write(nb, f)
print('wrote', outp)
