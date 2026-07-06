#!/usr/bin/env python
"""Build the new-panel SNP-vs-nonSNP notebook: raw-LFMM Manhattan + WZA Manhattan
(clq0.9 haploblocks), raw top hits, and gene-annotated class-specific blocks.

kMate env builds it; execute in the `basic` env (matplotlib + nbconvert).
"""
import os
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

C = []
C.append(new_markdown_cell(
    "# New-panel climate-GEA — SNP vs non-SNP (raw LFMM + WZA), clq0.9 haploblocks\n\n"
    "Raw **LFMM ridge K=16** (`calibrate='gif'`), env=bio1, on kMate gen9 pool Δp. "
    "**SNP** (1.99M records) vs **non-SNP = indel+SV pooled** (693k). Blocks = "
    "**clq0.9 haploblocks** (~82k genome-wide, median 188 bp), NOT the phase-1 LD "
    "blocks. Two views of each scan: the **raw per-record** LFMM, and the **WZA** "
    "(deg7-cap2000) block aggregation on the clq0.9 units."))

C.append(new_code_cell(
    "import sys, os\n"
    "PROJ='/global/scratch/users/tbellg/kmate'\n"
    "sys.path.insert(0, f'{PROJ}/analysis/grenenet_gea')            # lib.py\n"
    "sys.path.insert(0, f'{PROJ}/analysis/grenenet_gea/gea_newpanel')  # blocks_clq09\n"
    "import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl\n"
    "import scipy.stats as st, lib\n"
    "mpl.rcParams.update({'figure.dpi':110,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})\n"
    "GNP=f'{PROJ}/results/grenenet_gea/gea_newpanel'\n"
    "OUT=f'{GNP}/snp_vs_nonsnp'; os.makedirs(OUT, exist_ok=True)\n"
    "CHR=[f'Chr{i}' for i in range(1,6)]\n"
    "ACC={'snp':'#c1443c','nonsnp':'#2e7d5b'}; NAME={'snp':'SNP','nonsnp':'non-SNP (indel+SV)'}"))

C.append(new_markdown_cell("## Load per-record (raw LFMM) + per-block (WZA), build genome coords"))
C.append(new_code_cell(
    "rec, wza = {}, {}\n"
    "for cls in ['snp','nonsnp']:\n"
    "    d=pd.read_csv(f'{GNP}/lfmm_{cls}_gen9_bio1_clq09.csv')\n"
    "    d=d[np.isfinite(d.pval)].copy(); d['nlp']=-np.log10(d.pval.clip(1e-300))\n"
    "    rec[cls]=d\n"
    "    w=pd.read_csv(f'{GNP}/wza_{cls}_clq09.csv').rename(columns={'gene':'block','index':'block'})\n"
    "    w=w[np.isfinite(w.Z_pVal)].copy(); w['pval']=w.Z_pVal; w['nlp']=-np.log10(w.Z_pVal.clip(1e-300))\n"
    "    wza[cls]=w\n"
    "# genome offsets from max pos per chrom across everything\n"
    "chrlen={c:int(max(rec['snp'].loc[rec['snp'].chrom==c,'pos'].max(),\n"
    "                  rec['nonsnp'].loc[rec['nonsnp'].chrom==c,'pos'].max())) for c in CHR}\n"
    "off={}; run=0\n"
    "for c in CHR: off[c]=run; run+=chrlen[c]+int(2e6)\n"
    "ticks=[off[c]+chrlen[c]/2 for c in CHR]\n"
    "for D in list(rec.values())+list(wza.values()):\n"
    "    D['gx']=D.pos.astype(float)+D.chrom.map(off)\n"
    "for cls in ['snp','nonsnp']:\n"
    "    print(f'{cls}: {len(rec[cls]):,} records | {len(wza[cls]):,} WZA clq0.9 blocks '\n"
    "          f'| raw min p={rec[cls].pval.min():.1e} | WZA min p={wza[cls].Z_pVal.min():.1e}')"))

C.append(new_markdown_cell(
    "## (1) RAW LFMM Manhattan — per record, no aggregation\n"
    "SNP (top) vs non-SNP (bottom). p>0.1 thinned 20× for rendering; dashed = p=1e-5."))
C.append(new_code_cell(
    "cols={'Chr1':'#3b5b92','Chr2':'#8bb0d0','Chr3':'#3b5b92','Chr4':'#8bb0d0','Chr5':'#3b5b92'}\n"
    "def _bh(p):\n"
    "    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)\n"
    "    q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q,0,1)\n"
    "def manhattan(getD, title, fname, ptcol_by_chrom=True, cap=None):\n"
    "    fig,axes=plt.subplots(2,1,figsize=(13,6.5),sharex=True)\n"
    "    for ax,cls in zip(axes,['snp','nonsnp']):\n"
    "        d=getD(cls).copy()\n"
    "        # per-panel significance lines from THIS panel's p-values\n"
    "        pv=d.pval.to_numpy(float); pv=pv[np.isfinite(pv)]; nt=len(pv)\n"
    "        bonf=-np.log10(0.05/nt); pas=pv[_bh(pv)<0.05]\n"
    "        fdr=float(-np.log10(pas.max())) if pas.size else None\n"
    "        n_over=0\n"
    "        if cap is not None:\n"
    "            n_over=int((d.nlp>cap).sum()); d['nlp']=d.nlp.clip(upper=cap)\n"
    "        # seam-free thinning: UNIFORM random subsample (constant keep-prob across y,\n"
    "        # so no density discontinuity) + always keep the real hits (nlp>2)\n"
    "        base=d.sample(frac=min(1.0,250000/len(d)),random_state=0) if len(d)>250000 else d\n"
    "        t=pd.concat([base,d[d.nlp>2]]).drop_duplicates(subset=['gx','nlp'])\n"
    "        c=t.chrom.map(cols) if ptcol_by_chrom else ACC[cls]\n"
    "        ax.scatter(t.gx,t.nlp,s=3,c=c,rasterized=True,linewidths=0)\n"
    "        ax.axhline(bonf,ls='--',lw=1.0,color='#444',label=f'Bonferroni ({bonf:.1f})')\n"
    "        if fdr is not None: ax.axhline(fdr,ls=':',lw=1.3,color='#b8860b',label=f'FDR q<.05 ({fdr:.1f})')\n"
    "        else: ax.plot([],[],' ',label='FDR q<.05: none pass')\n"
    "        top=d.nlargest(1,'nlp').iloc[0]\n"
    "        ax.scatter([top.gx],[top.nlp],s=42,facecolors='none',edgecolors=ACC[cls],linewidths=1.6,zorder=5)\n"
    "        lab=f\"{top.chrom}:{int(top.pos):,} {top.block}\"\n"
    "        if cap is not None and n_over: lab+=f\"  (+{n_over} blocks capped at {cap}: Z_pVal underflow)\"\n"
    "        ax.annotate(lab,(top.gx,top.nlp),xytext=(6,-1),textcoords='offset points',fontsize=8,color=ACC[cls])\n"
    "        ax.set_ylabel('-log10 p' + (f' (capped {cap})' if cap else ''))\n"
    "        ax.text(0.995,0.90,NAME[cls],transform=ax.transAxes,ha='right',fontweight='bold',color=ACC[cls])\n"
    "        ytop=max(float(d.nlp.max()),bonf)*1.10\n"
    "        ax.set_ylim(-0.4, min(ytop,cap+0.6) if cap else ytop)\n"
    "        ax.legend(loc='upper left',fontsize=7.5,frameon=False,handlelength=1.6)\n"
    "    axes[1].set_xticks(ticks); axes[1].set_xticklabels(CHR); axes[1].set_xlabel('genome position')\n"
    "    fig.suptitle(title,y=0.98); fig.tight_layout()\n"
    "    fig.savefig(f'{OUT}/{fname}',dpi=150,bbox_inches='tight'); print('wrote',fname); plt.show()\n"
    "manhattan(lambda c: rec[c], 'RAW LFMM (K=16, GIF-calibrated) climate-GEA bio1 — per record',\n"
    "          'manhattan_raw_snp_vs_nonsnp.png')"))

C.append(new_markdown_cell(
    "## (2) WZA Manhattan — clq0.9 haploblock aggregation\n"
    "Same two scans, but each point is a **clq0.9 block** WZA `Z_pVal` (deg7-cap2000). "
    "This is the block-level test on the finer, more-trusted haploblocks."))
C.append(new_code_cell(
    "# WZA on the fine clq0.9 SNP blocks produces machine-zero (underflow) Z_pVal for a\n"
    "# handful of blocks -> cap display at 15 so both panels are readable/comparable.\n"
    "n_uf={c:int((wza[c].Z_pVal<=0).sum()) for c in ['snp','nonsnp']}\n"
    "print('WZA Z_pVal underflow (==0) blocks:', n_uf)\n"
    "manhattan(lambda c: wza[c], 'WZA (deg7-cap2000) on clq0.9 haploblocks — bio1',\n"
    "          'manhattan_wza_clq09_snp_vs_nonsnp.png', cap=15)"))

C.append(new_markdown_cell(
    "## (2b) Site-level PC1 raw LFMM Manhattan — the honest unit + composite axis\n"
    "The plot-level bio1 scan above has no genome-wide-significant hits (pseudoreplication: "
    "355 plots, 31 climate values). The **honest unit** is the 31 flower-weighted site means, "
    "one observation per climate value. Here the env is **PC1 of the 19 bioclim vars** (the "
    "best-calibrated composite axis, GIF≈1). p = deflate-only fair p; MAF = site-level MAF. "
    "This is where the non-SNP class produces real Bonferroni hits that SNPs miss."))
C.append(new_code_cell(
    "pcr={}\n"
    "for cls in ['snp','nonsnp']:\n"
    "    d=pd.read_csv(f'{GNP}/lfmm_{cls}_sitepc1_clq09.csv')\n"
    "    d=d[np.isfinite(d.pval)].copy(); d['nlp']=-np.log10(d.pval.clip(1e-300))\n"
    "    d['gx']=d.pos.astype(float)+d.chrom.map(off)\n"
    "    pcr[cls]=d\n"
    "# Bonferroni line per class (site-level, on the MAF>=0.05 set)\n"
    "for cls in ['snp','nonsnp']:\n"
    "    m=pcr[cls][pcr[cls].MAF>=0.05]; nb=len(m)\n"
    "    print(f'{cls}: PC1 site-level, MAF>=0.05 n={nb:,}  min p={m.pval.min():.1e}  '\n"
    "          f'Bonferroni-sig={(m.pval<0.05/nb).sum()}')\n"
    "manhattan(lambda c: pcr[c][pcr[c].MAF>=0.05], 'Site-level raw LFMM on PC1(bioclim) — honest 31-site unit',\n"
    "          'manhattan_sitepc1_snp_vs_nonsnp.png')"))

C.append(new_markdown_cell(
    "## (3) Raw per-record LFMM hits (no WZA) — annotated\n"
    "Top records straight from `lfmm_test`, with the overlapping TAIR10 gene (±2 kb)."))
C.append(new_code_cell(
    "genes=lib.load_genes()\n"
    "def annotate_records(t):\n"
    "    a=lib.annotate_svs(t.copy(), flank=2000, genes=genes)\n"
    "    t=t.copy(); t['gene_name']=a.gene_name.values; t['genes']=a.genes_all.values; return t\n"
    "def raw_top(cls,n=25):\n"
    "    t=rec[cls].nsmallest(n,'pval')[['chrom','pos','ref_len','alt_len','sv_size','MAF','block','pval']]\n"
    "    return annotate_records(t).reset_index(drop=True)\n"
    "raw_snp=raw_top('snp'); raw_non=raw_top('nonsnp')\n"
    "raw_snp.to_csv(f'{OUT}/raw_top_records_snp_annotated.csv',index=False)\n"
    "raw_non.to_csv(f'{OUT}/raw_top_records_nonsnp_annotated.csv',index=False)\n"
    "pd.set_option('display.width',220,'display.max_colwidth',38)\n"
    "print('records p<1e-5 — SNP %d | non-SNP %d'%((rec['snp'].pval<1e-5).sum(),(rec['nonsnp'].pval<1e-5).sum()))"))
C.append(new_markdown_cell("### Raw top SNP records")); C.append(new_code_cell("raw_snp"))
C.append(new_markdown_cell("### Raw top non-SNP records")); C.append(new_code_cell("raw_non"))

C.append(new_markdown_cell(
    "## (4) Class-specific clq0.9 blocks — do the scans hit the same blocks?\n"
    "Aggregate to clq0.9 blocks and find **SNP-specific** vs **non-SNP-specific** peaks, "
    "under BOTH the raw scan (block peak = min per-record p) and the WZA scan (block "
    "`Z_pVal`). A block is class-specific if it's a top-50 peak in one scan and rank "
    ">500 (or absent) in the other. Concordance = shared-block peak Spearman."))
C.append(new_code_cell(
    "def block_min_raw(cls):\n"
    "    g=rec[cls].groupby('block').agg(p=('pval','min'),n=('pval','size'),chrom=('chrom','first')).reset_index()\n"
    "    g['rank']=g.p.rank(method='min'); return g\n"
    "def block_wza(cls):\n"
    "    g=wza[cls][['block','Z_pVal','chrom','pos','SNPs']].rename(columns={'Z_pVal':'p','SNPs':'n'}).copy()\n"
    "    g['rank']=g.p.rank(method='min'); return g\n"
    "def specific(scoref, label):\n"
    "    bs=scoref('snp').rename(columns={'p':'p_snp','n':'n_snp','rank':'r_snp'})\n"
    "    bn=scoref('nonsnp').rename(columns={'p':'p_non','n':'n_non','rank':'r_non'})\n"
    "    m=bs[['block','p_snp','n_snp','r_snp','chrom']].merge(\n"
    "        bn[['block','p_non','n_non','r_non']],on='block',how='outer')\n"
    "    sh=m.dropna(subset=['p_snp','p_non'])\n"
    "    rho=st.spearmanr(-np.log10(sh.p_snp.clip(1e-300)),-np.log10(sh.p_non.clip(1e-300))).correlation\n"
    "    B=10**9; m['r_snp']=m.r_snp.fillna(B); m['r_non']=m.r_non.fillna(B)\n"
    "    snp_spec=m[(m.r_snp<=50)&(m.r_non>500)].sort_values('p_snp')\n"
    "    non_spec=m[(m.r_non<=50)&(m.r_snp>500)].sort_values('p_non')\n"
    "    print(f'[{label}] clq0.9 blocks snp={bs.block.nunique():,} non={bn.block.nunique():,} '\n"
    "          f'shared={len(sh):,} | peak Spearman rho={rho:.3f} | '\n"
    "          f'SNP-specific={len(snp_spec)} nonSNP-specific={len(non_spec)}')\n"
    "    return m, snp_spec, non_spec\n"
    "def annotate_blocks(spec, cls_lead):\n"
    "    rows=[]\n"
    "    for _,r in spec.iterrows():\n"
    "        sub=rec[cls_lead][rec[cls_lead].block==r.block]\n"
    "        if not len(sub):\n"
    "            continue\n"
    "        lead=sub.nsmallest(1,'pval').iloc[0]\n"
    "        rows.append(dict(block=r.block, chrom=lead.chrom, lead_pos=int(lead.pos),\n"
    "            sv_size=int(abs(lead.alt_len-lead.ref_len)), MAF=round(float(lead.MAF),3),\n"
    "            ref_len=int(lead.ref_len), p_snp=r.get('p_snp'), p_non=r.get('p_non')))\n"
    "    t=pd.DataFrame(rows)\n"
    "    if len(t):\n"
    "        a=lib.annotate_svs(t.rename(columns={'lead_pos':'pos'}),flank=2000,genes=genes)\n"
    "        t['gene_name']=a.gene_name.values; t['genes']=a.genes_all.values\n"
    "    return t\n"
    "mr,snp_spec_raw,non_spec_raw = specific(block_min_raw,'RAW')\n"
    "mw,snp_spec_wza,non_spec_wza = specific(block_wza,'WZA')\n"
    "# annotate the WZA-based class-specific blocks (block-level test) + save\n"
    "non_wza_ann=annotate_blocks(non_spec_wza,'nonsnp'); snp_wza_ann=annotate_blocks(snp_spec_wza,'snp')\n"
    "non_wza_ann.to_csv(f'{OUT}/wza_nonsnp_specific_blocks_annotated.csv',index=False)\n"
    "snp_wza_ann.to_csv(f'{OUT}/wza_snp_specific_blocks_annotated.csv',index=False)\n"
    "# also the raw-based ones\n"
    "annotate_blocks(non_spec_raw,'nonsnp').to_csv(f'{OUT}/raw_nonsnp_specific_blocks_annotated.csv',index=False)\n"
    "annotate_blocks(snp_spec_raw,'snp').to_csv(f'{OUT}/raw_snp_specific_blocks_annotated.csv',index=False)"))

C.append(new_markdown_cell("### non-SNP-specific clq0.9 blocks (WZA) — peaks SNPs miss"))
C.append(new_code_cell("non_wza_ann.head(20)"))
C.append(new_markdown_cell("### SNP-specific clq0.9 blocks (WZA) — peaks the non-SNP scan misses"))
C.append(new_code_cell("snp_wza_ann.head(20)"))

C.append(new_markdown_cell(
    "### Read-out\n"
    "- Two Manhattans (raw per-record + WZA on clq0.9) let you compare the scans at "
    "both resolutions; the finer clq0.9 blocks localize signal better than the phase-1 "
    "LD blocks.\n"
    "- The shared-block peak Spearman (printed above) quantifies overall SNP↔non-SNP "
    "concordance; the annotated tables name the genes under each class-specific block.\n"
    "- Non-SNP-specific, SV-led blocks are the candidates for 'signal SNPs miss'."))

nb=new_notebook(); nb['cells']=C
nb.metadata['kernelspec']={'name':'python3','display_name':'Python 3','language':'python'}
os.makedirs('notebooks',exist_ok=True)
with open('notebooks/snp_vs_nonsnp_lfmm.ipynb','w') as f: nbformat.write(nb,f)
print('wrote notebooks/snp_vs_nonsnp_lfmm.ipynb')
