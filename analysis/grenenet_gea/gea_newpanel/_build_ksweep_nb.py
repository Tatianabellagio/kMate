#!/usr/bin/env python
"""Notebook: Manhattan + QQ of the RAW LFMM (no GIF) at a given K, from the K-sweep
per-K p-values (ksweep_calib/raw_p/{cls}_{axis}_K{K}.f32, full-record float32).
Classes SNP/non-SNP/SV; axes bio5, pc1, bio12/13/16/19. Manhattans mark per-panel
Bonferroni + FDR; QQ grid annotates λ. Shows the slightly-inflated regime at high K.

Usage:  _build_ksweep_nb.py <K>   ->  notebooks/ksweep_K{K}_manhattan_qq.ipynb
"""
import os, sys
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

K = sys.argv[1] if len(sys.argv) > 1 else "4"
NB = f"ksweep_K{K}_manhattan_qq.ipynb"

C = []
C.append(new_markdown_cell(
    f"# RAW LFMM at **K={K}** — Manhattan + QQ (SNP / non-SNP / SV)\n\n"
    f"Uncalibrated LFMM (no GIF) at K={K} latent factors, honest 31-site unit, "
    f"site-MAF≥0.05. Axes: bio5, pc1, bio12/13/16/19. At higher K the temperature "
    f"axes move into the slightly-inflated λ≈1.1–1.4 range; the QQ λ and the "
    f"Bonferroni/FDR lines let you read whether hits appear where λ is well-behaved."))

C.append(new_code_cell(
    "import sys, os\n"
    "PROJ='/global/scratch/users/tbellg/kmate'\n"
    "sys.path.insert(0, f'{PROJ}/analysis/grenenet_gea')\n"
    "sys.path.insert(0, f'{PROJ}/analysis/grenenet_gea/gea_newpanel')\n"
    "import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl\n"
    "import scipy.stats as st, lib, blocks_clq09\n"
    "mpl.rcParams.update({'figure.dpi':110,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})\n"
    f"K='{K}'\n"
    "RAWP=f'{PROJ}/results/grenenet_gea/gea_newpanel/ksweep_calib/raw_p'\n"
    "SITE=f'{PROJ}/results/grenenet_gea/gea_newpanel/lfmm_site'\n"
    "CM=f'{PROJ}/results/grenenet_gea/phase1_replication/class_matrices'\n"
    "OUT=f'{PROJ}/results/grenenet_gea/gea_newpanel/snp_vs_nonsnp'; os.makedirs(OUT,exist_ok=True)\n"
    "CHR=[f'Chr{i}' for i in range(1,6)]\n"
    "CLS=['snp','nonsnp','sv']; NAME={'snp':'SNP','nonsnp':'non-SNP','sv':'SV'}\n"
    "ACC={'snp':'#c1443c','nonsnp':'#2e7d5b','sv':'#7048a8'}\n"
    "AXES=['bio5','pc1','bio12','bio13','bio16','bio19']\n"
    "_P0={'snp':'p0_snp.npy','nonsnp':'p0_nonsnp.npy','sv':'p0_nonsnp.npy'}\n"
    "def bh(p):\n"
    "    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)\n"
    "    q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q,0,1)\n"
    "def gif(p): p=np.clip(np.asarray(p,float),1e-300,1); return np.median(st.chi2.isf(p,1))/st.chi2.isf(0.5,1)"))

C.append(new_code_cell(
    "recs={}; mask={}\n"
    "for c in CLS:\n"
    "    r=pd.read_csv(f'{CM}/{c}_gen9.records.csv')\n"
    "    dims=np.loadtxt(f'{SITE}/lfmm_{c}_site_dims.txt',dtype=int); ns,nv=int(dims[0]),int(dims[1])\n"
    "    dp=np.fromfile(f'{SITE}/lfmm_{c}_site_Y.f64').reshape(ns,nv)\n"
    "    p0=np.load(f'{lib.AF_STORE}/{_P0[c]}')[r[\"col\"].to_numpy()]\n"
    "    mf=(dp+p0[None,:]).mean(0); sm=np.minimum(mf,1-mf); m=sm>=0.05; mask[c]=m\n"
    "    rr=r[m].copy(); rr['site_maf']=sm[m]\n"
    "    rr['block']=blocks_clq09.assign_clq09_blocks(rr.chrom.to_numpy(),rr.pos.to_numpy())\n"
    "    recs[c]=rr.reset_index(drop=True)\n"
    "chrlen={ch:int(max(recs[c].loc[recs[c].chrom==ch,'pos'].max() if (recs[c].chrom==ch).any() else 0 for c in CLS)) for ch in CHR}\n"
    "off={}; run=0\n"
    "for ch in CHR: off[ch]=run; run+=chrlen[ch]+int(2e6)\n"
    "ticks=[off[ch]+chrlen[ch]/2 for ch in CHR]\n"
    "for c in CLS: recs[c]['gx']=recs[c].pos.astype(float)+recs[c].chrom.map(off)\n"
    "def load_p(cls,ax):\n"
    "    p=np.fromfile(f'{RAWP}/{cls}_{ax}_K{K}.f32', dtype=np.float32)\n"
    "    return p[mask[cls]]\n"
    "print('loaded records + K=%s p-values'%K)"))

C.append(new_markdown_cell(
    f"## Peaks per bio var at K={K} — SNP / non-SNP / SV significant clq0.9 blocks (+ SNP∩SV)\n"
    "Computed from this notebook's own RAW p (site-MAF≥0.05, clq0.9 blocks), so it "
    "matches the Manhattans. `SNP∩SV` = blocks significant in BOTH the SNP and SV scans."))
C.append(new_code_cell(
    "def sig_blocks(cls, axname):\n"
    "    d=recs[cls].copy(); d['pval']=load_p(cls,axname); d=d[np.isfinite(d.pval)]\n"
    "    n=len(d); q=bh(d.pval.values)\n"
    "    bonf=set(d.loc[d.pval < 0.05/n,'block']); fdr=set(d.loc[q < 0.05,'block'])\n"
    "    bonf.discard(''); fdr.discard(''); return bonf, fdr\n"
    "SB={ax:{c:sig_blocks(c,ax) for c in CLS} for ax in AXES}\n"
    "rows=[]\n"
    "for ax in AXES:\n"
    "    r={'axis':ax}\n"
    "    for c in CLS:\n"
    "        b,f=SB[ax][c]; r[f'{NAME[c]}_Bonf']=len(b); r[f'{NAME[c]}_FDR']=len(f)\n"
    "    s,v=SB[ax]['snp'],SB[ax]['sv']\n"
    "    r['SNP∩SV_Bonf']=len(s[0]&v[0]); r['SNP∩SV_FDR']=len(s[1]&v[1])\n"
    "    rows.append(r)\n"
    "TBL=pd.DataFrame(rows).set_index('axis')\n"
    "TBL.to_csv(f'{OUT}/sig_blocks_ksweep_K{K}.csv')\n"
    "print(f'Peaks (significant clq0.9 blocks) at K={K}, RAW LFMM, site-MAF>=0.05 — matches Manhattans')\n"
    "TBL"))

C.append(new_markdown_cell(f"## QQ grid at K={K} — calibration per axis × class (λ annotated)"))
C.append(new_code_cell(
    "def qq(ax_, p, color, lab):\n"
    "    p=np.sort(np.asarray(p,float)); p=p[np.isfinite(p)]; n=len(p)\n"
    "    obs=-np.log10(np.clip(p,1e-300,1)); exp=-np.log10((np.arange(1,n+1)-0.5)/n)\n"
    "    keep=np.concatenate([np.where(obs>2)[0], np.random.RandomState(0).choice(np.where(obs<=2)[0], size=min(4000,int((obs<=2).sum())), replace=False)])\n"
    "    keep=np.sort(keep); ax_.scatter(exp[keep],obs[keep],s=3,c=color,rasterized=True,linewidths=0)\n"
    "    m=max(exp.max(),obs.max()); ax_.plot([0,m],[0,m],lw=0.7,color='0.5')\n"
    "    ax_.text(0.05,0.86,f'{lab}\\n$\\\\lambda$={gif(p):.2f}',transform=ax_.transAxes,fontsize=8,color=color)\n"
    "fig,axes=plt.subplots(len(AXES),3,figsize=(9,2.4*len(AXES)))\n"
    "for i,axname in enumerate(AXES):\n"
    "    for j,c in enumerate(CLS):\n"
    "        qq(axes[i,j], load_p(c,axname), ACC[c], f'{axname} {NAME[c]}')\n"
    "        if j==0: axes[i,j].set_ylabel('obs -log10 p')\n"
    "        if i==len(AXES)-1: axes[i,j].set_xlabel('exp -log10 p')\n"
    "fig.suptitle(f'QQ: RAW LFMM K={K} per axis x class',y=1.0); fig.tight_layout()\n"
    "fig.savefig(f'{OUT}/qq_ksweep_K{K}.png',dpi=140,bbox_inches='tight'); print('wrote qq'); plt.show()"))

C.append(new_markdown_cell(f"## Manhattans at K={K} — SNP / non-SNP / SV (Bonferroni + FDR marked)"))
C.append(new_code_cell(
    "cols={'Chr1':'#3b5b92','Chr2':'#8bb0d0','Chr3':'#3b5b92','Chr4':'#8bb0d0','Chr5':'#3b5b92'}\n"
    "def manhattan_axis(axname):\n"
    "    fig,axs=plt.subplots(3,1,figsize=(13,7.5),sharex=True)\n"
    "    for ax_,c in zip(axs,CLS):\n"
    "        d=recs[c].copy(); d['pval']=load_p(c,axname); d=d[np.isfinite(d.pval)]\n"
    "        d['nlp']=-np.log10(d.pval.clip(1e-300)); pv=d.pval.to_numpy(); nt=len(pv)\n"
    "        bonf=-np.log10(0.05/nt); pas=pv[bh(pv)<0.05]; fdr=(-np.log10(pas.max())) if pas.size else None\n"
    "        base=d.sample(frac=min(1.0,200000/len(d)),random_state=0) if len(d)>200000 else d\n"
    "        t=pd.concat([base,d[d.nlp>2]]).drop_duplicates(subset=['gx','nlp'])\n"
    "        ax_.scatter(t.gx,t.nlp,s=3,c=t.chrom.map(cols),rasterized=True,linewidths=0)\n"
    "        ax_.axhline(bonf,ls='--',lw=1.0,color='#444',label=f'Bonferroni ({bonf:.1f})')\n"
    "        if fdr is not None: ax_.axhline(fdr,ls=':',lw=1.3,color='#b8860b',label=f'FDR q<.05 ({fdr:.1f})')\n"
    "        else: ax_.plot([],[],' ',label='FDR q<.05: none')\n"
    "        top=d.nlargest(1,'nlp').iloc[0]\n"
    "        ax_.scatter([top.gx],[top.nlp],s=40,facecolors='none',edgecolors=ACC[c],linewidths=1.5,zorder=5)\n"
    "        ax_.annotate(f'{top.chrom}:{int(top.pos):,} {top.block}',(top.gx,top.nlp),xytext=(6,-1),\n"
    "                     textcoords='offset points',fontsize=7.5,color=ACC[c])\n"
    "        ax_.set_ylabel('-log10 p'); ax_.set_ylim(-0.4,max(float(d.nlp.max()),bonf)*1.10)\n"
    "        ax_.text(0.995,0.86,NAME[c],transform=ax_.transAxes,ha='right',fontweight='bold',color=ACC[c])\n"
    "        ax_.legend(loc='upper left',fontsize=7,frameon=False)\n"
    "    axs[2].set_xticks(ticks); axs[2].set_xticklabels(CHR); axs[2].set_xlabel('genome position')\n"
    "    lam=' | '.join(f'{NAME[c]} λ={gif(load_p(c,axname)):.2f}' for c in CLS)\n"
    "    fig.suptitle(f'{axname} — RAW LFMM K={K}   [{lam}]',y=0.995); fig.tight_layout()\n"
    "    fn=f'manhattan_{axname}_K{K}.png'; fig.savefig(f'{OUT}/{fn}',dpi=140,bbox_inches='tight')\n"
    "    print('wrote',fn); plt.show()"))
for ax in ["bio5", "pc1", "bio12", "bio13", "bio16", "bio19"]:
    C.append(new_markdown_cell(f"### {ax} (K={K})"))
    C.append(new_code_cell(f"manhattan_axis('{ax}')"))

nb = new_notebook(); nb['cells'] = C
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
os.makedirs('notebooks', exist_ok=True)
with open(f'notebooks/{NB}', 'w') as f:
    nbformat.write(nb, f)
print('wrote notebooks/' + NB)
