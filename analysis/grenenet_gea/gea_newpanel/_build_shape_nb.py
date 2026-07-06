#!/usr/bin/env python
"""Build the shape-GEA SNP-vs-SV notebook: does relaxing LFMM's single linear
term reveal non-linear climate responses (intermediate optimum / conditional
neutrality), and do SVs carry more of them than SNPs?

Consumes shape_gea.py outputs (shape_{cls}_{weight}.npz). kMate env builds it;
execute in `basic` (matplotlib + nbconvert).
"""
import os
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

C = []
C.append(new_markdown_cell(
    "# Response-function GEA — SNP vs SV, beyond the linear LFMM\n\n"
    "Standard LFMM tests **one** shape per locus: a monotone-linear climate/Δp "
    "relationship = **antagonistic pleiotropy / clinal** adaptation. It is blind to "
    "**intermediate optimum** (adaptive only mid-gradient → concave *hump*, linear "
    "term ≈0) and mis-specifies **conditional neutrality** (one-sided ramp). Here we "
    "basis-expand the environment (linear+quadratic, envs bio1 & bio12) at the honest "
    "**site level** (31 gardens), on the variance-stabilised **logit(p9)** scale, and "
    "calibrate every shape statistic with a **climate-permutation null** (no GIF "
    "fudge). Two weightings: `none` (OLS) and `flowers` (drift-aware ~ s-coef). "
    "Question: **are SVs enriched for curvature the linear scan misses, vs SNPs?**"))

C.append(new_code_cell(
    "import sys, os\n"
    "PROJ='/global/scratch/users/tbellg/kmate'\n"
    "sys.path.insert(0, f'{PROJ}/analysis/grenenet_gea')\n"
    "import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, lib\n"
    "mpl.rcParams.update({'figure.dpi':110,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})\n"
    "SG=f'{PROJ}/results/grenenet_gea/gea_newpanel/shape_gea'\n"
    "OUT=SG+'/fig'; os.makedirs(OUT, exist_ok=True)\n"
    "CLS=['snp','indel','sv']; WGT=['none','flowers']\n"
    "ACC={'snp':'#c1443c','indel':'#c9922e','sv':'#2e7d5b'}; NAME={'snp':'SNP','indel':'indel','sv':'SV'}\n"
    "def load(cls,w):\n"
    "    z=np.load(f'{SG}/shape_{cls}_{w}.npz', allow_pickle=True)\n"
    "    return pd.DataFrame({k:z[k] for k in z.files})\n"
    "D={(c,w):load(c,w) for c in CLS for w in WGT}\n"
    "for c in CLS: print(f'{c}: {len(D[(c,\"none\")]):,} loci')"))

C.append(new_markdown_cell(
    "## (1) Genome-wide shape rates — permutation p<0.05\n"
    "`lin` = clinal/AP (what LFMM sees). `quad` = curvature. **`quad_only`** = curved "
    "but NOT linear (p_lin>0.5) = the response the linear GEA is blind to. `concave` = "
    "hump (intermediate optimum) share among curved."))
C.append(new_code_cell(
    "rows=[]\n"
    "for w in WGT:\n"
    "  for c in CLS:\n"
    "    d=D[(c,w)]; sig_l=d.p_lin<0.05; sig_q=d.p_quad<0.05\n"
    "    qonly=sig_q & (d.p_lin>0.5)\n"
    "    conc=(d.qsign_bio1<0)\n"
    "    rows.append(dict(weight=w,cls=c,n=len(d),\n"
    "        lin_pct=100*sig_l.mean(), quad_pct=100*sig_q.mean(),\n"
    "        quad_only_pct=100*qonly.mean(),\n"
    "        concave_of_quad_pct=100*(sig_q&conc).sum()/max(sig_q.sum(),1)))\n"
    "tab=pd.DataFrame(rows); pd.set_option('display.float_format',lambda x:f'{x:.2f}')\n"
    "tab"))

C.append(new_markdown_cell(
    "## (2) Frequency-matched curvature enrichment — SV & indel vs SNP\n"
    "Curvature rate depends on the MAF/p0 spectrum, and SVs and SNPs differ there. "
    "So we **match on p0**: for each SV (indel) locus draw a SNP of nearest p0, and "
    "compare the `quad`-hit and `quad_only`-hit rates. Bootstrap over loci → 95% CI on "
    "the SV/SNP rate ratio. A ratio CI above 1 = SVs genuinely carry more curvature "
    "than frequency-matched SNPs."))
C.append(new_code_cell(
    "rng=np.random.default_rng(0)\n"
    "def matched_snp_idx(p0_target, p0_snp, k=1):\n"
    "    order=np.argsort(p0_snp); ps=p0_snp[order]\n"
    "    j=np.searchsorted(ps, p0_target)\n"
    "    j=np.clip(j,0,len(ps)-1); return order[j]\n"
    "def enrich(cls, w, nboot=500):\n"
    "    dv=D[(cls,w)]; ds=D[('snp',w)]\n"
    "    p0s=ds.p0.to_numpy()\n"
    "    mi=matched_snp_idx(dv.p0.to_numpy(), p0s)\n"
    "    def rates(mask_v, mask_s):\n"
    "        return mask_v.mean(), mask_s.mean()\n"
    "    out={}\n"
    "    for key,fn in [('quad',lambda d:(d.p_quad<0.05).to_numpy()),\n"
    "                   ('quad_only',lambda d:((d.p_quad<0.05)&(d.p_lin>0.5)).to_numpy())]:\n"
    "        mv=fn(dv); ms=fn(ds)[mi]\n"
    "        rv,rs=mv.mean(),ms.mean(); ratio=rv/max(rs,1e-9)\n"
    "        bs=[]\n"
    "        n=len(mv)\n"
    "        for _ in range(nboot):\n"
    "            b=rng.integers(0,n,n); bs.append(mv[b].mean()/max(ms[b].mean(),1e-9))\n"
    "        lo,hi=np.percentile(bs,[2.5,97.5])\n"
    "        out[key]=dict(rate_v=100*rv,rate_snp=100*rs,ratio=ratio,lo=lo,hi=hi)\n"
    "    return out\n"
    "erows=[]\n"
    "for w in WGT:\n"
    "  for cls in ['indel','sv']:\n"
    "    e=enrich(cls,w)\n"
    "    for key in ['quad','quad_only']:\n"
    "        r=e[key]; erows.append(dict(weight=w,cls=cls,stat=key,\n"
    "            rate_pct=r['rate_v'],matched_snp_pct=r['rate_snp'],\n"
    "            ratio=r['ratio'],ci=f\"[{r['lo']:.2f},{r['hi']:.2f}]\"))\n"
    "enr=pd.DataFrame(erows); enr"))

C.append(new_code_cell(
    "# plot the frequency-matched ratio + CI (SV & indel vs matched SNP), both weightings\n"
    "fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)\n"
    "for ax,key in zip(axes,['quad','quad_only']):\n"
    "    ys=[]; labs=[]; i=0\n"
    "    for w in WGT:\n"
    "      for cls in ['sv','indel']:\n"
    "        e=enrich(cls,w)[key]\n"
    "        ax.errorbar(i,e['ratio'],yerr=[[e['ratio']-e['lo']],[e['hi']-e['ratio']]],\n"
    "            fmt='o',ms=8,color=ACC[cls],capsize=4)\n"
    "        labs.append(f'{NAME[cls]}\\n{w}'); i+=1\n"
    "    ax.axhline(1,ls='--',lw=0.8,color='0.5')\n"
    "    ax.set_xticks(range(i)); ax.set_xticklabels(labs); ax.set_title(f'{key} rate ratio vs matched SNP')\n"
    "axes[0].set_ylabel('SV(indel) / matched-SNP rate ratio')\n"
    "fig.suptitle('Frequency-matched curvature enrichment (>1 = more than SNPs)',y=1.02)\n"
    "fig.tight_layout(); fig.savefig(f'{OUT}/curvature_enrichment.png',dpi=150,bbox_inches='tight'); plt.show()"))

C.append(new_markdown_cell(
    "## (3) Response-shape class mix — SNP vs indel vs SV (frequency-matched)\n"
    "Classify each locus: **none** (no shape), **linear/AP** (lin only), **hump/interm-"
    "opt** (quad sig, concave), **valley/disruptive** (quad sig, convex), **both**. "
    "SNP set is p0-matched to the SV set so the mix is comparable."))
C.append(new_code_cell(
    "def classify(d):\n"
    "    sl=d.p_lin<0.05; sq=d.p_quad<0.05; conc=d.qsign_bio1<0\n"
    "    cl=np.full(len(d),'none',object)\n"
    "    cl[sl&~sq]='linear/AP'\n"
    "    cl[sq&conc]='hump/interm-opt'\n"
    "    cl[sq&~conc]='valley/disrupt'\n"
    "    cl[sl&sq]='both'\n"
    "    return pd.Series(cl)\n"
    "ORD=['none','linear/AP','hump/interm-opt','valley/disrupt','both']\n"
    "COLC={'none':'#dddddd','linear/AP':'#c1443c','hump/interm-opt':'#2e7d5b','valley/disrupt':'#8bb0d0','both':'#6a3d9a'}\n"
    "fig,axes=plt.subplots(1,2,figsize=(11,4.2),sharey=True)\n"
    "for ax,w in zip(axes,WGT):\n"
    "    dsv=D[('sv',w)]; dsn=D[('snp',w)]\n"
    "    mi=matched_snp_idx(dsv.p0.to_numpy(), dsn.p0.to_numpy())\n"
    "    mats={'SV':classify(dsv),'SNP(matched)':classify(dsn.iloc[mi].reset_index(drop=True)),\n"
    "          'indel':classify(D[('indel',w)])}\n"
    "    labels=list(mats); bottom=np.zeros(len(labels))\n"
    "    for k in ORD:\n"
    "        vals=np.array([ (m==k).mean()*100 for m in mats.values()])\n"
    "        ax.bar(labels,vals,bottom=bottom,color=COLC[k],label=k); bottom+=vals\n"
    "    ax.set_title(f'weight={w}'); ax.set_ylabel('% of loci')\n"
    "axes[1].legend(bbox_to_anchor=(1.02,1),loc='upper left',fontsize=8,frameon=False)\n"
    "fig.suptitle('Response-shape class mix (freq-matched SNP)',y=1.02)\n"
    "fig.tight_layout(); fig.savefig(f'{OUT}/shape_class_mix.png',dpi=150,bbox_inches='tight'); plt.show()"))

C.append(new_markdown_cell(
    "## (4) Curvature Manhattan — where the intermediate-optimum hits land\n"
    "F_quad per locus, SNP (top) vs SV (bottom). Highlighted = **quad_only** (curved, "
    "linear-blind). These are candidate loci a standard LFMM would never flag."))
C.append(new_code_cell(
    "CHR=[f'Chr{i}' for i in range(1,6)]\n"
    "off={}; run=0\n"
    "allpos={c:max(D[(cl,'none')].loc[D[(cl,'none')].chrom==c,'pos'].max() for cl in CLS) for c in CHR}\n"
    "for c in CHR: off[c]=run; run+=int(allpos[c])+int(2e6)\n"
    "ticks=[off[c]+allpos[c]/2 for c in CHR]\n"
    "cols={'Chr1':'#3b5b92','Chr2':'#8bb0d0','Chr3':'#3b5b92','Chr4':'#8bb0d0','Chr5':'#3b5b92'}\n"
    "w='none'\n"
    "fig,axes=plt.subplots(2,1,figsize=(13,6),sharex=True)\n"
    "for ax,cls in zip(axes,['snp','sv']):\n"
    "    d=D[(cls,w)].copy(); d['gx']=d.pos.astype(float)+d.chrom.map(off)\n"
    "    d['nlq']=-np.log10(d.p_quad.clip(1e-4))\n"
    "    qonly=(d.p_quad<0.05)&(d.p_lin>0.5)\n"
    "    base=d[~qonly].sample(frac=min(1.0,200000/len(d)),random_state=0)\n"
    "    ax.scatter(base.gx,base.nlq,s=3,c=base.chrom.map(cols),rasterized=True,linewidths=0)\n"
    "    ax.scatter(d.gx[qonly],d.nlq[qonly],s=14,color=ACC[cls],linewidths=0,zorder=5)\n"
    "    ax.set_ylabel('-log10 perm-p (quad)')\n"
    "    ax.text(0.995,0.9,f'{NAME[cls]}  (quad_only n={int(qonly.sum()):,})',transform=ax.transAxes,\n"
    "            ha='right',fontweight='bold',color=ACC[cls])\n"
    "axes[1].set_xticks(ticks); axes[1].set_xticklabels(CHR); axes[1].set_xlabel('genome position')\n"
    "fig.suptitle('Curvature (intermediate-optimum) scan — highlighted = linear-blind',y=0.98)\n"
    "fig.tight_layout(); fig.savefig(f'{OUT}/curvature_manhattan.png',dpi=150,bbox_inches='tight'); plt.show()"))

C.append(new_markdown_cell(
    "## (5) Top linear-blind curved SV loci — annotated\n"
    "SVs with strong curvature (low perm-p_quad) but NO linear signal (p_lin>0.5): the "
    "SV response functions a standard climate-GEA cannot see. Gene = TAIR10 ±2 kb."))
C.append(new_code_cell(
    "genes=lib.load_genes(); w='none'\n"
    "d=D[('sv',w)].copy()\n"
    "cand=d[(d.p_lin>0.5)].nsmallest(30,'p_quad').copy()\n"
    "cand['shape']=np.where(cand.qsign_bio1<0,'hump(interm-opt)','valley(disrupt)')\n"
    "cand['ref_len']=1; cand['alt_len']=1+cand.sv_size  # for annotate_svs flank logic\n"
    "a=lib.annotate_svs(cand.rename(columns={}).copy(),flank=2000,genes=genes)\n"
    "cand['gene_name']=a.gene_name.values; cand['genes']=a.genes_all.values\n"
    "keep=['chrom','pos','sv_size','maf','p0','p_lin','p_quad','p_full','shape','gene_name','genes']\n"
    "cand=cand[keep].reset_index(drop=True)\n"
    "cand.to_csv(f'{SG}/fig/top_curved_sv_linear_blind_annotated.csv',index=False)\n"
    "pd.set_option('display.width',220,'display.max_colwidth',36); cand.head(25)"))

C.append(new_markdown_cell(
    "## (6) Weighting agreement — OLS vs drift-aware (flowers)\n"
    "Do the curved-SV calls survive drift-weighting? Overlap of `quad`-hit SV sets "
    "between `none` and `flowers`. Robust hits appear in both."))
C.append(new_code_cell(
    "for cls in CLS:\n"
    "    a=(D[(cls,'none')].p_quad<0.05).to_numpy(); b=(D[(cls,'flowers')].p_quad<0.05).to_numpy()\n"
    "    jac=(a&b).sum()/max((a|b).sum(),1)\n"
    "    print(f'{cls}: quad-hit none={a.sum():,} flowers={b.sum():,} both={(a&b).sum():,} Jaccard={jac:.3f}')"))

C.append(new_markdown_cell(
    "### Read-out\n"
    "- **(1)** quad-hit rate near ~5% genome-wide = permutation null is calibrated; the "
    "signal is in the SNP-vs-SV **difference**, not the absolute rate.\n"
    "- **(2)** is the headline: freq-matched SV/SNP curvature ratio with CI. Ratio CI >1 "
    "= SVs carry non-linear climate responses beyond what MAF explains — and `quad_only` "
    "isolates the part a linear LFMM is blind to.\n"
    "- **(3)** the shape-class mix shows whether SVs skew toward hump (intermediate "
    "optimum) vs the SNP linear/AP mode.\n"
    "- **(5)** names the candidate linear-blind SV loci; **(6)** checks they survive "
    "drift-weighting.\n"
    "- Caveat: 31 sites → low per-locus power; this is a **class-distribution** claim "
    "(are SVs *enriched* for curvature), not per-locus discovery."))

nb = new_notebook(); nb['cells'] = C
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
os.makedirs('notebooks', exist_ok=True)
with open('notebooks/shape_gea_snp_vs_sv.ipynb', 'w') as f:
    nbformat.write(nb, f)
print('wrote notebooks/shape_gea_snp_vs_sv.ipynb')
