#!/usr/bin/env python
"""Per-site parallelism ECDF-difference across the climate gradient (user 2026-07-03), analogous to
s_ecdf_difference_by_site. For each site, the ECDF-difference of the per-site drift-controlled
parallelism rho (mean^2/mean(slope^2) across the ~12 plots) — SV vs frequency-matched SNP. A positive
hump = SVs shifted to higher rho (more parallel/selected) at that site. Panels ordered + colour-coded
cold->hot (bio1); summary vs temperature and precipitation. Load-only (parallelism.npz); `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/notebooks/parallelism_by_site.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = """# Per-site parallelism vs climate — are SVs more parallel at hot / arid sites?

For each site, ρ = drift-controlled parallelism (mean²/mean(slope²) across the ~12 replicate plots)
per variant. We plot the **ECDF-difference of ρ, SV vs frequency-matched SNP**: a positive hump means
SVs are shifted to higher ρ (more parallel = more selected) than same-frequency SNPs at that site.
Panels ordered + coloured cold→hot (bio1); summary vs temperature (bio1) and dry-summer precip (bio18)."""

code = r"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
plt.rcParams.update({'figure.dpi':110,'font.size':8,'axes.linewidth':0.6})
G="/global/scratch/users/tbellg/kmate/results/grenenet_gea/sv_adaptive"
z=np.load(f"{G}/parallelism.npz"); rng=np.random.default_rng(0)
RSV=z["rhops_sv"]; RSN=z["rhops_snp"]; p0sv=z["p0_sv"]; p0sn=z["p0_snp"]
b1=z["site_bio1"]; b18=z["site_bio18"]; sites=z["used_sites"]
order=np.argsort(b1)                                  # cold -> hot
GRID=np.linspace(0,1,160)
def match_idx(tp0,sp0):
    e=np.quantile(np.concatenate([tp0,sp0]),np.linspace(0,1,26)); e[-1]+=1e-9
    tb=np.digitize(tp0,e[1:-1]); sb=np.digitize(sp0,e[1:-1]); idx=[]
    for k in range(25):
        pool=np.where(sb==k)[0]; nn=int(round((tb==k).mean()*15000))
        if nn and pool.size: idx.append(rng.choice(pool,nn,replace=True))
    return np.concatenate(idx)
MI=match_idx(p0sv,p0sn)                                # fixed matched-SNP index set (same each site)
cdf=lambda s:np.searchsorted(np.sort(s),GRID,side="right")/s.size

ncol,nrow=5,int(np.ceil(len(order)/5)); norm=mpl.colors.Normalize(b1.min(),b1.max()); cmap=mpl.cm.coolwarm
fig,axes=plt.subplots(nrow,ncol,figsize=(15,2.5*nrow),sharex=True,sharey=True); axes=axes.ravel()
for ax in axes[len(order):]: ax.axis("off")
exv=np.zeros(len(sites))                               # per-column SV mean-rho excess
for i,c in enumerate(order):
    ax=axes[i]; dd=cdf(RSN[MI,c])-cdf(RSV[:,c])       # pos = SV shifted to higher rho (more parallel)
    ax.axhline(0,color="k",lw=.6,ls=":")
    ax.fill_between(GRID,0,dd,color="#D55E00",alpha=.2,lw=0); ax.plot(GRID,dd,color="#D55E00",lw=1.6)
    exv[c]=RSV[:,c].mean()-RSN[MI,c].mean()
    col=cmap(norm(b1[c]))
    ax.set_title(f"#{i+1} site {int(sites[c])} bio1={b1[c]:.0f}°C KS={np.abs(dd).max():.02f}",fontsize=7,color=col,fontweight="bold")
    for sp in ax.spines.values(): sp.set_color(col); sp.set_linewidth(1.4)
    ax.tick_params(labelsize=6)
fig.supxlabel("parallelism ρ (right = more parallel across replicate plots)",y=0.005,fontsize=10)
fig.supylabel("ΔCDF: matched-SNP − SV   (positive hump = SVs MORE parallel at this site)",x=0.004,fontsize=9)
fig.suptitle("Per-site SV parallelism vs frequency-matched SNP, ordered cold→hot",fontsize=11,y=1.004)
fig.tight_layout(rect=[0.03,0.02,1,0.99]); fig.savefig(f"{G}/parallelism_by_site.png",dpi=130,bbox_inches="tight"); plt.show()

# ---- FIGURE A: excess-only (SV - matched-SNP) per site vs climate ----
figA,axA=plt.subplots(1,2,figsize=(12,4.4))
for a,(cv,lab) in zip(axA,[(b1,"bio1 mean annual temp"),(b18,"bio18 precip warmest qtr (low=arid)")]):
    a.axhline(0,color="k",lw=.6,ls=":"); a.scatter(cv,exv,c="#D55E00",s=42,zorder=3)
    m,c0=np.polyfit(cv,exv,1); xs=np.array([cv.min(),cv.max()]); a.plot(xs,c0+m*xs,color="#D55E00",ls="--")
    r=stats.spearmanr(cv,exv); a.set_xlabel(f"site {lab}"); a.set_ylabel("SV − matched-SNP mean ρ (per site)")
    a.set_title(f"{lab.split()[0]}: ρ={r.statistic:+.2f} p={r.pvalue:.3f}")
figA.suptitle("SV parallelism excess per site vs climate",y=1.02)
figA.tight_layout(); figA.savefig(f"{G}/parallelism_by_site_excess.png",dpi=130,bbox_inches="tight"); plt.show()

# ---- FIGURE B: both series (SV + matched-SNP baseline) + the gap ----
sv_mean=np.array([RSV[:,c].mean() for c in range(len(sites))])
sn_mean=np.array([RSN[MI,c].mean() for c in range(len(sites))])
fig2,ax2=plt.subplots(1,2,figsize=(12,4.6))
for a,(cv,lab) in zip(ax2,[(b1,"bio1 mean annual temp"),(b18,"bio18 precip warmest qtr (low=arid)")]):
    for x,lo,hi in zip(cv,sn_mean,sv_mean):
        a.plot([x,x],[lo,hi],color="#D55E00",lw=0.6,alpha=0.4,zorder=2)      # gap = SV excess
    a.scatter(cv,sn_mean,c="#888888",s=40,zorder=3,label="matched SNP (baseline)")
    a.scatter(cv,sv_mean,c="#D55E00",s=40,zorder=3,label="SV")
    for yy,col in [(sn_mean,"#888888"),(sv_mean,"#D55E00")]:
        m,c0=np.polyfit(cv,yy,1); xs=np.array([cv.min(),cv.max()]); a.plot(xs,c0+m*xs,color=col,ls="--",lw=1.2)
    rS=stats.spearmanr(cv,sv_mean); rN=stats.spearmanr(cv,sn_mean); rE=stats.spearmanr(cv,sv_mean-sn_mean)
    a.set_xlabel(f"site {lab}"); a.set_ylabel("mean parallelism ρ (per site)")
    a.set_title(f"{lab.split()[0]}: SV ρ={rS.statistic:+.2f} | SNP ρ={rN.statistic:+.2f} | excess ρ={rE.statistic:+.2f} (p={rE.pvalue:.3f})")
    a.legend(frameon=False,fontsize=8)
fig2.suptitle("Per-site parallelism ρ: SV vs matched-SNP baseline vs climate (gap = SV excess)",y=1.02)
fig2.tight_layout(); fig2.savefig(f"{G}/parallelism_by_site_summary.png",dpi=130,bbox_inches="tight"); plt.show()
print(f"SV mean-rho: corr bio1 {stats.spearmanr(b1,sv_mean).statistic:+.2f}; matched-SNP baseline corr bio1 {stats.spearmanr(b1,sn_mean).statistic:+.2f}")
print(f"SV mean-rho: corr bio18 {stats.spearmanr(b18,sv_mean).statistic:+.2f}; matched-SNP baseline corr bio18 {stats.spearmanr(b18,sn_mean).statistic:+.2f}")
print(f"EXCESS (SV-SNP): corr bio1 {stats.spearmanr(b1,sv_mean-sn_mean).statistic:+.2f} (p={stats.spearmanr(b1,sv_mean-sn_mean).pvalue:.3f}); "
      f"bio18 {stats.spearmanr(b18,sv_mean-sn_mean).statistic:+.2f} (p={stats.spearmanr(b18,sv_mean-sn_mean).pvalue:.3f})")
"""

md_take = """### Reading
- Each panel: **positive orange hump = SVs shifted to higher parallelism ρ than same-frequency SNPs**
  at that site (more consistently selected across the replicate plots). Ordered cold→hot.
- **Summary**: per-site SV parallelism excess vs temperature (bio1) and dry-summer precip (bio18) — a
  trend toward hot / arid would mirror the climate-slope β result; a flat trend would mean the
  parallelism enrichment is climate-agnostic (present everywhere).
- Same caveats: founder projection + insertion-polarity."""

nb=new_notebook(cells=[new_markdown_cell(md_title), new_code_cell(code), new_markdown_cell(md_take)])
ep=ExecutePreprocessor(timeout=900, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata":{"path":os.path.dirname(OUT)}})
with open(OUT,"w") as f: nbf.write(nb,f)
print(f"[built+executed] {OUT}")
