#!/usr/bin/env python
"""Build+execute `sv_main_figures_cr{pct}.ipynb`: the twin of `sv_main_figures.ipynb` with a RELAXED
founder call-rate floor (SV_CALL_MIN, default 0.7; the original uses 0.9).

The twin is derived from the CURRENT `sv_main_figures.ipynb` (cells copied, not re-typed), so it follows
whatever that notebook plots; only the inputs change -- every file is swapped for its `_cr{pct}` version
and `passes` is recomputed from the event catalog's `n_called`. The floor is applied to every variant
class alike (SV events, SNPs, short indels), as in the original. Spectra are projected to
M = min(200, floor(CALL_MIN * 231 / 10) * 10) because M cannot exceed the smallest k kept.
A last section puts the headline numbers of the two floors side by side.

Inputs (all with SV_CALL_MIN set, SV_UNIT=event):
  _sv_classic_sfs.py, _sv_sfs_context_origin.py, _sv_burden_build.py -> _sv_burden_lmm.py,
  _sv_grenenet_sfs.py, run_temporal_callrate.sbatch (temporal -> classify -> annotate),
  _sv_main_figures_checks.py
`basic` kernel, as the original.
"""
import os, re
import nbformat as nbf
from nbformat.v4 import new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.7))
PCT = int(round(CALL_MIN * 100)); CR = f"_cr{PCT}"
KMIN = int(-(-CALL_MIN * 231 // 1)); M = min(200, int(CALL_MIN * 231 // 10 * 10))
HERE = os.path.dirname(os.path.abspath(__file__))
NBDIR = os.path.join(os.path.dirname(HERE), "notebooks")
SRC = os.path.join(NBDIR, "sv_main_figures.ipynb"); OUT = os.path.join(NBDIR, f"sv_main_figures{CR}.ipynb")

# input / output files that carry the floor
FILES = ["classic_sfs.npz", "main_fig_checks_load.csv", "main_fig_checks_lmm.csv", "main_fig_checks_nosingleton.csv",
         "sv_sfs_context_origin_event.npz", "sv_burden_lmm_event.csv", "sv_burden_lmm_joint_event.csv",
         "grenenet_sfs_event.npz", "founder_burden_event.csv", "sv_purging_annotated_event2.csv",
         "sv_temporal_climate_matched_sites_event2.csv", "sv_temporal_climate_matched_persite_event2.npz",
         "main_fig1.png", "main_fig1.pdf", "main_fig2.png", "main_fig2.pdf", "main_fig3.png", "main_fig3.pdf"]


def tag(fn):
    stem, ext = fn.rsplit(".", 1); return f"{stem}{CR}.{ext}"


def fix_code(s):
    for fn in FILES:
        if fn in s:
            s = s.replace(fn, tag(fn))
    a = "ok=D.passes.to_numpy()"
    if a in s:
        s = s.replace(a, f"CALL_MIN={CALL_MIN}; ok=((D.n_carriers>=2)&(D.n_called>=CALL_MIN*231)).to_numpy()")
    return s


def fix_md(s):
    s = s.replace("# Selection on structural variants — main figures",
                  f"# Selection on structural variants — main figures (call rate ≥ {PCT}%)")
    s = s.replace("≥ 90%", f"≥ {PCT}%").replace("k ≥ 208", f"k ≥ {KMIN}")
    s = s.replace("M = 200", f"M = {M}").replace("j = 1…199", f"j = 1…{M - 1}").replace("a sample of 200", f"a sample of {M}")
    # sentences that quote numbers of the 90% run: point to this run's printed output instead
    s = re.sub(r" A per-garden Wilcoxon would give[^.]*\.[^.]*\.", " The per-garden counts and the null p are printed under Figure 2.", s)
    s = re.sub(r"extra covariates: .*?0-fold SNPs 0\.12\)\.",
               "extra covariates, and Spearman burden–fitness within origin tertiles (annotated in Fig. 2a–b); spectra "
               "without singletons (minor-allele count ≥ 2) printed under Figure 1.", s, flags=re.S)
    return s


md_note = f"""> **What is different from `sv_main_figures.ipynb`.** Same cells, same code; the founder call-rate
> floor is **{PCT}%** instead of 90% (variant called in ≥ {KMIN} of 231 founders), for every variant class
> (SV events, SNPs, short indels). All inputs are the `{CR}` reruns of the same scripts (`SV_CALL_MIN={CALL_MIN}`).
> Spectra are projected to M = {M} instead of 200 (M cannot exceed the smallest k kept), so Tajima's *D*
> values are not on the same scale as the 90% notebook — the last section re-projects the 90% set to
> M = {M} for a like-for-like comparison. Figure titles and legends below are copied from the original and
> describe what is plotted, not what this run shows. Built by `_build_sv_main_figures_callrate_nb.py`."""

md_cmp = f"""## Call-rate floor: 90% against {PCT}%, side by side

How many events each floor keeps, who the added events are, and the headline statistic of every panel
under both floors. Tajima's *D* of the 90% set is re-projected to M = {M} here so the two columns are
comparable."""

code_cmp = r'''
from scipy.stats import hypergeom
CR="__CR__"; N=231; Mc=int(np.load(f"{G}/classic_sfs{CR}.npz")["M"])
D0=pd.read_csv(f"{G}/sv_event_catalog.csv"); D0["cr"]=D0.n_called/N; d2=D0[D0.n_carriers>=2]
FLO={"90%":0.9,"__PCT__%":__CALL__}
# 1. what the floor keeps
t=pd.DataFrame({k:d2[d2.cr>=v].groupby("kind").size() for k,v in FLO.items()}); t.loc["ins, CDS"]=[((d2.cr>=v)&(d2.kind=="ins")&(d2.context=="CDS")).sum() for v in FLO.values()]
t.loc["del, CDS"]=[((d2.cr>=v)&(d2.kind=="del")&(d2.context=="CDS")).sum() for v in FLO.values()]; t["of all ≥2-carrier events"]=[f"{100*x/y:.0f}%" for x,y in zip(t.iloc[:,1],[ (d2.kind=='del').sum(),(d2.kind=='ins').sum(),((d2.kind=='ins')&(d2.context=='CDS')).sum(),((d2.kind=='del')&(d2.context=='CDS')).sum()])]
print("events with ≥ 2 carriers kept"); print(t.to_string())
add=d2[(d2.cr>=__CALL__)&(d2.cr<0.9)]; old=d2[d2.cr>=0.9]
print("\nwho is added (call rate __PCT__–90%) against who was already in (≥ 90%):")
print(pd.DataFrame({lab:x.groupby("kind").agg(n=("cr","size"),median_bp=("size_bp","median"),median_carriers=("n_carriers","median"),le3=("n_carriers",lambda s:(s<=3).mean()),
      alleles_per_event=("n_alleles","mean"),in_CDS=("context",lambda s:(s=="CDS").mean()),in_TE=("context",lambda s:(s=="TE").mean()),shared=("origin",lambda s:(s=="shared").mean())).stack() for lab,x in (("≥ 90%",old),("added",add))}).round(3).to_string())
# 2. spectra, both floors on M = Mc
def proj(c,k):
    u,cnt=np.unique(np.c_[c,k],axis=0,return_counts=True); j=np.arange(Mc+1); s=np.zeros(Mc+1)
    for (ci,ki),n in zip(u,cnt): s+=n*hypergeom.pmf(j,ki,ci,Mc)
    return s[1:Mc]
def tD(sfs):
    n=Mc; j=np.arange(1,n); S_=sfs.sum(); pi=(sfs*2*j*(n-j)/(n*(n-1))).sum(); a1=(1/j).sum(); a2=(1/j**2).sum()
    b1=(n+1)/(3*(n-1)); b2=2*(n*n+n+3)/(9*n*(n-1)); c1=b1-1/a1; c2=b2-(n+2)/(a1*n)+a2/a1**2; e1=c1/a1; e2=c2/(a1**2+a2)
    return (pi-S_/a1)/np.sqrt(e1*S_+e2*S_*(S_-1))
rows=[]
for lab,fn in (("90%","classic_sfs.npz"),("__PCT__%",f"classic_sfs{CR}.npz")):
    Z=np.load(f"{G}/{fn}"); ref=None
    for n_ in ("snp_4fold","snp_0fold","ins_CDS","ins_UTR","ins_intron","ins_TE","ins_promoter","ins_intergenic","del_CDS","del_intron"):
        c=Z[f"{n_}_c"].astype(int); k=Z[f"{n_}_k"].astype(int); f=c/k; v=np.minimum(f,1-f); ref=v if n_=="snp_4fold" else ref
        rows.append(dict(floor=lab,cls=n_,n=c.size,maf_lt5=(v<0.05).mean(),KS_D=stats.ks_2samp(v,ref).statistic,TajD=tD(proj(c,k))))
S2=pd.DataFrame(rows).pivot(index="cls",columns="floor",values=["n","maf_lt5","KS_D","TajD"]).reindex([r["cls"] for r in rows[:10]])
print(f"\nFigure 1 — spectra (Tajima's D of both floors projected to M = {Mc}):"); print(S2.round(3).to_string())
# 3. in-frame share of CDS insertions
print("\nFigure 1d — in-frame share:")
for lab,v in FLO.items():
    x=d2[(d2.cr>=v)&(d2.kind=="ins")]; c_=x[x.context=="CDS"]; o_=x[x.context!="CDS"]
    print(f"  {lab:<4} CDS insertions {(c_.size_bp%3==0).mean():.3f} (n = {len(c_):,});  other contexts {(o_.size_bp%3==0).mean():.3f} (n = {len(o_):,})")
# 4. burden, load, temporal slopes
print("\nFigure 2a — shared-variant burden β₁ per SD (p), CDS:")
for lab,sfx in (("90%",""),("__PCT__%",CR)):
    R_=pd.read_csv(f"{G}/sv_burden_lmm_event{sfx}.csv"); J_=pd.read_csv(f"{G}/sv_burden_lmm_joint_event{sfx}.csv"); K_=pd.read_csv(f"{G}/main_fig_checks_lmm{sfx}.csv").set_index("check")
    r=R_[(R_.freq=="all")&(R_.origin=="shared")&(R_.context=="CDS")].set_index("kind"); j=J_[(J_.model=="ins + SNP + indel")&(J_.term=="ins_shared_CDS_all")].iloc[0]; b=K_.loc["+ founder origin bio1"]
    print(f"  {lab:<4} "+"  ".join(f"{KL[k]} {r.loc[k,'beta']:+.2f} ({r.loc[k,'p']:.2g})" for k in ("ins","del","indel","snp"))+f"   | ins with SNP+indel load {j.beta:+.2f} ({j.p:.2g}); with origin bio1 {b.beta:+.2f} ({b.p:.2g})")
print("\nFigure 2c / 3a — load change vs matched shared CDS SNPs (mean over gardens; matched-SNP null p, one-sided lower; deletions two-sided):")
L0=pd.read_csv(f"{G}/main_fig_checks_load.csv").set_index("cls"); L1=pd.read_csv(f"{G}/main_fig_checks_load{CR}.csv").set_index("cls")
print(pd.DataFrame({("90%","n"):L0.n_shared,("90%","load"):L0.load,("90%","p"):np.where(L0.index=="all CDS deletions",L0.p_two_sided,L0.p_lower),
                    ("__PCT__%","n"):L1.n_shared,("__PCT__%","load"):L1.load,("__PCT__%","p"):np.where(L1.index=="all CDS deletions",L1.p_two_sided,L1.p_lower)}).round(4).to_string())
print("\nFigure 2d — garden-level slope of the median residual on garden bio1 (per °C):")
for lab,sfx in (("90%",""),("__PCT__%",CR)):
    V_=pd.read_csv(f"{G}/sv_purging_annotated_event2{sfx}.csv"); S_=pd.read_csv(f"{G}/sv_temporal_climate_matched_sites_event2{sfx}.csv").sort_values("bio1"); P_=np.load(f"{G}/sv_temporal_climate_matched_persite_event2{sfx}.npz")
    out=[]
    for o in ("cold-specific","warm-specific"):
        for cls,k in (("insertion","ins"),("deletion","del"),("indel","indel")):
            ri=V_[(V_.cls==cls)&(V_.origin==o)].ridx.to_numpy(); g=np.array([np.nanmedian(P_[f"{s}_res_{k}_matched"][ri]) for s in S_.site.astype(int)]); lr=stats.linregress(S_.bio1,g)
            out.append(f"{o[:4]} {k} {lr.slope:+.4f} (n={ri.size:,})")
    print(f"  {lab:<4} "+"; ".join(out))
    tt=V_.groupby("cls").agg(n=("frac_neg","size"),mean_frac_neg=("frac_neg","mean")); print("       tested variants, mean fraction of gardens declining: "+"; ".join(f"{c} n={int(r.n):,} {r.mean_frac_neg:.3f}" for c,r in tt.iterrows()))
'''.replace("__CR__", CR).replace("__PCT__", str(PCT)).replace("__CALL__", str(CALL_MIN))

nb = nbf.read(SRC, as_version=4)
cells = []
for c in nb.cells:
    if c.cell_type == "code":
        n = new_code_cell(fix_code(c.source))
    else:
        n = new_markdown_cell(fix_md(c.source))
    cells.append(n)
left = [fn for fn in FILES if any(fn in c.source for c in cells if c.cell_type == "code")]
assert not left, f"untagged inputs left in the twin: {left}"
cells.insert(1, new_markdown_cell(md_note))
cells += [new_markdown_cell(md_cmp), new_code_cell(code_cmp)]
nb.cells = cells
ep = ExecutePreprocessor(timeout=1800, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": NBDIR}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
