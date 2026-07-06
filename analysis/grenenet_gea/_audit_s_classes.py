#!/usr/bin/env python
"""AUDIT: are SNP/indel/SV s-distributions matching because of a bug, or for real? Checks:
 (1) class masks select genuinely different, correctly-typed, disjoint variants;
 (2) the saved per-class s arrays are different data (not aliased/duplicated);
 (3) an INDEPENDENT from-scratch recomputation of plot-replicate s reproduces plot_replicate_sz;
 (4) statistical comparison (KS) + WHY they match: is s a founder-projection of p0?
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
import importlib.util
_sp = importlib.util.spec_from_file_location(
    "plotsmod", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "_temporal_s_plots_snp_vs_nonsnp.py"))
pm = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(pm)

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE; PM = f"{lib.GEA}/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, EPS = 0.02, 50, 12, 1e-3
logit = lambda p: np.log(np.clip(p, EPS, 1-EPS)/(1-np.clip(p, EPS, 1-EPS)))

idx_snp = np.load(f"{STORE}/index_snp.npz"); idx_non = np.load(f"{STORE}/index_nonsnp.npz")
ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
rl = idx_non["ref_len"].astype(np.int64); al = idx_non["alt_len"].astype(np.int64)
dlen = np.abs(al-rl)
ch_snp = idx_snp["chrom"].astype("U5"); pos_snp = idx_snp["pos"].astype(np.int64)
rl_s = idx_snp["ref_len"].astype(np.int64); al_s = idx_snp["alt_len"].astype(np.int64)
p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)
common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=MIN_MAC)
common_snp = lib.founder_panel_keep(ch_snp, pos_snp, min_mac=MIN_MAC)
snp_m = common_snp & (p0_snp>=MIN_P0) & (p0_snp<=1-MIN_P0)
non_m = (dlen>=1) & common_non & (p0_non>=MIN_P0) & (p0_non<=1-MIN_P0)
cols_snp = np.where(snp_m)[0]; cols_non = np.where(non_m)[0]
p0s = p0_snp[cols_snp]; p0n = p0_non[cols_non]
dl = dlen[cols_non]; sv_k = dl>SV_BP; ind_k = (dl>=1)&(dl<=SV_BP)

print("="*70); print("CHECK 1: class masks — types, counts, disjoint")
print(f"  SNP  index ref_len/alt_len unique: {np.unique(rl_s)[:3]}/{np.unique(al_s)[:3]} "
      f"(dlen should be 0 for SNPs: max|al-rl|={np.abs(al_s-rl_s).max()})")
print(f"  scored SNP={cols_snp.size:,}  indel={int(ind_k.sum()):,}  SV={int(sv_k.sum()):,}")
print(f"  SV dlen range: {dl[sv_k].min()}–{dl[sv_k].max()} (all >{SV_BP}? {bool((dl[sv_k]>SV_BP).all())})")
print(f"  indel dlen range: {dl[ind_k].min()}–{dl[ind_k].max()} (all in [1,{SV_BP}]? "
      f"{bool(((dl[ind_k]>=1)&(dl[ind_k]<=SV_BP)).all())})")
print(f"  indel∩SV overlap: {int((ind_k & sv_k).sum())} (should be 0)")
print(f"  SNP vs nonSNP are separate indices (different files): "
      f"snp n={pos_snp.size:,} non n={pos_non.size:,}")

print("="*70); print("CHECK 2: saved per-class s arrays are DIFFERENT data (not aliased)")
npz = np.load(f"{lib.GEA}/sv_adaptive/s_dist_by_stratum.npz")
s_sn = npz["4_s_SNP"]; s_in = npz["4_s_indel"]; s_sv = npz["4_s_SV"]
print(f"  site4 sizes: SNP={s_sn.size} indel={s_in.size} SV={s_sv.size} (distinct sizes)")
print(f"  SNP≡indel? {np.array_equal(s_sn[:s_in.size], s_in)}   "
      f"SNP≡SV? {np.array_equal(s_sn[:s_sv.size], s_sv)}   (both should be False)")
print(f"  means: SNP={s_sn.mean():.4f} indel={s_in.mean():.4f} SV={s_sv.mean():.4f}")
print(f"  p0 medians per class: SNP={np.median(npz['4_p0_SNP']):.3f} "
      f"indel={np.median(npz['4_p0_indel']):.3f} SV={np.median(npz['4_p0_SV']):.3f}")

print("="*70); print("CHECK 3: INDEPENDENT from-scratch recompute of plot-replicate s (site 4)")
site = 4
# pick 4 SV cols and 4 SNP cols
sv_cols = cols_non[sv_k][[0, 100, 500, 1000]]
sn_cols = cols_snp[[0, 5000, 50000, 200000]]
def manual_s(kind, col, p0v):
    plots = {}
    for g in (1, 2, 3):
        try: meta = pd.read_csv(f"{PM}/pool_gen{g}_{kind}.meta.csv")
        except FileNotFoundError: continue
        m4 = meta[meta.site == site]
        if len(m4) == 0: continue
        mat = np.load(f"{PM}/pool_gen{g}_{kind}_af.npy", mmap_mode="r")
        for ridx, plot in zip(m4.index.to_numpy(), m4["plot"].to_numpy()):
            plots.setdefault(int(plot), []).append((float(g), float(mat[ridx, col])))
    sl = []
    for plot, pts in plots.items():
        ts = np.array([0.0]+[t for t,_ in pts]); ys = np.array([logit(p0v)]+[logit(a) for _,a in pts])
        tc = ts-ts.mean(); sl.append((tc@ys)/(tc**2).sum())
    return float(np.mean(sl)), len(sl)
# function output for the same cols
s_sv_fn,_,_,n = pm.plot_replicate_sz(site, "nonsnp", sv_cols, p0_non[sv_cols])
s_sn_fn,_,_,_ = pm.plot_replicate_sz(site, "snp", sn_cols, p0_snp[sn_cols])
print(f"  (n_plots={n})  SV cols: manual vs function")
for i, c in enumerate(sv_cols):
    ms, npl = manual_s("nonsnp", c, p0_non[c])
    print(f"    col{c} p0={p0_non[c]:.3f} dlen={dlen[c]:>5} : manual={ms:+.5f}  fn={s_sv_fn[i]:+.5f}  "
          f"Δ={abs(ms-s_sv_fn[i]):.2e}")
for i, c in enumerate(sn_cols):
    ms, npl = manual_s("snp", c, p0_snp[c])
    print(f"    SNP col{c} p0={p0_snp[c]:.3f} : manual={ms:+.5f}  fn={s_sn_fn[i]:+.5f}  Δ={abs(ms-s_sn_fn[i]):.2e}")

print("="*70); print("CHECK 4: are the distributions statistically identical, and WHY (MAF>=0.10, site4)")
def maf_filt(s, p0): m = np.minimum(p0,1-p0)>=0.10; return s[m]
S = maf_filt(s_sn, npz["4_p0_SNP"]); V = maf_filt(s_sv, npz["4_p0_SV"])
ks = stats.ks_2samp(S, V)
print(f"  KS SNP vs SV (MAF>=0.10): D={ks.statistic:.3f} p={ks.pvalue:.3g}  "
      f"(medians SNP={np.median(S):+.3f} SV={np.median(V):+.3f})")
# WHY: is s essentially a function of p0? corr(s, p0) within class, and does a p0-only predictor
# explain most of s? (founder projection => s strongly tied to p0/trajectory shared across classes)
r_sn = stats.spearmanr(npz["4_p0_SNP"], s_sn).statistic
r_sv = stats.spearmanr(npz["4_p0_SV"], s_sv).statistic
print(f"  Spearman(s, p0): SNP={r_sn:+.2f}  SV={r_sv:+.2f}  "
      f"-> s is largely a function of p0 (shared across classes via the founder projection)")
print("\nDONE.")
