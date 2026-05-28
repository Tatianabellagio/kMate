"""
Recreate AF_TRUTH_VS_ESTIMATE_v3qc_v3_mixedloose for Arch 3 chr1.

Builds three views on SNP-only records:
  1. RECIPE TRUTH = uniform-h projection through cn_var (h = 1/F for all founders)
     This is the AF you would observe in a pool where every founder contributes
     equally — a panel-intrinsic baseline that has no sequencing/projection error.
  2. hapFIRE AF (xwu's 2023 production HARP+CVXPY on the 231-panel greneNet_v1.1)
  3. cactus_em AF (h_v3 @ cn_var_arch3 from A6)

4-tuple join on (chrom, pos, ref, alt) so multi-allelic-split SNPs don't manufacture
off-diagonal scatter.

Output:  plots/AF_TRUTH_VS_ESTIMATE_arch3_chr1.png  (4-panel scatter + summary)
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.sparse import load_npz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
ARCH = ROOT / "arch3" / "chr1"
CN_VAR     = ARCH / "cn_var_231_arch3_chr1.cn_var.npz"
CN_CALLED  = ARCH / "cn_var_231_arch3_chr1.cn_var_called.npz"
META       = ARCH / "cn_var_231_arch3_chr1.meta.npz"
SEEDMIX    = ARCH / "SEEDMIX_S1_arch3_chr1.tsv"   # h_v3 @ cn_var_arch3
HAPFIRE    = "/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt"
HAPFIRE_VCF = "/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf"

OUT_PLOT = ROOT / "plots" / "AF_TRUTH_VS_ESTIMATE_arch3_chr1.png"
OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)

print("[load] cn_var + cn_var_called + meta")
cn_var       = load_npz(CN_VAR).tocsr()
cn_var_called = load_npz(CN_CALLED).tocsr()
meta = np.load(META, allow_pickle=True)
pos  = meta["pos"]
chrom = meta["chrom"]
ref_arr = meta["ref"]   # object array
alt_arr = meta["alt"]   # object array
ref_len = meta["ref_len"]
alt_len = meta["alt_len"]
F = cn_var.shape[0]
N = cn_var.shape[1]
print(f"  cn_var: {cn_var.shape}, {cn_var.nnz:,} nnz")

print("[build] recipe truth — uniform-h projection (h_i = 1/F)")
h_uniform = np.full(F, 1.0/F, dtype=np.float64)
numer = h_uniform @ cn_var          # shape (N,)
denom = h_uniform @ cn_var_called   # shape (N,)
recipe_truth = np.where(denom > 0, numer / denom, np.nan)
print(f"  recipe_truth: AF range [{np.nanmin(recipe_truth):.3f}, {np.nanmax(recipe_truth):.3f}], NaN={np.isnan(recipe_truth).sum():,}")

print("[load] SEEDMIX_S1 cactus_em arch3 projection")
sm = pd.read_csv(SEEDMIX, sep="\t")
print(f"  SEEDMIX TSV rows: {len(sm):,}  (matches meta records: {len(sm)==N})")
cem_af = sm["alt_freq"].to_numpy()  # row-aligned with meta

print("[subset] SNPs only (ref_len = alt_len = 1)")
snp_mask = (ref_len == 1) & (alt_len == 1)
print(f"  SNP records: {snp_mask.sum():,} / {N:,}")

# Convert ref/alt subset to plain Python strings (object array → list of str — keeps memory low)
ref_snp = [str(ref_arr[i]) for i in np.where(snp_mask)[0]]
alt_snp = [str(alt_arr[i]) for i in np.where(snp_mask)[0]]
chrom_snp = [str(chrom[i]) for i in np.where(snp_mask)[0]]
pos_snp = pos[snp_mask]

# Build DataFrame for SNP records
snp_df = pd.DataFrame({
    "chrom_num": [int(c.replace("Chr","")) for c in chrom_snp],
    "pos": pos_snp,
    "ref": ref_snp,
    "alt": alt_snp,
    "recipe_truth": recipe_truth[snp_mask],
    "cem_af": cem_af[snp_mask],
})
print(f"  SNP DataFrame: {len(snp_df):,} rows")

print("[load] hapFIRE AF + REF/ALT from greneNet_v1.1")
hf = pd.read_csv(HAPFIRE, sep="\t", header=None, names=["chrom_num","pos","af_hapfire"])
hf = hf[hf.chrom_num == 1].copy()
ra = pd.read_csv(HAPFIRE_VCF, sep="\t", comment="#", header=None,
                 usecols=[0,1,3,4], names=["chrom_num","pos","ref","alt"])
ra = ra[ra.chrom_num.astype(str) == "1"].copy()
ra["chrom_num"] = 1
hf_full = hf.merge(ra, on=["chrom_num","pos"], how="left").dropna(subset=["ref","alt"])
print(f"  hapFIRE chr1 SNPs (w/ REF/ALT): {len(hf_full):,}")

print("[join] 4-tuple (chrom, pos, ref, alt) — snp_df ∩ hapFIRE")
joined = snp_df.merge(hf_full[["chrom_num","pos","ref","alt","af_hapfire"]],
                      on=["chrom_num","pos","ref","alt"], how="inner")
print(f"  4-tuple intersect: {len(joined):,}")

def metrics(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    d = a[m] - b[m]
    return {
        "n": int(m.sum()),
        "MAE": float(np.abs(d).mean()),
        "RMSE": float(np.sqrt((d**2).mean())),
        "bias": float(d.mean()),
        "outlier_pct": float(100*(np.abs(d) > 0.10).sum() / max(1,m.sum())),
        "r": float(np.corrcoef(a[m], b[m])[0,1]),
    }

m_rec_vs_hf = metrics(joined["recipe_truth"].to_numpy(), joined["af_hapfire"].to_numpy())
m_rec_vs_cem = metrics(joined["recipe_truth"].to_numpy(), joined["cem_af"].to_numpy())
m_cem_vs_hf = metrics(joined["cem_af"].to_numpy(), joined["af_hapfire"].to_numpy())

print()
print("=== Summary metrics on 4-tuple intersect ===")
for name, m in [("recipe_truth vs hapFIRE", m_rec_vs_hf),
                ("recipe_truth vs cactus_em", m_rec_vs_cem),
                ("cactus_em vs hapFIRE", m_cem_vs_hf)]:
    print(f"  {name:<28} n={m['n']:>8,}  MAE={m['MAE']:.4f}  RMSE={m['RMSE']:.4f}  bias={m['bias']:+.4f}  r={m['r']:.4f}  |Δ|>0.10: {m['outlier_pct']:.2f}%")

# === Plot ===
print()
print("[plot] 4-panel scatter")
fig, axes = plt.subplots(2, 2, figsize=(14, 13))

def hexpanel(ax, x, y, xlabel, ylabel, title, m):
    hb = ax.hexbin(x, y, gridsize=80, bins="log", cmap="viridis", mincnt=1)
    ax.plot([0,1],[0,1],"r--",lw=1,alpha=0.6)
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nn={m['n']:,}  MAE={m['MAE']:.4f}  bias={m['bias']:+.4f}  r={m['r']:.4f}")
    plt.colorbar(hb, ax=ax, label="log10(count)")

x_rec = joined["recipe_truth"].to_numpy()
x_hf  = joined["af_hapfire"].to_numpy()
x_cem = joined["cem_af"].to_numpy()

hexpanel(axes[0,0], x_rec, x_hf,
         "recipe truth AF (uniform-h ∗ Arch 3 cn_var)", "hapFIRE AF",
         "Recipe truth vs hapFIRE\n(panel intrinsic vs 1001G short-read calls)", m_rec_vs_hf)
hexpanel(axes[0,1], x_rec, x_cem,
         "recipe truth AF", "cactus_em SEEDMIX_S1 AF (h_v3 ∗ Arch 3 cn_var)",
         "Recipe truth vs cactus_em\n(should be ~perfect if h ≈ uniform, slope deviation = sample composition)", m_rec_vs_cem)
hexpanel(axes[1,0], x_cem, x_hf,
         "cactus_em SEEDMIX_S1 AF", "hapFIRE AF",
         "cactus_em vs hapFIRE\n(pool-seq estimate vs canonical reference)", m_cem_vs_hf)

# Residual diagnostic — recipe truth - hapFIRE colored by something
ax = axes[1,1]
resid = x_rec - x_hf
hb = ax.hexbin(x_hf, resid, gridsize=80, bins="log", cmap="coolwarm", mincnt=1)
ax.axhline(0, color="black", lw=0.5, alpha=0.5)
ax.set_xlim(0,1); ax.set_ylim(-0.5, 0.5)
ax.set_xlabel("hapFIRE AF"); ax.set_ylabel("recipe_truth − hapFIRE")
ax.set_title("Recipe-truth residual vs hapFIRE\n(positive = our panel says higher AF than 1001G)")
plt.colorbar(hb, ax=ax, label="log10(count)")

plt.suptitle("Arch 3 chr1 — AF truth vs estimate (SEEDMIX_S1)", fontsize=14, y=1.00)
plt.tight_layout()
plt.savefig(OUT_PLOT, dpi=120, bbox_inches="tight")
print(f"  wrote {OUT_PLOT}")

# Also write the joined table for downstream notebook use
out_tsv = ARCH / "AF_TRUTH_VS_ESTIMATE_arch3_chr1.tsv"
joined.to_csv(out_tsv, sep="\t", index=False)
print(f"  wrote {out_tsv}  ({len(joined):,} rows)")
