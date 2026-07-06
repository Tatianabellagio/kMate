#!/usr/bin/env python3
"""Mechanism test for the pervasive SV-block JOINT excess (founder-GWAS, "selected at any
garden"): is it founder-STRUCTURE (clade) tagging, or a technical/estimation effect?

Established (sv_joint_diagnosis.py): the excess is ~1%, uniform across 30 gardens,
directionless, unconcentrated, locality at the noise floor. The founder-GWAS genotype is the
CLEAN panel (founder haplotype membership) and the trait is genome-wide ecotype h, so per-
block k-mer/estimation quality CANNOT drive per-block z (build_genotype docstring). What CAN:
the block's founder-genotype pattern. SVs often mark lineage/clade divergence -> SV-bearing
blocks may tag founder population structure that the genome-wide common-marker GRM doesn't
fully absorb, inflating a direction-agnostic 30-df omnibus (JOINT) but no 1-df climate axis.

This rebuilds the SAME founder genotype used by the GWAS (build_genotype('clq90'), same poly
filter), derives founder structure (PCs of the common-marker GRM), and asks:
  (a) do SV-bearing blocks load more on founder structure than size-matched SNP blocks?
  (b) is the SV chi2_joint excess concentrated in HIGH-structure-loading blocks?
  (c) mediation: does has_sv still predict chi2_joint after controlling structure + size?
  (d) genomic context: is the excess pericentromeric (low-recomb, structure-rich)?
  (e) Chr1 k-mer support (existing file): is the excess in low-k-mer-support blocks? (expect no)

Run in `basic` env from the repo root. Deterministic (seed=0). Reuses no new GWAS run.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import lib
from founder_genotype import build_genotype

OUT = Path("results/grenenet_gea/sv_adaptive")
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]
MAC_MIN, MAC_GRM = 3, 12
KPC = 20                                  # # founder structure axes (top GRM eigenvectors)
rng = np.random.default_rng(0)

# ---- rebuild the exact founder genotype behind the clq90 GWAS + same poly/grm filters ----
os.environ["MEMB_TAG"] = "clq90"
G, founders, reg = build_genotype("clq90"); nF = len(founders)
cnt = G.sum(0); blk = reg.block.to_numpy()
order = np.lexsort((-cnt, blk)); sbk = blk[order]
is_ref = np.zeros(len(cnt), bool)
is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
grm_set = (~is_ref) & (cnt >= MAC_GRM) & (cnt <= nF - MAC_GRM)
Gp = G[:, poly].astype(np.float64); regp = reg[poly].reset_index(drop=True)
print(f"{nF} founders | poly markers {poly.sum():,} | grm(common) {grm_set.sum():,}")

raw = lib.multisite_gwas_raw("clq90_pc1")
assert Gp.shape[1] == raw["M"], (Gp.shape[1], raw["M"])       # align with Z rows
regp_unit = (regp.chrom + ":" + regp.start.astype(str) + "-" + regp.end.astype(str)).to_numpy()
assert np.array_equal(regp_unit, raw["unit"]), "marker order drift vs GWAS npz"
chi2_joint = raw["chi2_joint"]; unit = raw["unit"]

# ---- founder structure: PCs of the common-marker GRM (same standardization as the GWAS) ----
sub = G[:, grm_set].astype(np.float64); p = sub.mean(0)
Zc = (sub - p) / np.sqrt(p * (1 - p) + 1e-9)
K = (Zc @ Zc.T) / grm_set.sum()
evals, evecs = np.linalg.eigh(K)
U = evecs[:, ::-1][:, :KPC]                                   # top-KPC founder axes (231 x KPC)
ve = evals[::-1] / evals.sum()
print(f"founder structure: top-{KPC} PCs explain {100*ve[:KPC].sum():.0f}% of GRM variance "
      f"(PC1 {100*ve[0]:.0f}%, PC1-5 {100*ve[:5].sum():.0f}%)")

# per-marker structure-loading = R^2 of its centered founder genotype on span(top-KPC PCs)
Gc = Gp - Gp.mean(0)
proj = U.T @ Gc                                               # (KPC, M)
ss_tot = (Gc ** 2).sum(0)
struct = (proj ** 2).sum(0) / np.clip(ss_tot, 1e-12, None)    # (M,) in [0,1]

# ---- block level: lead-JOINT marker per unit (self-consistent chi2_joint + structure) ----
# keep `unit` as a plain column throughout (join across mismatched index names drops it)
mk = pd.DataFrame({"unit": unit, "chi2_joint": chi2_joint, "struct": struct})
nmark = mk.groupby("unit").size().rename("n_markers").reset_index()   # tested hap-clusters/block
lead = mk.loc[mk.groupby("unit")["chi2_joint"].idxmax()].merge(nmark, on="unit")
L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv").rename(columns={"block_id": "unit"})
B = lead.merge(L[["unit", "n_kept", "has_sv", "dist_cen", "chrom"]], on="unit", how="inner")
B["bin"] = pd.cut(B.n_kept, EDGES, right=False, labels=False)
B = B[B.bin.notna()].reset_index(drop=True); B["bin"] = B["bin"].astype(int)
bins = B["bin"].to_numpy(); sv = np.where(B.has_sv.to_numpy() == 1)[0]
print(f"{len(B):,} blocks | {len(sv):,} SV-bearing | "
      f"n_markers/block: SNP-only median {B.loc[B.has_sv==0,'n_markers'].median():.0f}, "
      f"SV median {B.loc[B.has_sv==1,'n_markers'].median():.0f}")

# ---- (a) SV vs size-matched SNP: founder-structure loading ----
o, n, r, pp = lib.matched_perm_test(B.struct.to_numpy(), bins, sv, 5000)
print(f"\n(a) founder-structure loading: SV {o:.3f} vs size-matched {n:.3f}  x{r:.3f}  p={pp:.4f} "
      f"({'SV blocks tag structure MORE' if r > 1 and pp < 0.05 else 'n.s.'})")

# ---- (b) SV chi2_joint excess within structure-loading tertiles ----
print("\n(b) SV chi2_joint excess (size-matched) by founder-structure-loading tertile:")
qt = pd.qcut(B.struct.rank(method="first"), 3, labels=["low", "mid", "high"])
for lab in ["low", "mid", "high"]:
    idx = np.where((qt == lab).to_numpy())[0]
    sub_bins = B["bin"].to_numpy()
    svsub = idx[B.has_sv.to_numpy()[idx] == 1]
    if len(svsub) < 20:
        print(f"  {lab:4s}: too few SV blocks"); continue
    o, n, r, pp = lib.matched_perm_test(B.chi2_joint.to_numpy(), sub_bins, svsub, 5000)
    print(f"  {lab:4s} struct ({B.struct[idx].min():.2f}-{B.struct[idx].max():.2f}): "
          f"SV chi2_joint {o:.2f} vs matched {n:.2f}  x{r:.3f}  p={pp:.4f}  (n_sv={len(svsub)})")

# ---- (c) mediation: does has_sv survive controlling structure + size? ----
def z(v): return (v - v.mean()) / v.std()
y = z(np.log(B.chi2_joint + 1e-9))
X0 = np.column_stack([np.ones(len(B)), B.has_sv.to_numpy(float), z(np.log(B.n_kept))])
X1 = np.column_stack([X0, z(B.struct)])
b0 = np.linalg.lstsq(X0, y, rcond=None)[0]
b1 = np.linalg.lstsq(X1, y, rcond=None)[0]
print(f"\n(c) OLS z(log chi2_joint) ~ has_sv + size (+structure):")
print(f"  has_sv coef WITHOUT structure: {b0[1]:+.4f}")
print(f"  has_sv coef WITH   structure: {b1[1]:+.4f}   (structure coef {b1[3]:+.4f})")
print(f"  => structure explains {100*(1 - b1[1]/b0[1]):.0f}% of the has_sv->JOINT association"
      if b0[1] != 0 else "")

# ---- (d) genomic context: pericentromeric vs arm ----
print("\n(d) SV chi2_joint excess by genomic context (dist_cen split at median):")
med = B.dist_cen.median()
for lab, m in [("pericentromeric", B.dist_cen <= med), ("arm", B.dist_cen > med)]:
    idx = np.where(m.to_numpy())[0]; svsub = idx[B.has_sv.to_numpy()[idx] == 1]
    o, n, r, pp = lib.matched_perm_test(B.chi2_joint.to_numpy(), B["bin"].to_numpy(), svsub, 5000)
    print(f"  {lab:16s}: SV chi2_joint {o:.2f} vs matched {n:.2f}  x{r:.3f}  p={pp:.4f}  (n_sv={len(svsub)})")
rs = stats.spearmanr(B.struct, -np.log10(B.dist_cen + 1))
print(f"  structure-loading vs proximity-to-centromere: Spearman rho={rs.correlation:+.3f} p={rs.pvalue:.2g}")

# ---- (f) residual multiple-testing confound: do SV blocks resolve into more hap-clusters
# per block at matched n_kept? If so, block max-over-markers chi2_joint inflates for free.
# Test both: (i) n_markers excess at matched n_kept, (ii) re-match chi2_joint on n_markers. ----
print("\n(f) hap-cluster count per block as a residual confound:")
o, n, r, pp = lib.matched_perm_test(B.n_markers.to_numpy(float), bins, sv, 5000)
print(f"  (i) n_markers/block: SV {o:.2f} vs size(n_kept)-matched {n:.2f}  x{r:.3f}  p={pp:.4f}")
mbin = np.minimum(B.n_markers.to_numpy(), 12)                # integer bins 1..12+ on tested-marker count
o2, n2, r2, pp2 = lib.matched_perm_test(B.chi2_joint.to_numpy(), mbin, sv, 5000)
print(f"  (ii) SV chi2_joint re-matched on n_markers (not n_kept): x{r2:.3f}  p={pp2:.4f} "
      f"({'excess SURVIVES marker-count matching' if r2 > 1 and pp2 < 0.05 else 'excess GONE -> was a per-block test-count effect'})")
# combined double-matching (n_kept AND n_markers bin) for the headline
dbin = B["bin"].to_numpy() * 100 + mbin
o3, n3, r3, pp3 = lib.matched_perm_test(B.chi2_joint.to_numpy(), dbin, sv, 5000)
print(f"  (iii) SV chi2_joint double-matched (n_kept x n_markers): x{r3:.3f}  p={pp3:.4f}")

# ---- (e) Chr1 k-mer support (existing file): sanity that estimation quality is a non-factor ----
kf = OUT.parent / "blocks_mcf90/chr1_block_kmer_coverage.csv"
if kf.exists():
    kc = pd.read_csv(kf)
    kc["unit"] = "Chr1:" + kc.start_pos.astype(str) + "-" + kc.end_pos.astype(str)
    B1 = B[B.chrom == "Chr1"].merge(kc[["unit", "panel_kmers"]], on="unit", how="inner")
    print(f"\n(e) Chr1 k-mer support: {len(B1):,} blocks; SV chi2_joint excess by panel_kmers tertile:")
    B1 = B1.reset_index(drop=True); qk = pd.qcut(B1.panel_kmers.rank(method="first"), 3, labels=["low", "mid", "high"])
    for lab in ["low", "mid", "high"]:
        idx = np.where((qk == lab).to_numpy())[0]; svsub = idx[B1.has_sv.to_numpy()[idx] == 1]
        if len(svsub) < 15:
            print(f"  {lab:4s}: too few SV blocks (n={len(svsub)})"); continue
        o, n, r, pp = lib.matched_perm_test(B1.chi2_joint.to_numpy(), B1["bin"].to_numpy(), svsub, 5000)
        print(f"  {lab:4s} kmers ({int(B1.panel_kmers[idx].min())}-{int(B1.panel_kmers[idx].max())}): "
              f"SV chi2_joint x{r:.3f}  p={pp:.4f}  (n_sv={len(svsub)})")

B.to_csv(OUT / "sv_joint_mechanism_blocks.csv", index=False)
print(f"\n[done] -> {OUT}/sv_joint_mechanism_blocks.csv")
