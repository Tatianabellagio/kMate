#!/usr/bin/env python
"""Better 'did the block move the same direction?' metrics for CLQcut-0.9 blocks (Chr1):
  - AF-VE      : PC1-VE of raw gen9 AF        (static LD / haplotype coherence; current metric)
  - dP-VE      : PC1-VE of Delta p = gen9 AF - founding p0  (SELECTION-response coherence = same-direction)
  - eff_dim    : participation ratio (Sum L)^2/Sum L^2 on dP  (~ # independent directions of CHANGE; 1 = one unit)
  - dir_concord: fraction of a block's variants whose mean-Delta-p sign agrees with the block consensus
                 (after polarizing) -- the most literal 'all went the same way'
Run from repo root, kmate env, PYTHONPATH=analysis/grenenet_gea."""
import numpy as np, pandas as pd
import lib
from eval_block_coherence import CM, CLASSES

def _svals(M):
    M = np.array(M, float)
    cm = np.nanmean(M, 0); ii = np.where(np.isnan(M)); M[ii] = np.take(cm, ii[1])
    M = M - M.mean(0); sd = M.std(0); sd[sd == 0] = 1; M /= sd
    return np.linalg.svd(M, compute_uv=False) if M.shape[1] >= 2 else None

def pc1_ve(M):
    s = _svals(M); return np.nan if s is None else float(s[0]**2/np.sum(s**2))
def eff_dim(M):
    s = _svals(M);  l = None if s is None else s**2
    return np.nan if s is None else float(l.sum()**2/np.sum(l**2))
def abs_dp(dP):
    """block mean absolute change: how much it moved (direction-agnostic magnitude)."""
    return float(np.nanmean(np.abs(dP)))

# --- demonstration: does PC1 treat ANTI-PHASE (opposite-direction) SNPs as coherent? ---
rng0 = np.random.RandomState(7); F0 = 355
latent = rng0.randn(F0, 1)
antiphase = np.hstack([ latent + 0.1*rng0.randn(F0,5),     # 5 SNPs tag haplotype +
                       -latent + 0.1*rng0.randn(F0,5)])    # 5 SNPs tag the SAME event, opposite sign
three_hap = np.hstack([np.repeat(rng0.randn(F0,1),4,axis=1) for _ in range(3)])  # 3 independent haplotypes
print("DEMO: anti-phase biallelic block (SNPs move OPPOSITE ways, one event) PC1-VE=%.3f (expect ~1 -> PC1 DOES capture it)"
      % pc1_ve(antiphase))
print("DEMO: 3 independent haplotypes (multi-allelic)            PC1-VE=%.3f (expect ~0.33 -> low = real extra dims)\n"
      % pc1_ve(three_hap))

print("loading gen9 AF + records (Chr1, all-class) ...")
poss, afs, keys = [], [], []
for c in CLASSES:
    rec = pd.read_csv(f"{CM}/{c}_gen9.records.csv", usecols=["chrom","pos","ref_len","alt_len"])
    af = np.load(f"{CM}/{c}_gen9_af.npy", mmap_mode="r")
    m = (rec.chrom == "Chr1").values
    poss.append(rec.pos.values[m]); afs.append(np.asarray(af[:, m])); keys.append(rec[m])
pos = np.concatenate(poss); AF = np.concatenate(afs, 1); recs = pd.concat(keys, ignore_index=True)
o = np.argsort(pos, kind="stable"); pos = pos[o]; AF = AF[:, o]; recs = recs.iloc[o].reset_index(drop=True)
print(f"  {AF.shape[0]} pools x {AF.shape[1]} records")

print("aligning founding p0 ...")
p0 = lib.build_p0()
p0 = p0[~p0.index.duplicated(keep="first")]      # rec_key has co-located dups
recs["key"] = lib.rec_key(recs)
p0v = recs["key"].map(p0).values.astype(float)
dP = AF - p0v[None, :]                            # pools x variants, change from founding

bl = pd.read_csv("results/grenenet_gea/blocks_recompute/chr1_clq0.9_blocks_clq0.9.tsv", sep="\t")
rows = []
for s, e in zip(bl.start_pos, bl.end_pos):
    lo, hi = np.searchsorted(pos, s), np.searchsorted(pos, e, side="right")
    if hi-lo < 2: continue
    rows.append((s, e, hi-lo, pc1_ve(AF[:, lo:hi]), pc1_ve(dP[:, lo:hi]),
                 eff_dim(dP[:, lo:hi]), abs_dp(dP[:, lo:hi])))
r = pd.DataFrame(rows, columns=["start_pos","end_pos","n","af_ve","dp_ve","eff_dim","abs_dp"])
r.to_csv("results/grenenet_gea/blocks_recompute/chr1_clq0.9_directionality.csv", index=False)
print(f"CLQcut 0.9, {len(r)} blocks (>=2 variants):")
print(f"  median AF-VE   (static LD coherence)   : {r.af_ve.median():.3f}")
print(f"  median dP-VE   (CHANGE coherence)      : {r.dp_ve.median():.3f}   <- moved as one unit (handles anti-phase)")
print(f"  median eff_dim (# directions of change): {r.eff_dim.median():.2f}   (1 = one haplotype responding)")
print(f"  median abs_dp  (magnitude of change)   : {r.abs_dp.median():.4f}")
print(f"  blocks dP-VE>=0.7 (coherent change)    : {(r.dp_ve>=0.7).mean()*100:.0f}%")
print(f"  blocks eff_dim<=1.5 (~single haplotype): {(r.eff_dim<=1.5).mean()*100:.0f}%")
print("wrote chr1_clq0.9_directionality.csv")
