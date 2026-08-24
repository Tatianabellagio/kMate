#!/usr/bin/env python
"""Permutation-null check: is r2>=0.2 tagging real, or a max-over-many-
candidates artifact? For a sample of real testable SV anchors (MAC>=2 +
missingness floor, matching current production), shuffle the anchor's own
dose/call vector across founders (breaks true anchor<->SNP correspondence,
preserves the anchor's own MAF/missingness), recompute best r2 against the
SAME real candidate SNP window, and see how often the shuffled (null) best r2
also clears 0.2/0.5/0.8.
"""
import sys, numpy as np
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection")
import build_tagging_masked as btm

rng = np.random.default_rng(0)
CHROM = "Chr1"
MAC_FLOOR = 2
MIN_PAIR_FRAC = 0.5
N_SAMPLE = 500
N_PERM_PER_ANCHOR = 4

G, C, pos, rl, al, founders = btm._load_panel(CHROM)
NF = len(founders)
cls, size = btm.classify(rl, al)

ac = G.sum(0); an = C.sum(0)
mac = np.minimum(ac, an - ac)
mac_ok = mac >= MAC_FLOOR

anchor_mask = (cls == "SV") & mac_ok
anchor_idx = np.where(anchor_mask)[0]
anchor_pos = pos[anchor_idx]

snp_mask = (cls == "SNP") & mac_ok
snp_idx = np.where(snp_mask)[0]
snp_pos = pos[snp_idx]
Sfull = G[:, snp_idx].T
Scfull = C[:, snp_idx].T
order = np.argsort(snp_pos)
snp_pos, Sfull, Scfull = snp_pos[order], Sfull[order], Scfull[order]

min_pair = int(MIN_PAIR_FRAC * NF)
W = 50_000
lo_arr = np.searchsorted(snp_pos, anchor_pos - W, "left")
hi_arr = np.searchsorted(snp_pos, anchor_pos + W, "right")

A_dose_all = G[:, anchor_idx]
A_call_all = C[:, anchor_idx]

# restrict the diagnostic to anchors that actually have a non-trivial candidate
# window, sampled from those (matches what the production run would call "testable")
cand_counts = hi_arr - lo_arr
testable = np.where(cand_counts > 0)[0]
print(f"{len(anchor_idx):,} SV anchors (MAC>=2), {len(testable):,} with >=1 raw candidate in window", file=sys.stderr)

sample = rng.choice(testable, size=min(N_SAMPLE, len(testable)), replace=False)

real_best, null_best, n_window, n_at_best_real = [], [], [], []
for i in sample:
    lo, hi = lo_arr[i], hi_arr[i]
    S = Sfull[lo:hi]; Sc = Scfull[lo:hi]
    a_dose = A_dose_all[:, i]; a_call = A_call_all[:, i]

    r2, n = btm.masked_r2_block(a_dose, a_call, S, Sc, min_pair)
    valid = r2 >= 0
    n_window.append(int(valid.sum()))
    if valid.any():
        j = np.argmax(np.where(valid, r2, -1))
        real_best.append(r2[j])
        n_at_best_real.append(int(n[j]))
    else:
        real_best.append(np.nan)
        n_at_best_real.append(-1)

    # permutation null: shuffle (dose, call) jointly across the founder axis
    perm_bests = []
    for _ in range(N_PERM_PER_ANCHOR):
        perm = rng.permutation(NF)
        a_dose_p = a_dose[perm]; a_call_p = a_call[perm]
        r2p, np_ = btm.masked_r2_block(a_dose_p, a_call_p, S, Sc, min_pair)
        validp = r2p >= 0
        if validp.any():
            perm_bests.append(r2p[validp].max())
    null_best.append(np.nanmax(perm_bests) if perm_bests else np.nan)

real_best = np.array(real_best); null_best = np.array(null_best)
n_window = np.array(n_window); n_at_best_real = np.array(n_at_best_real)

rv = ~np.isnan(real_best); nv = ~np.isnan(null_best)
print(f"\nsampled {len(sample)} SV anchors, real testable {rv.sum()}, null-testable {nv.sum()}")
print(f"median # valid candidate SNPs in window (real): {np.median(n_window):.0f}")
print(f"median joint-called N at the REAL best-match pair: {np.median(n_at_best_real[n_at_best_real>=0]):.0f} (of {NF} founders)")
print("\nthreshold   real%    null%   (null = same anchor's own dose/call shuffled across founders,")
print("                              searched against the SAME real candidate window)")
for t in (0.2, 0.5, 0.8, 0.95):
    print(f"  r2>{t:<5}  {100*(real_best[rv]>t).mean():5.1f}%   {100*(null_best[nv]>t).mean():5.1f}%")
