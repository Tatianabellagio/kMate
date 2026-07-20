#!/usr/bin/env python
"""Genome-wide permutation-null sample for the SV/indel-SNP tagging survival
curves: for a sample of real (MAC>=2, min-pair-frac>=0.5, production-matching)
anchors per chrom x class x panel-mode, shuffle the anchor's own dose/call
vector across founders (destroys true anchor<->SNP correspondence, keeps the
anchor's own MAF/missingness) and re-search the SAME real candidate window,
recording the resulting "null best r2". Loads G/C once per chromosome and
reuses it across both anchor classes and both panel-modes to amortize the
expensive load.

Single sequential process by design (16GB interactive cgroup cap -- see
feedback_kmate_submit_dont_run_interactive.md); run via sbatch if this ever
needs to scale beyond a few thousand draws/combo.
"""
import sys, time
import numpy as np
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea")
import build_tagging_masked as btm

CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
MAC_FLOOR = 2
MIN_PAIR_FRAC = 0.5
N_SAMPLE_PER_CHROM = 400   # per (class, mode) per chrom -> ~2000 draws genome-wide per combo
OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/sv_snp_ld_v2/tagging_null.npz"

rng = np.random.default_rng(1)
records = []  # (cls, mode, chrom, real_r2, null_r2)

t_start = time.time()
for chrom in CHROMS:
    t0 = time.time()
    G, C, pos, rl, al, founders = btm._load_panel(chrom)
    NF = len(founders)
    cls_arr, size = btm.classify(rl, al)
    ac = G.sum(0); an = C.sum(0)
    mac = np.minimum(ac, an - ac)
    mac_ok = mac >= MAC_FLOOR
    print(f"[{chrom}] loaded in {time.time()-t0:.0f}s", file=sys.stderr)

    for anchor_class in ["sv", "indel"]:
        cls_label = "SV" if anchor_class == "sv" else "indel"
        anchor_mask = (cls_arr == cls_label) & mac_ok
        anchor_idx_all = np.where(anchor_mask)[0]
        anchor_pos_all = pos[anchor_idx_all]

        for panel_mode in ["panel", "shortread"]:
            tm0 = time.time()
            if panel_mode == "panel":
                snp_mask = (cls_arr == "SNP") & mac_ok
                snp_idx = np.where(snp_mask)[0]
                snp_pos = pos[snp_idx]
                Sfull = G[:, snp_idx].T
                Scfull = C[:, snp_idx].T
                n_shared = NF
                G_use, C_use = G, C
                anchor_idx, anchor_pos = anchor_idx_all, anchor_pos_all
            else:
                z = np.load(btm.SR_GENO.format(Cl=chrom), allow_pickle=True)
                sr_founders = [str(x) for x in z["founders"]]
                idx_map = {f: i for i, f in enumerate(sr_founders)}
                shared = [f for f in founders if f in idx_map]
                sr_col = [idx_map[f] for f in shared]
                panel_col = [founders.index(f) for f in shared]
                snp_pos = z["snp_pos"].astype(np.int64)
                Sfull = z["geno"][:, sr_col].astype(np.float32)
                Scfull = np.ones_like(Sfull)
                sr_ac = Sfull.sum(1) / 2.0
                sr_an = np.full(Sfull.shape[0], Sfull.shape[1], dtype=np.float64)
                sr_mac = np.minimum(sr_ac, sr_an - sr_ac)
                keep_snp = sr_mac >= MAC_FLOOR
                snp_pos, Sfull, Scfull = snp_pos[keep_snp], Sfull[keep_snp], Scfull[keep_snp]
                n_shared = len(shared)
                G_use, C_use = G[panel_col], C[panel_col]
                anchor_idx, anchor_pos = anchor_idx_all, anchor_pos_all

            order = np.argsort(snp_pos)
            snp_pos, Sfull, Scfull = snp_pos[order], Sfull[order], Scfull[order]
            min_pair = int(MIN_PAIR_FRAC * n_shared)
            W = 50_000
            lo_arr = np.searchsorted(snp_pos, anchor_pos - W, "left")
            hi_arr = np.searchsorted(snp_pos, anchor_pos + W, "right")
            cand_counts = hi_arr - lo_arr
            testable = np.where(cand_counts > 0)[0]
            if len(testable) == 0:
                continue

            n_sample = min(N_SAMPLE_PER_CHROM, len(testable))
            sample = rng.choice(testable, size=n_sample, replace=False)

            A_dose_all = G_use[:, anchor_idx]
            A_call_all = C_use[:, anchor_idx]

            for i in sample:
                lo, hi = lo_arr[i], hi_arr[i]
                S = Sfull[lo:hi]; Sc = Scfull[lo:hi]
                a_dose = A_dose_all[:, i]; a_call = A_call_all[:, i]

                r2, n = btm.masked_r2_block(a_dose, a_call, S, Sc, min_pair)
                valid = r2 >= 0
                real_r2 = r2[valid].max() if valid.any() else np.nan

                perm = rng.permutation(len(a_dose))
                a_dose_p = a_dose[perm]; a_call_p = a_call[perm]
                r2p, _ = btm.masked_r2_block(a_dose_p, a_call_p, S, Sc, min_pair)
                validp = r2p >= 0
                null_r2 = r2p[validp].max() if validp.any() else np.nan

                records.append((cls_label, panel_mode, chrom, real_r2, null_r2))
            print(f"  [{chrom}/{anchor_class}/{panel_mode}] {n_sample} draws in {time.time()-tm0:.0f}s "
                  f"(median candidates/window: {np.median(cand_counts[sample]):.0f})", file=sys.stderr)
    del G, C
    print(f"[{chrom}] total {time.time()-t0:.0f}s, cumulative {time.time()-t_start:.0f}s", file=sys.stderr)

cls_a = np.array([r[0] for r in records])
mode_a = np.array([r[1] for r in records])
chrom_a = np.array([r[2] for r in records])
real_a = np.array([r[3] for r in records], dtype=np.float32)
null_a = np.array([r[4] for r in records], dtype=np.float32)
np.savez(OUT, cls=cls_a, mode=mode_a, chrom=chrom_a, real_r2=real_a, null_r2=null_a)
print(f"\nwrote {len(records):,} draws to {OUT}; total time {time.time()-t_start:.0f}s", file=sys.stderr)

for cls_label in ["SV", "indel"]:
    for mode in ["panel", "shortread"]:
        m = (cls_a == cls_label) & (mode_a == mode)
        rv = real_a[m]; nv = null_a[m]
        rv, nv = rv[~np.isnan(rv)], nv[~np.isnan(nv)]
        print(f"{cls_label}/{mode}: n={m.sum()}  ", end="")
        for t in (0.2, 0.5, 0.8, 0.95):
            print(f"r2>{t}: real {100*(rv>t).mean():.1f}% null {100*(nv>t).mean():.1f}%  ", end="")
        print()
