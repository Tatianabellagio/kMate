#!/usr/bin/env python
"""Founding (gen-0) haplotype-frequency p0 from the 8 SEEDMIX window-mode reps.

Projects each SEEDMIX rep's per-block founder h onto the SAME dynld-K500 haplotype
units as the cohort matrix (results/grenenet_gea/hapfreq/) by reusing the canonical
load_mem + project_sample from build_hapfreq_matrix (STEP-1 deterministic membership).
This guarantees the p0 vector is column-aligned to hapfreq_matrix.npy.

  p0 = mean over the 8 reps of the per-haplotype frequency
  v0 = among-rep variance of logit(freq) / n_reps   (gen-0 point variance, Pipeline B)

Outputs (results/grenenet_gea/hapfreq/):
  hapfreq_p0_seedmix.npy   (n_haps,)  empirical founding haplotype freq
  hapfreq_v0_seedmix.npy   (n_haps,)  logit-variance of p0 / n_reps
  hapfreq_p0_seedmix_reps.npy  (n_reps, n_haps)  per-rep matrix (audit)

Env: kmate. Run after run_seedmix_window.sbatch (array 35169327) finishes.
"""
import glob, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_hapfreq_matrix as bhm

SM = "results/grenenet_kmate_window_seedmix"
OUT = os.environ.get("HF_DIR", "results/grenenet_gea/hapfreq")   # clq90: hapfreq_clq90
EPS = 1e-3


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    # point the canonical projector at the SEEDMIX window outputs; inherit membership tag +
    # h-source from env so p0 is column-aligned to the cohort matrix (clq90/global).
    bhm.WIN = SM
    bhm.MEMB_TAG = os.environ.get("MEMB_TAG", bhm.MEMB_TAG)
    bhm.H_SOURCE = os.environ.get("H_SOURCE", bhm.H_SOURCE)
    reps = sorted(os.path.basename(p).replace("_Chr1.h_blocks_per_chrom.npz", "")
                  for p in glob.glob(f"{SM}/*_Chr1.h_blocks_per_chrom.npz"))
    if not reps:
        sys.exit(f"no SEEDMIX h_blocks in {SM} — has the window array finished?")
    print(f"{len(reps)} SEEDMIX reps: {reps}")

    mem, reg, nh = bhm.load_mem(registry=True)

    # alignment check vs the cohort matrix columns
    cohort_reg = f"{OUT}/hapfreq_registry.csv"
    if os.path.exists(cohort_reg):
        import pandas as pd
        n_cohort = len(pd.read_csv(cohort_reg, usecols=["hap_id"]))
        assert n_cohort == nh == len(reg), \
            f"column mismatch: cohort {n_cohort} vs p0 {nh} (re-run STEP 1/2 with same membership)"
        print(f"aligned to cohort matrix: {nh:,} haplotypes")

    M = np.zeros((len(reps), nh), dtype=np.float32)
    for i, s in enumerate(reps):
        M[i] = bhm.project_sample(s, mem)
        print(f"  projected {s} (sum/unit≈{M[i].sum()/len(reg.unit_idx.unique()):.3f})", flush=True)

    p0 = M.mean(0)
    v0 = logit(M).var(0, ddof=1) / M.shape[0]
    np.save(f"{OUT}/hapfreq_p0_seedmix.npy", p0.astype(np.float32))
    np.save(f"{OUT}/hapfreq_v0_seedmix.npy", v0.astype(np.float32))
    np.save(f"{OUT}/hapfreq_p0_seedmix_reps.npy", M)
    print(f"[done] p0 mean {p0.mean():.4f}, v0 median {np.median(v0):.3f}; "
          f"{p0.size:,} haplotypes -> {OUT}/hapfreq_p0_seedmix.npy")


if __name__ == "__main__":
    main()
