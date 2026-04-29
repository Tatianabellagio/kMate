"""
Real-BAM end-to-end test on a visor_freqk simulation.

Setup:
- Pool: 231 GrENE-Net ecotypes uniformly at 1/231 each (truth from visor_freqk)
- Panel: 82 cactus founders, of which 80 overlap with GrENE-Net 231

Expected behavior:
- Founder freq vector h should sum to 1 (simplex constraint)
- 80 GrENE-overlap founders should each carry weight near 1/231 ≈ 0.0043,
  but inflated to absorb missing mass of the 151 unobserved GrENE founders
  (so each would be ~80×(1/231)/80 ≈ 0.0125 if perfectly equally absorbing,
  more typically a sparse subset gets the bulk)
- 2 non-GrENE founders in the panel should get near-zero weight

This is a *partial* test: the simulation pool covers founders not in our panel,
so we don't expect exact uniform recovery — we expect a sensible smear over the
panel members that are population-similar to the unobserved 151.

Pipeline tested:
- block_solver.solve_block_irls
- kmer_count.count_kmers_in_bam
- build_kmer_cn output (loaded from data/test_chr1_first200.*)
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.sparse import load_npz
from block_solver import solve_block_irls, solve_block_wls
from kmer_count import count_kmers_in_bam, estimate_coverage


def main():
    DATA = os.path.join(os.path.dirname(__file__), "..", "data")
    print("="*70)
    print("Real-BAM end-to-end test")
    print("="*70)

    # 1. Load cn matrix
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.cn.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    kmer_index = meta["kmer_index"]
    bubble_id = meta["bubble_id"]
    founders = meta["founders"]
    F, K = cn_sparse.shape
    cn = np.asarray(cn_sparse.todense()).astype(np.int8)
    print(f"\n[1] Loaded cn: F={F} founders × K={K:,} k-mers")
    print(f"    AC distribution: median={int(np.median(cn.sum(axis=0)))}  "
          f"max={cn.sum(axis=0).max()}  AC=1: {(cn.sum(axis=0)==1).sum()}")

    # 2. Pick a sim BAM (cov50, all 231 ecotypes uniformly, with a 1kb DEL at f=90%)
    bam = "/carnegie/nobackup/scratch/tbellagio/visor_freqk/data/reads_var/del/rep29/cov50/var_del_1kb_n231_f90_err001/sim.srt.bam"
    assert os.path.exists(bam), f"BAM not found: {bam}"
    print(f"\n[2] BAM: {bam}")

    # 3. Estimate per-base coverage and convert to k-mer coverage
    print(f"\n[3] Estimating coverage...")
    t0 = time.time()
    cov_per_base = estimate_coverage(bam)
    # k-mer coverage: each base contributes (read_length - k + 1)/read_length kmers
    # for 150bp reads, k=31: factor is 0.787
    cov_kmer = cov_per_base * (150 - 31 + 1) / 150
    print(f"    per-base coverage: {cov_per_base:.1f}x")
    print(f"    k-mer (k=31)  coverage: {cov_kmer:.1f}x")
    print(f"    [took {time.time()-t0:.1f}s]")

    # 4. Count k-mers from BAM
    print(f"\n[4] Counting {K:,} unique k-mers in BAM...")
    t0 = time.time()
    counts_dict = count_kmers_in_bam(bam, list(kmer_index), k=31, threads=4, hash_size="500M")
    counts = np.array([counts_dict[km] for km in kmer_index], dtype=np.int64)
    print(f"    [took {time.time()-t0:.1f}s]")
    print(f"    total k-mer hits: {counts.sum():,}")
    print(f"    nonzero k-mers: {(counts>0).sum():,}/{K:,}  ({(counts>0).mean():.1%})")
    print(f"    median count where nonzero: {np.median(counts[counts>0]):.1f}")
    print(f"    max count: {counts.max()}")

    # 5. Solve block: partition by bubble_id into "blocks" of contiguous bubbles
    # For this test, treat all 200 bubbles as a single "block" (Chr1 first 200)
    print(f"\n[5] Solving for h (block = all 200 bubbles)...")
    t0 = time.time()
    h_irls, obj_irls = solve_block_irls(counts, cn, coverage=cov_kmer,
                                         max_iter=5, tol=1e-4, verbose=True)
    print(f"    [took {time.time()-t0:.1f}s, obj={obj_irls:.2f}]")
    print(f"    h sum: {h_irls.sum():.4f}  (should be 1.0)")
    print(f"    h min: {h_irls.min():.4f}  max: {h_irls.max():.4f}")
    print(f"    h support (n founders >0.001): {(h_irls > 0.001).sum()}/{F}")

    # 6. Inspect distribution
    print(f"\n[6] Founder weight distribution:")
    sorted_idx = np.argsort(-h_irls)
    print(f"    Top 10 founders by inferred h:")
    for i, f in enumerate(sorted_idx[:10]):
        marker = ""
        # Mark non-GrENE-Net founders if known (from our cross-reference)
        # Simple: read founder list from the meta; everything starting with 'TAIR10' or
        # specific known non-GrENE names is non-GrENE
        print(f"      #{i+1:2d}: {founders[f]:>10s}   h={h_irls[f]:.4f}")

    # Distribution stats
    print(f"\n    Distribution stats:")
    print(f"      h_mean (over all 82): {h_irls.mean():.4f}  (uniform = 1/82 = {1/82:.4f})")
    print(f"      h_median:             {np.median(h_irls):.4f}")
    print(f"      effective n founders (1/Σh²): {1/np.sum(h_irls**2):.1f}")

    # 7. Sanity check: predicted vs observed counts
    pred = cov_kmer * (cn.T @ h_irls)
    nonzero = counts > 0
    if nonzero.any():
        rmse = np.sqrt(np.mean((counts[nonzero] - pred[nonzero])**2))
        # poisson noise expected: sqrt(mean count)
        expected_noise = np.sqrt(counts[nonzero].mean())
        print(f"\n[7] Goodness-of-fit:")
        print(f"    residual RMSE (nonzero kmers): {rmse:.2f}")
        print(f"    Poisson noise floor (sqrt(mean count)): {expected_noise:.2f}")
        print(f"    ratio: {rmse/expected_noise:.2f}  (1.0 = perfect Poisson fit)")


if __name__ == "__main__":
    main()
