"""
Comprehensive validation on the 2000-bubble kmer_pa matrix.

For each pool (uniform sim, skewed sim, SEEDMIX_S1):
  - Count k-mers against the 181K query set
  - Run baseline WLS solver
  - Compare recovered h to truth
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
import cvxpy as cp
from scipy.sparse import load_npz
from kmer_count import count_kmers_in_fasta


DATA = os.path.join(os.path.dirname(__file__), "..", "data")
_PROJ_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def solve_wls(counts, kmer_pa, coverage):
    K, F = counts.shape[0], kmer_pa.shape[0]
    weights = 1.0 / np.sqrt(counts + 1.0)
    h = cp.Variable(F, nonneg=True)
    mu = coverage * (kmer_pa.T @ h)
    residual = cp.multiply(weights, counts - mu)
    prob = cp.Problem(cp.Minimize(cp.sum_squares(residual)), [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def solve_kl(counts, kmer_pa, coverage):
    F = kmer_pa.shape[0]
    h = cp.Variable(F, nonneg=True)
    mu = coverage * (kmer_pa.T @ h) + 1e-3
    obj = cp.Minimize(cp.sum(cp.kl_div(counts, mu)))
    prob = cp.Problem(obj, [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def get_or_count(pool_name, fastqs, kmer_index):
    """Count k-mers for a pool, with cache."""
    cache_path = os.path.join(DATA, f"{pool_name}_counts_2k.npz")
    if os.path.exists(cache_path):
        return np.load(cache_path)["counts"]
    print(f"  Counting k-mers for {pool_name}...")
    t0 = time.time()
    cd = count_kmers_in_fasta(fastqs, list(kmer_index), k=31, threads=8, hash_size="3G")
    counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
    np.savez(cache_path, counts=counts)
    print(f"  [{time.time()-t0:.0f}s, nonzero {(counts>0).sum():,}/{len(counts):,}]")
    return counts


def evaluate(label, h_hat, h_true):
    h_hat = np.asarray(h_hat) / np.asarray(h_hat).sum()
    err = np.linalg.norm(h_hat - h_true)
    if h_true.std() > 1e-9:
        r = np.corrcoef(h_hat, h_true)[0, 1]
    else:
        r = float("nan")
    cv = h_hat.std() / max(1e-9, h_hat.mean())
    nz = (h_hat < 1e-4).sum()
    print(f"    {label:35s}  ||Δh||={err:.4f}  r={r:.3f}  CV={cv:.2f}  n_zero={nz}")
    return h_hat, err, r


def main():
    # Load 2000-bubble kmer_pa
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first2000.kmer_pa.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first2000.meta.npz"), allow_pickle=True)
    kmer_index = meta["kmer_index"]
    bubble_id = meta["bubble_id"]
    founders = meta["founders"]
    F, K = cn_sparse.shape
    kmer_pa = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = kmer_pa.sum(axis=0)
    print(f"="*72)
    print(f"2000-bubble kmer_pa validation: F={F}, K={K:,}")
    print(f"="*72)
    print(f"AC dist: AC=1: {(ac==1).sum():,}   AC 2-4: {((ac>=2)&(ac<=4)).sum():,}   "
          f"AC 5-10: {((ac>=5)&(ac<=10)).sum():,}   AC>10: {(ac>10).sum():,}")

    # Helpers for truth lookup
    panel_map = pd.read_csv(
        os.path.join(_PROJ_ROOT, "data/sv_panel_to_accession_id.tsv"),
        sep="\t",
    )
    asm_to_1001g = dict(zip(panel_map.Assembly_ID.astype(str), panel_map.Accession_ID.astype(str)))
    grenenet = set(open(os.path.join(_PROJ_ROOT, "data/vcf_samples_231.txt")).read().split())
    recipe = pd.read_csv(
        os.path.join(_PROJ_ROOT, "data/seedmix_recipe_normalized.tsv"), sep="\t"
    )
    recipe_dict = dict(zip(recipe.ID.astype(str), recipe.seed_prop))

    # ==================== Pool 1: UNIFORM SIM ====================
    print(f"\n{'-'*72}\n[1] UNIFORM SIM (82 founders at 1/82)\n{'-'*72}")
    truth = pd.read_csv(os.path.join(DATA, "sim_chr1", "uniform82_truth.tsv"), sep="\t")
    truth_dict = dict(zip(truth.founder.astype(str), truth.weight))
    h_true = np.array([truth_dict.get(str(f), 0.0) for f in founders])
    h_true = h_true / h_true.sum()

    counts = get_or_count("uniform82",
                           [os.path.join(DATA, "sim_chr1", "uniform82_pool_1.fq.gz"),
                            os.path.join(DATA, "sim_chr1", "uniform82_pool_2.fq.gz")],
                           kmer_index)
    cov = counts.sum() * F / ac.sum()
    print(f"  cov={cov:.1f}×, nonzero kmers: {(counts>0).sum():,}")
    h_wls = solve_wls(counts, kmer_pa, cov)
    evaluate("WLS", h_wls, h_true)
    # for KL on this big matrix, may take long - try
    print("  (KL solve may take a few minutes...)")
    t = time.time()
    h_kl = solve_kl(counts, kmer_pa, cov)
    print(f"  KL [{time.time()-t:.0f}s]")
    evaluate("KL  (Poisson)", h_kl, h_true)

    # ==================== Pool 2: SKEWED SIM ====================
    print(f"\n{'-'*72}\n[2] SKEWED SIM (5 founders at [0.40, 0.25, 0.15, 0.10, 0.10])\n{'-'*72}")
    truth = pd.read_csv(os.path.join(DATA, "sim_chr1_skewed", "skewed5_truth.tsv"), sep="\t")
    truth_dict = dict(zip(truth.founder.astype(str), truth.weight))
    h_true = np.array([truth_dict.get(str(f), 0.0) for f in founders])
    h_true = h_true / h_true.sum()

    counts = get_or_count("skewed5",
                           [os.path.join(DATA, "sim_chr1_skewed", "skewed5_pool_1.fq.gz"),
                            os.path.join(DATA, "sim_chr1_skewed", "skewed5_pool_2.fq.gz")],
                           kmer_index)
    cov = counts.sum() * F / ac.sum()
    print(f"  cov={cov:.1f}×, nonzero kmers: {(counts>0).sum():,}")
    h_wls = solve_wls(counts, kmer_pa, cov)
    evaluate("WLS", h_wls, h_true)
    h_kl = solve_kl(counts, kmer_pa, cov)
    evaluate("KL  (Poisson)", h_kl, h_true)
    print(f"  Top 5 by WLS:")
    for i in np.argsort(-h_wls/h_wls.sum())[:5]:
        h_norm = h_wls[i]/h_wls.sum()
        marker = "✓" if h_true[i] > 0 else " "
        print(f"    {marker} {founders[i]}: h={h_norm:.4f} truth={h_true[i]:.4f}")

    # ==================== Pool 3: SEEDMIX_S1 ====================
    print(f"\n{'-'*72}\n[3] SEEDMIX_S1 (real data, recipe-based truth)\n{'-'*72}")
    panel_1001g_ids = [asm_to_1001g.get(str(f), None) for f in founders]
    is_grenenet = np.array([fid in grenenet if fid else False for fid in panel_1001g_ids])
    expected_panel_prop = np.array([
        recipe_dict.get(fid, 0.0) if fid else 0.0 for fid in panel_1001g_ids
    ])
    panel_total = expected_panel_prop.sum()
    h_true = expected_panel_prop / panel_total if panel_total > 0 else expected_panel_prop
    print(f"  recipe panel mass: {panel_total*100:.1f}%, expected effective n: {1/np.sum(h_true**2):.1f}")

    counts = get_or_count("seedmix_S1",
                           ["/global/home/users/tbellg/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz",
                            "/global/home/users/tbellg/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz"],
                           kmer_index)
    cov = counts.sum() * F / ac.sum()
    print(f"  cov={cov:.1f}×, nonzero kmers: {(counts>0).sum():,}")
    h_wls = solve_wls(counts, kmer_pa, cov)
    evaluate("WLS (vs recipe)", h_wls, h_true)
    h_kl = solve_kl(counts, kmer_pa, cov)
    evaluate("KL  (vs recipe)", h_kl, h_true)
    # Mass on non-GrENE founders should be ~0
    h_wls_norm = h_wls / h_wls.sum()
    non_grene_mass = h_wls_norm[~is_grenenet].sum()
    print(f"  Mass on 2 non-GrENE-Net founders: {non_grene_mass:.4f} (expected ~0)")
    print(f"  Top 8 by WLS (✓ = in panel and in recipe):")
    for i in np.argsort(-h_wls_norm)[:8]:
        marker = "✓" if (is_grenenet[i] and h_true[i] > 0) else " "
        print(f"    {marker} {founders[i]}: h={h_wls_norm[i]:.4f} recipe={h_true[i]:.4f}")


if __name__ == "__main__":
    main()
