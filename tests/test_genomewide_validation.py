"""
Genome-wide validation: concatenate per-chromosome kmer_pa matrices, count k-mers
in each pool against the full set, run EM, compare to truth.

Pools:
  - uniform82 sim (all 82 founders at 1/82, 30× cov)
  - skewed5 sim (5 founders at [0.40, 0.25, 0.15, 0.10, 0.10], 30× cov)
  - SEEDMIX_S1 (real data, recipe-based truth on the 80 panel founders)

Solver: flat EM only. WLS via CVXPY/SCS does not scale to K~80M k-mers
(constructs an O(K) sparse linear system that exceeds 100 GB just for problem
data). EM is the production solver; per_sample_driver uses it for all 2,415
samples.

Reports R² and RMSE per the user's preference.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
from scipy.sparse import load_npz, hstack
from em_solver import solve_em
from kmer_count import count_kmers_in_fasta


DATA = os.path.join(os.path.dirname(__file__), "..", "data")
_PROJ_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def metrics(h_hat, h_true):
    h_hat = np.asarray(h_hat)
    if h_hat.sum() > 0:
        h_hat = h_hat / h_hat.sum()
    rmse = float(np.sqrt(np.mean((h_hat - h_true)**2)))
    if h_true.std() > 1e-9:
        ss_res = ((h_true - h_hat)**2).sum()
        ss_tot = ((h_true - h_true.mean())**2).sum()
        r2 = float(1 - ss_res / ss_tot)
    else:
        r2 = float("nan")
    return r2, rmse


def load_genomewide_cn():
    """Concatenate per-chromosome kmer_pa matrices (built by 56188)."""
    chroms = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
    cn_blocks = []
    kmer_indices = []
    bubble_ids = []
    bubble_offset = 0

    for c in chroms:
        cn_path = os.path.join(DATA, f"kmer_pa_{c}.kmer_pa.npz")
        meta_path = os.path.join(DATA, f"kmer_pa_{c}.meta.npz")
        if not os.path.exists(cn_path):
            print(f"  WARN: {c} kmer_pa missing, skipping")
            continue
        cn_b = load_npz(cn_path)
        meta = np.load(meta_path, allow_pickle=True)
        cn_blocks.append(cn_b)
        kmer_indices.append(meta["kmer_index"])
        bubble_ids.append(meta["bubble_id"] + bubble_offset)
        n_bubbles = int(meta["bubble_id"].max()) + 1
        bubble_offset += n_bubbles
        print(f"  {c}: kmer_pa shape={cn_b.shape}, bubbles={n_bubbles}, "
              f"kmers={cn_b.shape[1]:,}")

    if not cn_blocks:
        raise RuntimeError("No per-chrom kmer_pa files found. 56188 may not be done yet.")

    kmer_pa = hstack(cn_blocks, format="csr")
    kmer_index = np.concatenate(kmer_indices)
    bubble_id = np.concatenate(bubble_ids)
    print(f"\nFull kmer_pa shape: {kmer_pa.shape}, total bubbles: {int(bubble_id.max()) + 1}")
    return kmer_pa, kmer_index, bubble_id, np.load(
        os.path.join(DATA, "kmer_pa_Chr1.meta.npz"), allow_pickle=True)["founders"]


def get_or_count(pool_name, fastqs, kmer_index):
    cache_path = os.path.join(DATA, f"{pool_name}_counts_genomewide.npz")
    if os.path.exists(cache_path):
        cached = np.load(cache_path)
        # Verify dimensions match
        if len(cached["counts"]) == len(kmer_index):
            return cached["counts"]
        print(f"  cached counts for {pool_name} have different dimensions, re-counting")
    print(f"  Counting k-mers for {pool_name} ({len(kmer_index):,} queries, may take ~10 min)...")
    t = time.time()
    cd = count_kmers_in_fasta(fastqs, list(kmer_index), k=31, threads=8, hash_size="4G")
    counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
    np.savez(cache_path, counts=counts)
    print(f"    [{time.time()-t:.0f}s, nonzero {(counts>0).sum():,}/{len(counts):,}]")
    return counts


def evaluate_pool(label, h_true, fastqs, cn_f32, ac, kmer_index, F, K, pool_name):
    print(f"\n{'='*72}", flush=True)
    print(f"[{label}]", flush=True)
    print(f"{'='*72}", flush=True)
    counts = get_or_count(pool_name, fastqs, kmer_index)
    cov = counts.sum() * F / ac.sum()
    print(f"  coverage estimate: {cov:.1f}×, nonzero kmers: {(counts>0).sum():,}/{K:,}",
          flush=True)

    print(f"\n  {'method':<8} {'R²':>10} {'RMSE':>10} {'time':>8}", flush=True)
    print("  " + "-"*40, flush=True)

    t = time.time()
    h_em, _ = solve_em(counts, cn_f32, cov, max_iter=300, tol=1e-7)
    r2, rmse = metrics(h_em, h_true)
    print(f"  {'EM':<8} {r2:>10.4f} {rmse:>10.4f} {time.time()-t:>7.1f}s", flush=True)
    return h_em


def main():
    import gc
    print(f"="*72, flush=True)
    print(f"Genome-wide validation", flush=True)
    print(f"="*72, flush=True)

    kmer_pa, kmer_index, bubble_id, founders = load_genomewide_cn()
    F, K = kmer_pa.shape

    # Densify and cast to float32 ONCE; reused across all pools. Keeping the
    # int8 alive while EM internally casts to float32 was the source of the
    # 56202 segfault — peak ~2× memory per pool. Now: single 26 GB allocation.
    print(f"\n[{time.strftime('%H:%M:%S')}] Densifying + casting kmer_pa → float32 "
          f"({F*K*4/1e9:.1f} GB)...", flush=True)
    cn_f32 = np.asarray(kmer_pa.todense()).astype(np.float32)
    del kmer_pa
    gc.collect()
    ac = cn_f32.sum(axis=0)
    print(f"  AC dist: AC=1: {(ac==1).sum():,}  AC 2-4: {((ac>=2)&(ac<=4)).sum():,}  "
          f"AC 5-10: {((ac>=5)&(ac<=10)).sum():,}  AC>10: {(ac>10).sum():,}", flush=True)

    # 1) UNIFORM SIM
    truth = pd.read_csv(os.path.join(DATA, "sim_chr1", "uniform82_truth.tsv"), sep="\t")
    truth_dict = dict(zip(truth.founder.astype(str), truth.weight))
    h_true = np.array([truth_dict.get(str(f), 0.0) for f in founders])
    h_true = h_true / h_true.sum()
    evaluate_pool("UNIFORM 82", h_true,
                   [os.path.join(DATA, "sim_chr1", "uniform82_pool_1.fq.gz"),
                    os.path.join(DATA, "sim_chr1", "uniform82_pool_2.fq.gz")],
                   cn_f32, ac, kmer_index, F, K, "uniform82")

    # 2) SKEWED SIM
    truth = pd.read_csv(os.path.join(DATA, "sim_chr1_skewed", "skewed5_truth.tsv"), sep="\t")
    truth_dict = dict(zip(truth.founder.astype(str), truth.weight))
    h_true = np.array([truth_dict.get(str(f), 0.0) for f in founders])
    h_true = h_true / h_true.sum()
    evaluate_pool("SKEWED 5", h_true,
                   [os.path.join(DATA, "sim_chr1_skewed", "skewed5_pool_1.fq.gz"),
                    os.path.join(DATA, "sim_chr1_skewed", "skewed5_pool_2.fq.gz")],
                   cn_f32, ac, kmer_index, F, K, "skewed5")

    # 3) SEEDMIX_S1 (recipe-restricted truth)
    panel_map = pd.read_csv(
        os.path.join(_PROJ_ROOT, "data/sv_panel_to_accession_id.tsv"),
        sep="\t")
    asm_to_1001g = dict(zip(panel_map.Assembly_ID.astype(str),
                            panel_map.Accession_ID.astype(str)))
    recipe = pd.read_csv(
        os.path.join(_PROJ_ROOT, "data/seedmix_recipe_normalized.tsv"),
        sep="\t")
    recipe_dict = dict(zip(recipe.ID.astype(str), recipe.seed_prop))
    panel_1001g_ids = [asm_to_1001g.get(str(f), None) for f in founders]
    expected = np.array([
        recipe_dict.get(fid, 0.0) if fid else 0.0 for fid in panel_1001g_ids
    ])
    panel_total = expected.sum()
    h_true = expected / panel_total if panel_total > 0 else expected
    print(f"\n  recipe panel mass: {panel_total*100:.1f}%, "
          f"effective n on panel: {1/np.sum(h_true**2):.1f}", flush=True)
    evaluate_pool("SEEDMIX_S1", h_true,
                   ["/global/home/users/tbellg/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz",
                    "/global/home/users/tbellg/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz"],
                   cn_f32, ac, kmer_index, F, K, "seedmix_S1")


if __name__ == "__main__":
    main()
