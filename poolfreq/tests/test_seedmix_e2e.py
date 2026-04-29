"""
SEEDMIX end-to-end test.

SEEDMIX is the GrENE-Net founder pool, sequenced as 8 replicates. The
recipe gives expected proportions for each of the 231 founders.

For the 80 founders that are in our panel:
  expected h[f] = recipe_prop[f] / sum_over_panel_80(recipe_prop)
  (because our 82-founder simplex normalizes against panel mass only)

We expect:
- h sums to 1 (simplex constraint)
- Non-zero weight spread across the 80 GrENE-overlap panel founders
- Near-zero weight on the 2 non-GrENE-Net panel founders
- Order-of-magnitude correlation with recipe (recipe r > 0)

This is a meaningful real-data test, unlike the narrow visor_freqk simulation.
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from block_solver import solve_block_irls, solve_block_wls
from kmer_count import count_kmers_in_fasta


def main():
    DATA = os.path.join(os.path.dirname(__file__), "..", "data")
    print("="*70)
    print("SEEDMIX end-to-end test")
    print("="*70)

    # 1. Load cn matrix (200 bubbles, Chr1 first 100kb)
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.cn.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    kmer_index = meta["kmer_index"]
    founders = meta["founders"]
    F, K = cn_sparse.shape
    cn = np.asarray(cn_sparse.todense()).astype(np.int8)
    print(f"\n[1] Loaded cn: F={F}, K={K:,}")

    # 2. Map our 82 panel founders to 1001G IDs and check GrENE-overlap
    panel_map = pd.read_csv(
        "/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv",
        sep="\t",
    )
    asm_to_1001g = dict(zip(panel_map.Assembly_ID.astype(str), panel_map.Accession_ID.astype(str)))
    panel_1001g_ids = [asm_to_1001g.get(str(f), None) for f in founders]
    grenenet = set(open("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/vcf_samples_231.txt").read().split())
    is_grenenet = np.array([fid in grenenet if fid else False for fid in panel_1001g_ids])
    print(f"    {is_grenenet.sum()} of {F} panel founders are GrENE-Net 231 members")

    # 3. SEEDMIX recipe
    recipe = pd.read_csv(
        "/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/seedmix_recipe_normalized.tsv",
        sep="\t",
    )
    recipe_dict = dict(zip(recipe.ID.astype(str), recipe.seed_prop))
    expected_panel_prop = np.array([
        recipe_dict.get(fid, 0.0) if fid else 0.0 for fid in panel_1001g_ids
    ])
    panel_total_mass = expected_panel_prop.sum()
    if panel_total_mass > 0:
        expected_h = expected_panel_prop / panel_total_mass
    else:
        expected_h = expected_panel_prop
    print(f"    panel total seedmix mass: {panel_total_mass:.3f} (= {panel_total_mass*100:.1f}%)")
    print(f"    expected effective n founders (1/Σh²): {1/np.sum(expected_h**2):.1f}")

    # 4. Find SEEDMIX_S1 fastq pair
    SEED_DIR = "/home/tbellagio/scratch/pang/grenenet_reads/seed_mix"
    fq1 = os.path.join(SEED_DIR, "S1-1.1_P.fq.gz")
    fq2 = os.path.join(SEED_DIR, "S1-1.2_P.fq.gz")
    assert os.path.exists(fq1) and os.path.exists(fq2)
    print(f"\n[2] FASTQs: {os.path.basename(fq1)}, {os.path.basename(fq2)}")

    # 5. Estimate k-mer coverage from total bases / genome size
    # Use samtools-style estimate: for each fastq, count reads × read_length
    from subprocess import run
    bgzipped_size = sum(os.path.getsize(f) for f in [fq1, fq2])
    print(f"    fq.gz size: {bgzipped_size/1e9:.2f} GB")
    # Rough coverage estimate: gzip ratio ~3, ~150bp reads, genome 119Mb
    approx_cov_per_base = 3.0 * bgzipped_size / 119_146_348 / 4  # 4 lines per read
    cov_kmer = approx_cov_per_base * (150 - 31 + 1) / 150
    print(f"    rough k-mer coverage estimate: {cov_kmer:.1f}x")

    # 6. Count k-mers
    print(f"\n[3] Counting {K:,} k-mers in SEEDMIX_S1 fastqs...")
    t0 = time.time()
    counts_dict = count_kmers_in_fasta([fq1, fq2], list(kmer_index),
                                        k=31, threads=4, hash_size="2G")
    counts = np.array([counts_dict[km] for km in kmer_index], dtype=np.int64)
    print(f"    [took {time.time()-t0:.0f}s]")
    print(f"    nonzero k-mers: {(counts>0).sum():,}/{K:,} ({(counts>0).mean():.1%})")
    print(f"    median count where nonzero: {np.median(counts[counts>0]):.1f}")
    print(f"    total k-mer hits: {counts.sum():,}")

    # 7. WLS solve (more stable than IRLS for this system)
    print(f"\n[4] WLS solve...")
    t0 = time.time()
    h_wls, obj = solve_block_wls(counts, cn, coverage=cov_kmer)
    print(f"    [took {time.time()-t0:.0f}s, obj={obj:.0f}]")
    print(f"    h sum: {h_wls.sum():.4f}")
    print(f"    h support (>0.001): {(h_wls > 0.001).sum()}/{F}")
    print(f"    effective n founders (1/Σh²): {1/np.sum(h_wls**2):.1f}")

    # 8. Compare to expected (recipe-based)
    print(f"\n[5] Comparison to recipe-based expectation")
    df = pd.DataFrame({
        "founder": founders,
        "1001g_id": panel_1001g_ids,
        "is_grenenet": is_grenenet,
        "expected_h": expected_h,
        "h_wls": h_wls,
    })

    # Pearson r over the 80 GrENE-overlap founders
    r = np.corrcoef(df[df.is_grenenet]["expected_h"], df[df.is_grenenet]["h_wls"])[0, 1]
    print(f"    Pearson r (h_wls vs expected_h, over 80 panel founders): {r:.3f}")

    # Top panel founders by recipe weight
    top_recipe = df.sort_values("expected_h", ascending=False).head(10)
    print(f"\n    Top 10 by recipe (with h_wls):")
    for _, row in top_recipe.iterrows():
        marker = "✓" if row["is_grenenet"] else " "
        print(f"      {marker} {row['founder']:>8s} ({row['1001g_id']}): expected={row['expected_h']:.4f}  h_wls={row['h_wls']:.4f}")

    # Top inferred founders
    print(f"\n    Top 10 by h_wls:")
    top_h = df.sort_values("h_wls", ascending=False).head(10)
    for _, row in top_h.iterrows():
        marker = "✓" if row["is_grenenet"] else "✗"
        print(f"      {marker} {row['founder']:>8s} ({row['1001g_id']}): h_wls={row['h_wls']:.4f}  expected={row['expected_h']:.4f}")

    # Mass on non-GrENE panel founders (should be ~0 if pipeline is correct)
    non_grene_mass = df[~df.is_grenenet]["h_wls"].sum()
    print(f"\n    Mass on 2 non-GrENE-Net panel founders: {non_grene_mass:.4f} (expected ~0)")

    # save for inspection
    out_csv = os.path.join(DATA, "seedmix_S1_h_wls_chr1_first200.tsv")
    df.to_csv(out_csv, sep="\t", index=False)
    print(f"\n    Wrote {out_csv}")


if __name__ == "__main__":
    main()
