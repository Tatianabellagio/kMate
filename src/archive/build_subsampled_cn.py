"""
Per-(founder, ac-stratum) k-mer subsampling.

For each ac-stratum (bin of k-mer carrier count), subsample each founder's
carried k-mers down to a target count derived from the panel distribution
(default: median across founders for that stratum). Equalizes the rarity
distribution of carried k-mers across founders without amplifying per-cell
EM contribution (cn stays binary 0/1).

Post-subsample, k-mer ac_k drops (some doubletons become singletons). The
--refilt2 flag drops new ac<2 columns to keep the singleton-free invariant.

Reproducible via --seed.
"""
from __future__ import annotations
import argparse, os, time
import numpy as np
from pathlib import Path
from scipy.sparse import load_npz, save_npz, csr_matrix


def stratify_subsample(cn_csr, ac_k, bin_edges, target='median',
                        seed=42, verbose=True, protect_bins=0,
                        founder_class=None, match_ref=None):
    """
    cn_csr:    F × K sparse binary
    ac_k:      K-vec, panel carrier count per k-mer
    bin_edges: e.g. [0, 1, 2, 3, 5, 11, 26, 51, 101, 201, 232]
               → bins: (1,)=ac=1, (2,)=ac=2, (3,4)=ac=3-4, (5,...,10)=ac=5-10, ...
    target:    'median' | 'min' | 'p25'
    """
    rng = np.random.default_rng(seed)
    F, K = cn_csr.shape
    n_bins = len(bin_edges) - 1
    ac_bin = np.digitize(ac_k, bin_edges) - 1  # bin idx per k-mer (clamp to valid)
    ac_bin = np.clip(ac_bin, 0, n_bins - 1)
    if verbose:
        print(f'  bins: {n_bins}, edges: {bin_edges}')
        for b in range(n_bins):
            m = ac_bin == b
            print(f'    bin {b} ac∈[{bin_edges[b]},{bin_edges[b+1]-1}]: {m.sum():,} k-mers')

    # Per (founder, bin) count of carried k-mers
    if verbose: print(f'  computing per-(founder,bin) carriage...')
    t = time.time()
    counts = np.zeros((F, n_bins), dtype=np.int64)
    for f in range(F):
        s, e = cn_csr.indptr[f], cn_csr.indptr[f+1]
        cols = cn_csr.indices[s:e]
        if len(cols) == 0: continue
        b = ac_bin[cols]
        for bi in range(n_bins):
            counts[f, bi] = (b == bi).sum()
    if verbose: print(f'    [{time.time()-t:.0f}s]')

    # Target per bin
    if target == 'median':
        target_counts = np.median(counts, axis=0).astype(np.int64)
    elif target == 'min':
        target_counts = counts.min(axis=0)
    elif target == 'p25':
        target_counts = np.quantile(counts, 0.25, axis=0).astype(np.int64)
    else:
        raise ValueError(target)
    # Protect the lowest-AC bins (the discriminating/private k-mers): keep ALL
    # of them, subsample only the higher-AC (shared) bins. Capping the rare bins
    # to median is what costs founder-identifiability (→ leakage on absent
    # founders); the cactus class inflation lives in the shared/high-AC bins.
    if protect_bins > 0:
        target_counts[:protect_bins] = np.iinfo(np.int64).max
        if verbose:
            print(f'  PROTECTING bins 0..{protect_bins-1} (keep all, no subsample)')

    # Per-founder per-bin target. Default: broadcast the panel-wide target_counts.
    BIG = np.iinfo(np.int64).max
    per_founder_target = np.tile(target_counts, (F, 1))  # (F, n_bins)
    if match_ref is not None:
        # CLASS-MATCH: bring the non-reference class DOWN to the reference class's
        # per-bin median; keep every reference-class k-mer. This removes the
        # cactus-vs-PG data-source artifact (long-read founders carry excess
        # low-AC tags) without touching the reference (PG) founders — making the
        # panel balanced like a uniformly-genotyped SNP panel (hapFIRE-style).
        is_ref = np.array([c == match_ref for c in founder_class])
        ref_med = np.median(counts[is_ref], axis=0).astype(np.int64)  # per-bin PG median
        per_founder_target[is_ref, :] = BIG          # keep all ref-class tags
        per_founder_target[~is_ref, :] = ref_med[None, :]  # match others to ref median
        if verbose:
            print(f'  CLASS-MATCH to ref={match_ref}: ref n={is_ref.sum()}, '
                  f'other n={(~is_ref).sum()}; per-bin ref median target = {list(ref_med)}')
    if verbose:
        print(f'  per-bin counts: panel median per founder: {list(np.median(counts,axis=0).astype(int))}')
        print(f'  TARGET per bin: {list(target_counts)}')
        print(f'  panel max per bin:    {list(counts.max(axis=0))}')

    # Build new sparse cn: per founder, subsample within each bin
    if verbose: print(f'  subsampling per founder...')
    t = time.time()
    keep_rows = []
    keep_cols = []
    for f in range(F):
        s, e = cn_csr.indptr[f], cn_csr.indptr[f+1]
        cols = cn_csr.indices[s:e]
        if len(cols) == 0: continue
        b = ac_bin[cols]
        for bi in range(n_bins):
            in_bin = cols[b == bi]
            tgt = per_founder_target[f, bi]
            if len(in_bin) <= tgt:
                keep = in_bin
            else:
                keep = rng.choice(in_bin, size=tgt, replace=False)
            if len(keep) > 0:
                keep_rows.append(np.full(len(keep), f, dtype=np.int32))
                keep_cols.append(keep.astype(np.int32))
    if verbose: print(f'    [{time.time()-t:.0f}s]')

    rows = np.concatenate(keep_rows) if keep_rows else np.array([], dtype=np.int32)
    cols = np.concatenate(keep_cols) if keep_cols else np.array([], dtype=np.int32)
    new_cn = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                         shape=(F, K), dtype=np.int8)
    if verbose:
        print(f'  new cn: nnz={new_cn.nnz:,} (was {cn_csr.nnz:,}, drop {100*(1-new_cn.nnz/cn_csr.nnz):.1f}%)')
    return new_cn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-cn", required=True)
    ap.add_argument("--in-meta", required=True)
    ap.add_argument("--out-prefix", required=True,
                    help="writes {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--target", default="median", choices=["median","min","p25"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--refilt2", action="store_true",
                    help="drop columns with new ac_k < 2 after subsampling")
    ap.add_argument("--protect-bins", type=int, default=0,
                    help="keep ALL k-mers in the first N (lowest-AC) bins; only "
                         "subsample higher-AC bins to target. Preserves discriminating "
                         "rare k-mers while equalizing the shared-k-mer class inflation.")
    ap.add_argument("--class-match-json", default=None,
                    help="founder_split_cactus_pg.json; enables CLASS-MATCH mode")
    ap.add_argument("--match-ref", default="PG",
                    help="reference class to keep intact; the other class is "
                         "subsampled down to the reference class per-bin median "
                         "(default PG: cut cactus excess down to PG level).")
    args = ap.parse_args()

    print(f"[subsample] in_cn={args.in_cn}", flush=True)
    cn = load_npz(args.in_cn).tocsr()
    meta = np.load(args.in_meta, allow_pickle=True)
    F, K = cn.shape
    print(f"  shape F={F} K={K:,} nnz={cn.nnz:,}", flush=True)

    ac_k = np.asarray(cn.sum(axis=0)).flatten().astype(np.int32)
    print(f"  ac_k: min={ac_k.min()} max={ac_k.max()}, ac=1: {(ac_k==1).sum():,}, "
          f"ac=2: {(ac_k==2).sum():,}, ac>=10: {(ac_k>=10).sum():,}", flush=True)

    # Bins: handle ac=1 separately if present
    if ac_k.min() == 1:
        bin_edges = [1, 2, 3, 5, 11, 26, 51, 101, 201, 232]
    else:
        bin_edges = [2, 3, 5, 11, 26, 51, 101, 201, 232]

    founder_class = None
    if args.class_match_json:
        import json
        split = json.load(open(args.class_match_json))
        f2c = {}
        for cls, ids in split.items():
            for i in ids: f2c[str(i)] = cls
        founders_arr = np.asarray(meta["founders"]).astype(str)
        founder_class = [f2c.get(f, "UNK") for f in founders_arr]
        print(f"  class-match: {founders_arr.size} founders, "
              f"classes={ {c: founder_class.count(c) for c in set(founder_class)} }", flush=True)

    new_cn = stratify_subsample(cn, ac_k, bin_edges,
                                  target=args.target, seed=args.seed,
                                  protect_bins=args.protect_bins,
                                  founder_class=founder_class,
                                  match_ref=(args.match_ref if args.class_match_json else None))
    new_ac = np.asarray(new_cn.sum(axis=0)).flatten().astype(np.int32)
    print(f"  post-subsample ac_k: min={new_ac.min()} max={new_ac.max()}, "
          f"ac=0: {(new_ac==0).sum():,}, ac=1: {(new_ac==1).sum():,}, "
          f"ac>=2: {(new_ac>=2).sum():,}", flush=True)

    if args.refilt2:
        keep = (new_ac >= 2)
        print(f"  --refilt2: keeping {keep.sum():,}/{K:,} cols (drop {100*(1-keep.sum()/K):.1f}%)",
              flush=True)
        new_cn = new_cn[:, keep]
        kmer_index = np.asarray(meta['kmer_index'])[keep]
        bubble_id = np.asarray(meta['bubble_id'])[keep]
    else:
        # keep all cols but mark ac=0 as "empty"; downstream still uses cols
        keep_all = (new_ac > 0)
        new_cn = new_cn[:, keep_all]
        kmer_index = np.asarray(meta['kmer_index'])[keep_all]
        bubble_id = np.asarray(meta['bubble_id'])[keep_all]
        print(f"  dropping only ac=0 columns: kept {keep_all.sum():,}", flush=True)

    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)
    save_npz(args.out_prefix + '.cn.npz', new_cn)
    np.savez(args.out_prefix + '.meta.npz',
             kmer_index=kmer_index,
             bubble_id=bubble_id,
             bubble_chrom=meta['bubble_chrom'],
             bubble_start=meta['bubble_start'],
             bubble_end=meta['bubble_end'],
             founders=meta['founders'])
    print(f"[subsample] wrote {args.out_prefix}.cn.npz + .meta.npz "
          f"(F={F}, K_new={new_cn.shape[1]:,}, nnz={new_cn.nnz:,})", flush=True)


if __name__ == "__main__":
    main()
