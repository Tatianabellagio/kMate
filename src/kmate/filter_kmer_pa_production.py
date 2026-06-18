"""Production kmer_pa column filter.

Drops k-mer columns that carry no founder-discriminating signal:
  * ac == 0   : carried by no founder (dead; never matched in reads).
  * ac == 1   : private singletons (the filt2 rationale; sequencing-error /
                private-repeat prone, zero cross-founder discrimination).
  * ac == F   : invariant — carried by EVERY founder, i.e. monomorphic within
                the panel. mu_k = sum_f h_f * 1 = 1 regardless of h, so it adds
                the same constant to every founder's M-step term (a mild pull
                toward uniform). Harmless to the AF result but useless bloat.

Production keep rule: min_ac <= ac <= F - invariant_margin
  defaults min_ac=2, invariant_margin=1  ==>  keep 2 <= ac <= F-1.

This generalizes the old `build_kmer_pa_filt2_*.sh` heredoc, which kept ac>=2
only and thus left the ac==F invariants in. Reads any kmer_pa prefix and writes
the filtered matrix + correspondingly subset meta.

Usage:
  python src/filter_kmer_pa_production.py \
      --in-prefix  data/kmer_pa_231_arch3_raw/kmer_pa_Chr1 \
      --out-prefix data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1
  (Production normally applies this inline via build_kmer_pa.py --filter-production;
   this standalone path is for filtering a pre-built raw matrix.)
Optional:
  --min-ac N            lower keep bound (default 2).
  --invariant-margin M  drop ac > F-M (default 1 -> drop only ac==F;
                        M=2 also drops ac==F-1 near-invariants, etc.).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz, save_npz


def production_keep_mask(ac, F, min_ac=2, invariant_margin=1):
    """Boolean keep-mask over k-mer columns: min_ac <= ac <= F - invariant_margin.

    Drops ac<min_ac (dead + private) and ac>F-invariant_margin (invariant /
    near-invariant). `ac` is the per-column carrier count (sum over founders).
    This is the single source of truth for the production filter rule; both the
    standalone filter and build_kmer_pa's in-line filter call it.
    """
    ac = np.asarray(ac).ravel()
    hi = F - invariant_margin
    return (ac >= min_ac) & (ac <= hi)


def report_filter(ac, F, min_ac, invariant_margin, log=print):
    """Print the per-class drop breakdown for a production filter pass."""
    hi = F - invariant_margin
    keep = production_keep_mask(ac, F, min_ac, invariant_margin)
    K = len(ac)
    log(f"[filter] dropping:")
    log(f"           ac==0  (dead):        {int((ac == 0).sum()):,}")
    log(f"           ac==1  (private):     {int((ac == 1).sum()):,}")
    log(f"           ac> {hi} (invariant):  {int((ac > hi).sum()):,}"
        f"   [ac==F exactly: {int((ac == F).sum()):,}]")
    log(f"[filter] keep rule: {min_ac} <= ac <= {hi}")
    log(f"[filter] kept: {int(keep.sum()):,}/{K:,} ({keep.mean()*100:.2f}%)")
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-prefix", required=True,
                    help="input kmer_pa prefix (<prefix>.kmer_pa.npz + <prefix>.meta.npz)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-ac", type=int, default=2)
    ap.add_argument("--invariant-margin", type=int, default=1,
                    help="drop columns with ac > F - margin (default 1 = drop ac==F only)")
    args = ap.parse_args()

    kmer_pa = load_npz(args.in_prefix + ".kmer_pa.npz").tocsr()
    meta = np.load(args.in_prefix + ".meta.npz", allow_pickle=True)
    F, K = kmer_pa.shape
    ac = np.asarray(kmer_pa.sum(axis=0)).flatten()

    print(f"[filter] input  ({F}, {K:,}) nnz={kmer_pa.nnz:,}", flush=True)
    keep = report_filter(ac, F, args.min_ac, args.invariant_margin,
                         log=lambda m: print(m, flush=True))

    cn_f = kmer_pa.tocsc()[:, keep].tocsr()
    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    save_npz(str(out) + ".kmer_pa.npz", cn_f)

    new_meta = {}
    for k in meta.files:
        a = meta[k]
        new_meta[k] = a[keep] if (a.ndim == 1 and a.shape[0] == K) else a
    np.savez(str(out) + ".meta.npz", **new_meta)
    print(f"[filter] wrote {out}.kmer_pa.npz  shape={cn_f.shape} nnz={cn_f.nnz:,}", flush=True)


if __name__ == "__main__":
    main()
