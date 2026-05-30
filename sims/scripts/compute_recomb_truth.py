"""Compute per-VCF-record truth AF for a recombinant pool — MAR convention.

Inputs:
  --ancestry        ancestry.tsv (ind_id, chrom, start, end, founder)
  --weights         pool_weights.tsv (founder, count, weight) — realized pool composition
  --cn-var          cn_var .npz  (F × N_records sparse, rows = founders)
  --cn-var-meta     cn_var meta .npz with 'founders', 'chrom', 'pos', 'ref_len', 'alt_len'
  --source-weights  Optional: source_weights.tsv (founder, prob) — if absent, auto-detected
                    next to pool_weights.tsv. When present, also emits source_truth_af.
  --out             output TSV (chrom, pos, ref_len, alt_len, truth_af, info,
                                [source_truth_af, source_info])

Logic (MAR convention, 2026-05-21):
  Each individual contributes weight w_i. For each variant record at (chrom, pos),
  find the founder owning that position in each individual i. The realized-pool
  TRUTH AF is:

      truth_af[r] = (Σ_i w_i · cn_var[founder_i, r]) / (Σ_i w_i · cn_var_called[founder_i, r])

  When source_weights.tsv is present, ALSO compute the "infinite-pool" truth
  expected under the source population's founder probabilities p_f (Stage-1
  binomial sampling integrated out):

      source_truth_af[r] = (Σ_f p_f · cn_var[f, r]) / (Σ_f p_f · cn_var_called[f, r])

  The pool-vs-source gap is the Stage-1 sampling + drift contribution to error —
  a noise floor below which no method can recover the source-population AF.

  Previously the denominator was Σ_i w_i (treating missing GTs as REF/non-carrier).
  Switched to the MAR (missing-at-random) convention to match how cactus_em
  global-mode projects AF: only count called founders in the denominator. See
  MISSINGNESS_231PANEL.md for the production-data motivation — at SV records on
  the 231-panel, 75% of PG founders are `.`, and the "treat-./. -as-REF"
  convention systematically under-calls SV AF.

  Also emits an `info` column = h-mass observable at this record. Tells
  downstream consumers how much of the pool's haplotype mass we have direct
  evidence for at each record. info < 0.5 → low-information record.

Runtime: vectorized — for each individual, build a per-record founder index
via interval-bucket lookup, gather the cn_var + cn_var_called columns, sum
into running totals. O(n_indiv × n_records).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


def assign_founder_per_record(rec_chrom, rec_pos, ind_segs_for_chrom):
    """For one individual, return the founder owning each record's position.

    rec_chrom: array of record chroms (str)
    rec_pos:   array of record 1-based positions (int)
    ind_segs_for_chrom: dict {chrom: sorted list of (start, end, founder)}
        starts/ends are 1-based inclusive.

    Returns array of founder names (str) length len(rec_chrom). Unmapped → ''.
    """
    out = np.empty(len(rec_chrom), dtype=object)
    out[:] = ''
    for chrom, segs in ind_segs_for_chrom.items():
        if not segs:
            continue
        starts = np.array([s for s, _, _ in segs])
        ends = np.array([e for _, e, _ in segs])
        founders = [f for _, _, f in segs]
        mask = rec_chrom == chrom
        if not mask.any():
            continue
        pos = rec_pos[mask]
        # Binary search: for each pos, the segment whose start <= pos and end >= pos
        idx = np.searchsorted(starts, pos, side='right') - 1
        idx = np.clip(idx, 0, len(starts) - 1)
        # Verify pos is within [start[idx], end[idx]]
        in_range = (pos >= starts[idx]) & (pos <= ends[idx])
        founders_arr = np.array(founders, dtype=object)
        sub = np.where(in_range, founders_arr[idx], '')
        out[mask] = sub
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ancestry", required=True)
    ap.add_argument("--weights", required=True)
    ap.add_argument("--cn-var", required=True)
    ap.add_argument("--cn-var-meta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--source-weights", default=None,
                    help="Optional source_weights.tsv (founder, prob). If absent, "
                         "auto-detected next to --weights as 'source_weights.tsv'. "
                         "When found, emits source_truth_af + source_info columns.")
    args = ap.parse_args()

    print("loading inputs...", flush=True)
    cn = load_npz(args.cn_var).tocsr()  # F × N
    # Auto-detect cn_var_called next to cn_var. Required for MAR truth.
    called_path = args.cn_var.replace(".cn_var.npz", ".cn_var_called.npz")
    try:
        cn_called = load_npz(called_path).tocsr()
        print(f"  cn_var_called loaded from {called_path}", flush=True)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"MAR truth requires cn_var_called.npz next to cn_var.npz.\n"
            f"Expected: {called_path}\n"
            f"Rebuild cn_var with poolfreq/src/build_cn_var.py (it writes both)."
        )
    meta = np.load(args.cn_var_meta, allow_pickle=True)
    cn_founders = list(meta["founders"])
    founder_to_row = {f: i for i, f in enumerate(cn_founders)}
    rec_chrom = np.asarray(meta["chrom"]).astype(str)
    rec_pos = np.asarray(meta["pos"]).astype(np.int64)
    rec_ref = np.asarray(meta["ref_len"]).astype(np.int32)
    rec_alt = np.asarray(meta["alt_len"]).astype(np.int32)
    N_records = len(rec_chrom)
    F = len(cn_founders)
    print(f"  cn_var: {F} founders × {N_records:,} records", flush=True)
    print(f"  cn_var_called density: {cn_called.nnz/(F*N_records)*100:.2f}%", flush=True)

    # Auto-detect source_weights.tsv next to --weights (or use --source-weights).
    # When found, we compute the "infinite-pool" truth in addition to the realized
    # pool truth. The pool↔source gap = Stage-1 sampling noise + drift.
    src_path = args.source_weights
    if src_path is None:
        guess = str(Path(args.weights).parent / "source_weights.tsv")
        if Path(guess).exists():
            src_path = guess
    source_probs = None
    if src_path and Path(src_path).exists():
        sw = pd.read_csv(src_path, sep='\t')
        if {"founder", "prob"}.issubset(sw.columns):
            prob_lookup = dict(zip(sw["founder"].astype(str), sw["prob"].astype(float)))
            source_probs = np.array([prob_lookup.get(str(f), 0.0) for f in cn_founders],
                                    dtype=np.float64)
            s = source_probs.sum()
            if s > 0:
                source_probs = source_probs / s
                print(f"  source_weights loaded from {src_path}", flush=True)
                print(f"    nonzero source founders: "
                      f"{int((source_probs>0).sum())}/{F}  "
                      f"max p={source_probs.max():.4f}  "
                      f"min p>0={source_probs[source_probs>0].min():.6f}",
                      flush=True)
            else:
                print(f"  WARN: source_weights probs sum to 0 — skipping source_truth",
                      flush=True)
                source_probs = None
        else:
            print(f"  WARN: source_weights at {src_path} lacks (founder, prob) cols; "
                  f"skipping source_truth", flush=True)
    else:
        print(f"  no source_weights.tsv — only realized-pool truth will be written",
              flush=True)

    anc = pd.read_csv(args.ancestry, sep='\t')
    inds = sorted(anc.ind_id.unique())
    n_indiv = len(inds)
    print(f"  ancestry: {len(anc):,} segments across {n_indiv} individuals", flush=True)

    # Per-individual pool weights. Two supported file formats:
    #   1) `pool_weights.tsv` (founder, count, weight) — gen-0 founder weights.
    #      In this case each individual gets uniform 1/n_indiv weight in the pool
    #      (the original recomb-sim convention).
    #   2) `visor_pool_fractions.tsv` (ind_id, clone_dir, visor_pct) — per-mosaic
    #      pool fractions used by VISOR. Used by the SKEWED recomb sim.
    #      Each individual gets visor_pct/100 weight in the pool.
    weights_df = pd.read_csv(args.weights, sep='\t')
    if 'ind_id' in weights_df.columns and 'visor_pct' in weights_df.columns:
        # Skewed mode: per-individual weights
        ind_weight = dict(zip(weights_df['ind_id'].astype(str),
                              weights_df['visor_pct'].astype(float) / 100.0))
        ws = sum(ind_weight.values())
        if abs(ws - 1.0) > 1e-3:
            print(f"  WARN: pool weights sum to {ws:.4f}, not 1.0", flush=True)
        print(f"  per-individual pool weights from {args.weights} (skewed mode)", flush=True)
        print(f"    weight stats: min={min(ind_weight.values()):.4f} "
              f"max={max(ind_weight.values()):.4f} "
              f"sum={ws:.4f}", flush=True)
    else:
        # Uniform mode: each individual = 1/n_indiv
        ind_weight = {ind: 1.0 / n_indiv for ind in inds}
        print(f"  uniform per-individual pool weights (1/{n_indiv})", flush=True)

    # MAR truth accumulates BOTH a weighted carrier sum (numerator) and a
    # weighted called sum (denominator). Final truth_af = num / denom, and
    # info = denom (fraction of pool h-mass observable at this record).
    num = np.zeros(N_records, dtype=np.float64)
    denom = np.zeros(N_records, dtype=np.float64)
    n_assigned = np.zeros(N_records, dtype=np.int32)

    for ind_idx, ind_id in enumerate(inds):
        ind_anc = anc[anc.ind_id == ind_id]
        segs_by_chrom = {}
        for chrom, sub in ind_anc.groupby("chrom"):
            sub_sorted = sub.sort_values("start")
            segs_by_chrom[chrom] = list(zip(
                sub_sorted.start.astype(int).tolist(),
                sub_sorted.end.astype(int).tolist(),
                sub_sorted.founder.astype(str).tolist(),
            ))
        founder_per_rec = assign_founder_per_record(rec_chrom, rec_pos, segs_by_chrom)
        # For each record where the founder is known, pull both cn_var (carrier
        # mask) and cn_var_called (called mask). Both go into the running totals
        # weighted by this individual's pool weight.
        ind_carrier = np.zeros(N_records, dtype=np.float32)
        ind_called  = np.zeros(N_records, dtype=np.float32)
        for f in set(founder_per_rec):
            if f == '' or f not in founder_to_row:
                continue
            mask = founder_per_rec == f
            row_idx = founder_to_row[f]
            row_carrier = np.asarray(cn[row_idx, :].todense()).flatten()
            row_called  = np.asarray(cn_called[row_idx, :].todense()).flatten()
            ind_carrier[mask] = row_carrier[mask]
            ind_called[mask]  = row_called[mask]
            n_assigned[mask] += 1
        w = ind_weight.get(ind_id, 1.0 / n_indiv)
        num   += w * ind_carrier
        denom += w * ind_called
        if (ind_idx + 1) % 10 == 0:
            print(f"  processed {ind_idx+1}/{n_indiv} individuals", flush=True)

    # MAR projection. At records where denom is essentially 0 (no individual's
    # ancestry-founder was called), truth_af is undefined — emit NaN. Downstream
    # eval should filter or treat NaN explicitly. `info` is the denominator
    # itself (h-mass observable).
    INFO_EPS = 1e-9
    info = denom.astype(np.float32)
    truth_af = np.where(denom > INFO_EPS, num / np.maximum(denom, INFO_EPS),
                        np.nan).astype(np.float32)

    out_cols = {
        "chrom":    rec_chrom,
        "pos":      rec_pos,
        "ref_len":  rec_ref,
        "alt_len":  rec_alt,
        "truth_af": truth_af,
        "info":     info,
    }

    # If source_weights provided, compute the "infinite-pool" truth: project
    # source population probabilities through cn_var (MAR). This is the
    # population parameter the AF estimator is ultimately trying to recover —
    # the realized-pool truth contains additional Stage-1 binomial + drift noise.
    if source_probs is not None:
        # Use sparse matrix products: source_probs is (F,), cn is (F, N).
        src_num = np.asarray(source_probs @ cn).ravel()
        src_den = np.asarray(source_probs @ cn_called).ravel()
        source_info = src_den.astype(np.float32)
        source_truth_af = np.where(src_den > INFO_EPS,
                                   src_num / np.maximum(src_den, INFO_EPS),
                                   np.nan).astype(np.float32)
        out_cols["source_truth_af"] = source_truth_af
        out_cols["source_info"] = source_info

    out_df = pd.DataFrame(out_cols)
    out_df.to_csv(args.out, sep='\t', index=False, compression='gzip' if args.out.endswith('.gz') else None)
    print(f"\nwrote {args.out}", flush=True)
    print(f"  truth_af nonzero (>0):     {((truth_af > 0) & np.isfinite(truth_af)).sum():,} / {N_records:,}", flush=True)
    print(f"  truth_af in (0,1):         {((truth_af > 0) & (truth_af < 1)).sum():,}", flush=True)
    print(f"  truth_af = NaN (info=0):   {np.isnan(truth_af).sum():,}", flush=True)
    print(f"  info = 1 (fully observed): {(np.isclose(info, 1.0)).sum():,}", flush=True)
    print(f"  records with all individuals ancestry-assigned: {(n_assigned == n_indiv).sum():,}", flush=True)
    if source_probs is not None:
        gap = np.abs(source_truth_af - truth_af)
        gap = gap[np.isfinite(gap)]
        print(f"  source_truth_af nonzero (>0): {((source_truth_af > 0) & np.isfinite(source_truth_af)).sum():,}", flush=True)
        print(f"  |source - pool| truth gap:    mean={gap.mean():.4f}  median={np.median(gap):.4f}  "
              f"p99={np.percentile(gap, 99):.4f}", flush=True)
        print(f"    ↑ this gap is the Stage-1 binomial + drift noise floor", flush=True)


if __name__ == "__main__":
    main()
