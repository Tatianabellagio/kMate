#!/usr/bin/env python3
"""Per-block identifiability diagnostics for the local-only window mode.

Goal: replace / justify the scalar `--min-kmers-per-block` floor (currently 200
default, 50 in the bench) with a derivation. For each genomic block we record
BOTH the raw nonzero-k-mer count AND principled resolvability metrics derived
from the design (founder carriage matrix) and the observed Fisher information:

  nnz            : # nonzero-count k-mers in the block (the current proxy)
  n_present      : # founders carrying >=1 nonzero k-mer in the block
  effrank_design : effective rank of the block carriage matrix K_b (F x nnz,
                   binary). exp(spectral entropy of Gram eigenvalues). This is
                   the COVERAGE-FREE identifiability: how many independent
                   founder-distinguishing patterns the block contains.
  effrank_fisher : effective rank of the observed Fisher information J_b at the
                   fitted h_b, on the estimated support (reuses
                   h_uncertainty._resolvability_from_J). Resolvability AT THE
                   DATA (folds in coverage + collinearity).
  cond           : condition number of the support Fisher info (collinearity).

We fit the per-block EM with min_kmers_per_block=1 (fit EVERY non-empty block,
no global fallback) and NO HMM smoothing, so each block's projected AF is its
own raw local fit — exactly the quantity the floor is meant to gate. The
companion analysis joins the per-record est to sim truth, aggregates realized
AF error per block, and asks: which metric predicts block error, and where is
the elbow?

Outputs (per pool x unit), under <outdir>:
  <pool>_<unit>.blockdiag.tsv  one row per block (the metrics above + status)
  <pool>_<unit>.recest.tsv     per-record chrom,pos,ref_len,alt_len,alt_freq,block

Usage:
  block_floor_diag.py --pool <pool> --unit {w10kb,dynldK500} \
     --kmer-pa-prefix ... --var-pa ... --var-meta ... --var-called ... \
     --kmer-db <jf> --blocks-tsv <dynld.tsv> --chrom Chr1 \
     --threads 8 --outdir <dir>
"""
import argparse, time
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz

from kmate.per_sample_per_chrom import (_count_and_load_kmer_pa_dense,
                                        load_blocks_tsv, _resolvability_from_J)
from kmate.block_em import (define_windows, assign_kmers_to_blocks,
                            assign_records_to_blocks, solve_em_per_block,
                            project_blocks_to_records)
from kmate.h_uncertainty import fisher_information_h


def effrank_from_eigs(ev):
    """Effective rank = exp(spectral entropy) of nonneg eigenvalues (continuous
    count of dominant directions). Matches the convention used for the Fisher
    eff_rank so design and data metrics are comparable."""
    ev = np.clip(np.asarray(ev, float), 0, None)
    pos = ev[ev > ev.max() * 1e-12] if ev.max() > 0 else ev[:0]
    if pos.size == 0:
        return 0.0
    p = pos / pos.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--unit", required=True, choices=["w10kb", "dynldK500"])
    ap.add_argument("--kmer-pa-prefix", required=True)
    ap.add_argument("--var-pa", required=True)
    ap.add_argument("--var-meta", required=True)
    ap.add_argument("--var-called", required=True)
    ap.add_argument("--kmer-db", required=True)
    ap.add_argument("--blocks-tsv", default=None, help="dynld units TSV (unit=dynldK500)")
    ap.add_argument("--window-bp", type=int, default=10000)
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--variant", default="", help="optional tag suffix, e.g. filt2 / unfiltered")
    a = ap.parse_args()

    chrom = a.chrom
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    tag = f"{a.pool}_{a.unit}" + (f"_{a.variant}" if a.variant else "")

    # --- count k-mers + load this chrom's kmer_pa (queries the prebuilt JF) ---
    kmer_pa_dense, counts, meta, cov, F, K = _count_and_load_kmer_pa_dense(
        chrom, a.kmer_pa_prefix, None, a.threads, kmer_db=a.kmer_db)
    if kmer_pa_dense is None:
        raise SystemExit(f"kmer_pa missing for {chrom}")
    counts = counts.astype(np.float32)

    bubble_id = meta["bubble_id"]
    bubble_chrom = meta["bubble_chrom"]
    bubble_start = meta["bubble_start"]
    bubble_end = meta["bubble_end"]

    # --- blocks ---
    if a.unit == "dynldK500":
        assert a.blocks_tsv, "dynldK500 needs --blocks-tsv"
        blocks = load_blocks_tsv(a.blocks_tsv, chrom)
    else:
        blocks = define_windows(bubble_chrom, bubble_start, bubble_end,
                                window_bp=a.window_bp)
    n_blocks = len(blocks)
    print(f"[{tag}] {n_blocks} blocks", flush=True)

    kmer_block = assign_kmers_to_blocks(bubble_id, bubble_chrom,
                                        bubble_start, bubble_end, blocks)
    # inv_mb de-replication weight (matches production --kmer-weight inv_mb)
    m_b = np.bincount(bubble_id)[bubble_id].astype(np.float32)
    omega = (1.0 / m_b).astype(np.float32)

    # --- pure local fits: floor=1 (fit every non-empty block), no anchor, no smooth ---
    t = time.time()
    h_blocks, status, global_h = solve_em_per_block(
        counts, kmer_pa_dense, kmer_block, n_blocks, cov,
        em_max_iter=200, tol=1e-7, min_kmers_per_block=1,
        global_anchor_weight=0.0, omega=omega, local_only=True,
        n_workers=1)
    print(f"[{tag}] block-EM floor=1 {time.time()-t:.0f}s "
          f"({(status==0).sum()} fit / {(status==2).sum()} empty)", flush=True)

    # --- per-block diagnostics ---
    nz = counts > 0
    block_kmer_idx = [np.flatnonzero(kmer_block == b) for b in range(n_blocks)]
    rows = []
    t = time.time()
    for b in range(n_blocks):
        idxs = block_kmer_idx[b]
        idxs_nz = idxs[nz[idxs]]
        nnz = int(idxs_nz.size)
        blk = blocks[b]
        rec = dict(block=b, chrom=str(blk.chrom), start=int(blk.start),
                   end=int(blk.end), n_kmers=int(idxs.size), nnz=nnz,
                   status=int(status[b]))
        if nnz == 0:
            rec.update(n_present=0, effrank_design=0.0, effrank_fisher=0.0,
                       cond=np.inf, eff_n_h=np.nan)
            rows.append(rec); continue
        Kb = kmer_pa_dense[:, idxs_nz]                 # F x nnz (binary)
        present = np.flatnonzero(Kb.sum(axis=1) > 0)
        # design eff-rank: spectral entropy of Gram(Kb) restricted to present founders
        Kp = Kb[present, :].astype(np.float64)
        G = Kp @ Kp.T                                   # |present| x |present|
        ev_d = np.linalg.eigvalsh(G)
        effrank_design = effrank_from_eigs(ev_d)
        # fisher eff-rank at the fit, on the EM support
        hb = h_blocks[b]
        support = np.flatnonzero(hb > 1e-3)
        if support.size:
            Jb = fisher_information_h(hb, Kb, counts[idxs_nz], omega=omega[idxs_nz])
            effrank_fisher, cond = _resolvability_from_J(Jb, support)
        else:
            effrank_fisher, cond = 0.0, np.inf
        eff_n_h = float(1.0 / np.sum(hb.astype(np.float64) ** 2)) if hb.sum() > 0 else np.nan
        rec.update(n_present=int(present.size), effrank_design=effrank_design,
                   effrank_fisher=float(effrank_fisher), cond=float(cond),
                   eff_n_h=eff_n_h)
        rows.append(rec)
    print(f"[{tag}] diagnostics {time.time()-t:.0f}s", flush=True)

    import pandas as pd
    diag = pd.DataFrame(rows)
    diag_path = outdir / f"{tag}.blockdiag.tsv"
    diag.to_csv(diag_path, sep="\t", index=False)
    print(f"[{tag}] -> {diag_path} ({len(diag)} blocks)", flush=True)

    # --- project to records (NO smoothing, NaN fallback) + dump per-record block id ---
    var_pa = load_npz(a.var_pa)
    var_meta = np.load(a.var_meta, allow_pickle=True)
    var_called = load_npz(a.var_called)
    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    rec_pos = np.asarray(var_meta["pos"])
    ref_len = np.asarray(var_meta["ref_len"])
    alt_len = np.asarray(var_meta["alt_len"])
    idx = np.where(rec_chrom == str(chrom))[0]
    var_pa_chrom = var_pa[:, idx]
    var_called_chrom = var_called[:, idx]
    rec_block = assign_records_to_blocks(rec_chrom[idx], rec_pos[idx], blocks)
    nan_fb = np.full_like(global_h, np.nan)
    freqs, info = project_blocks_to_records(h_blocks, nan_fb, var_pa_chrom,
                                            rec_block, var_called=var_called_chrom)
    est = pd.DataFrame({
        "chrom": rec_chrom[idx], "pos": rec_pos[idx],
        "ref_len": ref_len[idx], "alt_len": alt_len[idx],
        "alt_freq": freqs, "block": rec_block,
    })
    est_path = outdir / f"{tag}.recest.tsv"
    est.to_csv(est_path, sep="\t", index=False)
    print(f"[{tag}] -> {est_path} ({len(est)} records, "
          f"{np.isfinite(freqs).mean()*100:.1f}% finite)", flush=True)


if __name__ == "__main__":
    main()
