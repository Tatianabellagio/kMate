#!/usr/bin/env python
"""Build the clq0.9 HAPLOBLOCK-frequency last-gen pool matrix for phase-1-style GEA.

The variant phase-1 replication tests per-record allele frequencies. This builds the
analogous input where the test UNIT is a clq0.9 haploblock's haplotype-cluster
frequency (a one-vs-rest allele), so the same model runners (run_kendall/lfmm/binomial)
run unchanged with `--class hap`.

Frequencies (as the user asked): project the per-sample GLOBAL founder h onto the
clq0.9 founder->haplotype membership (build_hap_membership clq90) -> per-sample
haplotype freq (= sum of global h over the founders carrying that haplotype), then
flower-weight-pool to the phase-1 LAST-GEN pools (each plot's last surviving
generation; 355 site_gen_plot pools = gen9). Pooling matches build_pool_matrix
exactly (flower-weighted mean over the pool's timepoint samples).

  CAVEAT (documented): global h is one 231-vector per chromosome, so EVERY clq0.9
  block on a chromosome is a linear function of the same founder vector. The blocks are
  LD-independent in the panel but NOT independent measurements -> genome-wide Kendall is
  redundant/inflated (the clade<->allele collinearity wall). This is inherent to
  global-h projection; window-mode local h would be needed for truly independent blocks.

Encoding: within each block keep the (k-1) NON-reference haplotypes (reference = the
most panel-common cluster; its column is the linear complement), then apply the phase-1
contemporary filter MAF>=--maf-min over the 355 pools + require variation. Each kept
haplotype also carries its block's SV content (does the block span overlap a panel SV?)
for the downstream SV-enrichment test.

Outputs (--out, default analysis/grenenet_gea/phase1_replication/results/class_matrices):
  hap_gen9_af.npy        float32 [355 x n_kept]        per-pool haplotype freq
  hap_gen9.records.csv   hap_id,chrom,pos,ref_len,alt_len,MAF,block(+rich cols)
                         row-aligned to the matrix columns; consumed by run_kendall.py
  hap_gen9.registry.csv  full per-haplotype metadata + per-block SV tags (for enrichment)
  (reuses the existing gen9.pools.csv; asserts pool order matches)

Env: kmate (pyclustering not needed here; membership already built). Run on a compute node.
  python build_hap_lastgen_matrix.py
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))   # analysis/grenenet_gea
import lib
import build_hapfreq_matrix as bhm          # load_mem + project_sample (global-h path)

WIN = "results/grenenet_kmate_window"
HM = "results/grenenet_gea/blocks_mcf90/hap_membership"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
SV_MIN_BP = 50


def build_registry(memb_tag):
    """Vectorized per-haplotype registry from the membership npz (avoids the slow
    per-hap Python loop in build_hapfreq_matrix.load_mem(registry=True)).

    Global hap id runs Chr1..Chr5 concatenated (== projection column order). Every
    per-unit field is already stored in the npz, so we just np.repeat by cluster count
    and recover per-hap founder counts by a single global bincount."""
    rows, gbase = [], 0
    for ch in CHROMS:
        M = np.load(f"{HM}/{ch.lower()}_hapmemb_{memb_tag}.npz", allow_pickle=True)
        labels = M["labels"]; off = M["hap_offset"].astype(np.int64)
        U = labels.shape[0]; nh = int(off[-1]); k = np.diff(off)          # (U,)
        gid = (off[:-1][:, None] + labels).ravel()                        # (U*231,)
        n_founders = np.bincount(gid, minlength=nh)                       # per-hap founder count
        cluster = np.concatenate([np.arange(kk) for kk in k])            # 0..k-1 within unit
        rep = lambda a: np.repeat(a, k)                                   # per-unit -> per-hap
        rows.append(pd.DataFrame({
            "hap_id": gbase + np.arange(nh), "chrom": ch,
            "unit_idx": rep(np.arange(U)), "unit_start": rep(M["unit_start"].astype(int)),
            "unit_end": rep(M["unit_end"].astype(int)), "cluster": cluster,
            "n_founders": n_founders, "panel_freq": n_founders / labels.shape[1],
            "unit_nvar": rep(M["unit_nvar"].astype(int)),
            "unit_n_eff": rep(M["unit_n_eff"].astype(float)),
            "covered": rep(M["unit_covered"]),
        }))
        gbase += nh
    return pd.concat(rows, ignore_index=True)


def build_pool_hapfreq(pools_df, mem, nh):
    """Flower-weighted per-pool haplotype freq [n_pools x nh], in pools_df order.

    For each last-gen pool, gather its timepoint samples from lib.pool_table and
    flower-weight-average their projected hap freq (matching build_pool_matrix)."""
    pt = lib.pool_table()
    pt["sampleid"] = pt.sampleid.astype(str)
    # restrict to samples that actually have a window-run h file
    pt = pt[pt.sampleid.apply(lambda s: os.path.exists(f"{WIN}/{s}_Chr1.h_blocks_per_chrom.npz"))]
    by_pool = {p: g for p, g in pt.groupby(pt.pool.astype(str))}

    cache = {}
    def proj(s):
        if s not in cache:
            cache[s] = bhm.project_sample(s, mem)
        return cache[s]

    P = np.full((len(pools_df), nh), np.nan, np.float64)
    n_missing = 0
    for i, pool in enumerate(pools_df.pool.astype(str).to_numpy()):
        g = by_pool.get(pool)
        if g is None or len(g) == 0:
            n_missing += 1
            continue
        w = g.flowerscollected.astype(float).to_numpy()
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if w.sum() <= 0:
            w = np.ones(len(g))                          # all-zero flowers -> unweighted
        acc = np.zeros(nh); wsum = 0.0
        for s, ws in zip(g.sampleid.to_numpy(), w):
            acc += ws * proj(s); wsum += ws
        P[i] = acc / wsum
        if (i + 1) % 50 == 0:
            print(f"  pooled {i+1}/{len(pools_df)} pools", flush=True)
    if n_missing:
        print(f"[warn] {n_missing} pools had no window-h member samples (left NaN)")
    print(f"[proj] projected {len(cache):,} unique member samples")
    return P


def sv_tags_per_block(reg_block):
    """For each clq0.9 block span, count overlapping panel SV records (|dlen|>SV_MIN_BP)
    and the largest |dlen|. reg_block: one row per block with chrom/block_start/block_end."""
    # panel record coords from one arch3 sample (all samples share the segregating panel)
    samp = lib.list_samples(lib.OUT)[0]
    meta = pd.read_csv(f"{lib.OUT}/{samp}.tsv", sep="\t",
                       usecols=["chrom", "pos", "ref_len", "alt_len"])
    meta["dlen"] = (meta.alt_len - meta.ref_len).abs()
    sv = meta[meta.dlen > SV_MIN_BP].copy()
    print(f"[sv] panel SV records (|dlen|>{SV_MIN_BP}): {len(sv):,}")
    n_sv = np.zeros(len(reg_block), int); sv_bp = np.zeros(len(reg_block), int)
    for ch, gb in reg_block.groupby("chrom"):
        s = sv[sv.chrom == ch]
        if s.empty:
            continue
        spos = s.pos.to_numpy(); send = (s.pos + s.ref_len).to_numpy(); sbp = s.dlen.to_numpy()
        order = np.argsort(spos); spos, send, sbp = spos[order], send[order], sbp[order]
        for j, r in zip(gb.index.to_numpy(), gb.itertuples()):
            bs, be = int(r.block_start), int(r.block_end)
            # SV span [spos, send) overlaps block [bs, be):  spos < be  and  send > bs
            lo = np.searchsorted(spos, be, "right")          # SVs with spos < be
            ov = np.where(send[:lo] > bs)[0]
            n_sv[reg_block.index.get_loc(j)] = ov.size
            sv_bp[reg_block.index.get_loc(j)] = int(sbp[ov].max()) if ov.size else 0
    return n_sv, sv_bp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--memb-tag", default="clq90")
    ap.add_argument("--maf-min", type=float, default=0.05)
    ap.add_argument("--out", default=f"{lib.GEA}/phase1_replication/class_matrices")
    ap.add_argument("--pools-csv", default=None, help="default: <out>/gen9.pools.csv")
    a = ap.parse_args()
    pools_csv = a.pools_csv or f"{a.out}/gen9.pools.csv"
    os.makedirs(a.out, exist_ok=True)

    # global-h projection onto the clq0.9 membership grid (registry=False -> fast; the full
    # per-hap registry is built vectorized by build_registry, not the slow load_mem loop)
    bhm.MEMB_TAG = a.memb_tag; bhm.H_SOURCE = "global"; bhm.WIN = WIN
    mem, _, nh = bhm.load_mem(registry=False)
    reg = build_registry(a.memb_tag)
    assert len(reg) == nh and (reg.hap_id.to_numpy() == np.arange(nh)).all(), "registry misaligned"
    reg["block"] = reg.chrom + "_" + reg.unit_idx.astype(str)
    print(f"[memb] {a.memb_tag}: {reg.block.nunique():,} clq0.9 blocks, {nh:,} haplotypes", flush=True)

    pools = pd.read_csv(pools_csv)
    pools["pool"] = pools.pool.astype(str)
    print(f"[pools] {len(pools)} last-gen pools, {pools.site.nunique()} sites from {pools_csv}")

    P = build_pool_hapfreq(pools, mem, nh)                # [355 x nh], flower-weighted

    # ---- contemporary filter over the pools + one-vs-rest (drop per-block reference) ----
    pbar = np.nanmean(P, 0)
    maf = np.minimum(pbar, 1 - pbar)
    finite_frac = np.isfinite(P).mean(0)
    ptp = np.nanmax(P, 0) - np.nanmin(P, 0)

    reg = reg.assign(block=reg.chrom + "_" + reg.unit_idx.astype(str),
                     block_start=reg.unit_start.astype(int), block_end=reg.unit_end.astype(int),
                     pool_pbar=pbar, maf=maf, finite_frac=finite_frac, ptp=ptp)
    # reference haplotype per block = the most panel-common cluster (its pool column is the complement)
    ref_idx = reg.groupby("block").panel_freq.idxmax()
    is_ref = np.zeros(nh, bool); is_ref[ref_idx.to_numpy()] = True
    reg["is_ref"] = is_ref

    keep = (~is_ref) & (maf >= a.maf_min) & (finite_frac >= 0.5) & (ptp > 0)
    print(f"[filter] keep {int(keep.sum()):,}/{nh:,} haplotypes "
          f"(non-ref & MAF>={a.maf_min} & finite>=50% & variable) "
          f"across {reg.loc[keep,'block'].nunique():,} blocks")

    regk = reg[keep].reset_index(drop=True)

    # per-block SV content (on kept blocks' spans; join back to every kept hap)
    blk = regk.drop_duplicates("block")[["block", "chrom", "block_start", "block_end"]].reset_index(drop=True)
    n_sv, sv_bp = sv_tags_per_block(blk)
    blk["n_sv"] = n_sv; blk["sv_bp_max"] = sv_bp; blk["has_sv"] = blk.n_sv > 0
    print(f"[sv] kept blocks overlapping >=1 SV: {int(blk.has_sv.sum()):,}/{len(blk):,} "
          f"({blk.has_sv.mean():.1%})")
    regk = regk.merge(blk[["block", "n_sv", "sv_bp_max", "has_sv"]], on="block", how="left")

    # ---- emit: matrix + run_kendall records contract + rich registry ----
    Ak = P[:, keep].astype(np.float32)
    np.save(f"{a.out}/hap_gen9_af.npy", Ak)

    recs = pd.DataFrame({
        "hap_id": regk.hap_id.to_numpy(),
        "chrom": regk.chrom.to_numpy(),
        "pos": regk.block_start.to_numpy(),          # block start as the locus anchor
        "ref_len": 1, "alt_len": 1,                  # dummy (runner passes through only)
        "maf": regk.maf.to_numpy(),
        "block": regk.block.to_numpy(),
    })
    recs.to_csv(f"{a.out}/hap_gen9.records.csv", index=False)

    regk[["hap_id", "chrom", "block", "block_start", "block_end", "cluster", "n_founders",
          "panel_freq", "unit_n_eff", "unit_nvar", "covered", "pool_pbar", "maf",
          "finite_frac", "n_sv", "sv_bp_max", "has_sv"]].to_csv(
        f"{a.out}/hap_gen9.registry.csv", index=False)

    print(f"[done] {Ak.shape[0]} pools x {Ak.shape[1]:,} haploblocks -> {a.out}/hap_gen9_af.npy")
    print(f"[done] records + registry written to {a.out}/hap_gen9.{{records,registry}}.csv")


if __name__ == "__main__":
    main()
