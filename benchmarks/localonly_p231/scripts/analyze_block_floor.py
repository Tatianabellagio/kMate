#!/usr/bin/env python3
"""Analyze per-block identifiability diagnostics: the elbow sweep + the derivation.

Consumes the floor=1 no-smooth outputs of block_floor_diag.py for each
(pool, unit) and asks two questions:

  1. ELBOW (empirical): if we abstain on blocks below an nnz floor T, how do
     call-rate and accuracy-on-called trade off as T rises? The single floor=1
     run gives the WHOLE curve in post (a record is "called at T" iff its
     block's nnz >= T). We do the same sweep on the identifiability metric
     effrank_design.

  2. DERIVATION (mechanism): per block, does realized AF error track the raw
     k-mer count (nnz) or the identifiability (effrank_design / effrank_fisher)?
     Spearman corr + scatter. If rank predicts error better than count, the
     floor should gate on rank; we then map the rank gate back to an equivalent
     nnz to compare with the current 50 / 200.

Outputs under results/:
  block_floor_sweep.png        call-rate & R2 vs nnz-floor and vs effrank-floor
  block_floor_derivation.png   per-block error vs predictor (+ effrank-vs-nnz)
  block_floor_summary.tsv      sweep table + predictor correlations
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KEYS = ["chrom", "pos", "ref_len", "alt_len"]
FD = Path("benchmarks/localonly_p231/floor_diag")
RES = Path("benchmarks/localonly_p231/results")
SIMS = Path("benchmarks/p231/sims")
UNITS = ["w10kb", "dynldK500"]
NNZ_GRID = [1, 10, 25, 50, 100, 200, 400, 800, 1600, 3200]
RANK_GRID = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
MIN_REC = 20            # min records/block for the per-block error regression
SHORT = {"cov10_n50_g0_s42_hotspots_p231_chr1": "n50 g0 outcross",
         "cov10_n231_g1_s42_self97_hotspots_p231_chr1": "n231 g1 SELFING",
         "cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1": "n50 g3 dom-sel"}


def load_joined(pool, unit):
    """Join per-record est to truth + attach the record's block diagnostics."""
    diag = pd.read_csv(FD / f"{pool}_{unit}.blockdiag.tsv", sep="\t")
    est = pd.read_csv(FD / f"{pool}_{unit}.recest.tsv", sep="\t")
    truth = pd.read_csv(SIMS / pool / "recomb_truth_raw.tsv.gz", sep="\t").dropna(subset=["truth_af"])
    truth["occ"] = truth.groupby(KEYS).cumcount()
    est["occ"] = est.groupby(KEYS).cumcount()
    m = truth.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz", "n_present", "effrank_design",
                      "effrank_fisher", "status"]], on="block", how="left")
    m["err"] = m["alt_freq"] - m["truth_af"]
    return m, diag


def r2(est, truth):
    if len(est) < 2:
        return np.nan
    ss = np.sum((truth - truth.mean()) ** 2)
    return 1 - np.sum((est - truth) ** 2) / ss if ss > 0 else np.nan


def sweep(m, col, grid):
    """For each floor T on `col`, call rate + R2 over records whose block col>=T."""
    out = []
    fin = np.isfinite(m["alt_freq"].values)
    n = len(m)
    for T in grid:
        keep = fin & (m[col].values >= T)
        nk = int(keep.sum())
        out.append(dict(floor=T, call_rate=nk / n,
                        R2=r2(m["alt_freq"].values[keep], m["truth_af"].values[keep]),
                        RMSE=float(np.sqrt(np.mean(m["err"].values[keep] ** 2))) if nk else np.nan))
    return pd.DataFrame(out)


def per_block_err(m):
    """Per-block realized RMSE (finite-est records only), merged with diagnostics."""
    mm = m[np.isfinite(m["alt_freq"].values)].copy()
    g = mm.groupby("block").agg(
        n_rec=("err", "size"),
        rmse=("err", lambda e: float(np.sqrt(np.mean(e ** 2)))),
        nnz=("nnz", "first"), n_present=("n_present", "first"),
        effrank_design=("effrank_design", "first"),
        effrank_fisher=("effrank_fisher", "first")).reset_index()
    return g[g.n_rec >= MIN_REC]


def main():
    pools = sorted({p.name.replace(f"_{u}.blockdiag.tsv", "")
                    for u in UNITS for p in FD.glob(f"*_{u}.blockdiag.tsv")})
    if not pools:
        sys.exit(f"no floor_diag outputs in {FD}")
    print("pools:", pools)

    summ_rows, corr_rows = [], []
    joined = {}   # (pool,unit) -> (m, diag, blkerr)
    for pool in pools:
        for unit in UNITS:
            if not (FD / f"{pool}_{unit}.blockdiag.tsv").exists():
                continue
            m, diag = load_joined(pool, unit)
            blk = per_block_err(m)
            joined[(pool, unit)] = (m, diag, blk)
            sw_n = sweep(m, "nnz", NNZ_GRID)
            for _, r in sw_n.iterrows():
                summ_rows.append(dict(pool=pool, unit=unit, gate="nnz", **r.to_dict()))
            sw_r = sweep(m, "effrank_design", RANK_GRID)
            for _, r in sw_r.iterrows():
                summ_rows.append(dict(pool=pool, unit=unit, gate="effrank_design", **r.to_dict()))
            # predictor correlations (Spearman |rho|, higher = better predicts error)
            if len(blk) > 5:
                c = {p: spearmanr(blk[p], blk.rmse).correlation
                     for p in ["nnz", "effrank_design", "effrank_fisher", "n_present"]}
            else:
                c = {p: np.nan for p in ["nnz", "effrank_design", "effrank_fisher", "n_present"]}
            corr_rows.append(dict(pool=pool, unit=unit, n_blocks=len(blk), **c))
            print(f"[{pool} {unit}] blocks(>= {MIN_REC} rec)={len(blk)}  "
                  f"spearman(rmse vs nnz)={c['nnz']:.3f}  vs effrank_design={c['effrank_design']:.3f}")

    summ = pd.DataFrame(summ_rows)
    corr = pd.DataFrame(corr_rows)
    summ.to_csv(RES / "block_floor_summary.tsv", sep="\t", index=False)
    corr.to_csv(RES / "block_floor_corr.tsv", sep="\t", index=False)
    print("\n=== predictor correlation (Spearman rmse vs metric) ===")
    print(corr.to_string(index=False))

    # ---- FIG 1: elbow sweep (rows=pool, col0=nnz gate, col1=effrank gate) ----
    fig, axes = plt.subplots(len(pools), 2, figsize=(11, 3.2 * len(pools)), squeeze=False)
    for i, pool in enumerate(pools):
        for j, (gate, grid) in enumerate([("nnz", NNZ_GRID), ("effrank_design", RANK_GRID)]):
            ax = axes[i][j]; ax2 = ax.twinx()
            for unit, c in zip(UNITS, ["C0", "C1"]):
                s = summ[(summ.pool == pool) & (summ.unit == unit) & (summ.gate == gate)]
                if s.empty:
                    continue
                ax.plot(s.floor, s.R2, "-o", color=c, ms=4, label=f"{unit} R²")
                ax2.plot(s.floor, s.call_rate, "--s", color=c, ms=3, alpha=.6)
            ax.set_xscale("log")
            if gate == "nnz":
                ax.axvline(50, color="gray", ls=":", lw=1)
                ax.axvline(200, color="k", ls=":", lw=1)
            ax.set_xlabel(f"{gate} floor"); ax.set_ylabel("R² on called (solid)")
            ax2.set_ylabel("call rate (dashed)"); ax2.set_ylim(0, 1.02)
            if i == 0 and j == 0:
                ax.legend(fontsize=7, loc="lower left")
            ax.set_title(f"{SHORT.get(pool, pool)} — gate on {gate}", fontsize=9)
    fig.suptitle("Local-only floor sweep: accuracy vs coverage as the floor rises\n"
                 "(grey ':'=current 50, black ':'=default 200)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(RES / "block_floor_sweep.png", dpi=130)
    print(f"-> {RES/'block_floor_sweep.png'}")

    # ---- FIG 2: derivation — per-block RMSE vs predictors ----
    fig, axes = plt.subplots(len(pools), 3, figsize=(13, 3.2 * len(pools)), squeeze=False)
    for i, pool in enumerate(pools):
        for j, unit in enumerate(["w10kb"]):   # use the cleaner uniform-window unit for the panels
            key = (pool, unit)
            if key not in joined:
                continue
            _, _, blk = joined[key]
            # col0: rmse vs nnz ; col1: rmse vs effrank_design ; col2: effrank vs nnz
            a0, a1, a2 = axes[i]
            a0.scatter(blk.nnz, blk.rmse, s=8, alpha=.4, color="C0")
            a0.set_xscale("log"); a0.set_xlabel("nnz k-mers"); a0.set_ylabel("block RMSE")
            rho_n = spearmanr(blk.nnz, blk.rmse).correlation if len(blk) > 5 else np.nan
            a0.set_title(f"{SHORT.get(pool,pool)}: RMSE vs nnz (ρ={rho_n:.2f})", fontsize=8.5)
            a1.scatter(blk.effrank_design, blk.rmse, s=8, alpha=.4, color="C2")
            a1.set_xlabel("effrank_design"); a1.set_ylabel("block RMSE")
            rho_r = spearmanr(blk.effrank_design, blk.rmse).correlation if len(blk) > 5 else np.nan
            a1.set_title(f"RMSE vs effrank (ρ={rho_r:.2f})", fontsize=8.5)
            a2.scatter(blk.nnz, blk.effrank_design, s=8, alpha=.4, color="C3")
            a2.set_xscale("log"); a2.set_xlabel("nnz k-mers"); a2.set_ylabel("effrank_design")
            a2.set_title("identifiability vs supply", fontsize=8.5)
    fig.suptitle("Derivation: does block error track k-mer COUNT or IDENTIFIABILITY? "
                 "(w10kb unit)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(RES / "block_floor_derivation.png", dpi=130)
    print(f"-> {RES/'block_floor_derivation.png'}")


if __name__ == "__main__":
    main()
