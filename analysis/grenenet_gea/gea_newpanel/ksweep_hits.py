#!/usr/bin/env python
"""K* selection + RAW-p hit calling for the LFMM K-sweep calibration.

Reads the per-class lambda_table_{cls}.csv (RAW GIF per axis x K) written by
ksweep_lfmm.R, and the RAW per-record p float32 dumps in ksweep_calib/raw_p/.

Per axis x class:
  * K* = the K whose RAW lambda is closest to 1 FROM ABOVE (min lambda s.t. >=1.0),
    preferring lambda in [1.0, 1.2]. If NO K reaches >=1.0 (deflated at all K down to
    K=1) -> flag 'deflated_all' and report the best-available (max) lambda.
  * At K* use the RAW p (no GIF). site-MAF>=0.05 filter (site_af = Y + p0[col],
    folded). Bonferroni (0.05 / n_tested) and BH-FDR q<0.05 -> clq0.9 blocks.
  * p-value histogram shape over the tested family.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea")
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel")
from blocks_clq09 import assign_clq09_blocks
from lib import annotate_svs, load_genes

ROOT = "/global/scratch/users/tbellg/kmate/results/grenenet_gea"
LFMM = f"{ROOT}/gea_newpanel/lfmm_site"
RECD = f"{ROOT}/phase1_replication/class_matrices"
OUTD = f"{ROOT}/gea_newpanel/ksweep_calib"
RAWD = f"{OUTD}/raw_p"

AXES = ["bio5", "pc1", "bio12", "bio13", "bio16", "bio19"]
CLASSES = ["snp", "nonsnp", "sv"]
NSITES = 31
MAF_MIN = 0.05


def bh_fdr(p):
    p = np.asarray(p, float); n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def pick_kstar(sub):
    """sub = rows for one axis x class with columns K, gif. Return dict."""
    sub = sub.sort_values("K")
    above = sub[sub["gif"] >= 1.0]
    if len(above):
        # closest to 1 from above = smallest gif among those >=1
        row = above.loc[above["gif"].idxmin()]
        kstar = int(row["K"]); lam = float(row["gif"])
        status = "calibrated" if lam <= 1.2 else "inflated_gt1.2"
        return dict(kstar=kstar, lam=lam, status=status)
    # all deflated
    row = sub.loc[sub["gif"].idxmax()]
    return dict(kstar=int(row["K"]), lam=float(row["gif"]), status="deflated_all")


def load_maf(cls):
    dims = open(f"{LFMM}/lfmm_{cls}_site_dims.txt").read().split()
    ns, nv = int(dims[0]), int(dims[1]); assert ns == NSITES
    Y = np.fromfile(f"{LFMM}/lfmm_{cls}_site_Y.f64", dtype=np.float64).reshape(ns, nv)
    rec = pd.read_csv(f"{RECD}/{cls}_gen9.records.csv"); assert len(rec) == nv
    p0src = "nonsnp" if cls == "sv" else cls
    p0 = np.load(f"{ROOT}/af_store/p0_{p0src}.npy").astype(np.float64)
    col = rec["col"].to_numpy()
    mean_af = (Y + p0[col][None, :]).mean(axis=0)
    maf = np.minimum(mean_af, 1.0 - mean_af)
    return rec, maf, nv


def hist_shape(p):
    """Qualitative label from a coarse histogram of raw p over tested loci."""
    b = np.histogram(p, bins=np.linspace(0, 1, 11))[0].astype(float)
    b /= b.sum()
    first = b[0]                      # fraction in [0,0.1)
    bulk = b[3:].mean()              # fraction density in [0.3,1)
    if first > 1.3 * bulk and first > 0.12:
        lab = "peak-near-0 (signal)"
    elif first < 0.85 * bulk:
        lab = "anti-conservative-depleted (deflated)"
    else:
        lab = "uniform (no signal)"
    return lab, b


def main():
    genes = load_genes()
    lam_tabs = {c: pd.read_csv(f"{OUTD}/lambda_table_{c}.csv") for c in CLASSES}
    full = pd.concat(lam_tabs.values(), ignore_index=True)
    full.to_csv(f"{OUTD}/lambda_table_all.csv", index=False)

    summary = []
    for cls in CLASSES:
        rec, maf, nv = load_maf(cls)
        maf_mask = maf >= MAF_MIN
        idx = np.where(maf_mask)[0]
        blocks = assign_clq09_blocks(rec["chrom"].to_numpy(), rec["pos"].to_numpy())
        rec = rec.assign(site_maf=maf, block_clq09=blocks)
        sub_base = rec.iloc[idx].reset_index(drop=True)
        ntest = len(idx)
        bonf_thr = 0.05 / ntest
        print(f"[{cls}] nv={nv:,} MAF>=.05 tested={ntest:,} Bonf_thr={bonf_thr:.2e}", flush=True)
        for axis in AXES:
            lt = lam_tabs[cls]
            ks = pick_kstar(lt[lt["axis"] == axis])
            K = ks["kstar"]
            p_all = np.fromfile(f"{RAWD}/{cls}_{axis}_K{K}.f32", dtype=np.float32).astype(np.float64)
            assert p_all.size == nv, (p_all.size, nv)
            p = p_all[idx]
            p = np.where(np.isfinite(p), p, 1.0)
            q = bh_fdr(p)
            bonf = p < bonf_thr
            fdr = q < 0.05
            lab, hbins = hist_shape(p)
            sub = sub_base.copy()
            sub["raw_p"] = p; sub["fdr_q"] = q
            sub["pass_bonf"] = bonf; sub["pass_fdr"] = fdr
            hits = sub[bonf | fdr].copy()
            if len(hits):
                hits = annotate_svs(hits, flank=2000, genes=genes)
            tag = f"{cls}_{axis}_K{K}"
            hits.sort_values("raw_p").to_csv(f"{OUTD}/hits_{tag}.csv", index=False)
            row = dict(cls=cls, axis=axis, kstar=K, lam=round(ks["lam"], 4),
                       status=ks["status"], n_tested=ntest,
                       n_bonf=int(bonf.sum()), n_fdr=int(fdr.sum()),
                       n_blocks_bonf=int(sub.loc[bonf, "block_clq09"].nunique()),
                       n_blocks_fdr=int(sub.loc[fdr, "block_clq09"].nunique()),
                       min_raw_p=float(p.min()), min_fdr_q=float(q.min()),
                       frac_p_lt05=float((p < 0.05).mean()), hist_shape=lab)
            summary.append(row)
            print(f"  {tag} lam={ks['lam']:.3f}[{ks['status']}] nBonf={bonf.sum()} "
                  f"nFDR={fdr.sum()} blkBonf={row['n_blocks_bonf']} blkFDR={row['n_blocks_fdr']} "
                  f"minP={p.min():.2e} minQ={q.min():.3f} frac<.05={row['frac_p_lt05']:.3f} "
                  f"hist={lab}", flush=True)
    sdf = pd.DataFrame(summary)
    sdf.to_csv(f"{OUTD}/kstar_hits_summary.csv", index=False)
    with open(f"{OUTD}/kstar_hits_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== K* HIT SUMMARY ===")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
