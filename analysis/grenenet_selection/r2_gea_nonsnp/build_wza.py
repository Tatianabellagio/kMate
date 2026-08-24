#!/usr/bin/env python
"""WZA (weighted-Z) on the hapFIRE LD blocks — block-level SV climate-GEA.

Aggregates the noisy per-SV two-stage permutation p-values into one block-level
score, as the GrENE-net phase-1 project did: drives a local copy of their WZA
(`wza_script.py` = Booker WZA + SNP-number spline correction) over our LD blocks.

Self-contained: derives each SV's LD block from the two-stage npz's own chrom/pos
(lib.assign_ld_blocks), caching to twostage_blocks.npz. Per the code review:
  - always passes --empiricalP (our p_perm is ALREADY a calibrated permutation p;
    re-ranking it, the script's default, just spreads boundary ties — skip it);
  - --maf-filter, --maf-source, --pcol, --tag are CLI (nothing load-bearing hardcoded);
  - --maf-source {founding,contemporary}: weight by MAF from p0 (founding) or the
    mean evolved-pool frequency (contemporary; Booker's intent) — sensitivity knob;
  - --pcol selects the per-SV stat: p_perm (two-sided assoc.) or p_perm_up
    (one-tailed up-in-warm -> a DIRECTIONAL/signed WZA).

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  # standard (two-sided, founding MAF, empiricalP):
  $PY analysis/grenenet_selection/r2_gea_nonsnp/build_wza.py --stat dp scoef
  # directional up-in-warm:        --pcol p_perm_up --tag up
  # contemporary-MAF sensitivity:  --maf-source contemporary --tag contempMAF
"""
from __future__ import annotations
import argparse, glob, os, sys, subprocess
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

GEA = lib.GEA; STORE = lib.AF_STORE; PM = f"{GEA}/pool_matrices"
WZA_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "wza_script.py")
BLKCACHE = f"{GEA}/gea/twostage_blocks.npz"
OUTDIR = f"{GEA}/gea/wza"
NC_MIN, SV_MIN_BP = 150, 50    # SV set (must match build_two_stage_gea)
FLOWER = {"AT4G00650": "FRI", "AT5G10140": "FLC", "AT1G65480": "FT",
          "AT2G45660": "SOC1", "AT5G61850": "LFY"}


def get_blocks(chrom, pos):
    """LD-block id per SV (cached). Derived from the SVs' own chrom/pos, so it is
    always aligned to the two-stage npz order — no external dependency."""
    if os.path.exists(BLKCACHE):
        b = np.load(BLKCACHE, allow_pickle=True)["block"].astype(str)
        if len(b) == len(pos):
            return b
    print("  computing SV->LD-block map (lib.assign_ld_blocks)...", flush=True)
    b = lib.assign_ld_blocks(chrom, pos)
    np.savez(BLKCACHE, block=np.asarray(b, dtype=object))
    return np.asarray(b, dtype=str)


def contemporary_maf():
    """MAF from the mean evolved-pool frequency (mean over all gen1/2/3 pools),
    aligned to the filtered SV order. Booker weights by contemporary heterozygosity."""
    size = None
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    size = np.abs(idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64))
    nc = np.asarray(np.load(sorted(glob.glob(f"{STORE}/nc_nonsnp/*.npy"))[0]))
    sv_idx = np.where((size > SV_MIN_BP) & (nc >= NC_MIN))[0]
    tot = np.zeros(len(sv_idx)); n = np.zeros(len(sv_idx))
    for g in (1, 2, 3):
        P = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy")[:, sv_idx].astype(np.float64)
        tot += np.nansum(P, axis=0); n += np.sum(np.isfinite(P), axis=0)
    pbar = tot / np.where(n > 0, n, np.nan)
    return np.minimum(pbar, 1 - pbar)


def run_stat(stat, climate, pcol, maf_source, maf_filter, tag, genes, cmaf=None):
    z = np.load(f"{GEA}/gea/twostage_{stat}_{climate}.npz", allow_pickle=True)
    if pcol not in z.files:
        raise SystemExit(f"{pcol} not in {stat} npz (rerun build_two_stage_gea.py)")
    chrom = z["chrom"].astype(str); pos = z["pos"]
    block = get_blocks(chrom, pos)
    maf = (np.minimum(z["p0"], 1 - z["p0"]) if maf_source == "founding" else cmaf)
    df = pd.DataFrame(dict(block=block, chrom=chrom, pos=pos, ref_len=z["ref_len"],
                           stat=z[pcol], beta=z["beta"], z_emp=z["z_emp"],
                           sv_size=z["sv_size"], MAF=maf))
    df = df[(df.block != "") & np.isfinite(df.stat) & np.isfinite(df.beta)
            & np.isfinite(df.MAF)].copy()

    in_csv = f"{OUTDIR}/wza_in_{stat}{tag}.csv"
    out_csv = f"{OUTDIR}/wza_out_{stat}{tag}.csv"
    df[["block", "stat", "MAF"]].to_csv(in_csv, index=False)
    cmd = [sys.executable, WZA_SCRIPT, "-c", in_csv, "-s", "stat", "-w", "block",
           "--output", out_csv, "--sep", ",", "--sample_snps", "0",
           "--maf_filter", str(maf_filter), "--empiricalP"]
    subprocess.run(cmd, cwd=OUTDIR, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    w = pd.read_csv(out_csv).rename(columns={"gene": "block"})

    sub = df[df.MAF > maf_filter]
    agg = sub.groupby("block").agg(
        chrom=("chrom", "first"), mid_pos=("pos", "median"),
        n_up=("beta", lambda b: int((b > 0).sum())),
        n_dn=("beta", lambda b: int((b < 0).sum())),
        lead_pos=("pos", lambda p: p.loc[sub.loc[p.index, "z_emp"].abs().idxmax()]),
        lead_ref=("ref_len", lambda r: r.loc[sub.loc[r.index, "z_emp"].abs().idxmax()]),
        mean_beta=("beta", "mean")).reset_index()
    w = w.merge(agg, on="block", how="left")
    w["dir"] = np.where(w.mean_beta > 0, "up-in-warm", "down-in-warm")
    w = w.sort_values("Z_pVal")
    w.to_csv(f"{OUTDIR}/wza_{stat}_{climate}{tag}.csv", index=False)

    top = w.head(80).copy()
    top["pos"] = top.lead_pos.fillna(top.mid_pos).astype(int)
    top["ref_len"] = top.lead_ref.fillna(1).astype(int)
    top = lib.annotate_svs(top, flank=2000, genes=genes)
    top["flower_locus"] = top.genes_all.apply(
        lambda s: ";".join(FLOWER[g] for g in str(s).split(";") if g in FLOWER))
    keep = ["block", "chrom", "mid_pos", "SNPs", "Z", "Z_pVal", "top_candidate_p",
            "dir", "n_up", "n_dn", "gene", "gene_name", "flower_locus"]
    top[keep].to_csv(f"{OUTDIR}/wza_{stat}_{climate}{tag}.top.csv", index=False)

    nb = w.Z_pVal.notna().sum(); bonf = 0.05 / nb
    print(f"[{stat}{tag}] pcol={pcol} maf={maf_source} | {nb:,} blocks | "
          f"Z_pVal<0.05: {(w.Z_pVal<0.05).sum():,} (exp {int(0.05*nb):,}) | "
          f"<Bonf: {(w.Z_pVal<bonf).sum()} | min {w.Z_pVal.min():.2e} "
          f"-> wza_{stat}_{climate}{tag}.csv")
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stat", nargs="+", default=["dp", "scoef"])
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--pcol", default="p_perm", help="per-SV stat column (p_perm | p_perm_up)")
    ap.add_argument("--maf-source", default="founding", choices=["founding", "contemporary"])
    ap.add_argument("--maf-filter", type=float, default=0.0)
    ap.add_argument("--tag", default="", help="output filename suffix, e.g. _up")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    genes = lib.load_genes()
    cmaf = contemporary_maf() if args.maf_source == "contemporary" else None
    for stat in args.stat:
        run_stat(stat, args.climate, args.pcol, args.maf_source, args.maf_filter,
                 args.tag, genes, cmaf=cmaf)


if __name__ == "__main__":
    main()
