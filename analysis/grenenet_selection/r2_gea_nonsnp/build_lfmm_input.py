#!/usr/bin/env python
"""Build LFMM inputs for the gen-3 SV climate GEA (mirrors the phase-1 pipeline).

Produces, in their format (loci x samples Δp matrix + per-sample standardized env):
  delta_p_gen{g}_sv.csv     rows = filtered SVs, columns = pool names (site_gen_plot),
                            values = Δp = (gen-g pool AF) - p0   [the LFMM Y, transposed in R]
  env_gen{g}_<bio>.csv      one standardized (z-scored) climate value per pool, same
                            column order as the Δp matrix       [the LFMM X]
  locus_index_gen{g}_sv.csv  chrom,pos,ref_len,alt_len,sv_size,p0,rec_index for each SV row

Filter = their MAF>=0.05 AND min-count (sum_pools af*flowers*2 > total_genomes*0.05),
last-generation samples — removes rare-SV artifacts (e.g. the Chr3 stripes).
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", type=int, default=3)
    ap.add_argument("--maf", type=float, default=0.05)
    ap.add_argument("--min-count-frac", type=float, default=0.05)
    ap.add_argument("--bios", nargs="+", default=["bio1"])
    ap.add_argument("--out", default=f"{lib.GEA}/r2_gea_nonsnp/results/lfmm")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    D = f"{lib.GEA}/common/results/pool_matrices"
    mt = pd.read_csv(f"{D}/pool_gen{args.gen}_nonsnp.meta.csv")
    af = np.load(f"{D}/pool_gen{args.gen}_nonsnp_af.npy")            # [pools x rec]
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    sv_size = np.abs(idx["alt_len"] - idx["ref_len"])
    sv = sv_size > 50

    maf = np.minimum(p0, 1 - p0)
    m_maf = (maf >= args.maf) & sv
    # min-count: sum_pools(af * flowers*2) > total_genomes * frac
    genomes = mt.total_flowers.to_numpy(float) * 2
    counts = (af[:, m_maf].astype(np.float64) * genomes[:, None]).sum(0)
    keep = counts > genomes.sum() * args.min_count_frac
    rec_idx = np.where(m_maf)[0][keep]
    print(f"SVs>50bp {int(sv.sum()):,} -> MAF>={args.maf} {int(m_maf.sum()):,} "
          f"-> +min-count {len(rec_idx):,}")

    pools = list(mt.pool)
    dp = (af[:, rec_idx].astype(np.float64) - p0[rec_idx]).T            # [loci x pools]
    dpdf = pd.DataFrame(dp, columns=pools)
    dpdf.to_csv(f"{args.out}/delta_p_gen{args.gen}_sv.csv", index=False)

    pd.DataFrame(dict(
        chrom=idx["chrom"][rec_idx], pos=idx["pos"][rec_idx],
        ref_len=idx["ref_len"][rec_idx], alt_len=idx["alt_len"][rec_idx],
        sv_size=sv_size[rec_idx], p0=p0[rec_idx], rec_index=rec_idx
    )).to_csv(f"{args.out}/locus_index_gen{args.gen}_sv.csv", index=False)

    for bio in args.bios:
        x = mt[bio].to_numpy(float)
        z = (x - x.mean()) / x.std()                                   # standardize
        pd.DataFrame({bio: z}).to_csv(
            f"{args.out}/env_gen{args.gen}_{bio}.csv", index=False)
    print(f"-> {args.out}/  (Δp {dp.shape[0]:,} loci x {dp.shape[1]} pools; "
          f"env: {', '.join(args.bios)})")


if __name__ == "__main__":
    main()
