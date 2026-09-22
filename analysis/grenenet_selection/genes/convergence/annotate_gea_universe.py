#!/usr/bin/env python
"""Region-annotate every record the non-SNP GEA tested, with the SAME classifier as the hits.

The GEA hits carry a region call from genes/dissection/screen_sig_blocks.classify; the tested
background did not (only an 80k random subsample, atac_enrich_background.csv.gz). Analysing the
GEA's own results by region therefore needs the full tested set annotated identically. This
runs classify() unchanged on every unique (chrom, pos, ref_len) in the non-SNP LFMM universe
(the region call depends on the REF span only), twice:

  tier_1kb   classify's own promoter window (PROM_BP = 1000) -- what the hits carry
  tier_2kb   PROM_BP = 2000 -- the promoter definition of the burden model
             (notebooks/sv_climate_maf.ipynb: promoter <= 2 kb upstream)

Writes results/gea_universe_regions.csv.gz. env: kmate. Compute node.
"""
import os, sys
from multiprocessing import Pool
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, GEA); sys.path.insert(0, os.path.join(GEA, "genes", "dissection"))
LFMM = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"


def run(args):
    chunk, prom = args
    import screen_sig_blocks as S
    S.PROM_BP = prom
    C, _ = S.classify(chunk.copy())
    return C[["chrom", "pos", "ref_len", "tier"]].rename(columns={"tier": f"tier_{prom // 1000}kb"})


def main(procs):
    U = pd.read_csv(f"{LFMM}/lfmm_nonsnp_gen9_bio1.csv", usecols=["chrom", "pos", "ref_len", "alt_len"])
    U = U.drop_duplicates(["chrom", "pos", "ref_len"]).reset_index(drop=True)
    U["alt_len"] = U.ref_len          # classify reads only the REF span
    edges = np.linspace(0, len(U), procs * 4 + 1).astype(int)
    chunks = [U.iloc[a:b] for a, b in zip(edges[:-1], edges[1:])]
    out = {}
    for prom in (1000, 2000):
        with Pool(procs) as p:
            out[prom] = pd.concat(p.map(run, [(c, prom) for c in chunks]), ignore_index=True)
        print(f"PROM_BP={prom}: {out[prom].iloc[:, -1].value_counts().to_dict()}", flush=True)
    R = out[1000].merge(out[2000], on=["chrom", "pos", "ref_len"])
    R.to_csv(f"{HERE}/results/gea_universe_regions.csv.gz", index=False)
    print(f"wrote {len(R):,} loci")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
