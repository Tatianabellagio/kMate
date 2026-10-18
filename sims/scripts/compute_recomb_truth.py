"""Per-record true allele frequency of a simulated pool.

For each record of a variant panel, the frequency of its alternative allele in the pool: individual i (share w_i,
shares.tsv) carries, at the record's position, the founder its ancestry track assigns there (ancestry.tsv), so

    truth_af = sum_i w_i * alt(founder_i) / sum_i w_i * called(founder_i)

over the individuals whose founder is called at the record (missing genotypes are left out, as kMate's projection does);
NaN where none is. info = sum_i w_i * called(founder_i), the share of the pool observed at the record.

Usage: python compute_recomb_truth.py --ancestry A --shares S --var-pa P.var_pa.npz --var-called P.var_called.npz
                                      --var-meta P.meta.npz --out truth.tsv.gz
Writes chrom, pos, ref_len, alt_len, truth_af, info, one row per panel record in panel order. Called by
benchmarks/simulations/simulate.py. Env: kmate.
"""
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


def founder_at(rec_chrom, rec_pos, segs):
    """Founder of each record's position on one individual's ancestry track ('' where no segment covers it).
    segs: ancestry rows of the individual (chrom, start, end, founder), 1-based inclusive."""
    out = np.full(len(rec_chrom), "", dtype=object)
    for chrom, s in segs.sort_values("start").groupby("chrom"):
        m = rec_chrom == chrom
        starts, ends, founders = s.start.to_numpy(), s.end.to_numpy(), s.founder.astype(str).to_numpy(dtype=object)
        k = np.clip(np.searchsorted(starts, rec_pos[m], side="right") - 1, 0, len(starts) - 1)
        out[m] = np.where((rec_pos[m] >= starts[k]) & (rec_pos[m] <= ends[k]), founders[k], "")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for arg in ("--ancestry", "--shares", "--var-pa", "--var-called", "--var-meta", "--out"):
        ap.add_argument(arg, required=True)
    a = ap.parse_args()
    alt, called = load_npz(a.var_pa).tocsr(), load_npz(a.var_called).tocsr()
    meta = np.load(a.var_meta, allow_pickle=True)
    row = {str(f): i for i, f in enumerate(meta["founders"])}
    rec_chrom, rec_pos = np.asarray(meta["chrom"]).astype(str), np.asarray(meta["pos"]).astype(np.int64)
    anc = pd.read_csv(a.ancestry, sep="\t", dtype={"founder": str})
    share = pd.read_csv(a.shares, sep="\t").set_index("ind_id").share
    assert set(share.index) == set(anc.ind_id) and abs(share.sum() - 1) < 1e-9, "shares.tsv does not match ancestry.tsv"
    assert set(anc.founder) <= set(row), "ancestry founders missing from the variant panel"

    num, den = np.zeros(len(rec_pos)), np.zeros(len(rec_pos))
    for ind, segs in anc.groupby("ind_id", sort=True):
        f_at = founder_at(rec_chrom, rec_pos, segs)
        carrier, obs = np.zeros(len(rec_pos), np.float32), np.zeros(len(rec_pos), np.float32)
        for f in set(f_at) - {""}:
            m = f_at == f
            carrier[m] = alt[row[f]].toarray().ravel()[m]
            obs[m] = called[row[f]].toarray().ravel()[m]
        num += share[ind] * carrier
        den += share[ind] * obs
    truth = np.where(den > 1e-9, num / np.maximum(den, 1e-9), np.nan).astype(np.float32)
    pd.DataFrame({"chrom": rec_chrom, "pos": rec_pos, "ref_len": np.asarray(meta["ref_len"]).astype(np.int32),
                  "alt_len": np.asarray(meta["alt_len"]).astype(np.int32), "truth_af": truth,
                  "info": den.astype(np.float32)}).to_csv(a.out, sep="\t", index=False)
    print(f"{a.out}: {len(truth):,} records, {np.isnan(truth).sum():,} without a called founder", flush=True)


if __name__ == "__main__":
    main()
