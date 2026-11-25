#!/usr/bin/env python
"""Attach the INSERTED-SEQUENCE class (insertion_tair10_class.csv) to every tested insertion in
`sv_purging_annotated{TAG}.csv`, plus log size, log own carriers and a TE-derived flag -- the input of
the size / TE section of sv_purging_fixed.ipynb.

Join: af_store non-SNP row (`afrow`; for events the representative allele, i.e. the member with most
carriers) -> arch3 panel record index -> key "Chr|rec".
Env: basic. SV_TAG selects the run (e.g. _rec2, _event2). Writes sv_purging_insertions{TAG}_seqclass.csv.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
TAG = os.environ.get("SV_TAG", "")


def main():
    V = pd.read_csv(f"{OUT}/sv_purging_annotated{TAG}.csv")
    I = V[(V.cls == "insertion") & V.origin.isin(["cold-specific", "shared", "warm-specific"])].copy()
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz"); ch = idx["chrom"].astype("U5")
    off, o = {}, 0
    for C in lib.CHROMS:
        off[C] = o; o += int((ch == C).sum())
    rec = np.zeros(len(I), np.int64)
    for C in lib.CHROMS:
        m = np.load(f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}.meta.npz", allow_pickle=True)
        nonsnp = np.where(~((m["ref_len"] == 1) & (m["alt_len"] == 1)))[0]
        s = (I.chrom == C).to_numpy(); rec[s] = nonsnp[I.afrow.to_numpy()[s] - off[C]]
    I["key"] = I.chrom + "|" + pd.Series(rec, index=I.index).astype(str)
    T = pd.read_csv(f"{OUT}/insertion_tair10_class.csv")[["key", "cls"]].rename(columns={"cls": "seq_class"})
    I = I.merge(T, on="key", how="left")
    print(f"sequence class missing: {int(I.seq_class.isna().sum())} of {len(I)}")
    I["lsize"] = np.log10(I.size_bp); I["lcar"] = np.log10(I.n_carriers)
    I["te"] = (I.seq_class == "te_derived").astype(int)
    I.to_csv(f"{OUT}/sv_purging_insertions{TAG}_seqclass.csv", index=False)
    print(pd.crosstab(I.origin, I.seq_class).to_string())
    print(f"[wrote] {OUT}/sv_purging_insertions{TAG}_seqclass.csv")


if __name__ == "__main__":
    main()
