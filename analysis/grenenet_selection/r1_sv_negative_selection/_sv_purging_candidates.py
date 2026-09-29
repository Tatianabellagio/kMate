#!/usr/bin/env python
"""Candidate purged insertions: cold-specific, declining globally in GrENE-Net.

1. Attach the INSERTED-SEQUENCE class to every tested insertion (insertion_tair10_class.csv:
   te_derived / novel / dispersed_dup / gene_dup / local_dup, plus TE family). Keyed Chr|rec where
   rec is the full arch3 panel record index (verified 2000/2000 on pos and size). NB `novel` is
   partly alignment power, not guaranteed-novel sequence.
2. Candidates = cold-specific insertions called GLOBAL (sign test p < 0.05 over 31 gardens),
   ranked by consistency (fraction of gardens declining) then mean residual. Expect roughly a third
   to be chance hits: restricted indels are called global at ~32% vs ~49% here.
3. Gene function for the overlapping / promoter gene via annotate_genes_tair_uniprot.py.

Env: basic for 1-2; step 3 shells out to the kmate python (needs outbound HTTPS).
Writes results/sv_adaptive/sv_purging_insertions_seqclass.csv, sv_purging_candidates.csv.
"""
import os, sys, subprocess
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
ANN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/annotate_genes_tair_uniprot.py"
KPY = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/python"


def main():
    A = pd.read_csv(f"{OUT}/sv_purging_annotated.csv")
    I = A[A.cls == "insertion"].copy()
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz"); ch = idx["chrom"].astype("U5")
    off, o = {}, 0
    for C in lib.CHROMS:
        off[C] = o; o += int((ch == C).sum())
    rec = np.zeros(len(I), np.int64)
    for C in lib.CHROMS:
        m = np.load(f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}.meta.npz", allow_pickle=True)
        prl = m["ref_len"].astype(np.int64); pal = m["alt_len"].astype(np.int64)
        nonsnp = np.where(~((prl == 1) & (pal == 1)))[0]
        sel = (I.chrom == C).to_numpy()
        rec[sel] = nonsnp[I.afrow.to_numpy()[sel] - off[C]]
    I["rec"] = rec; I["key"] = I.chrom + "|" + I.rec.astype(str)
    T = pd.read_csv(f"{OUT}/insertion_tair10_class.csv")[["key", "cls", "family"]]
    T = T.rename(columns={"cls": "seq_class", "family": "te_family"})
    I = I.merge(T, on="key", how="left")
    assert I.seq_class.notna().all(), "insertions without a sequence class"
    I.to_csv(f"{OUT}/sv_purging_insertions_seqclass.csv", index=False)
    print(pd.crosstab(I.origin, I.seq_class).to_string(), "\n")

    C = I[(I.origin == "cold-specific") & I.is_global].sort_values(["frac_neg", "mean_res"],
                                                                   ascending=[False, True])
    print(f"candidates: {len(C)} cold-specific insertions called global (of "
          f"{(I.origin=='cold-specific').sum()})")
    genes = sorted(g for g in C.gene.dropna().unique() if g)
    fn = pd.DataFrame(columns=["gene"])
    if genes:
        tmp = f"{OUT}/_cand_genes.csv"
        pd.DataFrame({"gene": genes}).to_csv(tmp, index=False)
        outg = f"{OUT}/sv_purging_candidate_genes.csv"
        r = subprocess.run([KPY, ANN, "--genes-from", tmp, "--out", outg], capture_output=True, text=True)
        print(r.stdout[-1500:], r.stderr[-1500:])
        if os.path.exists(outg):
            fn = pd.read_csv(outg)
        os.remove(tmp)
    keep = [c for c in ("gene", "name", "protein_name", "uniprot_function", "categories") if c in fn.columns]
    C = C.merge(fn[keep], on="gene", how="left") if "gene" in fn.columns and len(fn) else C
    cols = ["chrom", "pos", "size_bp", "seq_class", "te_family", "context", "gene"] + \
           [c for c in keep if c != "gene"] + ["n_cold", "n_carriers", "p0", "frac_neg", "mean_res", "p_global"]
    C[[c for c in cols if c in C.columns]].to_csv(f"{OUT}/sv_purging_candidates.csv", index=False)
    print(f"[wrote] {OUT}/sv_purging_candidates.csv ({len(C)} rows)")


if __name__ == "__main__":
    main()
