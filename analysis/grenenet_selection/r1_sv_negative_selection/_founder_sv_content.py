#!/usr/bin/env python
"""Per-founder SV content, measured several ways, plus divergence from the Col-0 reference.

WHY. The provenance analysis so far used ONE summary of a founder's SV content: `ins_frac`
= (number of common SV insertions carried) / (number of all common variants carried). Two
things are wrong with leaning on that alone:

  (1) It is one arbitrary choice among several. SV content can be counted as a NUMBER of
      SVs or as the SIZE of the genome affected (bp). A founder carrying one 20 kb
      insertion and a founder carrying twenty 1 kb insertions have the same bp load and a
      20x different count. These can behave differently and should both be reported.

  (2) An "insertion" in this panel means sequence present in the founder and absent from
      TAIR10/Col-0. The pangenome is Col-0-referenced, so insertion count is, definitionally,
      a measure of how much extra sequence a founder has relative to Col-0 -- i.e. a
      DIVERGENCE measure. If divergence from Col-0 is itself climate-structured (Col-0 is a
      central-European accession, and the most diverged A. thaliana lineages are relicts from
      the range margins), then "insertion-rich founders decline in heat" could be
      "Col-0-diverged founders decline in heat" with no SV content involved.

This script computes every measure on the same footing so those can be told apart.

Measures per founder (all on the SAME common-variant filter: MAC>=12, call rate>=0.9, to
match `_founder_load_test.py`):
    counts   n_snp, n_indel, n_ins, n_del, n_sv, tot_carry
    bp       bp_ins, bp_del, bp_sv          (sum of |alt_len - ref_len| over carried records)
    ratios   each of the above / tot_carry
    divergence  n_snp  (ALT alleles at SNPs = differences from Col-0)

Env: kmate. Reads panel/arch3 var_pa. Writes results/sv_adaptive/founder_sv_content.npz.
"""
import os, sys
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PROJ = lib.PROJ
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
MIN_MAC, CALL_MIN, SV_BP = 12, 0.9, 50


def main():
    acc = None
    for ch in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        cl = ch.lower()
        base = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        founders = meta["founders"].astype("U6")
        rl = meta["ref_len"].astype(np.int64)
        al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsr()
        vc = sp.load_npz(f"{base}.var_called.npz").tocsr()
        F = vp.shape[0]
        if acc is None:
            acc = {k: np.zeros(F) for k in
                   ("n_snp", "n_indel", "n_ins", "n_del", "tot_carry",
                    "bp_ins", "bp_del",
                    # per-founder CALL rate: the technical-quality covariate. A founder
                    # whose genotypes are better resolved can carry more of everything.
                    "called_all", "called_sv", "n_common", "n_sv_sites",
                    # UNFILTERED divergence: ALT over EVERY segregating record, no MAC
                    # floor, so rare + private variants (what actually marks divergent
                    # lineages) are included. n_snp is MAC>=12 and misses exactly those.
                    "n_alt_all", "n_snp_all")}
            acc["founders"] = founders
        assert list(founders) == list(acc["founders"]), f"{ch} founder order mismatch"

        n_alt = np.asarray(vp.sum(0)).ravel()
        n_cal = np.asarray(vc.sum(0)).ravel()
        common = (n_alt >= MIN_MAC) & (n_alt <= F - MIN_MAC) & (n_cal >= CALL_MIN * F)

        dlen = np.abs(al - rl)
        isins = al > rl
        snp = common & (rl == 1) & (al == 1)
        indel = common & (dlen >= 1) & (dlen <= SV_BP)
        sv = common & (dlen > SV_BP)
        vpc = vp.tocsc()

        def csum(mask):
            return np.asarray(vp[:, np.where(mask)[0]].sum(1)).ravel().astype(float)

        def bpsum(mask):
            # per founder: sum of |dlen| over the carried records in `mask`
            w = np.where(mask, dlen, 0).astype(float)
            return np.asarray(vpc @ w).ravel()

        acc["n_snp"] += csum(snp)
        acc["n_indel"] += csum(indel)
        acc["n_ins"] += csum(sv & isins)
        acc["n_del"] += csum(sv & ~isins)
        acc["tot_carry"] += csum(common)
        acc["bp_ins"] += bpsum(sv & isins)
        acc["bp_del"] += bpsum(sv & ~isins)
        acc["called_all"] += np.asarray(vc[:, np.where(common)[0]].sum(1)).ravel().astype(float)
        acc["called_sv"] += np.asarray(vc[:, np.where(sv)[0]].sum(1)).ravel().astype(float)
        acc["n_alt_all"] += np.asarray(vp.sum(1)).ravel().astype(float)
        acc["n_snp_all"] += csum((rl == 1) & (al == 1))
        acc["n_common"] += float(common.sum())
        acc["n_sv_sites"] += float(sv.sum())
        print(f"[{ch}] common={int(common.sum()):,}  snp={int(snp.sum()):,} "
              f"indel={int(indel.sum()):,} SV ins={int((sv & isins).sum()):,} "
              f"del={int((sv & ~isins).sum()):,}", flush=True)

    acc["n_sv"] = acc["n_ins"] + acc["n_del"]
    acc["bp_sv"] = acc["bp_ins"] + acc["bp_del"]
    # n_common / n_sv_sites were accumulated into per-founder arrays of identical values
    acc["call_rate_all"] = acc["called_all"] / np.maximum(acc["n_common"], 1)
    acc["call_rate_sv"] = acc["called_sv"] / np.maximum(acc["n_sv_sites"], 1)
    tot = np.maximum(acc["tot_carry"], 1)
    for k in ("n_snp", "n_indel", "n_ins", "n_del", "n_sv", "bp_ins", "bp_del", "bp_sv"):
        acc[k + "_frac"] = acc[k] / tot

    np.savez_compressed(f"{OUT}/founder_sv_content.npz", **acc)
    f = acc["founders"]
    print(f"\nfounders={len(f)}")
    print(f"  n_snp      median {np.median(acc['n_snp']):,.0f}   (ALT at common SNPs, MAC>=12)")
    print(f"  n_snp_all  median {np.median(acc['n_snp_all']):,.0f}   (ALT at ALL SNPs -- better divergence proxy)")
    print(f"  n_alt_all  median {np.median(acc['n_alt_all']):,.0f}   (ALT over every segregating record)")
    print(f"  n_ins      median {np.median(acc['n_ins']):,.0f}    bp_ins median {np.median(acc['bp_ins']):,.0f} bp")
    print(f"  n_del      median {np.median(acc['n_del']):,.0f}    bp_del median {np.median(acc['bp_del']):,.0f} bp")
    print(f"  tot_carry  median {np.median(acc['tot_carry']):,.0f}")
    for nm in ("6909", "6911"):
        if nm in list(f):
            i = list(f).index(nm)
            print(f"  founder {nm}: n_snp={acc['n_snp'][i]:,.0f} n_ins={acc['n_ins'][i]:,.0f} "
                  f"(Col-0 would be near-zero divergence)")
    print(f"\n[wrote] {OUT}/founder_sv_content.npz")


if __name__ == "__main__":
    main()
