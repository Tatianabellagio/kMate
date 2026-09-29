#!/usr/bin/env python
"""Map every SV record of the arch3 panel to the graph BUBBLE (snarl) it came from.

WHY. arch3 is the HPRC/PanGenie representation: the graph VCF is multi-allelic, one record per
bubble, and each biallelic panel record is an atomic variant nested in a bubble, propagated by
symbolic-ID matching (panel/arch3/README.md). Counting biallelic records therefore counts ALLELES,
and the number of alleles at a site scales with how many assemblies are in the graph — one
insertion can appear as up to ~78 records, one per cactus assembly. For presence/absence questions
(is this insertion in cold founders and absent from warm ones?) the unit is the bubble.

This script recovers that unit exactly, instead of the chrom:pos proxy: symbolic IDs are read from
INFO/ID of the annotated multi-allelic VCF (one bubble per record, `ID` column = snarl id), and
matched to the ID column of the merged biallelic panel, whose record order is the var_pa column
order.

Env: kmate (only stdlib + numpy needed). Writes results/sv_adaptive/arch3_bubble_map.npz:
per chromosome, for every SV record (net length change > 50 bp) its column index in var_pa and an
integer bubble code (-1 = no bubble found; those fall back to chrom:pos:type downstream).
"""
import os, sys, gzip
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PANEL = "/global/scratch/users/tbellg/kmate/panel/arch3"
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP = 50


def panel_sv_records(ch):
    """Stream the merged biallelic panel VCF: return (col index, symbolic ID, pos, dlen) for SV records."""
    cl = ch.lower()
    path = f"{PANEL}/{cl}/merged_231_{cl}_final.vcf.gz"
    cols, ids, poss, dlens = [], [], [], []
    i = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 6)
            dlen = len(f[4]) - len(f[3])
            if abs(dlen) > SV_BP:
                cols.append(i); ids.append(f[2]); poss.append(int(f[1])); dlens.append(dlen)
            i += 1
    return np.array(cols), ids, np.array(poss), np.array(dlens), i


def bubble_of(ch, wanted):
    """Stream the annotated multi-allelic VCF; return {symbolic ID -> bubble id} for ids in `wanted`."""
    cl = ch.lower()
    path = f"{PANEL}/{cl}/{cl}_135_annotated.sorted.vcf.gz"
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            bub = f[2]
            j = f[7].find("ID=")
            if j < 0:
                continue
            info = f[7][j + 3:].split(";", 1)[0]
            for alt in info.split(","):
                for sid in alt.split(":"):
                    if sid in wanted:
                        out[sid] = bub
    return out


def main():
    acc = {}
    for ch in lib.CHROMS:
        cols, ids, poss, dlens, nrec = panel_sv_records(ch)
        meta = np.load(f"{lib.PROJ}/panel/arch3/{ch.lower()}/var_pa_231_arch3_{ch.lower()}.meta.npz", allow_pickle=True)
        mpos = meta["pos"].astype(np.int64); mdl = meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64)
        assert nrec == mpos.size, f"{ch}: panel VCF has {nrec:,} records, meta has {mpos.size:,}"
        assert np.array_equal(mpos[cols], poss) and np.array_equal(mdl[cols], dlens), f"{ch}: record order differs from meta"
        m = bubble_of(ch, set(ids))
        keys = [m.get(s, "") for s in ids]
        uniq = sorted({k for k in keys if k})
        code = {k: i for i, k in enumerate(uniq)}
        bub = np.array([code.get(k, -1) for k in keys], np.int64)
        acc[f"{ch}_rec"] = cols; acc[f"{ch}_bubble"] = bub
        nb = len(uniq); miss = int((bub < 0).sum())
        # how much does the bubble unit differ from the chrom:pos proxy?
        pos_units = len({(p, d > 0) for p, d in zip(poss, dlens)})
        print(f"[{ch}] {cols.size:,} SV records -> {nb:,} bubbles ({miss:,} unmatched); chrom:pos+type proxy gives {pos_units:,} units", flush=True)
    np.savez_compressed(f"{OUT}/arch3_bubble_map.npz", chroms=np.array(lib.CHROMS, dtype=str), **acc)
    print(f"[wrote] {OUT}/arch3_bubble_map.npz")


if __name__ == "__main__":
    main()
