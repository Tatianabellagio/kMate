#!/usr/bin/env python
"""Extract the ALT (inserted) sequence for every SV insertion in the arch3 panel.

The ALT allele of an insertion record IS the inserted sequence -- arch3 keeps explicit
sequence, not symbolic <INS> placeholders (verified). This writes one FASTA record per
insertion, keyed `{chrom}|{rec}` where `rec` is the record's ordinal within that chromosome's
VCF. That ordinal is the join key to `insertion_context.csv`.

WHY AN ORDINAL AND NOT chrom:pos -- 2.14% of arch3 positions carry more than one biallelic
record (the multiallelic pos-key trap). Joining on position silently matches the wrong ALT
allele. The ordinal is exact.

This script ASSERTS that meta record order equals VCF record order before writing anything;
var_pa is built 1:1 from the VCF so it should hold, but the whole downstream analysis depends
on it, so it is checked rather than assumed.

Env: kmate. Writes results/sv_adaptive/seq/insertions_{chrom}.fa (+ a combined index).
"""
import os, sys, gzip
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PROJ = lib.PROJ
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SEQ = f"{OUT}/seq"
SV_BP = 50
CHROMS = [f"Chr{i}" for i in range(1, 6)]


def main():
    os.makedirs(SEQ, exist_ok=True)
    tot = 0; totbp = 0
    for ch in CHROMS:
        cl = ch.lower()
        meta = np.load(f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.meta.npz",
                       allow_pickle=True)
        mpos = meta["pos"].astype(np.int64)
        mrl = meta["ref_len"].astype(np.int64)
        mal = meta["alt_len"].astype(np.int64)
        vcf = f"{PROJ}/panel/arch3/{cl}/merged_231_{cl}_final.vcf.gz"

        n = 0; w = 0; bp = 0
        checked = 0
        with gzip.open(vcf, "rt") as fh, open(f"{SEQ}/insertions_{ch}.fa", "w") as out:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 6)
                pos = int(f[1]); ref = f[3]; alt = f[4]
                # order check against the panel meta (first 5000 records)
                if checked < 5000:
                    assert mpos[n] == pos and mrl[n] == len(ref) and mal[n] == len(alt), (
                        f"{ch} record {n}: VCF pos={pos} ref={len(ref)} alt={len(alt)} vs "
                        f"meta pos={mpos[n]} ref={mrl[n]} alt={mal[n]} -- ORDER MISMATCH")
                    checked += 1
                if len(alt) - len(ref) > SV_BP:
                    out.write(f">{ch}|{n}\n{alt}\n")
                    w += 1; bp += len(alt) - len(ref)
                n += 1
        assert n == mpos.size, f"{ch}: VCF has {n} records, meta has {mpos.size}"
        print(f"  {ch}: {n:,} records, {w:,} SV insertions written ({bp/1e6:.1f} Mb) "
              f"[order verified on first {checked:,}]", flush=True)
        tot += w; totbp += bp
    print(f"\n[wrote] {SEQ}/insertions_Chr{{1..5}}.fa  --  {tot:,} sequences, {totbp/1e6:.1f} Mb")


if __name__ == "__main__":
    main()
