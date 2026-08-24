#!/usr/bin/env python
"""SNP overlap: production arch3 panel vs the old GrENE-Net 231 SNP catalog
(`greneNet_final_v1.1.recode.vcf`, the panel hapFIRE consumes natively).

This is the arch3 re-run of the v3-era analysis in
`old_docs/panel_overlap_v3_grenenet_results_summary.md` (archived at
`archive/panel_overlap_v3_grenenet/`, 2026-05-14: v3 covered only 55.05% of the
GrENE-Net catalog). arch3's decomposition (annotate_vcf + convert-to-biallelic,
avoiding `bcftools norm -m -any`) was expected to raise SNP coverage -- this
re-run gives the actual number.

Reports THREE overlap definitions:
  - position-only  : same (chrom, pos), any allele. Overcounts at the ~210K
    arch3 positions carrying >1 SNP record (multi-allelic decomposition) where
    the specific ALT differs from GrENE-Net's.
  - exact (chrom, pos, ref, alt) : the 4-tuple key `jobA7_compare_vs_hapfire.sh`
    already established as the correct join (its header warns "joining on pos
    alone manufactures off-diagonal scatter at multi-allelic split records").
    This is the number to trust over position-only.
  - exact, arch3 singletons (AC=1) dropped : same 4-tuple join, but arch3's SNP
    set is first restricted to founder AC>=2 (dropping the ~30% of arch3 SNP
    records carried by exactly 1 of 231 founders -- the least-reliable, most
    idiosyncratic calls; same AC definition PANEL_STATS.md's allele-frequency
    section and the K_pa filt2inv filter use). GrENE-Net's side is left as-is
    (its own singleton status isn't the question being asked).

Inputs:
  - arch3 SNP records: panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.meta.npz
    (ref_len==1 & alt_len==1) -- carries actual `ref`/`alt` base strings, not
    just lengths. Requires numpy>=2 to unpickle (run under the `kmate` env).
  - GrENE-Net 231 SNP records: read directly from the local mirror of
    greneNet_final_v1.1.recode.vcf at
    /global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/
    (verified: per-chrom record counts 827765/536387/617020/531361/722947
    match the archived old-test counts exactly; confirmed 0 multi-allelic
    rows and 0 duplicate (chrom,pos) in this VCF -- so for GrENE-Net,
    position-only and exact-4-tuple keying give identical set sizes; all the
    position-only/exact discrepancy comes from the arch3 side).

Output: analysis/grenenet_gea/panel_overlap_grenenet_summary.csv
  per-chrom and genome-wide, both overlap definitions.
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
import lib

PROJ = lib.PROJ
GEA = lib.GEA
GRENENET_VCF = ("/global/scratch/projects/fc_moilab/projects/grenenet-phase1/"
                "vcf/greneNet_final_v1.1.recode.vcf")
OUT_CSV = f"{GEA}/panel_overlap_grenenet_summary.csv"


def arch3_snp_records(chrom_n: int):
    """(pos, ref, alt, ac) for every SNP record on this chrom (arch3 keeps duplicate
    positions at multi-allelic SNP sites -- decomposed into separate biallelic rows).
    ac = founder allele count (carriers among the 231), from var_pa (CSR: rows=founders,
    cols=records) -- bincount of column indices over all nonzeros, same formula as
    PANEL_STATS.md's AC/AF section."""
    cl = f"chr{chrom_n}"
    base = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    m = np.load(f"{base}.meta.npz", allow_pickle=True)
    vp = np.load(f"{base}.var_pa.npz", allow_pickle=True)
    n_records = int(vp["shape"][1])
    ac_all = np.bincount(vp["indices"], minlength=n_records)
    snp = (m["ref_len"] == 1) & (m["alt_len"] == 1)
    pos = m["pos"][snp].astype(np.int64)
    ref = m["ref"][snp].astype(str)
    alt = m["alt"][snp].astype(str)
    ac = ac_all[snp]
    return pos, ref, alt, ac


def grenenet_snp_records(chrom_n: int):
    """(pos, ref, alt) for this chrom, parsed directly from the local GrENE-Net VCF
    mirror (plain-text, no bcftools needed -- confirmed biallelic-only, no dedup)."""
    pos, ref, alt = [], [], []
    with open(GRENENET_VCF) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 5)
            if f[0] != str(chrom_n):
                continue
            pos.append(int(f[1])); ref.append(f[3]); alt.append(f[4])
    return (np.array(pos, dtype=np.int64), np.array(ref, dtype=str), np.array(alt, dtype=str))


def main():
    rows = []
    for c in range(1, 6):
        ap, ar, aa, ac = arch3_snp_records(c)
        gp, gr, ga = grenenet_snp_records(c)

        # position-only sets
        a_pos, g_pos = set(ap.tolist()), set(gp.tolist())
        shared_pos = a_pos & g_pos

        # exact (pos, ref, alt) sets
        a_tup = set(zip(ap.tolist(), ar.tolist(), aa.tolist()))
        g_tup = set(zip(gp.tolist(), gr.tolist(), ga.tolist()))
        shared_tup = a_tup & g_tup

        # exact, arch3 singletons (AC=1) dropped
        keep = ac >= 2
        a_tup_ns = set(zip(ap[keep].tolist(), ar[keep].tolist(), aa[keep].tolist()))
        shared_tup_ns = a_tup_ns & g_tup

        rows.append(dict(
            chrom=f"Chr{c}",
            n_arch3_pos=len(a_pos), n_grenenet_pos=len(g_pos), n_shared_pos=len(shared_pos),
            pct_grenenet_in_arch3_pos=100 * len(shared_pos) / len(g_pos),
            n_arch3_exact=len(a_tup), n_grenenet_exact=len(g_tup), n_shared_exact=len(shared_tup),
            pct_grenenet_in_arch3_exact=100 * len(shared_tup) / len(g_tup),
            n_arch3_exact_nosingleton=len(a_tup_ns), n_shared_exact_nosingleton=len(shared_tup_ns),
            pct_grenenet_in_arch3_exact_nosingleton=100 * len(shared_tup_ns) / len(g_tup),
            n_arch3_singleton_dropped=int((~keep).sum()),
        ))
        print(f"Chr{c}: position-only shared={len(shared_pos):,} "
              f"({100*len(shared_pos)/len(g_pos):.2f}% of GrENE-Net)  |  "
              f"exact(pos,ref,alt) shared={len(shared_tup):,} "
              f"({100*len(shared_tup)/len(g_tup):.2f}% of GrENE-Net)  |  "
              f"exact, no arch3 singletons shared={len(shared_tup_ns):,} "
              f"({100*len(shared_tup_ns)/len(g_tup):.2f}% of GrENE-Net; "
              f"dropped {int((~keep).sum()):,} arch3 SNP singletons)")

    df = pd.DataFrame(rows)
    tot = {"chrom": "TOTAL"}
    for col in df.columns:
        if col != "chrom":
            tot[col] = df[col].sum()
    tot["pct_grenenet_in_arch3_pos"] = 100 * tot["n_shared_pos"] / tot["n_grenenet_pos"]
    tot["pct_grenenet_in_arch3_exact"] = 100 * tot["n_shared_exact"] / tot["n_grenenet_exact"]
    tot["pct_grenenet_in_arch3_exact_nosingleton"] = (
        100 * tot["n_shared_exact_nosingleton"] / tot["n_grenenet_exact"])
    df = pd.concat([df, pd.DataFrame([tot])], ignore_index=True)
    df.to_csv(OUT_CSV, index=False)
    print(df.to_string(index=False))
    print(f"[saved] {OUT_CSV}")


if __name__ == "__main__":
    main()
