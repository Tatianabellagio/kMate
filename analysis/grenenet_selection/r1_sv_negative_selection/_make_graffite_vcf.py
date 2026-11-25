#!/usr/bin/env python
"""Build a sites-only SV VCF for GraffiTE, with IDs that join back to our analysis.

GraffiTE requires variant IDs of at most 50 characters; the arch3 panel's IDs are full bubble paths
running to thousands. IDs are rewritten as `{chrom}_{rec}` where `rec` is the record's ordinal within
that chromosome's VCF -- the same key `_extract_insertion_seqs.py` uses, so every GraffiTE annotation
joins directly to `insertion_context.csv`, `sv_event_catalog.csv` and the burden tables.

Keeps biallelic indels with |len(ALT) - len(REF)| >= 50 -- insertions AND deletions, because GraffiTE
annotates the REF sequence of a deletion the same way it annotates the ALT of an insertion. Genotypes are
dropped: GraffiTE's annotation path reads only REF and ALT, and `--genotype false` skips the step that
would need them.
"""
import gzip, os

PROJ = "/global/scratch/projects/fc_moilab/tbellg/kmate"
OUT = f"{os.path.dirname(os.path.abspath(__file__))}/results/sv_adaptive/graffite"
os.makedirs(OUT, exist_ok=True)
MINSV = 50
dst = f"{OUT}/sv_sites.vcf"

n_ins = n_del = n_tot = 0
with open(dst, "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    for i in range(1, 6):
        ch, cl = f"Chr{i}", f"chr{i}"
        src = f"{PROJ}/panel/arch3/{cl}/merged_231_{cl}_final.vcf.gz"
        rec_in_chrom = 0
        with gzip.open(src, "rt") as fh:
            for line in fh:
                if line.startswith("##"):
                    if i == 1 and line.startswith("##contig"):
                        o.write(line)
                    continue
                if line.startswith("#"):
                    if i == 1:
                        o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                    continue
                f = line.rstrip("\n").split("\t")
                ref, alt = f[3], f[4]
                if "," not in alt and "<" not in alt and len(ref) > 0 and len(alt) > 0:
                    d = len(alt) - len(ref)
                    if abs(d) >= MINSV:
                        o.write(f"{f[0]}\t{f[1]}\t{f[0]}_{rec_in_chrom}\t{ref}\t{alt}\t.\t.\t.\n")
                        n_ins += d > 0
                        n_del += d < 0
                rec_in_chrom += 1
        n_tot += rec_in_chrom
        print(f"  {ch}: {rec_in_chrom:,} records scanned", flush=True)
print(f"\nwrote {dst}")
print(f"  {n_ins:,} insertions + {n_del:,} deletions = {n_ins+n_del:,} SV records "
      f"(of {n_tot:,} panel records)")
