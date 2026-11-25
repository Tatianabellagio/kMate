#!/usr/bin/env python
"""Collect GraffiTE's raw per-chromosome masking output, and write the SV records it skipped as FASTA.

WHY. GraffiTE builds its RepeatMasker query with `bcftools view --types indels`, which drops every record
whose alleles are not a pure gain/loss after trimming shared bases. In the arch3 panel those are the
HPRC `annotate_vcf` COMPLEX records (ALT path replaces reference nodes) plus INS/DEL-labelled records
that also swap a reference base: 71,865 of the 227,851 SV records in `sv_sites.vcf`. GraffiTE's authors
state complex SVs are outside its design (github.com/cgroza/GraffiTE/issues/20), so these are masked here
with GraffiTE's own RepeatMasker + OneCode steps (`_submit_graffite_skipped.sh`) but get no TSD search.

Sequence convention is GraffiTE's (`bin/repmask_vcf.sh`): ALT for insertions, REF for deletions, whole
allele including any anchor / replaced bases. IDs are `{chrom}_{rec}` as in `_make_graffite_vcf.py`.

The raw outputs (unfiltered RepeatMasker `.out`, OneCode `.onecode.out`, the query FASTA) exist only in
Nextflow's work directory; `pangenome.vcf` keeps only `total_match_span > 0.80`. They are copied to
results/sv_adaptive/graffite/raw/{chrom}/ so the work directory can be cleaned.

Env: basic. Writes results/sv_adaptive/graffite/{raw,skipped}/{chrom}/.
"""
import glob, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
G = f"{HERE}/results/sv_adaptive/graffite"
WORK = "/global/scratch/projects/fc_moilab/tbellg/kmate/telib/nxf_work"
KEEP = ["indels.fa", "indels.fa.fai", "indels.length", "span", "genotypes_repmasked.vcf.gz",
        "repeatmasker_dir/indels.fa.out", "repeatmasker_dir/indels.fa.onecode.out",
        "repeatmasker_dir/indels.fa.tbl", "repeatmasker_dir/OneCode_LTR.dic", "repeatmasker_dir/onecode.log"]

# the five repeatmask_VCF task directories, one per chromosome
tasks = {}
for f in glob.glob(f"{WORK}/*/*/.command.sh"):
    if "repmask_vcf.sh" in open(f).read():
        d = os.path.dirname(f)
        assert open(f"{d}/.exitcode").read().strip() == "0", d
        with open(f"{d}/genotypes.sorted.vcf") as fh:
            ch = next(l for l in fh if l[0] != "#").split("\t")[0]
        assert ch not in tasks, f"two task dirs for {ch} (stale -resume copy?)"
        tasks[ch] = d
assert sorted(tasks) == [f"Chr{i}" for i in range(1, 6)], tasks

n_in = n_masked = n_skip = 0
for ch, d in sorted(tasks.items()):
    raw = f"{G}/raw/{ch}"
    os.makedirs(f"{raw}/repeatmasker_dir", exist_ok=True)
    for k in KEEP:
        shutil.copy2(f"{d}/{k}", f"{raw}/{k}")
    masked = {l[1:].strip() for l in open(f"{d}/indels.fa") if l[0] == ">"}
    sk = f"{G}/skipped/{ch}"
    os.makedirs(sk, exist_ok=True)
    ni = ns = 0
    with open(f"{d}/genotypes.sorted.vcf") as fh, open(f"{sk}/indels.fa", "w") as o:
        for l in fh:
            if l[0] == "#":
                continue
            f = l.split("\t", 6)
            ni += 1
            if f[2] in masked:
                continue
            seq = f[4] if len(f[4]) > len(f[3]) else f[3]
            o.write(f">{f[2]}\n{seq}\n")
            ns += 1
    assert ni == len(masked) + ns, (ch, ni, len(masked), ns)
    print(f"  {ch}: {ni:,} records = {len(masked):,} masked by GraffiTE + {ns:,} skipped  ({d})", flush=True)
    n_in += ni; n_masked += len(masked); n_skip += ns
print(f"\n{n_in:,} SV records = {n_masked:,} masked by GraffiTE + {n_skip:,} skipped")
print(f"raw outputs -> {G}/raw/  skipped FASTA -> {G}/skipped/")
