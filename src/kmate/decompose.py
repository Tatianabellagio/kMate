#!/usr/bin/env python
"""Decompose a multi-allelic pangenome VCF to the biallelic form kMate requires.

`kmate build-var-pa` needs one ALT per record. Getting there by pairwise
realignment (`bcftools norm -m -any`, `--atomize`, `vcfwave`) is *wrong on a
pangenome graph*: it scatters carriers across shifted positions and silently drops
them where ALT paths converge on the same atomic variant. Measured on this
project's panel: up to ~99% carrier loss at a SNP co-located with a multi-allelic
indel.

The correct approach is **symbolic-ID propagation**, developed for the Human
Pangenome Reference Consortium: every atomic variant nested inside a graph bubble
carries a symbolic ID, and per-sample genotypes are transferred from the
multi-allelic record to the biallelic catalog by *matching those IDs* rather than
by alignment, so convergence and position shifts are handled by construction.

TWO WAYS IN
-----------
  * `--gfa` (founders ARE the graph): pass the Minigraph-Cactus VCF as
    `--genotyped-vcf` and its GFA. The VCF is annotated against the graph and its
    own genotypes are decomposed. This is the usual case: founder assemblies ->
    `cactus-pangenome` -> `kmate decompose --gfa`.
  * `--annotated-vcf` + `--biallelic-catalog` (founders were GENOTYPED on a graph,
    e.g. with PanGenie): pass the annotation from the graph's own VCF, and the
    genotyped VCF separately. INFO/ID is copied across first (`transfer-id`).

ATTRIBUTION
-----------
The method is not kMate's. This command runs two scripts from PanGenie
(Jana Ebler, MIT), bundled unmodified under `kmate/_vendor/pangenie/`:

  * `annotate_vcf.py` — decomposes each bubble into its nested variants and emits
    the annotated multi-allelic VCF and the matching biallelic catalog, both
    carrying INFO/ID.
  * `convert-to-biallelic.py` — propagates each sample's GT onto every catalog
    record whose symbolic ID is on the called ALT path.

Cite Ebler et al. (2022) Nature Genetics 54:518-525 and Liao et al. (2023)
Nature 617:312-324 if the decomposition matters to your results.

The only step original to kMate is `transfer-id` (`kmate.transfer_id`), which
copies INFO/ID from the annotated catalog onto a genotyped VCF produced from the
same graph. It exists because `bcftools annotate -c INFO/ID` corrupts the
angle-bracketed graph-node IDs.

KNOWN LIMITATION (inherited)
----------------------------
`annotate_vcf.py` builds its atomic catalog with `vcfwave` internally, so a small
fraction of atomic variants (~0.8% on this project's panel) are absent from the
catalog and cannot be recovered by ID matching.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

CITATION = """\
The decomposition method here is not kMate's: it is HPRC symbolic-ID propagation,
run via annotate_vcf.py and convert-to-biallelic.py from PanGenie
(https://github.com/eblerjana/pangenie, MIT; bundled in kmate/_vendor/pangenie/).
  Ebler et al. (2022) Nature Genetics 54:518-525
  Liao et al. (2023) Nature 617:312-324
"""

VENDOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_vendor", "pangenie")


def _need(tool: str, hint: str) -> str:
    p = shutil.which(tool)
    if p is None:
        sys.exit(f"ERROR: `{tool}` not found on PATH.\n  {hint}")
    return p


def _script(path: str | None, name: str) -> str:
    """An explicit path wins; otherwise the bundled copy."""
    if path:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {path} does not exist")
        return path
    bundled = os.path.join(VENDOR, name)
    if not os.path.isfile(bundled):
        sys.exit(f"ERROR: bundled {name} missing from {VENDOR}; reinstall kmate")
    return bundled


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        prog="kmate decompose",
        description="Multi-allelic -> biallelic VCF by HPRC symbolic-ID propagation. "
                    "Give either --gfa, or --annotated-vcf with --biallelic-catalog. "
                    "Runs PanGenie's scripts (see --citation).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=CITATION,
    )
    ap.add_argument("--genotyped-vcf", required=True,
                    help="the VCF to decompose. With --gfa: the Minigraph-Cactus VCF of "
                         "that graph (phased or haploid GTs). Otherwise: a VCF genotyped "
                         "on the same graph, which may lack INFO/ID")
    ap.add_argument("--gfa", default=None,
                    help="the graph's GFA (gzipped, as cactus-pangenome --gfa writes it). "
                         "Annotates --genotyped-vcf against it")
    ap.add_argument("--annotated-vcf", default=None,
                    help="annotated MULTI-allelic VCF carrying INFO/ID (from annotate_vcf.py)")
    ap.add_argument("--biallelic-catalog", default=None,
                    help="biallelic catalog VCF with the same symbolic IDs (from annotate_vcf.py)")
    ap.add_argument("--out", required=True, help="output VCF(.gz)")
    ap.add_argument("--annotate-vcf-script", default=None,
                    help="use this annotate_vcf.py instead of the bundled copy")
    ap.add_argument("--convert-to-biallelic", default=None,
                    help="use this convert-to-biallelic.py instead of the bundled copy")
    ap.add_argument("--haploidize", action="store_true",
                    help="emit haploid GTs, as kMate's panel builders require "
                         "(see --het for how heterozygotes are handled)")
    ap.add_argument("--het", choices=("missing", "split"), default="missing",
                    help="with --haploidize: 'missing' (default) for INBRED founders -- a "
                         "het call is treated as an artefact and set missing; 'split' for "
                         "OUTBRED, PHASED founders -- each sample becomes two haplotype "
                         "columns (sample.h1/.h2) and the founder axis doubles.")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--keep-temp", action="store_true")
    ap.add_argument("--citation", action="store_true", help="print attribution and exit")

    # --citation must work on its own, before argparse enforces the required inputs
    if "--citation" in (sys.argv[1:] if argv is None else argv):
        print(CITATION)
        return

    args = ap.parse_args(argv)
    explicit = args.annotated_vcf is not None or args.biallelic_catalog is not None
    if args.gfa and explicit:
        ap.error("give --gfa, or --annotated-vcf with --biallelic-catalog, not both")
    if not args.gfa and not (args.annotated_vcf and args.biallelic_catalog):
        ap.error("give --gfa, or both --annotated-vcf and --biallelic-catalog")

    sys.stderr.write("[kmate decompose] method: HPRC symbolic-ID propagation "
                     "(see --citation)\n")

    bcftools = _need("bcftools", "Install with: mamba install -c bioconda bcftools")
    convert = _script(args.convert_to_biallelic, "convert-to-biallelic.py")

    tmp = tempfile.mkdtemp(prefix="kmate_decompose_")
    try:
        if args.gfa:
            # --- step 1: annotate the graph VCF against its GFA (third-party) ---
            annotate = _script(args.annotate_vcf_script, "annotate_vcf.py")
            sys.stderr.write("[kmate decompose] step 1/4 annotate_vcf (PanGenie; slow, "
                             "reads the whole GFA)\n")
            plain = os.path.join(tmp, "graph.vcf")       # annotate_vcf reads plain text, twice
            subprocess.run([bcftools, "view", "-Ov", "-o", plain, args.genotyped_vcf],
                           check=True)
            prefix = os.path.join(tmp, "annotated")
            subprocess.run([sys.executable, annotate, "-vcf", plain, "-gfa", args.gfa,
                            "-o", prefix], check=True, stdout=sys.stderr)
            annotated_gt = prefix + ".vcf"               # already carries the GTs
            catalog = prefix + "_biallelic.vcf.gz"
            subprocess.run([bcftools, "view", "-Oz", "-o", catalog,
                            prefix + "_biallelic.vcf"], check=True)
        else:
            # --- step 1: transfer INFO/ID onto the genotyped VCF (kMate's own step) ---
            annotated_gt = os.path.join(tmp, "genotyped.annotated.vcf")
            catalog = args.biallelic_catalog
            sys.stderr.write("[kmate decompose] step 1/4 transfer INFO/ID (kMate)\n")
            from .transfer_id import main as transfer_main
            argv_save = sys.argv
            sys.argv = ["transfer-id", "--cactus", args.annotated_vcf,
                        "--pg", args.genotyped_vcf, "--out", annotated_gt]
            try:
                transfer_main()
            finally:
                sys.argv = argv_save

        # --- step 2: symbolic-ID propagation (third-party) ---
        sys.stderr.write("[kmate decompose] step 2/4 convert-to-biallelic (PanGenie; slow)\n")
        bial = os.path.join(tmp, "biallelic.vcf")
        with open(bial, "w") as out_fh:
            p1 = subprocess.Popen([bcftools, "view", annotated_gt], stdout=subprocess.PIPE)
            p2 = subprocess.Popen([sys.executable, convert, catalog],
                                  stdin=p1.stdout, stdout=out_fh)
            p1.stdout.close()
            if p2.wait() != 0 or p1.wait() != 0:
                sys.exit("ERROR: convert-to-biallelic failed")

        # --- step 3: sort. Nested variants land at shifted positions, so the output
        # is out of order, and build-var-pa / build-kmer-pa fetch by region (index). ---
        sys.stderr.write("[kmate decompose] step 3/4 sort\n")
        bial_sorted = os.path.join(tmp, "biallelic.sorted.vcf.gz")
        subprocess.run([bcftools, "sort", "-T", os.path.join(tmp, "sort"),
                        "-Oz", "-o", bial_sorted, bial], check=True)

        # --- step 4: fill tags, optionally haploidize ---
        sys.stderr.write("[kmate decompose] step 4/4 fill-tags%s\n"
                         % (" + haploidize" if args.haploidize else ""))
        filled = os.path.join(tmp, "filled.vcf.gz")
        subprocess.run([bcftools, "+fill-tags", bial_sorted, "--threads", str(args.threads),
                        "-Oz", "-o", filled, "--", "-t", "AC,AN,F_MISSING"], check=True)

        if args.haploidize:
            from .haploidize import haploidize_vcf
            haploidize_vcf(filled, args.out, bcftools=bcftools,
                           threads=args.threads, het=args.het)
        else:
            shutil.move(filled, args.out)

        if args.out.endswith(".gz"):
            subprocess.run([bcftools, "index", "-t", "-f", args.out], check=True)
        sys.stderr.write(f"[kmate decompose] wrote {args.out}\n")
    finally:
        if not args.keep_temp:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
