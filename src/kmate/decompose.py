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

ATTRIBUTION
-----------
This command is mostly ORCHESTRATION of third-party tools. It does not
reimplement them, and they are not bundled — you must obtain them separately.

  * `annotate_vcf.py` — Human Pangenome Reference Consortium, from the
    `prepare-vcf-MC` workflow in HPRC/genotyping-pipelines. Produces the annotated
    multi-allelic VCF and the matching biallelic catalog, both carrying INFO/ID.
      Liao et al. (2023) "A draft human pangenome reference", Nature 617:312-324.

  * `convert-to-biallelic.py` — Jana Ebler (eblerjana), from `pangenie-tools`
    / the PanGenie run-from-callset pipeline. Propagates each sample's GT onto
    every catalog record whose symbolic ID appears in the called ALT path.
      Ebler et al. (2022) "Pangenome-based genome inference allows efficient and
      accurate genotyping across a wide spectrum of variant classes",
      Nature Genetics 54:518-525.

  * `bcftools` for tag filling.
      Danecek et al. (2021) "Twelve years of SAMtools and BCFtools", GigaScience 10:giab008.

The only step original to kMate is `transfer-id` (`kmate.transfer_id`), which
copies INFO/ID from the annotated catalog onto a genotyped VCF produced from the
same graph. It exists because `bcftools annotate -c INFO/ID` corrupts the
angle-bracketed graph-node IDs.

If you publish results that used this command, cite the works above — the method
is theirs.

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
kmate decompose orchestrates third-party tools; the method is not ours:
  annotate_vcf.py         HPRC prepare-vcf-MC   — Liao et al. 2023, Nature 617:312-324
  convert-to-biallelic.py eblerjana/pangenie-tools — Ebler et al. 2022, Nat Genet 54:518-525
  bcftools                                      — Danecek et al. 2021, GigaScience 10:giab008
Only the INFO/ID transfer step is kMate's own. Please cite the above.
"""


def _need(tool: str, hint: str) -> str:
    p = shutil.which(tool)
    if p is None:
        sys.exit(f"ERROR: `{tool}` not found on PATH.\n  {hint}")
    return p


def _need_script(path: str | None, name: str, where: str) -> str:
    if path and os.path.isfile(path):
        return path
    found = shutil.which(name)
    if found:
        return found
    sys.exit(
        f"ERROR: `{name}` not found.\n"
        f"  It is third-party and is NOT bundled with kMate.\n"
        f"  Obtain it from: {where}\n"
        f"  Then pass its path explicitly (see `kmate decompose --help`)."
    )


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        prog="kmate decompose",
        description="Multi-allelic -> biallelic VCF by HPRC symbolic-ID propagation. "
                    "Orchestrates third-party tools (see --citation).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=CITATION,
    )
    ap.add_argument("--annotated-vcf", required=True,
                    help="annotated MULTI-allelic VCF carrying INFO/ID (from annotate_vcf.py)")
    ap.add_argument("--biallelic-catalog", required=True,
                    help="biallelic catalog VCF with the same symbolic IDs (from annotate_vcf.py)")
    ap.add_argument("--genotyped-vcf", required=True,
                    help="the genotyped VCF to decompose (same graph, may lack INFO/ID)")
    ap.add_argument("--out", required=True, help="output VCF(.gz)")
    ap.add_argument("--convert-to-biallelic", default=None,
                    help="path to eblerjana's convert-to-biallelic.py (third-party, not bundled)")
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

    sys.stderr.write(CITATION + "\n")

    bcftools = _need("bcftools", "Install with: mamba install -c bioconda bcftools")
    convert = _need_script(
        args.convert_to_biallelic, "convert-to-biallelic.py",
        "https://github.com/eblerjana/pangenie-tools "
        "(pipelines/run-from-callset/scripts/convert-to-biallelic.py)",
    )

    tmp = tempfile.mkdtemp(prefix="kmate_decompose_")
    try:
        # --- step 1: transfer INFO/ID onto the genotyped VCF (kMate's own step) ---
        annotated_gt = os.path.join(tmp, "genotyped.annotated.vcf")
        sys.stderr.write("[kmate decompose] step 1/3 transfer INFO/ID (kMate)\n")
        from .transfer_id import main as transfer_main
        argv_save = sys.argv
        sys.argv = ["transfer-id", "--cactus", args.annotated_vcf,
                    "--pg", args.genotyped_vcf, "--out", annotated_gt]
        try:
            transfer_main()
        finally:
            sys.argv = argv_save

        # --- step 2: symbolic-ID propagation (third-party) ---
        sys.stderr.write("[kmate decompose] step 2/3 convert-to-biallelic (eblerjana; slow)\n")
        bial = os.path.join(tmp, "biallelic.vcf")
        with open(bial, "w") as out_fh:
            p1 = subprocess.Popen([bcftools, "view", annotated_gt], stdout=subprocess.PIPE)
            p2 = subprocess.Popen([sys.executable, convert, args.biallelic_catalog],
                                  stdin=p1.stdout, stdout=out_fh)
            p1.stdout.close()
            if p2.wait() != 0 or p1.wait() != 0:
                sys.exit("ERROR: convert-to-biallelic failed")

        # --- step 3: fill tags, optionally haploidize ---
        sys.stderr.write("[kmate decompose] step 3/3 fill-tags%s\n"
                         % (" + haploidize" if args.haploidize else ""))
        filled = os.path.join(tmp, "filled.vcf.gz")
        subprocess.run([bcftools, "+fill-tags", bial, "--threads", str(args.threads),
                        "-Oz", "-o", filled, "--", "-t", "AC,AN,F_MISSING"], check=True)

        if args.haploidize:
            from .haploidize import haploidize_vcf
            haploidize_vcf(filled, args.out, bcftools=bcftools,
                           threads=args.threads, het=args.het)
        else:
            shutil.move(filled, args.out)

        subprocess.run([bcftools, "index", "-t", "-f", args.out], check=False)
        sys.stderr.write(f"[kmate decompose] wrote {args.out}\n")
    finally:
        if not args.keep_temp:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
