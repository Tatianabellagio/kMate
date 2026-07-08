#!/usr/bin/env python
"""Completeness check for a per_sample_per_chrom cohort (guards Finding 1: a
preempted task can leave a truncated per-chrom TSV that resumes as 'done' and
concatenates into a chromosome-short genome-wide TSV).

For every sample in MANIFEST, verify the final <sample>.tsv has EXACTLY the
expected record count (1 header + Σ_chrom n_records), and that each per-chrom
<sample>_<Chr>.tsv has its chrom's expected count. Reports: complete / short
(which chrom) / missing. With --fix, removes the bad final + short per-chrom
TSVs (so a rerun regenerates them) and writes a redo-manifest.

Usage:
  python grenenet/validate_cohort.py --manifest data/sample_manifest_usesample.tsv \
      --out-dir results/grenenet_gea/rerun_kfw_hb/evolved [--fix]
"""
import argparse, os, subprocess
import numpy as np

ROOT = "/global/scratch/users/tbellg/kmate"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
VAR_PA_DIR = f"{ROOT}/panel/arch3"
VAR_PA_TAG = "var_pa_231_arch3"


def expected_counts():
    exp = {}
    for c in CHROMS:
        cl = c.lower()
        m = np.load(f"{VAR_PA_DIR}/{cl}/{VAR_PA_TAG}_{cl}.meta.npz", allow_pickle=True)
        exp[c] = int(len(m["pos"]))
    return exp, sum(exp.values())


def wc(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return -1
    return int(subprocess.run(["wc", "-l", path], capture_output=True, text=True).stdout.split()[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--fix", action="store_true", help="remove bad final + short per-chrom TSVs and write a redo-manifest")
    args = ap.parse_args()

    exp, total = expected_counts()
    print(f"expected records/chrom: {exp}")
    print(f"expected genome-wide TSV lines = 1 header + {total:,} records = {total+1:,}\n")

    with open(args.manifest) as fh:
        header = fh.readline()
        rows = [ln for ln in fh if ln.strip()]
    sids = [ln.split("\t")[0] for ln in rows]

    complete = missing = short = 0
    redo = []
    for sid, ln in zip(sids, rows):
        final = f"{args.out_dir}/{sid}.tsv"
        n = wc(final)
        if n == -1:
            missing += 1
            redo.append(ln)  # not yet run (or removed) — normal if cohort still draining
            continue
        if n == total + 1:
            complete += 1
            continue
        # short/over — pinpoint the offending chrom(s)
        short += 1
        bad = []
        for c in CHROMS:
            pc = f"{args.out_dir}/{sid}_{c}.tsv"
            nc = wc(pc)
            if nc != exp[c] + 1:
                bad.append(f"{c}({nc-1 if nc>0 else 'missing'}/{exp[c]})")
                if args.fix and os.path.exists(pc):
                    os.remove(pc)
        print(f"  SHORT {sid}: final {n-1:,}/{total:,} records; bad chroms: {', '.join(bad) or '(concat stale?)'}")
        if args.fix and os.path.exists(final):
            os.remove(final)
        redo.append(ln)

    print(f"\ncomplete={complete}  short/corrupt={short}  missing/not-yet-run={missing}  (of {len(sids)})")
    if args.fix and (short or missing):
        redo_path = f"{args.out_dir}/_redo_manifest.tsv"
        with open(redo_path, "w") as fh:
            fh.write(header)
            fh.writelines(redo)
        print(f"[--fix] wrote redo-manifest ({len(redo)} samples) -> {redo_path}")
        print("  requeue with: sbatch --array=1-N%120 ... --export=ALL,MANIFEST=%s,... run_site_array_perchrom.sh" % redo_path)


if __name__ == "__main__":
    main()
