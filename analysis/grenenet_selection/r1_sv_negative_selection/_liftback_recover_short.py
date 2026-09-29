#!/usr/bin/env python
"""Recover the SHORT insertions that `-x asm5` failed to lift back.

Measured dropout bias: of the 157,300 insertions assigned a carrier assembly, 126,576 lifted
and ~30,700 did not -- and the failures are strongly length-biased (median 97 bp vs 1,017 bp
for successes). That is a minimap2 preset artefact, not biology: `-x asm5` uses k=19,w=19
minimizers and minimum chain/alignment scores tuned for assembly-scale contigs, so a ~100 bp
query often produces no reportable alignment.

Left uncorrected this would silently condition every class proportion on longer insertions.
This re-aligns only the missing keys with short-query settings (k=11, w=5, low min-score) and
appends the recovered hits.

Env: kmate. Writes seq/byasm_short/*.paf; `_annotate_liftback.py` picks them up on re-run.
"""
import os, sys, glob, subprocess
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
BY = f"{OUT}/seq/byasm"; SH = f"{OUT}/seq/byasm_short"
ASM_DIR = "/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only"
MM2 = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/minimap2"
THREADS = int(os.environ.get("KMATE_THREADS", "16"))


def main():
    os.makedirs(SH, exist_ok=True)
    A = pd.read_csv(f"{OUT}/liftback_annotation.csv")
    done = set(A.key)
    # asm IDs are numeric-looking; read as str or pandas makes them float ("100306.0")
    S = pd.read_csv(f"{OUT}/liftback_assignment.csv",
                    dtype={"asm": str}, keep_default_na=False)
    S = S[S.asm != ""]
    miss = S[~S.key.isin(done)]
    print(f"[recover] {len(miss):,} assigned insertions have no alignment yet")
    by = {a: set(g.key) for a, g in miss.groupby(miss.asm)}
    print(f"[recover] spread over {len(by)} assemblies")

    for i, (a, keys) in enumerate(sorted(by.items()), 1):
        paf = f"{SH}/{a}.paf"
        if os.path.exists(paf) and os.path.getsize(paf) > 0:
            continue
        src = f"{BY}/{a}.fa"
        if not os.path.exists(src):
            continue
        q = f"{SH}/{a}.fa"
        with open(src) as fh, open(q, "w") as o:
            k = None
            for line in fh:
                if line[0] == ">":
                    k = line[1:].strip()
                elif k in keys:
                    o.write(f">{k}\n{line}")
        if os.path.getsize(q) == 0:
            continue
        cmd = [MM2, "-k", "11", "-w", "5", "-m", "40", "-s", "40", "-c",
               "--secondary=no", "-t", str(THREADS), f"{ASM_DIR}/{a}.chr.fa", q]
        with open(paf, "w") as o, open(f"{SH}/{a}.log", "w") as e:
            subprocess.run(cmd, stdout=o, stderr=e, check=True)
        if i % 10 == 0 or i == len(by):
            print(f"  [{i}/{len(by)}] {a}", flush=True)
    n = sum(sum(1 for _ in open(p)) for p in glob.glob(f"{SH}/*.paf"))
    print(f"\n[done] {n:,} recovery alignments in {SH}")


if __name__ == "__main__":
    main()
