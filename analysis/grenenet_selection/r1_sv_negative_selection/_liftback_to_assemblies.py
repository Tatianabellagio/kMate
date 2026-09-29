#!/usr/bin/env python
"""STEP 0 -- lift each inserted sequence back into a carrier's assembly, then read off the
annotations that were computed there WITH full genomic context.

The move. An inserted sequence is, by construction, an exact substring of the assembly of any
founder that carries it. So aligning it back is a ~100%-identity lookup, not a fuzzy search.
Once we have assembly coordinates we can intersect the annotation tracks that already exist
for 82 assemblies:

    13_annotation_helixer_v0.3.5   de novo genes, including genes absent from Col-0
    12_annotation_liftoff_TAIR10   duplications of KNOWN Col-0 genes (blind to novel ones)
    02_annotation_RepeatMasker     TE family
    03_annotation_TRASH_v2         tandem / satellite repeats

Why this rather than annotating the fragments directly: Helixer's minimum record length is
25 kbp and its land-plant window is 21-107 kbp. Our median insertion is 754 bp. Running any
ab initio finder on the fragments is out-of-domain, not merely less accurate. The assemblies
were annotated with the context the models expect; this borrows that work.

Coverage limit, measured rather than assumed: only 80 of the 231 founders have an assembly,
and 57.5% of insertions are private to a single founder. Insertions private to a PanGenie
founder cannot be lifted at all. The point of this step is to quantify that residue.

Env: kmate (minimap2). Resumable -- existing PAFs are skipped.
Writes results/sv_adaptive/seq/byasm/*.fa, *.paf and liftback_assignment.csv.
"""
import os, sys, json, glob, subprocess
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PROJ = lib.PROJ
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SEQ = f"{OUT}/seq"; BY = f"{SEQ}/byasm"
ASM_DIR = "/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only"
MM2 = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/minimap2"
CHROMS = [f"Chr{i}" for i in range(1, 6)]
THREADS = int(os.environ.get("KMATE_THREADS", "16"))


def main():
    os.makedirs(BY, exist_ok=True)

    # ---- accession -> assembly id (only assemblies actually present in the graph dir)
    have = {os.path.basename(f).split(".")[0] for f in glob.glob(f"{ASM_DIR}/*.chr.fa")}
    req = pd.read_csv(f"{PROJ}/data/request_assemblies_for_Moi.csv")
    req["Assembly_ID"] = req["Assembly_ID"].astype(str)
    req["Accession_ID"] = req["Accession_ID"].astype(str)
    acc2asm = {}
    for src in (req[req.Best_for_Moi == "yes"], req):
        for _, r in src.iterrows():
            if r["Accession_ID"] not in acc2asm and r["Assembly_ID"] in have:
                acc2asm[r["Accession_ID"]] = r["Assembly_ID"]
    print(f"[asm] {len(acc2asm)} accessions map to an assembly present on disk")

    # ---- for each insertion pick ONE carrier that has an assembly
    assign = []
    for ch in CHROMS:
        cl = ch.lower()
        base = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        founders = meta["founders"].astype("U6")
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        ins = np.where((al - rl) > 50)[0]
        # founder rows that have an assembly, ordered so the choice is deterministic
        asm_rows = np.array([i for i, f in enumerate(founders) if f in acc2asm])
        asm_ids = np.array([acc2asm[founders[i]] for i in asm_rows])
        sub = vp[asm_rows][:, ins].tocsc()
        indptr, indices = sub.indptr, sub.indices
        for k, col in enumerate(ins):
            lo, hi = indptr[k], indptr[k + 1]
            if hi > lo:
                assign.append((ch, int(col), asm_ids[indices[lo]], int(hi - lo)))
            else:
                assign.append((ch, int(col), "", 0))
        print(f"  {ch}: {len(ins):,} insertions, "
              f"{sum(1 for a in assign[-len(ins):] if a[2]):,} liftable", flush=True)

    A = pd.DataFrame(assign, columns=["chrom", "rec", "asm", "n_asm_carriers"])
    A["key"] = A.chrom + "|" + A.rec.astype(str)
    A.to_csv(f"{OUT}/liftback_assignment.csv", index=False)
    lift = A[A.asm != ""]
    print(f"\n[assign] {len(lift):,}/{len(A):,} ({100*len(lift)/len(A):.1f}%) insertions have "
          f"at least one carrier with an assembly")
    print(f"[assign] spread over {lift.asm.nunique()} assemblies "
          f"(median {int(lift.asm.value_counts().median()):,} per assembly)")

    # ---- write one query FASTA per assembly, single streaming pass per chrom
    want = dict(zip(lift.key, lift.asm))
    handles = {}
    try:
        for ch in CHROMS:
            with open(f"{SEQ}/insertions_{ch}.fa") as fh:
                key = None
                for line in fh:
                    if line[0] == ">":
                        key = line[1:].strip()
                    elif key is not None:
                        a = want.get(key)
                        if a:
                            if a not in handles:
                                handles[a] = open(f"{BY}/{a}.fa", "w")
                            handles[a].write(f">{key}\n{line}")
            print(f"  routed {ch}", flush=True)
    finally:
        for h in handles.values():
            h.close()
    print(f"[fasta] wrote {len(handles)} per-assembly query files")

    # ---- minimap2 each query set against its own assembly
    todo = sorted(handles)
    for i, a in enumerate(todo, 1):
        paf = f"{BY}/{a}.paf"
        if os.path.exists(paf) and os.path.getsize(paf) > 0:
            continue
        ref = f"{ASM_DIR}/{a}.chr.fa"
        cmd = [MM2, "-x", "asm5", "-c", "--secondary=no", "-t", str(THREADS), ref, f"{BY}/{a}.fa"]
        with open(paf, "w") as o, open(f"{BY}/{a}.log", "w") as e:
            subprocess.run(cmd, stdout=o, stderr=e, check=True)
        n = sum(1 for _ in open(paf))
        print(f"  [{i}/{len(todo)}] {a}: {n:,} alignments", flush=True)

    print(f"\n[done] PAFs in {BY}")


if __name__ == "__main__":
    main()
