#!/usr/bin/env python
"""Read off what each lifted-back insertion overlaps in its source assembly's annotation.

Consumes the PAFs from `_liftback_to_assemblies.py` and intersects each insertion's assembly
coordinates against the four precomputed tracks. Those tracks were built on whole assemblies,
so they carry the genomic context a 754 bp fragment does not -- which is the entire point of
lifting back rather than annotating fragments.

WHAT EACH TRACK ACTUALLY CONTAINS (checked, not assumed):
  13_annotation_helixer_v0.3.5   de novo gene models. Finds genes ABSENT from Col-0.
  12_annotation_liftoff_TAIR10   known Col-0 genes projected onto the assembly; carries
                                 `extra_copy_number`, so it identifies DUPLICATIONS of known
                                 genes. Blind to novel genes by construction.
  02_annotation_RepeatMasker     NOT a TE-family track despite the name -- it holds only
                                 centromere / telomere / 45S_rDNA / 5S_rDNA / chloroplast /
                                 mitochondria / N_stretch. Useful for (a) genomic-compartment
                                 context and (b) the ORGANELLAR-INSERTION screen, which
                                 matters because singleton insertions in a pangenome graph are
                                 enriched for organellar and contamination artefacts.
                                 TE families still require our own RepeatMasker run.
  03_annotation_TRASH_v2         satellite / tandem repeat arrays.

Env: kmate. Writes results/sv_adaptive/liftback_annotation.csv.
"""
import os, sys, glob, gzip
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
BY = f"{OUT}/seq/byasm"
REL = "/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso"
TRACKS = {
    "helixer": (f"{REL}/13_annotation_helixer_v0.3.5", "{a}_helixer.gff3.gz", {"gene"}),
    "liftoff": (f"{REL}/12_annotation_liftoff_TAIR10_sc0.90", "{a}.TAIR10.genes.gff", {"gene"}),
    "repeat":  (f"{REL}/02_annotation_RepeatMasker", "{a}.Repeats_merged.gff", None),
    "trash":   (f"{REL}/03_annotation_TRASH_v2", "{a}.scaffolds_contigs.reformatted.gff", None),
}
MIN_COV = 0.80          # query coverage required to trust a lift-back
MIN_ID = 0.95           # identity: the ALT is a substring of the assembly, so this is high


def load_gff(path, feats):
    """seqid -> (starts, ends, labels) sorted by start."""
    out = {}
    op = gzip.open if path.endswith(".gz") else open
    try:
        with op(path, "rt") as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if feats is not None and f[2] not in feats:
                    continue
                out.setdefault(f[0], []).append((int(f[3]), int(f[4]), f[2], f[8]))
    except FileNotFoundError:
        return {}
    for k, v in out.items():
        v.sort()
        st = np.array([x[0] for x in v]); en = np.array([x[1] for x in v])
        out[k] = (st, en, [x[2] for x in v], [x[3] for x in v])
    return out


def overlaps(qs, qe, iv):
    """Indices of intervals overlapping [qs,qe]. iv sorted by start."""
    st, en, _, _ = iv
    j = np.searchsorted(st, qe, side="right")          # all starts <= qe
    if j == 0:
        return []
    lo = max(0, j - 3000)
    idx = np.arange(lo, j)
    hit = idx[en[idx] >= qs]
    return hit.tolist()


def main():
    # pass 1 (asm5, long) + pass 2 (k11/w5 recovery for short queries)
    groups = {}
    for d in (BY, f"{OUT}/seq/byasm_short"):
        for f in glob.glob(f"{d}/*.paf"):
            groups.setdefault(os.path.basename(f)[:-4], []).append(f)
    items = sorted(groups.items())
    print(f"[paf] {len(items)} assemblies, "
          f"{sum(len(v) for v in groups.values())} PAF files")
    rows = []
    for n, (a, pafs_a) in enumerate(items, 1):
        best = {}
        for line in (l for f in pafs_a for l in open(f)):
            f = line.split("\t")
            if len(f) < 12:
                continue
            q = f[0]; ql = int(f[1]); tn = f[5]
            ts = int(f[7]); te = int(f[8]); nm = int(f[9]); al = int(f[10])
            cur = best.get(q)
            if cur is None or nm > cur[4]:
                best[q] = (tn, ts, te, ql, nm, al)
        if not best:
            continue
        tr = {k: load_gff(os.path.join(d, p.format(a=a)), ft)
              for k, (d, p, ft) in TRACKS.items()}
        for q, (tn, ts, te, ql, nm, al) in best.items():
            cov = al / max(ql, 1); ident = nm / max(al, 1)
            rec = dict(key=q, asm=a, tname=tn, tstart=ts, tend=te, qlen=ql,
                       cov=cov, ident=ident, ok=(cov >= MIN_COV and ident >= MIN_ID))
            for k in TRACKS:
                iv = tr[k].get(tn)
                if iv is None:
                    rec[k] = ""
                    continue
                hits = overlaps(ts, te, iv)
                if not hits:
                    rec[k] = ""
                elif k in ("helixer", "liftoff"):
                    ids = []
                    for i in hits:
                        at = iv[3][i]
                        gid = ""
                        for part in at.split(";"):
                            if part.startswith("ID="):
                                gid = part[3:]; break
                        ids.append(gid)
                    rec[k] = ",".join(sorted(set(ids))[:5])
                else:
                    rec[k] = ",".join(sorted({iv[2][i] for i in hits})[:5])
            rows.append(rec)
        if n % 10 == 0 or n == len(items):
            print(f"  [{n}/{len(items)}] {a}: {len(best):,} lifted "
                  f"(running total {len(rows):,})", flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/liftback_annotation.csv", index=False)
    print(f"\n[lifted] {len(D):,} insertions aligned back")
    print(f"  passing cov>={MIN_COV} & id>={MIN_ID}: {int(D.ok.sum()):,} ({100*D.ok.mean():.1f}%)")
    G = D[D.ok]
    print(f"\n  carries a Helixer de novo gene : {int((G.helixer != '').sum()):,} "
          f"({100*(G.helixer != '').mean():.1f}%)")
    print(f"  overlaps a Liftoff known gene  : {int((G.liftoff != '').sum()):,} "
          f"({100*(G.liftoff != '').mean():.1f}%)")
    print(f"  in a TRASH satellite array     : {int((G.trash != '').sum()):,} "
          f"({100*(G.trash != '').mean():.1f}%)")
    print(f"  in a curated repeat compartment: {int((G.repeat != '').sum()):,} "
          f"({100*(G.repeat != '').mean():.1f}%)")
    if (G.repeat != "").any():
        print("\n  repeat-compartment breakdown:")
        print(G.loc[G.repeat != "", "repeat"].value_counts().head(10).to_string())
    print(f"\n[wrote] {OUT}/liftback_annotation.csv")


if __name__ == "__main__":
    main()
