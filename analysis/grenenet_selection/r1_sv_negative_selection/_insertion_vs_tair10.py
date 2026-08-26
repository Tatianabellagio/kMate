#!/usr/bin/env python
"""What IS the inserted sequence? Align it back to TAIR10 and classify the cargo.

One alignment answers both open questions at once:

  DUPLICATION vs NOVEL   an insertion is absent from Col-0 *at that locus*, but the sequence
                         often exists elsewhere in Col-0. Where it matches says where it came
                         from -- a local tandem duplication, a dispersed copy, or nothing.
  TE-DERIVED             if the TAIR10 match lands inside a `transposable_element`, the cargo
                         is TE-derived -- and TAIR10's TE records carry the FAMILY in their
                         `Alias` attribute (ATCOPIA24, ATREP4, ATHILA, VANDAL...), so we get
                         family-level calls without RepeatMasker or any external library.

WHY dc-megablast AND NOT minimap2. minimap2 `-x asm5/asm10` uses k=19,w=19 minimizers and
minimum chain scores tuned for assembly-scale contigs. Our median query is 754 bp and copies
may be diverged, so minimap2 silently misses short and diverged hits -- the same preset trap
that length-biased the first lift-back pass. Discontiguous megablast is sensitive to ~75%
identity and is the right instrument at this query length.

Classification per insertion (priority order):
  te_<FAMILY>      best hit overlaps a TAIR10 transposable_element -> TE-derived
  gene_dup         best hit overlaps a protein-coding gene
  local_dup        best hit within LOCAL_BP of the insertion's own site, same chromosome
  dispersed_dup    hit elsewhere, neither TE nor gene
  novel            no hit passing threshold -> absent from Col-0 entirely

CAVEAT to carry with the `novel` class: "no homology hit" is weak evidence of novelty and gets
weaker the shorter the query, because a short query has little statistical power in an
alignment search. The unclassified rate MUST be read length-stratified, which the summary does.

Env: novelseq (blast). Writes results/sv_adaptive/insertion_tair10_class.csv.
"""
import os, sys, glob, subprocess, gzip
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SEQ = f"{OUT}/seq"; BL = f"{OUT}/blast"
REF = "/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa"
ENV = "/global/home/users/tbellg/miniforge3/envs/novelseq/bin"
CHROMS = [f"Chr{i}" for i in range(1, 6)]
THREADS = int(os.environ.get("KMATE_THREADS", "16"))
LOCAL_BP = 10_000
MIN_COV = 0.50          # fraction of the inserted sequence covered by TAIR10 hits
FMT = "6 qseqid qlen sseqid sstart send pident length evalue bitscore"


def load_intervals():
    """TAIR10 TE intervals (with family alias) and gene intervals, per chrom."""
    te, gn = {c: [] for c in CHROMS}, {c: [] for c in CHROMS}
    with open(lib.TAIR10_GENES_TE) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in CHROMS:
                continue
            s, e = int(f[3]), int(f[4])
            if f[2] == "transposable_element":
                fam = ""
                for p in f[8].split(";"):
                    if p.startswith("Alias="):
                        fam = p[6:]; break
                te[f[0]].append((s, e, fam))
            elif f[2] == "gene":
                gn[f[0]].append((s, e))
    out = {}
    for nm, d in (("te", te), ("gene", gn)):
        out[nm] = {}
        for c, v in d.items():
            v.sort()
            out[nm][c] = (np.array([x[0] for x in v]), np.array([x[1] for x in v]),
                          [x[2] for x in v] if nm == "te" else None)
    print(f"  [tair10] TE {sum(len(out['te'][c][0]) for c in CHROMS):,}  "
          f"genes {sum(len(out['gene'][c][0]) for c in CHROMS):,}")
    return out


def hit_iv(pos, iv):
    st, en, fam = iv
    if st.size == 0:
        return None
    j = np.searchsorted(st, pos, side="right")
    lo = max(0, j - 500)
    for i in range(j - 1, lo - 1, -1):
        if en[i] >= pos:
            return fam[i] if fam is not None else ""
    return None


def main():
    os.makedirs(BL, exist_ok=True)
    if not os.path.exists(f"{REF}.nin") and not os.path.exists(f"{BL}/tair10.nin"):
        print("[db] makeblastdb")
        subprocess.run([f"{ENV}/makeblastdb", "-in", REF, "-dbtype", "nucl",
                        "-out", f"{BL}/tair10"], check=True,
                       stdout=subprocess.DEVNULL)
    db = f"{BL}/tair10"

    for ch in CHROMS:
        o = f"{BL}/{ch}.tsv"
        if os.path.exists(o) and os.path.getsize(o) > 0:
            continue
        print(f"[blast] {ch}", flush=True)
        subprocess.run([f"{ENV}/blastn", "-task", "dc-megablast",
                        "-query", f"{SEQ}/insertions_{ch}.fa", "-db", db,
                        "-outfmt", FMT, "-evalue", "1e-10", "-perc_identity", "70",
                        "-max_target_seqs", "20", "-num_threads", str(THREADS),
                        "-out", o], check=True)

    print("[parse] loading TAIR10 annotations")
    IV = load_intervals()
    D = pd.read_csv(f"{OUT}/insertion_context.csv")
    D["key"] = D.chrom + "|" + D.rec.astype(str)
    site = dict(zip(D.key, zip(D.chrom, D.pos)))

    best, covsum = {}, {}
    for ch in CHROMS:
        f = f"{BL}/{ch}.tsv"
        if not os.path.exists(f):
            continue
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            q, ql, sc = p[0], int(p[1]), p[2]
            ss, se, ident, alen, bits = int(p[3]), int(p[4]), float(p[5]), int(p[6]), float(p[8])
            covsum[q] = covsum.get(q, 0) + alen
            cur = best.get(q)
            if cur is None or bits > cur[0]:
                best[q] = (bits, sc, min(ss, se), max(ss, se), ql, ident, alen)
        print(f"  parsed {ch}", flush=True)

    rows = []
    for k, (ch0, pos0) in site.items():
        b = best.get(k)
        if b is None:
            rows.append(dict(key=k, cls="novel", family="", hit_chrom="", hit_pos=-1,
                             pident=np.nan, qcov=0.0))
            continue
        bits, sc, ss, se, ql, ident, alen = b
        qcov = min(covsum.get(k, 0) / max(ql, 1), 1.0)
        if qcov < MIN_COV:
            rows.append(dict(key=k, cls="novel", family="", hit_chrom=sc, hit_pos=ss,
                             pident=ident, qcov=qcov))
            continue
        mid = (ss + se) // 2
        fam = hit_iv(mid, IV["te"].get(sc, (np.array([]), np.array([]), []))) if sc in CHROMS else None
        if fam is not None:
            cls, family = "te_derived", fam
        elif sc in CHROMS and hit_iv(mid, IV["gene"][sc]) is not None:
            cls, family = "gene_dup", ""
        elif sc == ch0 and abs(mid - pos0) <= LOCAL_BP:
            cls, family = "local_dup", ""
        else:
            cls, family = "dispersed_dup", ""
        rows.append(dict(key=k, cls=cls, family=family, hit_chrom=sc, hit_pos=mid,
                         pident=ident, qcov=qcov))

    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/insertion_tair10_class.csv", index=False)
    M = D.merge(R, on="key", how="left")
    print(f"\n[class] {len(R):,} insertions")
    print((M.cls.value_counts(normalize=True) * 100).round(2).to_string())
    print("\nlength-stratified 'novel' rate (the key caveat):")
    for lo, hi in [(51, 100), (101, 300), (301, 1000), (1001, 5000), (5001, 10**9)]:
        m = (M["size"] >= lo) & (M["size"] <= hi)
        if m.sum():
            print(f"  {lo:>6}-{hi if hi < 10**9 else 'max':<6} n={int(m.sum()):>7,}  "
                  f"novel {100*(M.loc[m, 'cls'] == 'novel').mean():>5.1f}%  "
                  f"te_derived {100*(M.loc[m, 'cls'] == 'te_derived').mean():>5.1f}%")
    print("\ntop TE families:")
    print(M.loc[M.cls == "te_derived", "family"].value_counts().head(15).to_string())
    print(f"\n[wrote] {OUT}/insertion_tair10_class.csv")


if __name__ == "__main__":
    main()
