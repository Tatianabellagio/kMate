#!/usr/bin/env python
"""Is the dinucleotide-shuffle null calibrated for LONG sequence? A length-matched control.

`cargo_sites_null.py` (audit 2026-09-20) rebuilt the cargo enrichment on merged SITES
instead of motif hits and left 25 of 193 loci at p <= 0.05. Before re-enabling L_cargo for
those 25, one pattern has to be explained: the significance rate is almost entirely a
function of insertion LENGTH, not of biology.

    <200 bp    4/81   4.9%    <- nominal
    200-500    2/54   3.7%    <- nominal
    500-2 kb   7/39  17.9%
    >2 kb     12/19  63.2%    <- enrichments of only 1.14-1.19x, p at the 1/101 floor

That is the signature of a MIS-SPECIFIED NULL, not of signal. A dinucleotide shuffle
preserves mononucleotide and dinucleotide composition and nothing else; real genomic
sequence also carries higher-order structure (k-mer reuse, degenerate internal repeats, TE
internal architecture) that creates motif matches a dinucleotide shuffle cannot. That bias
is per-base and roughly constant, while the shuffle's sampling noise shrinks as sequence
grows -- so a fixed ~15% excess becomes unbeatable at 5 kb and invisible at 150 bp. Exactly
the observed pattern, and 15 of the 25 survivors are content_class == TE.

The test: run the IDENTICAL statistic on sequence that carries no candidate signal at all.

  arm 'genomic'  random reference windows, length-matched to each real insertion
  arm 'te'       random windows length-matched AND >=50% covered by a TAIR10
                 annotated transposable_element -- isolates "is it the TE-ness?"

Neither arm can contain a cargo signal. So any p <= 0.05 there is false positive, and the
false-positive rate per length bin IS the null's calibration curve. If long control windows
come out "enriched" at the rate the long insertions do, the 25 survivors are the null
breaking down and L_cargo cannot be re-enabled at any length.

Everything downstream of sequence choice is imported from `cargo_sites_null.py` unchanged
(same 100 dinucleotide shuffles, same seed, same FIMO threshold, same complexity / repeat /
telomere filters, same family-merge), so the arms differ only in where the DNA came from.

Writes results/cargo_null_control.csv. env: kmate + meme. Compute node.
"""
import os, sys, subprocess, tempfile
from multiprocessing import Pool
import numpy as np, pandas as pd, pysam

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, f"{GEA}/cargo")
import tfbs_turnover as TF                                              # noqa: E402
from tfbs_turnover import MOTIFS, FIMO, REF_FA, low_complexity, repeat_tract   # noqa: E402
from tfbs_insertion_payload import long_repeat_tract, TELO              # noqa: E402
from cargo_sites_null import counts, FAMOF, SHUF, N                     # noqa: E402
import lib                                                             # noqa: E402

OUT = f"{HERE}/results"
CARGO = f"{GEA}/cargo/results/tfbs_payload"
REPS = 3                      # control windows per locus per arm
SEED = 20260920
ACGT_MIN = 0.99               # reject windows with N / IUPAC ambiguity


def chrom_sizes() -> dict:
    fa = pysam.FastaFile(REF_FA)
    return {c: l for c, l in zip(fa.references, fa.lengths) if c.startswith("Chr")
            and c not in ("ChrM", "ChrC")}


def te_intervals() -> pd.DataFrame:
    """TAIR10 annotated transposable_element features, as a plain interval table."""
    rows = []
    for line in open(lib.TAIR10_GENES_TE):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 4 and f[2] == "transposable_element" and f[0].startswith("Chr"):
            rows.append((f[0], int(f[3]), int(f[4])))
    T = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    T["len"] = T.end - T.start + 1
    return T


def draw(arm: str, want: int, rng, sizes: dict, TE: pd.DataFrame, fa) -> tuple | None:
    """A length-`want` window: anywhere (genomic) or >=50% inside a TE (te)."""
    for _ in range(200):
        if arm == "genomic":
            ch = rng.choice(list(sizes)); lo = int(rng.integers(1, sizes[ch] - want))
        else:
            t = TE.iloc[int(rng.integers(0, len(TE)))]
            need = int(np.ceil(want / 2))                 # >=50% of the window inside the TE
            if t["len"] < need:
                continue
            # window start such that overlap with [t.start, t.end] >= need
            lo_min = max(1, int(t.start) - (want - need))
            lo_max = min(int(t.end) - need + 1, sizes[t.chrom] - want)
            if lo_max < lo_min:
                continue
            ch = t.chrom; lo = int(rng.integers(lo_min, lo_max + 1))
        seq = fa.fetch(ch, lo - 1, lo - 1 + want).upper()
        if len(seq) != want:
            continue
        if sum(seq.count(b) for b in "ACGT") / len(seq) < ACGT_MIN:
            continue
        return ch, lo, seq
    return None


def test_seq(seq: str) -> dict:
    """The cargo statistic, verbatim: 100 dinucleotide shuffles, FIMO, filters, merge."""
    with tempfile.TemporaryDirectory() as td:
        open(f"{td}/r.fa", "w").write(f">real\n{seq}\n")
        with open(f"{td}/s.fa", "w") as o:
            subprocess.run([SHUF, "-kmer", "2", "-copies", str(N), "-seed", "7",
                            f"{td}/r.fa"], stdout=o, check=True)
        open(f"{td}/all.fa", "w").write(open(f"{td}/r.fa").read() + open(f"{td}/s.fa").read())
        subprocess.run([FIMO, "--oc", f"{td}/o", "--thresh", "1e-4", "--verbosity", "1",
                        MOTIFS, f"{td}/all.fa"], capture_output=True, check=True)
        d = pd.read_csv(f"{td}/o/fimo.tsv", sep="\t", comment="#")
    d = d.dropna(subset=["motif_id", "sequence_name"]).rename(columns={"p-value": "pval"})
    d = d[pd.to_numeric(d.pval, errors="coerce") < 1e-4]
    s = d.matched_sequence.astype(str)
    d = d[~s.map(low_complexity) & ~s.map(repeat_tract) & ~s.map(long_repeat_tract)
          & ~s.map(lambda x: bool(TELO.search(x)))].copy()
    d["tf_family"] = d.motif_id.astype(str).map(FAMOF).fillna("?")
    real = d[d.sequence_name == "real"]
    o_any, o_fam = counts(real)
    h = np.zeros(N); a = np.zeros(N); f = np.zeros(N)
    for name, g in d[d.sequence_name != "real"].groupby("sequence_name"):
        i = int(str(name).rsplit("_shuf_", 1)[1]) - 1
        h[i] = len(g); a[i], f[i] = counts(g)
    m = lambda obs, arr: ((obs / arr.mean()) if arr.mean() else np.nan,
                          (1 + int((arr >= obs).sum())) / (N + 1))
    he, hp = m(len(real), h); ae, ap = m(o_any, a); fe, fp = m(o_fam, f)
    return dict(hits_obs=len(real), hits_exp=round(h.mean(), 2),
                hits_enrich=round(he, 2) if he == he else np.nan, hits_p=round(hp, 4),
                any_sites_obs=o_any, any_sites_exp=round(a.mean(), 2),
                any_sites_enrich=round(ae, 2) if ae == ae else np.nan, any_sites_p=round(ap, 4),
                fam_sites_obs=o_fam, fam_sites_exp=round(f.mean(), 2),
                fam_sites_enrich=round(fe, 2) if fe == fe else np.nan, fam_sites_p=round(fp, 4))


def one(job) -> dict:
    lid, arm, rep, ch, lo, want = job
    base = dict(locus_id=lid, arm=arm, rep=rep, chrom=ch, pos=lo, ins_bp=want)
    try:
        fa = pysam.FastaFile(REF_FA)
        seq = fa.fetch(ch, lo - 1, lo - 1 + want).upper()
        return {**base, **test_seq(seq)}
    except Exception as e:                                              # noqa: BLE001
        return {**base, "error": str(e)[:120]}


def main():
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else len(os.sched_getaffinity(0))
    R = pd.read_csv(f"{OUT}/cargo_sites_null.csv")
    R = R[R.fam_sites_p.notna() & (R.ins_bp >= 20)]
    sizes = chrom_sizes(); TE = te_intervals(); fa = pysam.FastaFile(REF_FA)
    rng = np.random.default_rng(SEED)
    print(f"{len(R)} real loci to match; {len(TE)} TAIR10 TEs; {len(sizes)} chroms", flush=True)

    jobs, misses = [], 0
    for r in R.itertuples():
        for arm in ("genomic", "te"):
            for rep in range(1, REPS + 1):
                w = draw(arm, int(r.ins_bp), rng, sizes, TE, fa)
                if w is None:
                    misses += 1; continue
                jobs.append((r.locus_id, arm, rep, w[0], w[1], int(r.ins_bp)))
    print(f"{len(jobs)} control windows ({misses} unplaceable), {procs} procs", flush=True)

    with Pool(procs) as p:
        rows = p.map(one, jobs)
    C = pd.DataFrame(rows)
    C.to_csv(f"{OUT}/cargo_null_control.csv", index=False)

    ok = C[C.fam_sites_p.notna()].copy()
    ok["szbin"] = pd.cut(ok.ins_bp, [0, 200, 500, 2000, 1e9],
                         labels=["<200", "200-500", "500-2k", ">2k"])
    print("\nfalse-positive rate of the cargo statistic on signal-free sequence")
    print("(any p <= 0.05 here is a false positive by construction)\n")
    for stat in ("hits_p", "fam_sites_p", "any_sites_p"):
        print(f"  {stat}")
        t = ok.groupby(["arm", "szbin"], observed=True).agg(
            n=("locus_id", "size"), fpr=(stat, lambda x: (x <= 0.05).mean()),
            med_enrich=(stat.replace("_p", "_enrich"), "median"))
        print(t.to_string(), "\n")
    print(f"wrote {OUT}/cargo_null_control.csv")


if __name__ == "__main__":
    main()
