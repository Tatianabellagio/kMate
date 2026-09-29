#!/usr/bin/env python
"""Is the candidates' cargo different from what an ORDINARY insertion carries?

The question `cargo/tfbs_insertion_payload.py` asked -- and the one `cargo_sites_null.py`
and `cargo_null_control.py` between them destroyed -- was "does this inserted sequence carry
more TF motif sites than its OWN base composition predicts?", scored against a
dinucleotide shuffle. That null is not usable: on 1,158 signal-free length-matched windows
its false-positive rate runs 9.5% -> 70.2% from <200 bp to >2 kb, because the shuffle's bias
is per-base and constant while its sampling noise shrinks with length.

This asks a different question, and it never touches a shuffle:

    do the 193 GEA/GWAS-hit insertions carry different cargo from
    length-, class- and frequency-matched insertions that are NOT hits?

Both arms are real genomic insertion sequence, so the length artefact cannot arise -- it
would have to act equally on both. `r1_sv_negative_selection` characterised every arch3
insertion (172,220 records), which is the comparison set; 2,337 records sitting at a pooled
GEA/GWAS hit position are excluded from it.

MATCHING, per candidate: same BLAST class `cls` (te_derived / gene_dup / local_dup /
dispersed_dup / novel -- sequence class is the strongest confounder for motif content), same
allele-frequency stratum (candidates ran 0.004-0.78; rare insertions are younger and
compositionally different), and inserted length within +/-25%, widened to +/-50% then
+/-100% only if fewer than K matches exist. Sampling is without replacement across the
whole draw, so no background record serves two candidates.

SCORED AS PRESENCE, not enrichment. A canonical stress cis-element is identified by its
CORE sequence inside a family-restricted motif hit, not by "a PWM matched": ABRE/G-box =
ACGTG in a bZIP hit, W-box = TTGAC in a WRKY hit, DRE/CRT = [AG]CCGAC in an ERF hit.
Overlapping hits are merged, so one element is one site however many motifs name it -- the
mistake that produced the original result. The cargo tree's complexity / repeat / telomere
filters are kept, identically on both arms.

TEST. Each candidate and its K matches form an exchangeable stratum under H0, so the test is
a within-stratum permutation: relabel one of the K+1 sequences in each stratum as the
pseudo-candidate, 20,000 times, and read the p-value off that null. This is used rather than
a rank statistic because the outcome is heavily zero-inflated -- 103 of 193 strata have no
ABRE anywhere in all 11 sequences -- and a mid-rank percentile then has null expectation
(K+2)/2/(K+1) = 0.545, not 0.5. An earlier version of this script compared the mean
percentile to 0.5 and reported p = 0.0009 for W-box cluster size; the permutation gives
p = 0.066 for the same data. Permutation handles the ties by construction.

Note the p-value FLOOR: with K = 10 matches, any statistic driven by a single extreme
candidate cannot go below 1/(K+1) = 0.09, because under H0 that sequence is the labelled
candidate 1 time in 11. To sharpen one locus, deepen its own stratum instead (see the SKS3
follow-up in the README).

Writes results/cargo_vs_background_{sites,loci,summary}.csv. env: kmate + meme. Compute node.
"""
import os, sys, re, subprocess, tempfile
from multiprocessing import Pool
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, GEA); sys.path.insert(0, f"{GEA}/cargo")
import lib                                                              # noqa: E402
import sv_content as SC                                                 # noqa: E402
import build_sv_cargo as BC                                             # noqa: E402
from tfbs_turnover import MOTIFS, FIMO, low_complexity, repeat_tract    # noqa: E402
from tfbs_insertion_payload import long_repeat_tract, TELO              # noqa: E402
from cargo_sites_null import FAMOF                                      # noqa: E402

OUT = f"{HERE}/results"
PAY = f"{GEA}/cargo/results/tfbs_payload"
R1 = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
K = 10                        # matched background per candidate
SEED = 20260920
MIN_BP = 20
CHUNK = 60                    # sequences per FIMO call

# canonical stress cis-elements: (TF family, core regex). Family-restricted so a core that
# merely occurs in the sequence does not count unless a motif of the right family calls it.
ELEMENTS = {"ABRE_Gbox": ("bZIP", r"ACGTG"),
            "Wbox":      ("WRKY", r"TTGAC"),
            "DRE_CRT":   ("ERF",  r"[AG]CCGAC")}
FREQ_BINS = [0, 0.02, 0.05, 0.15, 0.30, 1.01]


def merged(iv):
    if not len(iv):
        return []
    iv = sorted(iv); out = []; cs, ce = iv[0]
    for a, b in iv[1:]:
        if a > ce:
            out.append((cs, ce)); cs, ce = a, b
        else:
            ce = max(ce, b)
    out.append((cs, ce)); return out


def fimo_chunk(args):
    """FIMO over one chunk of sequences; returns the filtered hit table."""
    idx, seqs = args
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/q.fa", "w") as f:
            for n, s in seqs:
                f.write(f">{n}\n{s}\n")
        r = subprocess.run([FIMO, "--text", "--thresh", "1e-4", "--verbosity", "1",
                            "--max-stored-scores", "10000000", MOTIFS, f"{td}/q.fa"],
                           capture_output=True, text=True, check=True)
    if not r.stdout.strip():
        return pd.DataFrame()
    from io import StringIO
    d = pd.read_csv(StringIO(r.stdout), sep="\t", comment="#")
    d = d.dropna(subset=["motif_id", "sequence_name"]).rename(columns={"p-value": "pval"})
    d = d[pd.to_numeric(d.pval, errors="coerce") < 1e-4]
    s = d.matched_sequence.astype(str)
    d = d[~s.map(low_complexity) & ~s.map(repeat_tract) & ~s.map(long_repeat_tract)
          & ~s.map(lambda x: bool(TELO.search(x)))].copy()
    d["tf_family"] = d.motif_id.astype(str).map(FAMOF).fillna("?")
    return d


def score(d: pd.DataFrame, names: list) -> pd.DataFrame:
    """Per sequence: presence and largest merged cluster of each canonical element."""
    rows = {n: {"seq_name": n} for n in names}
    for n in names:
        for e in ELEMENTS:
            rows[n][f"{e}_n"] = 0; rows[n][f"{e}_bp"] = 0
    if len(d):
        up = d.matched_sequence.astype(str).str.upper()
        # Some PlantTFDB matrices are defined on the opposite strand, so FIMO's
        # matched_sequence carries the core reverse-complemented -- at Chr5:19,636,028 the
        # same 44 bp ABRE shows up both as GATGATGACGTGGCA (ACGTG) and as TGCCACGTCATCATC
        # (CACGT). Testing one orientation only would drop those hits, and drop them
        # unevenly across elements. So test the core against both orientations.
        rc = up.str.translate(str.maketrans("ACGT", "TGCA")).str[::-1]
        for e, (fam, core) in ELEMENTS.items():
            has = up.str.contains(core, regex=True) | rc.str.contains(core, regex=True)
            h = d[(d.tf_family == fam) & has]
            for n, g in h.groupby("sequence_name"):
                m = merged(list(zip(g.start, g.stop)))
                if str(n) in rows:
                    rows[str(n)][f"{e}_n"] = len(m)
                    rows[str(n)][f"{e}_bp"] = max((b - a + 1 for a, b in m), default=0)
    return pd.DataFrame(rows.values())


def build_sets():
    """Candidates joined to their panel record, plus K matched non-hit background each."""
    D = BC.insertion_records()
    V = pd.read_csv(f"{OUT}/variants_classified.csv")
    hits = set(zip(V.chrom.astype(str), V.pos.astype(int)))
    D["is_hit"] = [(c, int(p)) in hits for c, p in zip(D.chrom.astype(str), D.pos)]
    D["fbin"] = pd.cut(D.freq, FREQ_BINS, labels=False, include_lowest=True)

    L = pd.read_csv(f"{PAY}/payload_loci.csv", keep_default_na=False)
    DS = pd.read_csv(f"{OUT}/sv_hit_dossier.csv", keep_default_na=False)
    DS["locus_id"] = DS.chrom + "_" + DS.pos.astype(str)
    C = (L[["locus_id", "ins_bp", "site_cat", "content_class", "atac_peak"]]
         .merge(DS[["locus_id", "dom_key", "host_gene", "host_gene_basis"]], on="locus_id", how="left")
         .merge(D[["key", "cls", "freq", "fbin", "ref_len", "alt_len", "chrom", "pos"]]
                .rename(columns={"key": "dom_key"}), on="dom_key", how="left"))
    assert C.cls.notna().all(), "a candidate did not join to a panel record"

    BG = D[~D.is_hit].copy()
    rng = np.random.default_rng(SEED)
    used, picks = set(), []
    for r in C.itertuples():
        pool = BG[(BG.cls == r.cls) & (BG.fbin == r.fbin) & ~BG.key.isin(used)]
        for tol in (0.25, 0.50, 1.00):
            cand = pool[(pool["size"] >= (1 - tol) * r.ins_bp) & (pool["size"] <= (1 + tol) * r.ins_bp)]
            if len(cand) >= K:
                break
        take = cand if len(cand) <= K else cand.iloc[rng.choice(len(cand), K, replace=False)]
        used |= set(take.key)
        for t in take.itertuples():
            picks.append(dict(locus_id=r.locus_id, key=t.key, chrom=t.chrom, pos=t.pos,
                              ref_len=t.ref_len, alt_len=t.alt_len, cls=t.cls, freq=t.freq,
                              match_tol=tol, genic=t.genic, te_overlap=t.te_overlap))
    B = pd.DataFrame(picks)
    print(f"candidates {len(C)}; background {len(B)} "
          f"({B.groupby('locus_id').size().min()}-{B.groupby('locus_id').size().max()} per candidate); "
          f"tolerance used: {B.match_tol.value_counts().sort_index().to_dict()}", flush=True)
    return C, B


def sequences(C, B):
    """inserted_part for every candidate and background record, via the panel FASTA."""
    keys = set(C.dom_key) | set(B.key)
    ref = SC.load_ref(); alt = SC.load_alt(f"{R1}/seq", keys)
    out = {}
    for k, ch, p, rl in list(zip(C.dom_key, C.chrom, C.pos, C.ref_len)) + \
                        list(zip(B.key, B.chrom, B.pos, B.ref_len)):
        ins = SC.inserted_part(ref[ch][int(p) - 1:int(p) - 1 + int(rl)], alt[k])
        if len(ins) >= MIN_BP:
            out[k] = ins
    return out


NPERM = 20000


def permute(M, rng):
    """Within-stratum permutation. M is (n_strata, K+1) with column 0 the real candidate.
    Returns observed and null for presence count, mean cluster bp and max cluster bp."""
    n, kp1 = M.shape
    obs = dict(present=float((M[:, 0] > 0).sum()), mean_bp=float(M[:, 0].mean()),
               max_bp=float(M[:, 0].max()))
    idx = rng.integers(0, kp1, size=(NPERM, n))
    P = M[np.arange(n)[None, :], idx]
    null = dict(present=(P > 0).sum(1), mean_bp=P.mean(1), max_bp=P.max(1))
    out = {}
    for k, v in obs.items():
        out[k] = (v, float(null[k].mean()), (1 + int((null[k] >= v).sum())) / (NPERM + 1))
    return out


def main():
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else len(os.sched_getaffinity(0))
    C, B = build_sets()
    seqs = sequences(C, B)
    print(f"{len(seqs)} sequences with >= {MIN_BP} bp of new DNA; {procs} procs", flush=True)

    items = [(k.replace("|", "_"), s) for k, s in seqs.items()]
    chunks = [(i, items[i:i + CHUNK]) for i in range(0, len(items), CHUNK)]
    with Pool(procs) as p:
        parts = p.map(fimo_chunk, chunks)
    d = pd.concat([x for x in parts if len(x)], ignore_index=True) if any(len(x) for x in parts) else pd.DataFrame()
    print(f"filtered motif hits: {len(d)}", flush=True)
    Sc = score(d, [n for n, _ in items])
    Sc["key"] = Sc.seq_name.str.replace("_", "|", regex=False)
    Sc.to_csv(f"{OUT}/cargo_vs_background_sites.csv", index=False)

    Cs = C.merge(Sc, left_on="dom_key", right_on="key", how="left", suffixes=("", "_s"))
    Bs = B.merge(Sc, on="key", how="left")
    Cs["arm"] = "candidate"; Bs["arm"] = "background"
    keep = ["locus_id", "arm", "key", "cls", "freq"] + [f"{e}_{s}" for e in ELEMENTS for s in ("n", "bp")]
    allrows = pd.concat([Cs.assign(key=Cs.dom_key)[keep], Bs[keep]], ignore_index=True)
    allrows.to_csv(f"{OUT}/cargo_vs_background_loci.csv", index=False)

    print("\n=== candidate vs length/class/frequency-matched NON-HIT insertions ===")
    print("(presence of a canonical stress cis-element in the inserted DNA; no shuffle null)")
    print(f"within-stratum permutation, {NPERM} draws\n")
    rng = np.random.default_rng(SEED)
    rows = []
    for e in ELEMENTS:
        mat = []
        for lid, g in allrows.groupby("locus_id"):
            c = g[g.arm == "candidate"][f"{e}_bp"].to_numpy(float)
            b = g[g.arm == "background"][f"{e}_bp"].to_numpy(float)
            if len(c) == 1 and len(b) == K:
                mat.append(np.concatenate([c, b]))
        M = np.array(mat)
        r = permute(M, rng)
        rows.append(dict(element=e, n_strata=len(M),
                         cand_present=int(r["present"][0]), null_present=round(r["present"][1], 1),
                         p_presence=round(r["present"][2], 4),
                         cand_mean_bp=round(r["mean_bp"][0], 2), null_mean_bp=round(r["mean_bp"][1], 2),
                         p_mean_bp=round(r["mean_bp"][2], 4),
                         cand_max_bp=int(r["max_bp"][0]), null_max_bp=round(r["max_bp"][1], 1),
                         p_max_bp=round(r["max_bp"][2], 4)))
    T = pd.DataFrame(rows)
    print(T.to_string(index=False))
    alpha = 0.05 / (3 * len(ELEMENTS))
    sig = [(r.element, c) for r in T.itertuples()
           for c in ("p_presence", "p_mean_bp", "p_max_bp") if getattr(r, c) < alpha]
    print(f"\nBonferroni over {len(ELEMENTS)} elements x 3 statistics: alpha = {alpha:.4f}")
    print(f"  surviving: {sig if sig else 'NONE'}")
    T.to_csv(f"{OUT}/cargo_vs_background_summary.csv", index=False)

    print("\n=== how unusual is the SKS3 configuration? ===")
    bg = allrows[allrows.arm == "background"]
    for thr in (20, 30, 44):
        n = int((bg.ABRE_Gbox_bp >= thr).sum())
        print(f"  background insertions with an ABRE cluster >= {thr} bp: {n}/{len(bg)} = {100*n/len(bg):.2f}%")
    c = allrows[(allrows.arm == "candidate")]
    print(f"  candidates with an ABRE cluster >= 44 bp: {int((c.ABRE_Gbox_bp >= 44).sum())}/{len(c)}")
    print(f"\nwrote {OUT}/cargo_vs_background_{{sites,loci,summary}}.csv")


if __name__ == "__main__":
    main()
