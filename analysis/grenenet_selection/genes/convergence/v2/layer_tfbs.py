#!/usr/bin/env python
"""v2 layer: REF-vs-ALT TF-motif turnover per UNIT (same call as tfbs_pool.py --scope extra: tfbs_turnover.run_one,
pad 400, FIMO p<1e-4, q<0.005, low-complexity dropped, repeat-tract matches separated).
Sequences come from the unit's representative record (units.rep_alleles); no size cap (as in the old code).
  --worker i --nworkers n   scan a slice -> tfbs/_summary_<i>.csv
  --merge                   -> v2_layer_tfbs.csv
Scope: ftier in F3_UTR/F4_promoter, plus in_atac units among ftier F1..F6.
"""
import os, sys, glob, argparse, time
import numpy as np, pandas as pd, pysam
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, CONV)
import units as UN                                              # noqa: E402
import tfbs_turnover as TF                                      # noqa: E402
import repeat_context as RC                                     # noqa: E402
R = f"{CONV}/results/candidate_list_v2_events_cr50"
TF.OUT = f"{R}/tfbs"
COLS = ["unit", "ftier", "target_gene", "tfbs_n_ref", "tfbs_lost", "tfbs_gained", "tfbs_lost_families", "tfbs_repeat_hits",
        "tfbs_lost_motifs", "tfbs_gained_motifs", "tfbs_n_families", "tfbs_sites_per_motif", "tfbs_tandem",
        "tfbs_redundant", "tfbs_microsat", "tfbs_repeat"]


def scope():
    F = pd.read_csv(f"{R}/v2_functional_units.csv"); A = pd.read_csv(f"{R}/v2_layer_atac.csv")
    F = F.merge(A[["unit", "in_atac"]], on="unit", how="left")
    F["in_atac"] = F.in_atac.fillna(False).astype(bool)
    reg = F.ftier.isin(["F3_UTR", "F4_promoter"])
    atac = F.in_atac & F.ftier.str[:2].isin([f"F{i}" for i in range(1, 7)])
    return F[reg | atac][["unit", "ftier", "target_gene", "vclass"]].reset_index(drop=True)


def clean_counts(tag):
    H = pd.read_csv(f"{TF.OUT}/{tag}_hits.csv"); f = f"{TF.OUT}/{tag}_turnover.csv"
    T = pd.read_csv(f) if os.path.exists(f) and os.path.getsize(f) > 1 else pd.DataFrame()
    rep = set(map(tuple, H[H.kept & H.repeat_tract][["motif", "strand", "gstart"]].values))
    if not len(T):
        return 0, 0, "", T, H
    raw = T                                                     # distinct-motif columns use the unfiltered turnover (as tfbs_pool.resummarize)
    T = T[[(m, s, g) not in rep for m, s, g in zip(T.motif, T.strand, T.gstart)]]
    L, G = T[T.effect == "LOST"], T[T.effect == "GAINED"]
    return len(L), len(G), ",".join(sorted(set(L.family.astype(str)))), raw, H


def main(worker, nw, resum=False):
    S = scope().iloc[worker::nw]
    A = UN.rep_alleles(S.unit.values).set_index("unit")
    fa = pysam.FastaFile(TF.REF_FA); rows = []
    for d in S.itertuples():
        a = A.loc[d.unit]; tag = f"unit_{d.unit}"
        row = dict(unit=d.unit, ftier=d.ftier, target_gene=d.target_gene, skipped="")
        ref_a, alt_a = a.ref.upper(), a.alt.upper()
        pos = int(a.pos)
        if set(ref_a + alt_a) - set("ACGT"):
            row["skipped"] = "non-ACGT allele (symbolic/ambiguous)"; rows.append(row); continue
        g = fa.fetch(a.chrom, pos - 1, pos - 1 + len(ref_a)).upper()
        if g != ref_a:
            row["skipped"] = "REF does not match genome"; rows.append(row); continue
        TF.fetch_allele = lambda *x, _r=ref_a, _a=alt_a: (_r, _a)   # representative record, not a (pos,len) lookup
        rec = pd.Series(dict(symbol=tag, target_gene=d.target_gene, chrom=a.chrom, pos=pos, ref_len=len(ref_a),
                             alt_len=len(alt_a), set="v2", ftier=d.ftier))
        if not (resum and os.path.exists(f"{TF.OUT}/{tag}_hits.csv")):
            try:
                TF.run_one(rec, 400, 1e-4, TF.FIMO, fa)
            except BaseException as e:                          # noqa: BLE001
                row["skipped"] = f"error: {type(e).__name__} {str(e)[:100]}"; rows.append(row); continue
        nl, ng, fam, T, H = clean_counts(tag)
        r = dict(tfbs_n_ref=int((H.kept & (H.seq == "ref")).sum()), tfbs_lost=nl, tfbs_gained=ng, tfbs_lost_families=fam,
                 tfbs_repeat_hits=int((H.kept & H.repeat_tract).sum()), tfbs_lost_motifs=0, tfbs_gained_motifs=0,
                 tfbs_n_families=0, tfbs_sites_per_motif=float("nan"), tfbs_tandem=float("nan"))
        if len(T):
            T = T[T.effect.isin(["LOST", "GAINED"])]
        if len(T):
            rl = len(ref_a)
            fl = fa.fetch(a.chrom, pos - 26, pos - 1).upper(); fr = fa.fetch(a.chrom, pos - 1 + rl, pos - 1 + rl + 25).upper()
            r.update(tfbs_lost_motifs=T[T.effect == "LOST"].motif.nunique(), tfbs_gained_motifs=T[T.effect == "GAINED"].motif.nunique(),
                     tfbs_n_families=T.family.nunique(), tfbs_sites_per_motif=round(len(T) / T.motif.nunique(), 2),
                     tfbs_tandem=round(max(RC.tandem_fraction(fl + ref_a + fr), RC.tandem_fraction(fl + alt_a + fr)), 2))
        row.update(r); rows.append(row)
    pd.DataFrame(rows).to_csv(f"{TF.OUT}/_summary_{worker}.csv", index=False)
    print(f"worker {worker}: {len(rows)} units", flush=True)


def merge():
    S = pd.concat([pd.read_csv(f) for f in glob.glob(f"{TF.OUT}/_summary_*.csv")]).sort_values("unit")
    S["tfbs_redundant"] = S.tfbs_sites_per_motif >= RC.SPM_MAX
    S["tfbs_microsat"] = S.tfbs_tandem >= RC.TANDEM_MIN
    S["tfbs_repeat"] = S.tfbs_redundant | S.tfbs_microsat
    sc = scope().set_index("unit").vclass
    S["vclass"] = S.unit.map(sc)
    S = S.rename(columns={"skipped": "skipped_reason"})
    S[COLS + ["vclass", "skipped_reason"]].to_csv(f"{R}/v2_layer_tfbs.csv", index=False)
    ok = S.skipped_reason.fillna("") == ""
    print(len(S), "units;", int(ok.sum()), "scanned;", int((~ok).sum()), "skipped")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--worker", type=int, default=0); ap.add_argument("--nworkers", type=int, default=1)
    ap.add_argument("--merge", action="store_true"); ap.add_argument("--n", action="store_true"); ap.add_argument("--resum", action="store_true")
    a = ap.parse_args()
    if a.n:
        s = scope(); print(len(s), s.vclass.value_counts().to_dict())
    else:
        merge() if a.merge else main(a.worker, a.nworkers, a.resum)
