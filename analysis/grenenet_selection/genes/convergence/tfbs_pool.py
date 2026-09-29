#!/usr/bin/env python
"""REF-vs-ALT motif turnover for every promoter / 5'UTR variant in the evidence matrix.

Until now tfbs_turnover.py had been run on ~17 hand-picked variants, so the mechanism line
for the rest of the pool was position only. Same call (tfbs_turnover.run_one: pad 400,
FIMO p<1e-4, q<0.005, low-complexity dropped), per variant, with per-variant file names
(<gene>_<pos>_<ref>_<alt>) written under results/tfbs/pool/ so genes with several variants
do not overwrite each other. Repeat-tract matches (the AT4G13200 lesson) are separated from
clean turnover.

Writes results/tfbs/pool/_summary_<worker>.csv; merged by --merge into
results/tfbs_pool_summary.csv.
env: kmate + meme (fimo). Compute node.
"""
import os, sys, glob, argparse
import pandas as pd, pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tfbs_turnover as TF                                       # noqa: E402
import repeat_context as RC                                     # noqa: E402

OUT = f"{HERE}/results"
TF.OUT = f"{OUT}/tfbs/pool"


def clean_counts(tag):
    H = pd.read_csv(f"{TF.OUT}/{tag}_hits.csv")
    f = f"{TF.OUT}/{tag}_turnover.csv"
    T = pd.read_csv(f) if os.path.exists(f) and os.path.getsize(f) > 1 else pd.DataFrame()
    rep = set(map(tuple, H[H.kept & H.repeat_tract][["motif", "strand", "gstart"]].values))
    if not len(T):
        return 0, 0, ""
    ok = [(m, s, g) not in rep for m, s, g in zip(T.motif, T.strand, T.gstart)]
    T = T[ok]
    L, G = T[T.effect == "LOST"], T[T.effect == "GAINED"]
    return len(L), len(G), ",".join(sorted(set(L.family.astype(str))))


def main(worker, nworkers):
    os.makedirs(TF.OUT, exist_ok=True)
    D = pd.read_csv(f"{OUT}/evidence_matrix.csv")
    D = D[D["mode"].isin(["promoter", "5'UTR"]) & (D.vclass != "mnp")].reset_index(drop=True)
    D = D.iloc[worker::nworkers]
    fa = pysam.FastaFile(TF.REF_FA)
    rows = []
    for d in D.itertuples():
        tag = f"{d.target_gene}_{d.pos}_{d.ref_len}_{d.alt_len}"
        rec = pd.Series(dict(symbol=tag, target_gene=d.target_gene, chrom=d.chrom, pos=int(d.pos),
                             ref_len=int(d.ref_len), alt_len=int(d.alt_len), set="pool", ftier=""))
        try:
            res = TF.run_one(rec, 400, 1e-4, TF.FIMO, fa)
        except Exception as e:                                     # noqa: BLE001
            print(f"  {tag}: {e}"); res = None
        if not res:
            continue
        nl, ng, fam = clean_counts(tag)
        rows.append(dict(store_row=d.store_row, tfbs_n_ref=res["n_ref_hits"],
                         tfbs_lost=nl, tfbs_gained=ng, tfbs_lost_families=fam,
                         tfbs_repeat_hits=res["n_repeat_tract_hits"]))
    pd.DataFrame(rows).to_csv(f"{TF.OUT}/_summary_{worker}.csv", index=False)
    print(f"worker {worker}: {len(rows)} variants", flush=True)


def resummarize():
    """Rebuild the summary from the saved hits/turnover files, no FIMO rerun.

    Adds distinct-motif counts and the two repeat flags described in repeat_context.py.
    Nothing is dropped: raw site counts stay beside the new columns, and the flag is a
    column the scorer reads.
    """
    D = pd.read_csv(f"{OUT}/evidence_matrix.csv")
    D = D[D["mode"].isin(["promoter", "5'UTR"]) & (D.vclass != "mnp")]
    fa = pysam.FastaFile(TF.REF_FA)
    rows = []
    for d in D.itertuples():
        tag = f"{d.target_gene}_{d.pos}_{d.ref_len}_{d.alt_len}"
        if not os.path.exists(f"{TF.OUT}/{tag}_hits.csv"):
            continue
        nl, ng, fam = clean_counts(tag)
        H = pd.read_csv(f"{TF.OUT}/{tag}_hits.csv")
        f = f"{TF.OUT}/{tag}_turnover.csv"
        T = pd.read_csv(f) if os.path.exists(f) and os.path.getsize(f) > 1 else pd.DataFrame()
        r = dict(store_row=d.store_row, tfbs_n_ref=int((H.kept & (H.seq == "ref")).sum()),
                 tfbs_lost=nl, tfbs_gained=ng, tfbs_lost_families=fam,
                 tfbs_repeat_hits=int((H.kept & H.repeat_tract).sum()),
                 tfbs_lost_motifs=0, tfbs_gained_motifs=0, tfbs_n_families=0,
                 tfbs_sites_per_motif=float("nan"), tfbs_tandem=float("nan"))
        if len(T):
            T = T[T.effect.isin(["LOST", "GAINED"])]
        if len(T):
            ch, pos, rl = d.chrom, int(d.pos), int(d.ref_len)
            ref_a, alt_a = TF.fetch_allele(ch, pos, rl, int(d.alt_len))
            fl = fa.fetch(ch, pos - 26, pos - 1).upper()
            fr = fa.fetch(ch, pos - 1 + rl, pos - 1 + rl + 25).upper()
            r.update(tfbs_lost_motifs=T[T.effect == "LOST"].motif.nunique(),
                     tfbs_gained_motifs=T[T.effect == "GAINED"].motif.nunique(),
                     tfbs_n_families=T.family.nunique(),
                     tfbs_sites_per_motif=round(len(T) / T.motif.nunique(), 2),
                     tfbs_tandem=round(max(RC.tandem_fraction(fl + ref_a + fr),
                                           RC.tandem_fraction(fl + alt_a + fr)), 2))
        rows.append(r)
    S = pd.DataFrame(rows)
    S["tfbs_redundant"] = S.tfbs_sites_per_motif >= RC.SPM_MAX
    S["tfbs_microsat"] = S.tfbs_tandem >= RC.TANDEM_MIN
    S["tfbs_repeat"] = S.tfbs_redundant | S.tfbs_microsat
    S.to_csv(f"{OUT}/tfbs_pool_summary.csv", index=False)
    turn = (S.tfbs_lost + S.tfbs_gained) > 0
    print(f"{len(S)} variants; turnover {int(turn.sum())}; flagged redundant "
          f"{int(S.tfbs_redundant.sum())}, microsatellite {int(S.tfbs_microsat.sum())}, "
          f"turnover and not flagged {int((turn & ~S.tfbs_repeat).sum())}")


def merge():
    S = pd.concat([pd.read_csv(f) for f in glob.glob(f"{TF.OUT}/_summary_*.csv")])
    S.to_csv(f"{OUT}/tfbs_pool_summary.csv", index=False)
    print(f"merged {len(S)} variants; with >=1 clean site lost or gained: "
          f"{int(((S.tfbs_lost + S.tfbs_gained) > 0).sum())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--nworkers", type=int, default=1); ap.add_argument("--merge", action="store_true")
    ap.add_argument("--resummarize", action="store_true")
    a = ap.parse_args()
    if a.resummarize:
        resummarize()
    else:
        merge() if a.merge else main(a.worker, a.nworkers)
