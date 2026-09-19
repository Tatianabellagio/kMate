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


def merge():
    S = pd.concat([pd.read_csv(f) for f in glob.glob(f"{TF.OUT}/_summary_*.csv")])
    S.to_csv(f"{OUT}/tfbs_pool_summary.csv", index=False)
    print(f"merged {len(S)} variants; with >=1 clean site lost or gained: "
          f"{int(((S.tfbs_lost + S.tfbs_gained) > 0).sum())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--nworkers", type=int, default=1); ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()
    merge() if a.merge else main(a.worker, a.nworkers)
