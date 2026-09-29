#!/usr/bin/env python
"""Cargo enrichment recomputed on MERGED SITES, not motif hits.

Audit finding (2026-09-20). `cargo/tfbs_insertion_payload.py` counts FIMO hits on the real
and shuffled sides alike. The ratio is a fair comparison of hit counts, but hits are not
independent: one element is matched by every motif in its family at every offset, so a
single real ABRE becomes dozens of counts while a dinucleotide shuffle, which destroys
clustered elements, rarely produces such a pile. Recomputing six loci with overlapping hits
merged showed the result does not survive -- Chr5:19,636,028 (SKS3) goes from 2.10x at
p = 0.03 to 1.07x at p = 0.40 (family-merged) or 0.68x at p = 0.99 (all-motif merged).

This rebuilds the statistic for all 193 loci with two merge conventions:

  fam_sites  overlapping hits of the SAME TF family merged -- one element of one family is
             one site. Matches how the dossier's panel F already displays them. Primary.
  any_sites  overlapping hits merged regardless of family. A harsher bound: two different
             TFs whose sites overlap by 1 bp collapse into one. Sensitivity check.

Everything else is kept identical to theirs: the same inserted sequence (sv_content.
inserted_part, which strips the shared prefix AND suffix -- a naive alt[len(ref):] is wrong
for records that are not left-anchored), the same 100 dinucleotide shuffles, the same FIMO
threshold and the same complexity/repeat/telomere filters. Real and shuffles go in ONE FIMO
input per locus so the background model is identical on both sides.

Writes results/cargo_sites_null.csv. env: kmate + meme. Compute node.
"""
import os, sys, subprocess, tempfile
from multiprocessing import Pool
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, f"{GEA}/cargo")
import tfbs_turnover as TF, sv_content as SC                      # noqa: E402
from tfbs_turnover import MOTIFS, FIMO, FAMILY, low_complexity, repeat_tract   # noqa: E402
from tfbs_insertion_payload import long_repeat_tract, TELO        # noqa: E402

OUT = f"{HERE}/results"
CARGO = f"{GEA}/cargo/results/tfbs_payload"
SHUF = "/global/home/users/tbellg/miniforge3/envs/meme/bin/fasta-shuffle-letters"
N = 100
_fam = pd.read_csv(FAMILY, sep="\t").rename(columns={"Gene_id": "motif_id", "Family": "tf_family"})
FAMOF = dict(zip(_fam.motif_id.astype(str), _fam.tf_family.astype(str)))


def merged(iv) -> int:
    if not len(iv):
        return 0
    iv = sorted(iv); n = 1; end = iv[0][1]
    for s, e in iv[1:]:
        if s > end:
            n += 1; end = e
        else:
            end = max(end, e)
    return n


def counts(g) -> tuple:
    return (merged(list(zip(g.start, g.stop))),
            sum(merged(list(zip(x.start, x.stop))) for _, x in g.groupby("tf_family")))


def one(args) -> dict:
    lid, chrom, pos, rl, al = args
    try:
        ref, alt = TF.fetch_allele(chrom, int(pos), int(rl), int(al))
        ins = SC.inserted_part(ref, alt) if ref else ""
        if len(ins) < 20:
            return dict(locus_id=lid, ins_bp=len(ins))
        with tempfile.TemporaryDirectory() as td:
            open(f"{td}/r.fa", "w").write(f">real\n{ins}\n")
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
        return dict(locus_id=lid, ins_bp=len(ins),
                    hits_obs=len(real), hits_exp=round(h.mean(), 2),
                    hits_enrich=round(he, 2) if he == he else np.nan, hits_p=round(hp, 4),
                    any_sites_obs=o_any, any_sites_exp=round(a.mean(), 2),
                    any_sites_enrich=round(ae, 2) if ae == ae else np.nan, any_sites_p=round(ap, 4),
                    fam_sites_obs=o_fam, fam_sites_exp=round(f.mean(), 2),
                    fam_sites_enrich=round(fe, 2) if fe == fe else np.nan, fam_sites_p=round(fp, 4),
                    shuf_zero_reps=int((h == 0).sum()))
    except Exception as e:                                        # noqa: BLE001
        return dict(locus_id=lid, error=str(e)[:120])


def main():
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else len(os.sched_getaffinity(0))
    L = pd.read_csv(f"{CARGO}/payload_loci.csv")
    jobs = list(zip(L.locus_id, L.chrom, L.pos, L.ref_len, L.alt_len))
    print(f"{len(jobs)} loci, {procs} procs", flush=True)
    with Pool(procs) as p:
        rows = p.map(one, jobs)
    R = pd.DataFrame(rows).merge(
        L[["locus_id", "site_cat", "nearest_gene", "atac_peak", "content_class",
           "enrich", "emp_p"]].rename(columns={"enrich": "their_enrich", "emp_p": "their_p"}),
        on="locus_id", how="left")
    R.to_csv(f"{OUT}/cargo_sites_null.csv", index=False)
    ok = R[R.fam_sites_p.notna()]
    print(f"\nloci recomputed: {len(ok)}")
    for c, lab in [("hits_p", "hits (their statistic)"), ("fam_sites_p", "family-merged sites"),
                   ("any_sites_p", "all-motif-merged sites")]:
        print(f"  p <= 0.05 on {lab}: {int((ok[c] <= 0.05).sum())}/{len(ok)}"
              f"   median enrichment {ok[c.replace('_p','_enrich')].median():.2f}")
    print(f"\nwrote {OUT}/cargo_sites_null.csv")


if __name__ == "__main__":
    main()
