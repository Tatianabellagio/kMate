#!/usr/bin/env python
"""Recover WHICH allele each per-garden GWAS hit is, at multiallelic positions.

`persite_gwas_{cls}.npz` stores chrom/pos only. At the 2.14% of arch3 positions carrying more
than one biallelic record (memory `panel-multiallelic-pos-key-trap`) that key is ambiguous:
`build_gwas_pool.py` assigns the LARGEST record at the position, and
`plot_variant_garden_grid.sig_gardens` reads the FIRST npz row at the position. Neither is
necessarily the allele that was significant. Found 2026-09-16 when the CYP71B4 / CML25 /
EPFL5 grids showed 0/30 significant gardens for alleles the pool listed as 3-4 garden hits.

The rows can be recovered exactly: `gemma_persite_gwas.scan` tests `pos[idx]` with
idx = where(MAC >= MIN_MAC & call >= CALL_MIN & class), in panel order, and `assemble`
concatenates the per-chrom parts in CHROMS order, then drops markers untestable everywhere.
Replaying the same filter on the panel meta gives ref_len/alt_len per npz row; the positions
are checked against the stored ones before anything is written.

Outputs -> results/
  gwas_hits_allele_resolved.csv          one row per (class, chrom, pos, ref_len, alt_len)
      significant in >= 1 garden, with n_gardens / gardens / best_nlp for THAT allele,
      n_alleles_tested_at_pos and n_sig_alleles_at_pos (a position can carry more than
      one significant allele).
  gwas_hits_allele_resolved_records.csv  one row per (allele x significant garden) --
      the allele-resolved replacement for `persite_bonferroni_hits.csv`, which is keyed
      on position only. Carries garden, bio1, nlp, and the panel MAC/MAF of THAT allele.

env: kmate. Compute node (loads the per-chrom panel presence matrices). ~5 min.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402
R3 = f"{lib.GEA}/r3_persite_gwas"
sys.path.insert(0, R3)
from class_split_gwas import _panel_order, _load_chrom_212, _classes    # noqa: E402
import gemma_persite_gwas as gp                                         # noqa: E402

GWASD = f"{R3}/results/gemma_gwas"
BONF = {"nonsnp": 7.021, "sv": 5.408}
OUT = f"{HERE}/results/gwas_hits_allele_resolved.csv"


def main():
    order, _, _, _ = _panel_order()
    N = len(order)
    per_chrom = {}
    for cl in gp.CHROMS:
        pos, rl, al, _, _, n_alt, n_cal = _load_chrom_212(cl, order)
        mac = np.minimum(n_alt, N - n_alt)
        keep = (mac >= gp.MIN_MAC) & (n_cal >= gp.CALL_MIN * N)
        per_chrom[cl] = (pos, rl, al, _classes(rl, al), keep, mac, n_cal)
        print(f"{cl}: {keep.sum():,} testable records", flush=True)

    rows, rec = [], []
    for cls in ("nonsnp", "sv"):
        z = np.load(f"{GWASD}/persite_gwas_{cls}.npz", allow_pickle=True)
        P_, RL, AL, KR, MC, NC = [], [], [], [], [], []
        for cl in gp.CHROMS:
            pos, rl, al, vc, keep, mac, n_cal = per_chrom[cl]
            want = {"nonsnp": vc != 0, "sv": vc == 2}[cls]
            idx = np.where(keep & want)[0]
            part = np.load(f"{GWASD}/parts/{cl}_{cls}.npz", allow_pickle=True)
            if not np.array_equal(part["pos"], pos[idx]):
                sys.exit(f"{cl}/{cls}: replayed marker set does not match the parts file")
            fin = np.isfinite(part["Z"]).any(1)
            P_.append(pos[idx][fin]); RL.append(rl[idx][fin]); AL.append(al[idx][fin])
            # rank of each record among ALL panel records sharing (pos, ref_len, alt_len):
            # same-length ALT alleles collide on that key; the af_store lists them in the
            # same order (checked 2026-09-17 via p0 vs panel carrier fraction), so the rank
            # picks the store row
            kk = pd.Series(list(zip(pos, rl, al)))
            rank = kk.groupby(kk).cumcount().to_numpy()
            KR.append(rank[idx][fin]); MC.append(mac[idx][fin]); NC.append(n_cal[idx][fin])
        pos_all, rl_all, al_all, kr_all, mac_all, ncal_all = map(
            np.concatenate, (P_, RL, AL, KR, MC, NC))
        if not np.array_equal(pos_all, z["pos"]):
            sys.exit(f"{cls}: replayed positions do not match persite_gwas_{cls}.npz")
        nlp = -np.log10(np.clip(z["P"], 1e-300, None))
        sig = nlp >= BONF[cls]
        df = pd.DataFrame({"cls": cls, "chrom": z["chrom"], "pos": z["pos"],
                           "ref_len": rl_all, "alt_len": al_all, "key_rank": kr_all,
                           "mac": mac_all.astype(int),
                           # MAF over ALL founders, as persite_bonferroni_hits.csv does
                           # (mac / N, not mac / n_called)
                           "maf": np.round(mac_all / N, 4)})
        key = df.chrom + ":" + df.pos.astype(str)
        df["n_alleles_tested_at_pos"] = key.map(key.value_counts())
        hit = np.where(sig.any(1))[0]
        sites = z["sites"].astype(int)
        bio1 = z["bio1"] if "bio1" in z.files else np.full(len(sites), np.nan)
        for i in hit:
            s = sites[sig[i]]
            j = int(np.nanargmax(nlp[i]))
            rows.append({**df.iloc[i].to_dict(), "n_gardens": int(sig[i].sum()),
                         "gardens": ",".join(map(str, sorted(s))),
                         "best_nlp": round(float(nlp[i, j]), 2), "best_garden": int(sites[j]),
                         "best_Z": round(float(z["Z"][i, j]), 2),
                         "Z_sig_gardens": ",".join(f"{v:.1f}" for v in z["Z"][i][sig[i]])})
            for si in np.where(sig[i])[0]:
                rec.append({"cls": cls, "chrom": str(z["chrom"][i]), "pos": int(z["pos"][i]),
                            "ref_len": int(rl_all[i]), "alt_len": int(al_all[i]),
                            "key_rank": int(kr_all[i]), "garden": int(sites[si]),
                            "bio1": float(bio1[si]), "mac": int(mac_all[i]),
                            "maf": round(float(mac_all[i] / N), 4),
                            "nlp": round(float(nlp[i, si]), 2),
                            "Z": round(float(z["Z"][i, si]), 2)})
        print(f"{cls}: {len(hit)} significant allele records", flush=True)

    H = pd.DataFrame(rows)
    k = H.cls + H.chrom + ":" + H.pos.astype(str)
    H["key_rank"] = H.key_rank.astype(int)
    H["n_sig_alleles_at_pos"] = k.map(k.value_counts())
    H.to_csv(OUT, index=False)
    pd.DataFrame(rec).to_csv(OUT.replace(".csv", "_records.csv"), index=False)
    print(f"wrote {OUT}: {len(H)} rows ({len(rec)} allele x garden records); "
          f"multiallelic positions: {(H.n_alleles_tested_at_pos > 1).sum()}")


if __name__ == "__main__":
    main()
