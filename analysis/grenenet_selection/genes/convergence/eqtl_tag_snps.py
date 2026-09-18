#!/usr/bin/env python
"""Which SNP does the cis-eQTL readout actually stand on, for each candidate?

The eQTL is run on SNPs, so a candidate's expression evidence is read at whichever SNP
carries its haplotype. Two ways to find that SNP, in order of preference:

  founder_r2   `best_r2_snp_pos` from build_snp_tagging.py -- real LD, measured in the 231
               founder genotypes. Preferred when available.
  pool_af      the SNP within +-WIN whose ALT-frequency vector across the 352 gen9 pools
               correlates best with the variant's. Used where founder r2 is UNTESTABLE --
               GPX6 and FUS3, whose SVs are called in only ~34% of founders, so no r2 can
               be computed at all.

The pool-AF route is not LD and must not be called that: two variants can track each other
because both respond to climate. But it is measured on dense data (the AF matrices have no
missingness) rather than on the 34% of founders that happen to be callable, and where it
returns r ~ 1 with a matching MAF the two are the same haplotype for practical purposes.
GPX6 is the clear case: Chr4:7,013,143 has pool MAF 0.117746 against the deletion's
0.117746 and an identical garden-level climate gradient.

env: kmate. Compute node, ~2 min. Writes results/eqtl_tag_snps.csv.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
import lib                                                       # noqa: E402

OUT = f"{HERE}/results"
WIN = 25_000

CAND = {  # sym -> (gene, chrom, pos, ref_len, alt_len)
    "GPX6": ("AT4G11600", "Chr4", 7011705, 1225, 61),
    "FUS3": ("AT3G26790", "Chr3", 9856500, 3, 1),
    "AT2G30000": ("AT2G30000", "Chr2", 12805666, 1, 2),
    "AT4G13200": ("AT4G13200", "Chr4", 7669175, 1, 3),
    "CRK14": ("AT4G23220", "Chr4", 12157244, 64, 1),
    "SSL7": ("AT3G51450", "Chr3", 19091417, 2, 13),
    "SCPL34": ("AT5G23210", "Chr5", 7812576, 3, 1),
    "GLR1.3": ("AT5G48410", "Chr5", 19619649, 2, 1),
    "LPP1": ("AT2G01180", "Chr2", 108729, 1, 2),
}


def main():
    F = pd.read_csv(f"{OUT}/functional_variants.csv")
    GP = pd.read_csv(f"{OUT}/gea_pool.csv")        # MAF, to resolve colliding keys
    Mn = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy",
                 mmap_mode="r")
    Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy",
                 mmap_mode="r")
    idx_n = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    idx_s = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
    p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
    kn = pd.DataFrame({"chrom": idx_n["chrom"].astype(str), "pos": idx_n["pos"],
                       "ref_len": idx_n["ref_len"], "alt_len": idx_n["alt_len"]})
    kn["row"] = np.arange(len(kn))
    cs, ps_ = idx_s["chrom"].astype(str), idx_s["pos"]

    rows = []
    for sym, (gene, ch, pos, rl, al) in CAND.items():
        m = F[(F.chrom == ch) & (F.pos == pos) & (F.ref_len == rl) & (F.alt_len == al)]
        tag_pos, r2, method, note = None, np.nan, None, ""
        if len(m):
            r = m.iloc[0]
            if not bool(r.tag_untestable) and np.isfinite(r.best_r2_snp_pos) \
                    and r.best_r2_snp_pos > 0:
                tag_pos, r2, method = int(r.best_r2_snp_pos), float(r.best_r2_snp), "founder_r2"

        if tag_pos is None:                       # untestable -> pool-AF proxy
            hit = kn[(kn.chrom == ch) & (kn.pos == pos) & (kn.ref_len == rl)
                     & (kn.alt_len == al)]
            if not len(hit):
                rows.append(dict(symbol=sym, gene=gene, tag_pos=np.nan, method="none"))
                continue
            # NEVER .iloc[0] here. GPX6's key matches TWO af_store rows -- the GEA record
            # (gen9 MAF 0.117746) and a collider (0.034201) -- and taking the first one
            # silently correlated the WRONG allele, which reported a best proxy of r=0.69
            # at 14 kb when the true answer is r=1.000 at 1.4 kb. Resolve by gen9 MAF
            # against the GEA pool, the same rule as audit_gea_colliders.py.
            row = int(hit.row.iloc[0])
            if len(hit) > 1:
                want = GP[(GP.chrom == ch) & (GP.pos == pos) & (GP.ref_len == rl)
                          & (GP.alt_len == al)].MAF
                if not len(want):
                    raise SystemExit(f"{sym}: ambiguous key and no GEA MAF to resolve it")
                target = float(want.iloc[0])
                cands = []
                for rr_ in hit.row:
                    mu = float(np.nanmean(np.asarray(Mn[:, int(rr_)], dtype=np.float64)))
                    cands.append((abs(min(mu, 1 - mu) - target), int(rr_)))
                d, row = min(cands)
                if d > 1e-6:
                    raise SystemExit(f"{sym}: no af_store row matches the GEA MAF {target}")
            v = np.asarray(Mn[:, row], dtype=np.float64)
            sel = np.where((cs == ch) & (ps_ >= pos - WIN) & (ps_ <= pos + WIN))[0]
            sel = sel[np.minimum(p0s[sel], 1 - p0s[sel]) > 0.05]
            S = np.asarray(Ms[:, sel], dtype=np.float64)
            vz = v - np.nanmean(v)
            Sz = S - np.nanmean(S, 0)
            rr = (np.nansum(Sz * vz[:, None], 0)
                  / np.sqrt(np.nansum(Sz ** 2, 0) * np.nansum(vz ** 2)))
            j = int(np.nanargmax(np.abs(rr)))
            tag_pos, r2, method = int(ps_[sel][j]), float(rr[j]), "pool_af"
            note = (f"founder r2 untestable; best pool-AF match over 352 pools "
                    f"(r={rr[j]:+.4f})")
        rows.append(dict(symbol=sym, gene=gene, chrom=ch, var_pos=pos,
                         tag_pos=tag_pos, tag_stat=round(r2, 4), method=method,
                         dist_bp=tag_pos - pos, note=note))

    T = pd.DataFrame(rows)
    # what the eQTL says at that SNP
    if os.path.exists(f"{OUT}/eqtl/cis_eqtl_all.csv"):
        Q = pd.read_csv(f"{OUT}/eqtl/cis_eqtl_all.csv")
        nlp, pct, exact = [], [], []
        for t in T.itertuples():
            q = Q[Q.symbol == t.symbol]
            if not len(q) or not np.isfinite(t.tag_pos):
                nlp.append(np.nan); pct.append(np.nan); exact.append(False); continue
            h = q[q.ps == t.tag_pos]
            if len(h):
                v = h.nlp.iloc[0]; exact.append(True)
            else:
                d = (q.ps - t.tag_pos).abs(); v = q.loc[d.idxmin(), "nlp"]; exact.append(False)
            nlp.append(round(v, 3)); pct.append(round(100 * (q.nlp < v).mean(), 1))
        T["eqtl_nlp"], T["eqtl_pct"], T["tag_tested"] = nlp, pct, exact
    T.to_csv(f"{OUT}/eqtl_tag_snps.csv", index=False)
    print(T.to_string(index=False))
    print(f"\nwrote {OUT}/eqtl_tag_snps.csv")


if __name__ == "__main__":
    main()
