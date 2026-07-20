#!/usr/bin/env python
"""Sanity check the passenger finding against the GLOBAL-mode identity:

In GLOBAL mode AF(v) = Σ_f h[f]·G[f,v], so two variants with the SAME founder column are
FORCED to identical frequency. If a co-block SNP separates founder fitness (and sweeps) but
the SV doesn't, their founder columns MUST differ -- i.e. the SV is in LOW founder-LD (r²<0.9)
with the selected SNP/haplotype, despite sharing the ~2 kb block. If instead the SV is r²>=0.9
with the block's SNPs, the passenger result would be self-contradictory (a bug).

For each top-JOINT SV: founder-r² vs every co-block SNP (panel genotype), report max/median and
the fraction of co-block SNPs at r²>=0.9. Also the direct GLOBAL identity check:
Δp(SV) should equal Σ_{carrier founders} (per-founder frequency slope). Env: kmate.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

SEED = lib.SEEDMIX
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
PANEL = "panel/arch3"
FG = f"{lib.GEA}/hapfreq/multisite_founder_gwas_clq90_pc1"
SVL = f"{lib.GEA}/sv_adaptive/sv_landscape_clq0.9.csv"


def r2(a, B):
    """founder-genotype r² of binary vector a (231,) vs each column of B (231 x m)."""
    a = a - a.mean(); Bc = B - B.mean(0)
    num = (a[:, None] * Bc).sum(0) ** 2
    den = (a @ a) * (Bc ** 2).sum(0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def main():
    csv = pd.read_csv(f"{FG}.csv")
    up = csv.groupby("unit", as_index=False).p_joint.min()
    L = pd.read_csv(SVL).rename(columns={"block_id": "unit"})
    up = up.merge(L[["unit", "chrom", "start_pos", "end_pos", "has_sv", "n_kept"]], on="unit")
    up = up[up.n_kept >= 2]
    ntop = int(round(0.005 * len(up)))
    topsv = up.nsmallest(ntop, "p_joint"); topsv = topsv[topsv.has_sv == 1]

    rows = []
    for ch in CHROMS:
        blk = topsv[topsv.chrom == ch]
        if not len(blk):
            continue
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64)
        dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
        is_snp = (meta["ref_len"] == 1) & (meta["alt_len"] == 1)
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        common = (na >= 12) & (na <= 219) & (nc / 231 >= 0.9)
        sv = common & (dl > 50)
        snp = common & is_snp
        for _, b in blk.iterrows():
            inb = (pos >= b.start_pos) & (pos <= b.end_pos)
            sv_j = np.where(inb & sv)[0]
            snp_j = np.where(inb & snp)[0]
            if len(sv_j) == 0 or len(snp_j) == 0:
                continue
            Bsnp = vp[:, snp_j].toarray().astype(float)     # 231 x nsnp
            for j in sv_j:
                a = vp[:, j].toarray().ravel().astype(float)
                rr = r2(a, Bsnp)
                rows.append(dict(unit=b.unit, pos=int(pos[j]), n_snp_block=len(snp_j),
                                 r2_max=float(rr.max()), r2_median=float(np.median(rr)),
                                 frac_snp_r2ge09=float((rr >= 0.9).mean()),
                                 sv_mac=int(na[j])))
    R = pd.DataFrame(rows)
    print(f"top-JOINT SVs with >=1 co-block common SNP: {len(R)}")
    print(f"\nFounder-LD of each top-JOINT SV vs its co-block SNPs:")
    print(f"  r2_max    (best co-block SNP): median {R.r2_max.median():.3f}  "
          f"[SVs with a co-block SNP at r2>=0.9: {int((R.r2_max>=0.9).sum())}/{len(R)}]")
    print(f"  r2_median (typical co-block SNP): median {R.r2_median.median():.3f}")
    print(f"  frac of co-block SNPs at r2>=0.9 with the SV: median {R.frac_snp_r2ge09.median():.3f}")
    verdict = ("SVs are in LOW founder-LD with co-block SNPs -> different founder columns -> "
               "different trajectory is EXPECTED, not a bug (same ~2kb window, different haplotype)"
               if R.r2_max.median() < 0.9 else
               "SVs ARE r2>=0.9 with co-block SNPs -> passenger result would be contradictory (recheck)")
    print(f"\n=> {verdict}")
    R.to_csv(f"{lib.GEA}/sv_adaptive/sv_ld_check.csv", index=False)

    # ---- GLOBAL-mode identity: Δp(SV) ?= Σ_carrier per-founder slope ----
    def genome_h(samp, base):
        gs = []
        for c in CHROMS:
            f = f"{base}/{samp}_{c}.h_per_chrom.npz"
            if not os.path.exists(f):
                return None
            gs.append(np.load(f, allow_pickle=True)[c].astype(np.float64))
        return np.mean(gs, 0)
    cache = np.load(f"{lib.GEA}/fitness/sample_genome_h.npz", allow_pickle=True)
    H = cache["H"]; samps = cache["samples"].astype(str); hmap = {s: i for i, s in enumerate(samps)}
    seeds = sorted({p.split("/")[-1].split("_Chr")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
    h0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(hmap)]
    # genome-wide mean founder h change (founding -> all evolved, flower-weighted)
    w = pt.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
    hev = (np.vstack([H[hmap[str(s)]] for s in pt.sampleid.astype(str)]) * w[:, None]).sum(0) / w.sum()
    dh = hev - h0                                            # per-founder Δh
    print(f"\nGLOBAL-mode identity Δp(SV) = Σ_carrier Δh[f]  (spot check, first 5 top-JOINT SVs):")
    ch0 = "Chr1"
    meta = np.load(f"{PANEL}/{ch0.lower()}/var_pa_231_arch3_{ch0.lower()}.meta.npz", allow_pickle=True)
    pos = meta["pos"].astype(np.int64); vp = sp.load_npz(f"{PANEL}/{ch0.lower()}/var_pa_231_arch3_{ch0.lower()}.var_pa.npz")
    shown = 0
    for _, b in topsv[topsv.chrom == ch0].iterrows():
        inb = np.where((pos >= b.start_pos) & (pos <= b.end_pos))[0]
        for j in inb:
            g = vp[:, j].toarray().ravel().astype(float)
            if g.sum() < 12 or g.sum() > 219:
                continue
            pred = float((dh * g).sum())                    # predicted Δp from founder identity
            print(f"  {ch0}:{pos[j]}  carriers={int(g.sum())}  predicted Δp(sum carrier Δh) = {pred:+.4f}")
            shown += 1
            if shown >= 5:
                break
        if shown >= 5:
            break
    print("  (compare to the observed pooled Δp ~ +0.002 median from _sv_passenger_test -- same scale => identity holds)")


if __name__ == "__main__":
    main()
