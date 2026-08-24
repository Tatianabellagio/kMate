#!/usr/bin/env python
"""DECISIVE test: in the top-JOINT SV blocks, is it the SV or the HAPLOTYPE that carries the
founder-fitness signal?

_sv_founder_direction.py found the SVs' own founder-carrier partition shows NO fitness gap
(0.91x baseline, balanced +/- across sites). But that's only meaningful if SOMETHING in those
blocks DOES separate founders by fitness — otherwise the reconstruction is just underpowered.

Positive control: for each top-JOINT SV block, compute the founder-fitness gap of its best
HAP-CLUSTER (the founder-GWAS's actual unit) and compare to the block's SV gap. If hap-cluster
gap >> SV gap -> the SV is a PASSENGER; the block scores on its haplotype structure, the SV
just co-occurs. Also: the SV's own Δp in evolved pools (indexed to just the SV columns).
Env: kmate.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from founder_genotype import build_genotype

SEED = lib.SEEDMIX
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
PANEL = "panel/arch3"
FG = f"{lib.GEA}/hapfreq/multisite_founder_gwas_clq90_pc1"
SVL = f"{lib.GEA}/sv_adaptive/sv_landscape_clq0.9.csv"


def genome_h(samp, base):
    gs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_per_chrom.npz"
        if not os.path.exists(f):
            return None
        gs.append(np.load(f, allow_pickle=True)[ch].astype(np.float64))
    return np.mean(gs, 0)


def per_founder_fitness():
    cache = np.load(f"{lib.GEA}/fitness/sample_genome_h.npz", allow_pickle=True)
    H = cache["H"]; samples = cache["samples"].astype(str); founders = cache["founders"].astype(str)
    hmap = {s: i for i, s in enumerate(samples)}
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(hmap)]
    Sg = []
    for site, sd in pt.groupby("site"):
        cell = {}
        for gen, g in sd.groupby("generation"):
            if int(gen) not in (1, 2, 3):
                continue
            hs, ws = [], []
            for _, r in g.iterrows():
                hs.append(H[hmap[str(r.sampleid)]])
                wv = r.flowerscollected if np.isfinite(r.flowerscollected) and r.flowerscollected > 0 else 1.0
                ws.append(wv)
            ws = np.asarray(ws); cell[int(gen)] = (np.vstack(hs) * ws[:, None]).sum(0) / ws.sum()
        present = sorted(cell)
        if 1 not in present:
            continue
        t = np.array([0.0] + [float(g) for g in present]); tc = t - t.mean()
        Y = np.vstack([p0] + [np.clip(cell[g], 0, 1) for g in present])
        Sg.append((tc[:, None] * Y).sum(0) / (tc @ tc))
    S = np.vstack(Sg)                                       # (nsite x 231)
    return founders, S, S.mean(0)


def main():
    founders, S, fit = per_founder_fitness()
    # top-JOINT SV blocks
    csv = pd.read_csv(f"{FG}.csv")
    up = csv.groupby("unit", as_index=False).p_joint.min()
    L = pd.read_csv(SVL).rename(columns={"block_id": "unit"})
    up = up.merge(L[["unit", "chrom", "start_pos", "end_pos", "has_sv", "n_kept"]], on="unit")
    up = up[up.n_kept >= 2]
    ntop = int(round(0.005 * len(up)))
    topsv = up.nsmallest(ntop, "p_joint")
    topsv = topsv[topsv.has_sv == 1]
    units = set(topsv.unit)

    # hap-cluster genotype for these blocks
    G, gf, reg = build_genotype()
    assert list(gf) == list(founders), "founder order mismatch"
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    col_by_unit = {u: g.index.to_numpy() for u, g in reg.groupby("unit") if u in units}

    # per unit: best hap-cluster |global fitness gap| (carrier vs non-carrier)
    hap_gap, hap_absmax = [], []
    for u, cols in col_by_unit.items():
        gaps = []
        pmax = []
        for c in cols:
            car = G[:, c] == 1; nc = int(car.sum())
            if nc < 2 or nc > 229:
                continue
            gaps.append(fit[car].mean() - fit[~car].mean())
            pmax.append(np.abs(S[:, car].mean(1) - S[:, ~car].mean(1)).max())
        if gaps:
            hap_gap.append(max(gaps, key=abs)); hap_absmax.append(max(pmax))
    hap_gap = np.array(hap_gap); hap_absmax = np.array(hap_absmax)

    # SV carriers in these blocks
    sv_gap, sv_absmax, sv_cols_by_ch = [], [], {}
    for ch in CHROMS:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64)
        dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); ncl = np.asarray(vc.sum(0)).ravel()
        sv = (dl > 50) & (na >= 12) & (na <= 219) & (ncl / 231 >= 0.9)
        blk = topsv[topsv.chrom == ch]
        keep_pos = []
        for j in np.where(sv)[0]:
            if ((pos[j] >= blk.start_pos.to_numpy()) & (pos[j] <= blk.end_pos.to_numpy())).any():
                car = vp[:, j].toarray().ravel() == 1; nc = int(car.sum())
                if nc < 2 or nc > 229:
                    continue
                sv_gap.append(fit[car].mean() - fit[~car].mean())
                sv_absmax.append(np.abs(S[:, car].mean(1) - S[:, ~car].mean(1)).max())
                keep_pos.append(int(pos[j]))
        sv_cols_by_ch[ch] = keep_pos
    sv_gap = np.array(sv_gap); sv_absmax = np.array(sv_absmax)

    print(f"top-JOINT SV blocks: {len(units)} | hap-clusters scored {len(hap_gap)} | SVs scored {len(sv_gap)}")
    print("\nPOSITIVE CONTROL — founder-fitness gap (carrier vs non-carrier), same blocks:")
    print(f"  block best HAP-CLUSTER : |global gap| median {np.median(np.abs(hap_gap)):.4f} | "
          f"per-site absmax median {np.median(hap_absmax):.4f}")
    print(f"  the SVs in those blocks: |global gap| median {np.median(np.abs(sv_gap)):.4f} | "
          f"per-site absmax median {np.median(sv_absmax):.4f}")
    ratio = np.median(np.abs(hap_gap)) / max(np.median(np.abs(sv_gap)), 1e-9)
    print(f"  -> hap-cluster fitness gap is {ratio:.1f}x the SV gap "
          f"=> {'SV is a PASSENGER (block scores on its haplotype, not the SV)' if ratio > 1.8 else 'SV tracks the haplotype fitness'}")

    # SV OWN Δp -- index ONLY the SV columns (no full-matrix load)
    idx_non = np.load(f"{lib.AF_STORE}/index_nonsnp.npz")
    nchrom = idx_non["chrom"].astype(str); npos = idx_non["pos"].astype(np.int64)
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy").astype(np.float64)
    kmap = {f"{c}:{p}": i for c, p, i in zip(nchrom, npos, range(len(npos)))}
    want = [kmap[f"{ch}:{p}"] for ch in CHROMS for p in sv_cols_by_ch.get(ch, []) if f"{ch}:{p}" in kmap]
    want = np.array(sorted(set(want)))
    PM = f"{lib.GEA}/pool_matrices"
    num = np.zeros(len(want)); den = np.zeros(len(want))
    for g in (1, 2, 3):
        m = pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv")
        mat = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy", mmap_mode="r")
        sub = np.asarray(mat[:, want]).astype(np.float64)   # only the SV columns
        wt = m["total_flowers"].to_numpy(float); wt = np.where(np.isfinite(wt) & (wt > 0), wt, 1.0)
        ok = np.isfinite(sub)
        num += np.nansum(np.where(ok, sub, 0) * wt[:, None], 0); den += np.nansum(np.where(ok, wt[:, None], 0), 0)
    ev = num / np.where(den > 0, den, np.nan)
    dp = ev - p0[want]
    print(f"\nSV OWN Δfrequency (evolved pooled mean - founding p0), n={np.isfinite(dp).sum()}:")
    print(f"  mean {np.nanmean(dp):+.4f} | median {np.nanmedian(dp):+.4f} | rose in {100*np.nanmean(dp>0):.0f}% | "
          f"|Δp| median {np.nanmedian(np.abs(dp)):.4f}  (founding p0 median {np.median(p0[want]):.3f})")


if __name__ == "__main__":
    main()
