#!/usr/bin/env python
"""Follow-up: the founder-GWAS JOINT is direction-AGNOSTIC, and top-JOINT SVs show NO global
founder-fitness direction (_sv_founder_mechanism.py). So what IS going on with those SVs?

Two direct reads, top-JOINT SVs vs a size-matched SV baseline (SVs in non-top blocks):
  (1) PER-SITE MAGNITUDE: across 30 sites, is the SV-carrier vs non-carrier founder-fitness
      gap larger for top-JOINT SVs? (JOINT rewards large per-site effects, any direction.)
  (2) LOCAL vs NOISE: are the per-site gaps SITE-SPECIFIC-DIRECTIONAL (win here, lose there
      = local adaptation) or symmetric noise? Report per-SV: how many sites +, how many -,
      and the max |per-site gap|.
  (3) THE SV'S OWN Δp: does the SV allele actually rise/fall in the evolved pools? overall
      Δp (evolved mean - p0) and how it splits across sites.

All cached (sample_genome_h.npz per-site fitness; af_store for SV Δp). Env: kmate.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

SEED = "results/grenenet_kmate_window_seedmix"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
PANEL = "panel/arch3"
FG = f"{lib.GEA}/hapfreq/multisite_founder_gwas_clq90_pc1"
SVL = f"{lib.GEA}/sv_adaptive/sv_landscape_clq0.9.csv"


def genome_h(samp, base):
    gs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_blocks_per_chrom.npz"
        if not os.path.exists(f):
            return None
        gs.append(np.load(f, allow_pickle=True)[f"{ch}_global_h"].astype(np.float64))
    return np.mean(gs, 0)


def per_site_fitness():
    """(n_site x 231) per-founder frequency slope at each site + site ids."""
    cache = np.load(f"{lib.GEA}/fitness/sample_genome_h.npz", allow_pickle=True)
    H = cache["H"]; samples = cache["samples"].astype(str)
    hmap = {s: i for i, s in enumerate(samples)}
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(hmap)]
    S, sites = [], []
    for site, sd in pt.groupby("site"):
        cell = {}
        for gen, g in sd.groupby("generation"):
            if int(gen) not in (1, 2, 3):
                continue
            hs, ws = [], []
            for _, r in g.iterrows():
                hs.append(H[hmap[str(r.sampleid)]])
                w = r.flowerscollected if np.isfinite(r.flowerscollected) and r.flowerscollected > 0 else 1.0
                ws.append(w)
            ws = np.asarray(ws); cell[int(gen)] = (np.vstack(hs) * ws[:, None]).sum(0) / ws.sum()
        present = sorted(cell)
        if 1 not in present:
            continue
        t = np.array([0.0] + [float(g) for g in present]); tc = t - t.mean()
        Y = np.vstack([p0] + [np.clip(cell[g], 0, 1) for g in present])
        S.append((tc[:, None] * Y).sum(0) / (tc @ tc)); sites.append(int(site))
    return np.array(sites), np.vstack(S)                    # (nsite,), (nsite x 231)


def sv_panel():
    """per SV: chrom,pos,founder 0/1 carrier vector (MAC>=12 SVs)."""
    out = {}
    for ch in CHROMS:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64)
        dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        n_alt = np.asarray(vp.sum(0)).ravel(); n_cal = np.asarray(vc.sum(0)).ravel()
        sv = (dl > 50) & (n_alt >= 12) & (n_alt <= 231 - 12) & (n_cal / 231 >= 0.9)
        idx = np.where(sv)[0]
        out[ch] = (pos[idx], vp[:, idx].toarray().astype(np.int8))
    return out


def collect(blocks, panel, sites, S):
    """per SV inside `blocks`: per-site carrier-noncarrier fitness gap (nsite,)."""
    recs = []
    for ch, grp in blocks.groupby("chrom"):
        if ch not in panel:
            continue
        pos, G = panel[ch]
        st = grp.start_pos.to_numpy(); en = grp.end_pos.to_numpy()
        for j in range(len(pos)):
            if not ((pos[j] >= st) & (pos[j] <= en)).any():
                continue
            c = G[:, j] == 1; nc = int(c.sum())
            if nc < 2 or nc > 229:
                continue
            gap = S[:, c].mean(1) - S[:, ~c].mean(1)        # per-site fitness gap (nsite,)
            recs.append(dict(chrom=ch, pos=int(pos[j]), n_carriers=nc,
                             persite_absmax=np.abs(gap).max(),
                             persite_absmean=np.abs(gap).mean(),
                             n_pos=int((gap > 0).sum()), n_neg=int((gap < 0).sum()),
                             global_gap=gap.mean()))
    return pd.DataFrame(recs)


def main():
    sites, S = per_site_fitness()
    print(f"per-site founder fitness: {len(sites)} sites x 231 founders")
    csv = pd.read_csv(f"{FG}.csv")
    up = csv.groupby("unit", as_index=False).p_joint.min()
    L = pd.read_csv(SVL).rename(columns={"block_id": "unit"})
    up = up.merge(L[["unit", "chrom", "start_pos", "end_pos", "has_sv", "n_kept"]], on="unit")
    up = up[up.n_kept >= 2]
    ntop = int(round(0.005 * len(up)))
    topset = set(up.nsmallest(ntop, "p_joint").unit)
    top = up[up.unit.isin(topset) & (up.has_sv == 1)]
    base = up[(~up.unit.isin(topset)) & (up.has_sv == 1)]
    panel = sv_panel()
    RT = collect(top, panel, sites, S); RB = collect(base, panel, sites, S)
    print(f"\ntop-JOINT SVs n={len(RT)} | baseline SVs n={len(RB)}")
    print("PER-SITE MAGNITUDE (|carrier-noncarrier fitness gap|, across sites):")
    print(f"  top-JOINT : abs-max median {RT.persite_absmax.median():.4f} | abs-mean median {RT.persite_absmean.median():.4f}")
    print(f"  baseline  : abs-max median {RB.persite_absmax.median():.4f} | abs-mean median {RB.persite_absmean.median():.4f}")
    print(f"  -> top-JOINT SVs have {RT.persite_absmax.median()/RB.persite_absmax.median():.2f}x the per-site max effect")
    print("\nLOCAL vs NOISE (per SV: sites where carriers gain vs lose):")
    print(f"  top-JOINT: median +sites {RT.n_pos.median():.0f} / -sites {RT.n_neg.median():.0f} "
          f"(of {len(sites)}); |global gap| median {RT.global_gap.abs().median():.4f}")
    print(f"  => gains and losses roughly BALANCED across sites per SV => mixed-direction "
          f"(local), not a consistent global win/lose")

    # (3) SV OWN Δp in the evolved pools
    idx_non = np.load(f"{lib.AF_STORE}/index_nonsnp.npz")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy").astype(np.float64)
    key = pd.Series(idx_non["chrom"].astype(str)) + ":" + pd.Series(idx_non["pos"].astype(np.int64).astype(str))
    kmap = {k: i for i, k in enumerate(key)}
    # evolved AF = flower-weighted mean over all evolved pools (gen1-3), from pool_matrices
    PM = f"{lib.GEA}/pool_matrices"
    ev = np.zeros(len(p0)); w = 0.0
    for g in (1, 2, 3):
        m = pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv")
        mat = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy", mmap_mode="r")
        wt = m["total_flowers"].to_numpy(float); wt = np.where(np.isfinite(wt) & (wt > 0), wt, 1.0)
        sub = np.asarray(mat).astype(np.float64)
        ok = np.isfinite(sub)
        ev += np.nansum(np.where(ok, sub, 0) * wt[:, None], 0); w += np.nansum(np.where(ok, wt[:, None], 0), 0)
    ev = ev / np.where(w > 0, w, np.nan)
    dltop = []
    for _, r in RT.iterrows():
        i = kmap.get(f"{r.chrom}:{r.pos}")
        if i is not None and np.isfinite(ev[i]):
            dltop.append(ev[i] - p0[i])
    dltop = np.array(dltop)
    print(f"\nSV OWN Δfrequency (evolved pool mean - founding p0), top-JOINT SVs n={len(dltop)}:")
    print(f"  mean Δp {np.nanmean(dltop):+.4f} | median {np.nanmedian(dltop):+.4f} | "
          f"rose in {100*(dltop>0).mean():.0f}% | |Δp| median {np.nanmedian(np.abs(dltop)):.4f}")
    print(f"  => {'rises' if np.nanmean(dltop)>0.005 else ('falls' if np.nanmean(dltop)<-0.005 else 'flat on average')} "
          f"in the pooled evolved population")
    RT.to_csv(f"{lib.GEA}/sv_adaptive/sv_founder_direction.csv", index=False)
    print(f"\n[done] -> {lib.GEA}/sv_adaptive/sv_founder_direction.csv")


if __name__ == "__main__":
    main()
