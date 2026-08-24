#!/usr/bin/env python
"""Mechanism check for the founder-GWAS SV signal (×2.07 top-0.5% JOINT sv_frac):

Q1. Do SV-bearing top-JOINT blocks actually SEPARATE high- vs low-fitness founders?
Q2. Which DIRECTION — do the SV-carrying founder haplotypes rise or fall in frequency?

DATA (all cached / no new EM):
  - per-founder fitness = genome-wide founder-frequency SLOPE, per site, from the cached
    per-sample h (fitness/sample_genome_h.npz) + seedmix p0. Averaged across sites = GLOBAL
    founder fitness (how much each of the 231 ecotypes won/lost overall). This is exactly the
    founder-GWAS trait (build_trait_site), just reconstructed here per founder.
  - which founders carry each SV = panel var_pa at the SV position (231-founder 0/1).
  - top-JOINT SV-bearing blocks = sv_landscape (MAC>=12) ∩ top-0.5% by founder-GWAS p_joint.

For every SV inside a top-JOINT block: split the 231 founders into SV-carriers vs non-carriers,
compare their GLOBAL fitness slope (Q1 = separation, Q2 = sign). Also report the SV allele's own
Δfrequency in the evolved pools (does the SV go up or down?). Env: kmate (numpy2/scipy.sparse).
"""
import os, sys, glob
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

WIN = lib.OUT
SEED = lib.SEEDMIX
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
Tg = np.array([0.0, 1.0, 2.0, 3.0])
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


def per_founder_global_fitness():
    """Per-founder frequency SLOPE (p0->gen1->2->3), flower-weighted per site-gen, then
    averaged over sites = each ecotype's overall win/loss rate. Uses cached sample h."""
    cache = np.load(f"{lib.GEA}/fitness/sample_genome_h.npz", allow_pickle=True)
    H = cache["H"]; samples = cache["samples"].astype(str); founders = cache["founders"].astype(str)
    hmap = {s: i for i, s in enumerate(samples)}
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).isin(hmap)]
    slopes = []
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
        slopes.append((tc[:, None] * Y).sum(0) / (tc @ tc))
    S = np.vstack(slopes)                                  # (n_site, 231)
    return founders, S.mean(0), p0                         # global fitness slope per founder


def panel_carriers():
    """Founder 0/1 alt-presence + positions for every panel SV (|dlen|>50, MAC>=12)."""
    rows = []
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
        rows.append((ch, pos[idx], vp[:, idx].toarray().astype(np.int8)))   # (231 x nsv)
    return rows


def main():
    founders, fit, p0 = per_founder_global_fitness()
    print(f"per-founder GLOBAL fitness slope: mean {fit.mean():+.4f}, sd {fit.std():.4f} "
          f"(winners rise, losers fall over gens)")

    # top-JOINT blocks (founder-GWAS), aggregate hap-clusters -> unit min-p
    csv = pd.read_csv(f"{FG}.csv")
    up = csv.groupby("unit", as_index=False).p_joint.min()
    L = pd.read_csv(SVL); L = L.rename(columns={"block_id": "unit"})
    up = up.merge(L[["unit", "chrom", "start_pos", "end_pos", "n_sv", "has_sv", "n_kept"]], on="unit")
    up = up[up.n_kept >= 2]
    ntop = int(round(0.005 * len(up)))
    top = up.nsmallest(ntop, "p_joint")
    topsv = top[top.has_sv == 1]
    print(f"top-0.5% JOINT blocks: {len(top)} | carrying an SV: {len(topsv)} "
          f"({100*len(topsv)/len(top):.0f}%)")

    # founder carriers for every panel SV, restricted to SVs inside a top-JOINT SV block
    carr = panel_carriers()
    fmean = fit.mean()
    recs = []
    for ch, pos, G in carr:
        blk = topsv[topsv.chrom == ch]
        if not len(blk):
            continue
        for _, b in blk.iterrows():
            inb = (pos >= b.start_pos) & (pos <= b.end_pos)
            for j in np.where(inb)[0]:
                carriers = G[:, j] == 1
                nc = int(carriers.sum())
                if nc < 2 or nc > 229:
                    continue
                fit_c = fit[carriers].mean(); fit_n = fit[~carriers].mean()
                recs.append(dict(unit=b.unit, pos=int(pos[j]), n_carriers=nc,
                                 fit_carrier=fit_c, fit_noncarrier=fit_n,
                                 delta=fit_c - fit_n, carrier_p0=float(p0[carriers].mean() if False else np.nan)))
    R = pd.DataFrame(recs)
    print(f"\nSVs inside top-JOINT SV blocks: {len(R)}  (each split into founder carriers vs not)")
    print(f"  Q1 SEPARATION  |Δfitness| median = {R.delta.abs().median():+.4f}  "
          f"vs genome-wide founder fitness sd {fit.std():.4f}")
    up_frac = (R.delta > 0).mean()
    print(f"  Q2 DIRECTION   SV-carrier founders FITTER (delta>0) in {100*up_frac:.0f}% of SVs "
          f"(median Δ = {R.delta.median():+.4f})")
    # signed test: are SV-carriers systematically fitter or less fit than non-carriers?
    w = stats.wilcoxon(R.delta) if len(R) > 10 else None
    print(f"  Wilcoxon signed-rank on Δ (carrier-noncarrier fitness): "
          f"p={w.pvalue:.2e}" if w else "  (too few)")

    # BASELINE: same for SVs in NON-selected blocks (should be ~0 / balanced)
    nonsv_blocks = up[(up.has_sv == 1) & (~up.unit.isin(set(top.unit)))]
    base = []
    nb = set(nonsv_blocks.unit)
    for ch, pos, G in carr:
        blk = nonsv_blocks[nonsv_blocks.chrom == ch]
        if not len(blk):
            continue
        st = blk.start_pos.to_numpy(); en = blk.end_pos.to_numpy()
        for j in range(len(pos)):
            if not ((pos[j] >= st) & (pos[j] <= en)).any():
                continue
            carriers = G[:, j] == 1; ncj = int(carriers.sum())
            if ncj < 2 or ncj > 229:
                continue
            base.append(fit[carriers].mean() - fit[~carriers].mean())
    base = np.array(base)
    print(f"\nBASELINE (SVs in NON-top-JOINT blocks, n={len(base)}): "
          f"carrier fitter in {100*(base>0).mean():.0f}%, median Δ {np.median(base):+.4f}")

    R.to_csv(f"{lib.GEA}/sv_adaptive/sv_founder_mechanism.csv", index=False)
    print(f"\n[done] -> {lib.GEA}/sv_adaptive/sv_founder_mechanism.csv")


if __name__ == "__main__":
    main()
