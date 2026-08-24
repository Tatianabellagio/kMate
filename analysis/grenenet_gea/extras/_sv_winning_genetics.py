#!/usr/bin/env python
"""Are SVs part of the SHARED WINNING genetics across GrENE-net sites?

Operationalizes the multiple-winners identifiability argument (user, 2026-07-02): selection
at each site elevates MANY founders (median effective winners = 15; all 31 sites >=3), so the
genetics SHARED by the winners (vs the losers), and REPLICATED across sites that have different
winners, is identifiable in GLOBAL mode -- even though a single SV cannot be separated from a
co-segregating SNP on the same winner-carrier background (that residual = the causal tier).

Per common variant v (founder MAC>=12, called>=0.9), per site s:
    beta[v,s] = mean(sel[f,s] | founder f carries v) - mean(sel[f,s] | f does not)
a winner-vs-loser contrast computed at the variant's OWN founder-carrier pattern (the correct
unit -- lesson from the block-vs-haplotype retraction), using the well-estimated per-founder
selection coefficient sel[f,s]. beta>0 => carriers of v were the winners at site s.

sel[f,s] under BOTH definitions (report only conclusions holding under both):
    DH    : h_lastgen - p0            (flower-weighted last available generation)
    SLOPE : OLS logit-slope of h over t=[0, available gens]   (uses intermediate timepoints)

TEST -- SVs vs founder-MAC + carrier-genome-divergence matched SNPs (kills the rare-variant
swing artifact AND the 'SV sits on a divergent haplotype' artifact; the confound check found
corr(mean-selection, SV-load) = -0.23, i.e. no global SV-rich-genome bias, and SV-load ~ SNP-load
per founder so divergence-matching removes the residual):
    (1) directional : is mean beta over SVs > matched-SNP null (SVs on the winning side)?
    (2) magnitude   : is mean |beta| over SVs elevated (SVs more selection-aligned)?
    significance = CROSS-SITE replication (site-level permutation; Stouffer over sites), the
    parallelism test the old founder-GWAS JOINT axis never did. Plus leave-top-winner-out.

Causal tier (--causal-only): restrict to SVs in low founder-LD (max r2 < 0.2 to any common SNP
in its cis window) -- the only SVs where no co-segregating SNP can be the driver.

Env: kmate.  Writes analysis/grenenet_gea/sv_adaptive/sv_winning_genetics{_suffix}.{csv,json}.
"""
import os, sys, json, glob
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from lib import genome_h

os.chdir("/global/scratch/users/tbellg/kmate")
CH = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
PANEL = "panel/arch3"
SEED = lib.SEEDMIX
nF = 231
MAC_MIN = 12                      # common founder floor (MAF>=5%), matches sv_landscape / LD blocks
SNP_SAMPLE = 300_000             # random common-SNP pool for the matched null (plenty for binning)
RNG = np.random.default_rng(0)
EPS = 1e-3


def per_founder_selection():
    """sel_DH[site,f] and sel_SLOPE[site,f] (nsite x 231), + site ids and p0."""
    c = np.load("analysis/grenenet_gea/fitness/sample_genome_h.npz", allow_pickle=True)
    H = c["H"]; samples = c["samples"].astype(str); founders = c["founders"].astype(str)
    hmap = {s: i for i, s in enumerate(samples)}
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(hmap)]
    dh_rows, slope_rows, g1_rows, sites = [], [], [], []
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
        if not present:
            continue
        # DH: last gen - p0
        dh_rows.append(np.clip(cell[present[-1]], 0, 1) - p0)
        # SLOPE: OLS logit-slope over [0]+present
        t = np.array([0.0] + [float(g) for g in present]); tc = t - t.mean()
        Y = np.vstack([p0] + [np.clip(cell[g], 0, 1) for g in present])
        L = np.log(np.clip(Y, EPS, 1 - EPS) / (1 - np.clip(Y, EPS, 1 - EPS)))
        slope_rows.append((tc[:, None] * L).sum(0) / (tc @ tc))
        # G1: uniform p0->gen1 window (window-length-immune control); NaN if no gen1
        g1_rows.append(np.clip(cell[1], 0, 1) - p0 if 1 in cell else np.full_like(p0, np.nan))
        sites.append(int(site))
    return (np.vstack(dh_rows), np.vstack(slope_rows), np.vstack(g1_rows),
            np.array(sites), founders, p0)


def per_founder_snp_load():
    """genome-wide common-SNP allele load per founder (for carrier-divergence matching)."""
    load = np.zeros(nF)
    for ch in CH:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        snp = (rl == 1) & (al == 1) & (na >= MAC_MIN) & (na <= nF - MAC_MIN) & (nc / nF >= 0.9)
        load += np.asarray(vp[:, snp].sum(1)).ravel()
    return load


def collect_variants(snp_load):
    """Common variants across the genome: all SV + smallindel, sampled SNPs. Returns
    G (231 x M csc, 0/1), class array, mac, pos, chrom, carrier-divergence."""
    # how many common SNPs per chrom to hit ~SNP_SAMPLE genome-wide
    snp_counts = {}
    for ch in CH:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        common = (na >= MAC_MIN) & (na <= nF - MAC_MIN) & (nc / nF >= 0.9)
        snp_counts[ch] = int((common & (rl == 1) & (al == 1)).sum())
    tot_snp = sum(snp_counts.values())
    Gs, cls, mac, pos, chrom = [], [], [], [], []
    for ch in CH:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        dl = np.abs(al - rl); p = meta["pos"].astype(np.int64)
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz").tocsc()
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        common = (na >= MAC_MIN) & (na <= nF - MAC_MIN) & (nc / nF >= 0.9)
        is_sv = common & (dl > 50); is_snp = common & (rl == 1) & (al == 1)
        is_ind = common & ~is_sv & ~is_snp
        # subsample SNPs
        snp_idx = np.where(is_snp)[0]
        k = int(round(len(snp_idx) * SNP_SAMPLE / max(tot_snp, 1)))
        keep_snp = RNG.choice(snp_idx, size=min(k, len(snp_idx)), replace=False) if len(snp_idx) else snp_idx
        for idx, lab in [(np.where(is_sv)[0], "sv"), (np.where(is_ind)[0], "indel"),
                         (np.sort(keep_snp), "snp")]:
            if len(idx) == 0:
                continue
            Gs.append(vp[:, idx]); cls.append(np.full(len(idx), lab)); mac.append(na[idx])
            pos.append(p[idx]); chrom.append(np.full(len(idx), ch))
    G = sp.hstack(Gs).tocsc().astype(np.float64)
    cls = np.concatenate(cls); mac = np.concatenate(mac).astype(int)
    pos = np.concatenate(pos); chrom = np.concatenate(chrom)
    car = np.asarray(G.sum(0)).ravel()
    div = (G.T @ snp_load) / np.where(car > 0, car, 1)      # mean carrier SNP-load
    return G, cls, mac, pos, chrom, car, div


def beta_matrix(G, car, sel):
    """beta[v,s] = mean(sel|carrier) - mean(sel|non). G:(231xM), sel:(nsite x 231)."""
    Gt_sel = G.T @ sel.T                                   # (M x nsite)
    tot = sel.sum(1)                                       # (nsite,)
    car = car[:, None]; ncar = (nF - car)
    mean_car = Gt_sel / np.where(car > 0, car, np.nan)
    mean_non = (tot[None, :] - Gt_sel) / np.where(ncar > 0, ncar, np.nan)
    return mean_car - mean_non                             # (M x nsite)


def matched_test(beta, cls, mac, div, label, n_null=1000, target="sv"):
    """Per site: mean beta and mean|beta| over TARGET-class variants vs MAC+div-matched SNP
    draws; then combine across sites (Stouffer on per-site z + sign consistency).
    target='sv' (default) or 'indel' -- indel is the SV-SPECIFICITY control (also a
    non-reference allele; if it depletes equally the effect is not SV-specific)."""
    is_sv = cls == target; is_snp = cls == "snp"
    # bins on mac (deciles among all) and div (deciles among SNPs)
    mb = np.clip(np.digitize(mac, np.quantile(mac, np.linspace(0, 1, 11))[1:-1]), 0, 9)
    dq = np.quantile(div[is_snp], np.linspace(0, 1, 11))[1:-1]
    db = np.clip(np.digitize(div, dq), 0, 9)
    key = mb * 10 + db
    snp_by_key = {k: np.where(is_snp & (key == k))[0] for k in np.unique(key[is_sv])}
    sv_idx = np.where(is_sv)[0]; sv_keys = key[sv_idx]
    usable = np.array([k in snp_by_key and len(snp_by_key[k]) > 0 for k in sv_keys])
    sv_idx = sv_idx[usable]; sv_keys = sv_keys[usable]
    n_sv = len(sv_idx)
    # per-key sv counts, for vectorized batch sampling of the matched null
    ukeys, kcounts = np.unique(sv_keys, return_counts=True)
    nS = beta.shape[1]
    zdir = np.full(nS, np.nan); zmag = np.full(nS, np.nan)
    fold_mag = np.full(nS, np.nan); obs_dir = np.full(nS, np.nan)
    for s in range(nS):
        b = beta[:, s]; ab = np.abs(b)
        obs_d = np.nanmean(b[sv_idx]); obs_m = np.nanmean(ab[sv_idx])
        nd = np.zeros(n_null); nm = np.zeros(n_null)          # sums over sv, /n_sv at end
        for k, c in zip(ukeys, kcounts):                     # draw c matched SNPs per key, all perms at once
            mem = snp_by_key[k]
            idx = RNG.integers(0, len(mem), size=(n_null, int(c)))
            nd += b[mem[idx]].sum(1); nm += ab[mem[idx]].sum(1)
        nd /= n_sv; nm /= n_sv
        zdir[s] = (obs_d - nd.mean()) / (nd.std() + 1e-12)
        zmag[s] = (obs_m - nm.mean()) / (nm.std() + 1e-12)
        fold_mag[s] = obs_m / (np.median(nm) + 1e-12); obs_dir[s] = obs_d
    stouffer_dir = np.nansum(zdir) / np.sqrt(np.isfinite(zdir).sum())
    stouffer_mag = np.nansum(zmag) / np.sqrt(np.isfinite(zmag).sum())
    p_dir = 2 * stats.norm.sf(abs(stouffer_dir)); p_mag = stats.norm.sf(stouffer_mag)
    frac_pos = np.mean(obs_dir[np.isfinite(obs_dir)] > 0)
    frac_mag_up = np.mean(zmag[np.isfinite(zmag)] > 0)
    summ = dict(label=label, n_sv=int(len(sv_idx)), n_sites=int(np.isfinite(zdir).sum()),
                stouffer_dir=float(stouffer_dir), p_dir=float(p_dir), frac_sites_beta_pos=float(frac_pos),
                stouffer_mag=float(stouffer_mag), p_mag=float(p_mag), frac_sites_mag_up=float(frac_mag_up),
                median_fold_mag=float(np.nanmedian(fold_mag)),
                mean_z_dir=float(np.nanmean(zdir)), mean_z_mag=float(np.nanmean(zmag)))
    persite = dict(zdir=zdir, zmag=zmag, obs_dir=obs_dir, fold_mag=fold_mag)
    return summ, persite


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--drop-top-winner", action="store_true",
                    help="leave-one-out: drop the single most recurrently-winning founder")
    ap.add_argument("--n-null", type=int, default=1000)
    args = ap.parse_args()
    suffix = "_droptop" if args.drop_top_winner else ""

    print("[1/4] per-founder selection (DH + SLOPE + G1) ...")
    dh, slope, g1, sites, founders, p0 = per_founder_selection()
    print(f"      {len(sites)} sites x {nF} founders")
    if args.drop_top_winner:
        top = int(np.argmax(dh.mean(0)))
        print(f"      dropping top recurrent winner: founder {founders[top]}")
        keep = np.ones(nF, bool); keep[top] = False
    else:
        keep = np.ones(nF, bool)

    print("[2/4] per-founder SNP load + common variants ...")
    snp_load = per_founder_snp_load()
    G, cls, mac, pos, chrom, car, div = collect_variants(snp_load)
    print(f"      {G.shape[1]:,} common variants: "
          f"sv={int((cls=='sv').sum()):,} indel={int((cls=='indel').sum()):,} snp={int((cls=='snp').sum()):,}")

    if args.drop_top_winner:
        G = G[keep]; dh = dh[:, keep]; slope = slope[:, keep]; g1 = g1[:, keep]
        car = np.asarray(G.sum(0)).ravel()

    print("[3/4] beta matrices + matched test, CROSS-SITE and PER-SITE (both sel defs) ...")
    clim = lib.load_climate()
    bio1 = np.array([float(clim["bio1"].reindex([s]).iloc[0]) if s in clim.index else np.nan
                     for s in sites])
    rows = []; ps = pd.DataFrame(dict(site=sites, bio1=bio1))
    for name, sel in [("DH", dh), ("SLOPE", slope), ("G1", g1)]:
        beta = beta_matrix(G, car, sel)
        for tgt in ("sv", "indel"):                        # indel = SV-specificity control
            res, persite = matched_test(beta, cls, mac, div, name, n_null=args.n_null, target=tgt)
            res["target"] = tgt; res["sel_def"] = name; rows.append(res)
            if tgt == "sv":
                ps[f"zdir_{name}"] = persite["zdir"]; ps[f"zmag_{name}"] = persite["zmag"]
                ps[f"obsdir_{name}"] = persite["obs_dir"]; ps[f"foldmag_{name}"] = persite["fold_mag"]
            print(f"  [{name}/{tgt:>5}] CROSS-SITE dir Stouffer={res['stouffer_dir']:+.2f} "
                  f"({res['frac_sites_beta_pos']:.0%} sites beta>0) | mag Stouffer={res['stouffer_mag']:+.2f}")

    # PER-SITE: which sites have SVs on the winning side, and does it track climate?
    ps["zdir_mean"] = ps[["zdir_DH", "zdir_SLOPE"]].mean(1)      # both-def average
    ps["both_pos"] = (ps.zdir_DH > 0) & (ps.zdir_SLOPE > 0)
    ps_sorted = ps.sort_values("zdir_mean", ascending=False)
    print("\nPER-SITE (SVs on the winning side; z>0 => yes; both_pos = agrees under DH & SLOPE):")
    print(ps_sorted[["site", "bio1", "zdir_DH", "zdir_SLOPE", "zmag_DH", "zmag_SLOPE", "both_pos"]]
          .round(2).to_string(index=False))
    n_sig = int(((ps.zdir_DH > 1.64) & (ps.zdir_SLOPE > 1.64)).sum())
    print(f"\nsites with SVs significantly on winning side (both z>1.64): {n_sig}/{len(ps)}")
    for name in ("DH", "SLOPE", "G1"):
        m = np.isfinite(ps[f"zdir_{name}"]) & np.isfinite(ps.bio1)
        if m.sum() > 3:
            r = stats.spearmanr(ps.bio1[m], ps[f"zdir_{name}"][m])
            print(f"  climate dependence [{name:>5}]: corr(bio1, zdir) = {r.statistic:+.2f} (p={r.pvalue:.3f}, "
                  f"n={int(m.sum())}) -> -ve = SVs purged more at HOT sites  [G1 = window-immune control]")

    print("[4/4] writing ...")
    out = f"{lib.GEA}/sv_adaptive/sv_winning_genetics{suffix}"
    pd.DataFrame(rows).to_csv(f"{out}.csv", index=False)
    ps.to_csv(f"{out}_persite.csv", index=False)
    json.dump(dict(mac_min=MAC_MIN, sites=sites.tolist(), n_null=args.n_null,
                   drop_top_winner=args.drop_top_winner, results=rows,
                   n_sites_sig_both=n_sig), open(f"{out}.json", "w"), indent=2)
    print(f"[done] {out}.csv + _persite.csv + .json")
    print("\nCONCLUSION holds only if BOTH DH and SLOPE agree. Cross-site p_dir<0.05 w/ frac>>0.5 "
          "=> SVs recurrently winning; per-site z>1.64 => SVs matter AT THAT SITE (may be climate-patterned).")


if __name__ == "__main__":
    main()
