#!/usr/bin/env python
"""Which SV events are present in cold-origin founders but absent from warm-origin ones?

Founders are split into origin-bio1 TERTILES (cold / temperate / warm, 77 each). SV records
(>50 bp net length change; insertion = alt longer, deletion = ref longer) are collapsed to Truvari
EVENTS. arch3 is the HPRC/PanGenie biallelic representation of a multi-allelic graph VCF, so
alternative alleles of one event are separate records -- up to ~78 of them, one per cactus
assembly -- and counting records counts alleles, with a count that scales with how many assemblies
are in the graph. `_build_event_map.py` collapses alleles of one bubble by Truvari-collapse
sequence/size similarity (`arch3_event_map.npz`, MODE=truvari defaults: pctsize/pctseq 0.95,
refdist 500bp) -- a stricter, more conservative grouping than the graph bubble alone, which can hold
genuinely different insertions at the same site; carriers = union over the event's records. For
every event the number of carriers in each group is counted; an event is "present" in a group if it
has >= k carriers there (k = 1, 2, 3), giving the 7 Venn regions plus "below k everywhere", by
carrier-count stratum.

NULL (called-founder): for each event, its n carriers are redrawn NDRAW times uniformly from the
founders CALLED there (multivariate hypergeometric on the per-group called counts). This respects
that carriers can only come from called founders, and the cold tertile has more cactus (long-read)
founders, which are called at more insertion records. Draws are independent across events, so
pooled null intervals ignore linkage.

Env: kmate. Writes results/sv_adaptive/sv_origin_partition_{founders.csv,variants_event.csv,null_event.npz}
(variants_event.csv is one row per EVENT). The older graph-bubble version
(sv_origin_partition_{variants,null}.csv/.npz, kept for the record) grouped alleles by pangenome
topology only, not sequence similarity.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
KEY = "/global/scratch/users/tbellg/gea_grene-net/key_files"
SV_BP, NDRAW, SEED = 50, 200, 0
GROUPS = ["cold", "temperate", "warm"]
KS = (1, 2, 3)
STRATA = [(1, 1), (2, 3), (4, 7), (8, 20), (21, 10 ** 6)]
FRAC_MIN_CARRIERS, FRAC_BINS = 5, 10          # for the "fraction of carriers that are cold" histogram


def stratum_of(n_alt):
    s = np.full(n_alt.size, -1)
    for i, (a, b) in enumerate(STRATA):
        s[(n_alt >= a) & (n_alt <= b)] = i
    return s


def pattern_counts(C, k, strat):
    """C [n x 3] group carrier counts -> counts [n_strata x 8] of presence patterns (bit0 cold, bit1 temperate, bit2 warm)."""
    pat = ((C[:, 0] >= k) * 1 + (C[:, 1] >= k) * 2 + (C[:, 2] >= k) * 4).astype(int)
    return np.bincount(strat * 8 + pat, minlength=len(STRATA) * 8).reshape(len(STRATA), 8)


def frac_hist(C, n_alt):
    m = n_alt >= FRAC_MIN_CARRIERS
    out = np.zeros((3, FRAC_BINS))
    for g in range(3):
        # integer-safe bin index (linspace edges carry float noise); fraction 1.0 goes to the last bin
        b = np.minimum((C[m, g] * FRAC_BINS) // n_alt[m], FRAC_BINS - 1).astype(int)
        out[g] = np.bincount(b, minlength=FRAC_BINS)
    return out


def called_null(n, cg, rng):
    """n [L] carriers, cg [L x 3] called founders per group -> [NDRAW x L x 3] carriers per group,
    drawn without replacement from each locus's called founders."""
    out = np.zeros((NDRAW, n.size, 3), np.int32)
    key = np.column_stack([n, cg]).astype(np.int64)
    uk, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    order = np.argsort(inv, kind="stable"); bounds = np.searchsorted(inv[order], np.arange(len(uk) + 1))
    for u, (nn, c0, c1, c2) in enumerate(uk):
        idx = order[bounds[u]:bounds[u + 1]]
        d = rng.multivariate_hypergeometric([c0, c1, c2], nn, size=(NDRAW, idx.size))
        out[:, idx, :] = d
    return out


def main():
    rng = np.random.default_rng(SEED)
    meta0 = np.load(f"{lib.PROJ}/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz", allow_pickle=True)
    founders = meta0["founders"].astype("U6")
    eco = pd.read_csv(f"{KEY}/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
    eco["ecotypeid"] = eco.ecotypeid.astype(str); em = eco.set_index("ecotypeid")
    bio1 = np.array([float(em.loc[f, "bio1"]) for f in founders])
    grp = pd.qcut(bio1, 3, labels=GROUPS).astype(str)
    half = pd.read_csv(f"{OUT}/founder_sv_content.csv", dtype={"founder": str}).set_index("founder")
    F = pd.DataFrame(dict(founder=founders, bio1=bio1, group=grp,
                          is_cactus=half.loc[founders, "is_cactus"].astype(bool).to_numpy()))
    F.to_csv(f"{OUT}/sv_origin_partition_founders.csv", index=False)
    edges = np.quantile(bio1, [1 / 3, 2 / 3])
    print(f"[1] {len(founders)} founders -> tertiles by origin bio1: cold <= {edges[0]:.2f}, warm >= {edges[1]:.2f} C; "
          f"sizes {F.group.value_counts().to_dict()}; cactus share {F.groupby('group').is_cactus.mean().round(3).to_dict()}")
    G = np.column_stack([(grp == g).astype(float) for g in GROUPS])            # [231 x 3]

    rows = []
    obs_pat = {t: {k: np.zeros((len(STRATA), 8), int) for k in KS} for t in ("ins", "del")}
    null_pat = {t: {k: np.zeros((NDRAW, len(STRATA), 8), int) for k in KS} for t in ("ins", "del")}
    obs_frac = {t: np.zeros((3, FRAC_BINS)) for t in ("ins", "del")}
    null_frac = {t: np.zeros((NDRAW, 3, FRAC_BINS)) for t in ("ins", "del")}
    nrec = {t: 0 for t in ("ins", "del")}

    EM = np.load(f"{OUT}/arch3_event_map.npz")
    print("[2] per-chromosome events, carrier counts per group, called-founder null")
    for ch in lib.CHROMS:
        cl = ch.lower(); base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        assert np.array_equal(meta["founders"].astype("U6"), founders), f"{ch}: founder order differs"
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64); pos = meta["pos"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        n_alt_all = np.asarray(vp.sum(0)).ravel()
        for t, sel in (("ins", (al - rl) > SV_BP), ("del", (rl - al) > SV_BP)):
            cols = np.where(sel & (n_alt_all >= 1))[0]; nrec[t] += cols.size
            bcode = pd.Series(EM[f"{ch}_bubble"], index=EM[f"{ch}_rec"]).reindex(cols)
            assert bcode.notna().all() and (bcode >= 0).all(), f"{ch} {t}: records without an event"
            ub, li = np.unique(bcode.to_numpy().astype(np.int64), return_inverse=True); li = li.ravel()
            M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
            bpos = np.zeros(ub.size, np.int64); np.maximum.at(bpos, li, pos[cols])
            P = ((M @ vp[:, cols].T.tocsc()) > 0).astype(np.float64)               # [loci x 231] carriers (union)
            Cd = ((M @ vc[:, cols].T.tocsc()) > 0).astype(np.float64)              # [loci x 231] called (union)
            n = np.asarray(P.sum(1)).ravel().astype(int)
            C = np.asarray(P @ G).astype(int); CG = np.asarray(Cd @ G).astype(int)
            assert np.array_equal(C.sum(1), n) and (C <= CG).all()
            size = np.abs(al - rl)[cols]
            smax = np.zeros(ub.size, np.int64); np.maximum.at(smax, li, size)
            smin = np.full(ub.size, np.iinfo(np.int64).max); np.minimum.at(smin, li, size)
            npos = pd.Series(pos[cols]).groupby(li).nunique().to_numpy()      # distinct positions inside the event
            strat = stratum_of(n)
            rows.append(pd.DataFrame(dict(chrom=ch, event=ub, pos=bpos, n_positions=npos, type=t, n_records=np.bincount(li), size_min=smin, size_max=smax,
                                          n_alt=n, n_called=CG.sum(1),
                                          n_cold=C[:, 0], n_temperate=C[:, 1], n_warm=C[:, 2],
                                          called_cold=CG[:, 0], called_temperate=CG[:, 1], called_warm=CG[:, 2])))
            for k in KS:
                obs_pat[t][k] += pattern_counts(C, k, strat)
            obs_frac[t] += frac_hist(C, n)
            Nd = called_null(n, CG, rng)
            for b in range(NDRAW):
                for k in KS:
                    null_pat[t][k][b] += pattern_counts(Nd[b], k, strat)
                null_frac[t][b] += frac_hist(Nd[b], n)
            print(f"  {ch} {t}: {cols.size:,} records -> {ub.size:,} events", flush=True)

    D = pd.concat(rows, ignore_index=True)
    D.to_csv(f"{OUT}/sv_origin_partition_variants_event.csv", index=False)
    np.savez_compressed(f"{OUT}/sv_origin_partition_null_event.npz",
                        groups=np.array(GROUPS), ks=np.array(KS), strata=np.array(STRATA), ndraw=NDRAW,
                        frac_min_carriers=FRAC_MIN_CARRIERS, tertile_edges=edges,
                        nrec_ins=nrec["ins"], nrec_del=nrec["del"],
                        **{f"obs_pat_{t}_k{k}": obs_pat[t][k] for t in obs_pat for k in KS},
                        **{f"null_pat_{t}_k{k}": null_pat[t][k] for t in null_pat for k in KS},
                        **{f"obs_frac_{t}": obs_frac[t] for t in obs_frac},
                        **{f"null_frac_{t}": null_frac[t] for t in null_frac})

    names = ["below k everywhere", "cold only", "temperate only", "cold+temperate", "warm only",
             "cold+warm", "temperate+warm", "all three"]
    for t in ("ins", "del"):
        print(f"\n[3] {t}: {nrec[t]:,} records -> {int((D.type == t).sum()):,} events -- observed / called-null expected, strata >= 2 carriers")
        for k in KS:
            o = obs_pat[t][k][1:].sum(0); e = null_pat[t][k][:, 1:].sum(1).mean(0)
            print(f"  k={k}: " + " | ".join(f"{names[i]} {o[i]:,}/{e[i]:,.0f}={o[i]/e[i]:.2f}x" for i in (1, 4, 2, 7)))
    print(f"\n[wrote] {OUT}/sv_origin_partition_{{founders.csv,variants_event.csv,null_event.npz}}")


if __name__ == "__main__":
    main()
