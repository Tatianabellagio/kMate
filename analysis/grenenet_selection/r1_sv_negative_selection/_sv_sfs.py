#!/usr/bin/env python
"""Site-frequency spectra of SV insertions / deletions, counted as ALLELES and as EVENTS.

WHY TWO UNITS. arch3 is the biallelic representation of a multi-allelic graph VCF, so one event
appears as one record per distinct inserted sequence (up to ~78, one per cactus assembly).
Counting records gives the ALLELE spectrum; merging a Truvari event's records (union of carriers,
`_build_event_map.py`, `arch3_event_map.npz`, truvari-collapse defaults) gives the EVENT spectrum,
i.e. presence/absence of the event. The two can only differ where an event has several alleles, so
the contrast is itself the diagnostic.

PROJECTION. Call rate varies per record, so raw carrier counts are not comparable. Every variant is
down-projected hypergeometrically to a fixed sample of NPROJ called founders (records called in
fewer are dropped) -- the standard fix, and the same one that makes the two founder sets below
comparable.

FOUNDER SETS (three, because neither half of the panel is neutral).
  `all`      = 231 founders.
  `pangenie` = the 151 genotyped founders only: none contributed assemblies to the graph, so every
               one is scored against the same allele panel and the assembly-richness imbalance
               between climate groups cannot act. Tertiles 45 / 54 / 52.
  `cactus`   = the 80 long-read assembly founders only: every one is assembly-resolved, so private
               insertions are discoverable symmetrically and allele-splitting is at its maximum --
               the upper bound on the allele/bubble gap. Tertiles 32 / 23 / 25, so WITHIN-GROUP
               frequency resolution is coarse (1/23 = 4.3%): a "<= 5%" group bin is singletons only.
               Use this set for the pooled class contrast, not for fine within-group spectra.
Group spectra (cold / temperate / warm origin tertiles) are computed within each set.

REFERENCE CLASSES. SNPs and short indels (1-50 bp) from the same panel and the same founders,
subsampled, as the neutral-ish comparison: demography moves every class, selection does not.

Env: kmate. Writes results/sv_adaptive/sv_sfs_event.npz.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.special import gammaln

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP, NPROJ, NPROJ_G, SEED, SNP_KEEP, CALL_MIN = 50, 60, 30, 0, 0.05, 0.9   # NPROJ_G: per-group projection (each tertile is ~50 founders in the PanGenie-only set)      # NPROJ called founders; SNP_KEEP = fraction of SNP/indel records kept
CLASSES = ("ins_allele", "ins_event", "del_allele", "del_event", "snp", "indel")
GROUPS = ("all", "cold", "temperate", "warm")


def project(c, m, k=NPROJ):
    """Hypergeometric down-projection: rows = variants (c carriers of m called), cols = 0..k."""
    j = np.arange(k + 1)
    lc = gammaln(c[:, None] + 1) - gammaln(j[None, :] + 1) - gammaln(np.maximum(c[:, None] - j[None, :], 0) + 1)
    lo = gammaln((m - c)[:, None] + 1) - gammaln(k - j[None, :] + 1) - gammaln(np.maximum((m - c)[:, None] - (k - j[None, :]), 0) + 1)
    ln = gammaln(m[:, None] + 1) - gammaln(k + 1) - gammaln(m[:, None] - k + 1)
    ok = (j[None, :] <= c[:, None]) & ((k - j[None, :]) <= (m - c)[:, None])
    return np.where(ok, np.exp(lc + lo - ln), 0.0)


def main():
    rng = np.random.default_rng(SEED)
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    BM = np.load(f"{OUT}/arch3_event_map.npz")
    is_cac = F.is_cactus.to_numpy().astype(bool)
    sets = {"all": np.ones(len(F), bool), "pangenie": ~is_cac, "cactus": is_cac}
    grp = F.group.to_numpy()

    def n_in(s, g):
        """founders available in this (set, group)."""
        m = sets[s]
        return int(m.sum()) if g == "all" else int((m & (grp == g)).sum())

    def kproj(s, g):
        """Projection size. NPROJ / NPROJ_G as before, but never more than 75% of the founders that
        (set, group) actually has -- the cactus tertiles are 23-32 founders, so a fixed 30 would
        empty them. Unchanged for the `all` and `pangenie` sets."""
        return min(NPROJ if g == "all" else NPROJ_G, int(0.75 * n_in(s, g)))

    for s in sets:
        print(f"[set] {s}: {n_in(s,'all')} founders "
              + ", ".join(f"{g} {n_in(s,g)} (k={kproj(s,g)})" for g in GROUPS if g != "all"))
    acc = {f"{s}_{g}_{c}": np.zeros(kproj(s, g) + 1) for s in sets for g in GROUPS for c in CLASSES}
    freq = {f"{s}_{g}_{c}": [] for s in sets for g in GROUPS for c in CLASSES}     # raw ALT frequency among called founders
    freqw = {f"{s}_{g}_{c}": [] for s in sets for g in GROUPS for c in CLASSES}    # same, WELL-CALLED only (>= 90% of the set called)
    # freqw2: as freqw but PRIVATE units dropped -- panel-wide carrier count >= 2. An SV carried by a
    # single ecotype cannot be told apart from an assembly artefact, so any COUNT of variants should
    # exclude them; the spectra themselves keep them (a frequency spectrum without its rarest class
    # is not one), which is why both are stored.
    freqw2 = {f"{s}_{g}_{c}": [] for s in sets for g in GROUPS for c in CLASSES}
    # freqw12: only COMMON units -- panel-wide carrier count >= 12 (MAF ~5%, the same floor the rest
    # of this project uses for "common"). Same construction as freqw2, stricter threshold.
    freqw12 = {f"{s}_{g}_{c}": [] for s in sets for g in GROUPS for c in CLASSES}
    nvar = {f"{s}_{g}_{c}": 0 for s in sets for g in GROUPS for c in CLASSES}

    def add(sname, cname, carriers, called, skip_freqw=False):
        """carriers/called: [n_units x n_founders_in_set] presence and called masks -> projected SFS.
        skip_freqw: the well-called spectra of insertion/deletion alleles and bubbles are built
        together by add_freqw_paired below, on a shared event set."""
        # PANEL-LEVEL call-rate filter (2026-09-18). The well-called set is decided once on the whole
        # founder set, not per climate tertile. Filtering per tertile made the three group curves
        # different variant sets, and since call rate correlates with origin climate among the
        # PanGenie founders (rho = -0.40) that difference was itself climate-structured. A tertile now
        # contributes only "does this variant segregate here" (>= 1 carrier).
        nall = carriers.shape[1]
        mac_all = carriers.sum(1)
        panel_ok = (called.sum(1) >= CALL_MIN * nall) & (mac_all > 0)
        panel_ok2 = panel_ok & (mac_all >= 2)                 # private (panel-wide MAC 1) dropped
        panel_ok12 = panel_ok & (mac_all >= 12)               # common only
        for g in GROUPS:
            gm = np.ones(carriers.shape[1], bool) if g == "all" else (grp[sets[sname]] == g)
            k = kproj(sname, g)
            c = carriers[:, gm].sum(1).astype(np.int64); m = called[:, gm].sum(1).astype(np.int64)
            keep = (m >= k) & (c > 0) & (c < m)
            if not keep.any():
                continue
            acc[f"{sname}_{g}_{cname}"] += project(c[keep], m[keep], k).sum(0)
            nvar[f"{sname}_{g}_{cname}"] += int(keep.sum())
            # raw frequency (carriers / called), no projection -- the simple view, as in
            # sv_insertion_landing_vs_purging.ipynb; keep every variant with >=1 carrier and >=30 called
            nset = int(gm.sum())
            r = (m >= min(30, int(0.8 * nset))) & (c > 0)
            freq[f"{sname}_{g}_{cname}"].append((c[r] / m[r]).astype(np.float32))
            if not skip_freqw:
                rw = panel_ok & (c > 0) & (m > 0)
                freqw[f"{sname}_{g}_{cname}"].append((c[rw] / m[rw]).astype(np.float32))
                rw2 = panel_ok2 & (c > 0) & (m > 0)
                freqw2[f"{sname}_{g}_{cname}"].append((c[rw2] / m[rw2]).astype(np.float32))
                rw12 = panel_ok12 & (c > 0) & (m > 0)
                freqw12[f"{sname}_{g}_{cname}"].append((c[rw12] / m[rw12]).astype(np.float32))

    def add_freqw_paired(sname, acls, bcls, P, li, Pb, Cb):
        """WELL-CALLED spectra for an allele/bubble pair, on the SAME set of events.

        The call-rate filter is applied at the BUBBLE level and that bubble's alleles are carried
        through with it. Filtering each unit separately (what this script did before 2026-09-17)
        kept 29.6% of insertion alleles but 64.2% of bubbles in the 231-founder set, so the two
        curves were drawn from different event sets and were not a like-for-like contrast.

        Both units use the BUBBLE's called set as the denominator: a founder genotyped at the bubble
        is genotyped at all of its alleles (missingness in this graph VCF is per bubble, not per
        biallelic record), so the allele frequencies of a bubble then sum to the bubble's frequency
        -- the property that makes the two units comparable at all.

        The call-rate filter is PANEL-LEVEL (2026-09-18): decided once over the whole founder set, so
        every climate tertile is drawn from the same set of well-called bubbles and only has to
        supply >= 1 carrier. See `add` for why."""
        nall = Pb.shape[1]
        mac_b = Pb.sum(1); mac_a = P.sum(1)
        panel_b = (Cb.sum(1) >= CALL_MIN * nall) & (mac_b > 0)
        panel_b2 = panel_b & (mac_b >= 2)                     # bubble carried by >= 2 ecotypes
        priv_a = mac_a >= 2                                   # allele carried by >= 2 ecotypes
        panel_b12 = panel_b & (mac_b >= 12)                   # common bubble
        com_a = mac_a >= 12                                   # common allele
        for g in GROUPS:
            gm = np.ones(P.shape[1], bool) if g == "all" else (grp[sets[sname]] == g)
            cb = Pb[:, gm].sum(1).astype(np.int64); mb = Cb[:, gm].sum(1).astype(np.int64)
            keep_b = panel_b & (cb > 0) & (mb > 0)
            freqw[f"{sname}_{g}_{bcls}"].append((cb[keep_b] / mb[keep_b]).astype(np.float32))
            ca = P[:, gm].sum(1).astype(np.int64); den = mb[li]
            keep_a = keep_b[li] & (ca > 0)
            freqw[f"{sname}_{g}_{acls}"].append((ca[keep_a] / den[keep_a]).astype(np.float32))
            keep_b2 = panel_b2 & (cb > 0) & (mb > 0)
            freqw2[f"{sname}_{g}_{bcls}"].append((cb[keep_b2] / mb[keep_b2]).astype(np.float32))
            keep_a2 = keep_b2[li] & priv_a & (ca > 0)
            freqw2[f"{sname}_{g}_{acls}"].append((ca[keep_a2] / den[keep_a2]).astype(np.float32))
            keep_b12 = panel_b12 & (cb > 0) & (mb > 0)
            freqw12[f"{sname}_{g}_{bcls}"].append((cb[keep_b12] / mb[keep_b12]).astype(np.float32))
            keep_a12 = keep_b12[li] & com_a & (ca > 0)
            freqw12[f"{sname}_{g}_{acls}"].append((ca[keep_a12] / den[keep_a12]).astype(np.float32))

    for ch in lib.CHROMS:
        cl = ch.lower(); base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        assert np.array_equal(meta["founders"].astype("U6"), F.founder.to_numpy()), f"{ch}: founder order"
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64); dl = al - rl
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        bmap = pd.Series(BM[f"{ch}_bubble"], index=BM[f"{ch}_rec"])
        sub = rng.random(dl.size) < SNP_KEEP
        sel = {"snp": (rl == 1) & (al == 1) & sub, "indel": (np.abs(dl) >= 1) & (np.abs(dl) <= SV_BP) & sub,
               "ins_allele": dl > SV_BP, "del_allele": -dl > SV_BP}
        for sname, fmask in sets.items():
            for cname, s in sel.items():
                cols = np.where(s)[0]
                if cols.size == 0:
                    continue
                P = np.asarray(vp[fmask][:, cols].todense()).T.astype(bool)     # [nrec x nfounders]
                C = np.asarray(vc[fmask][:, cols].todense()).T.astype(bool)
                paired = cname in ("ins_allele", "del_allele")
                add(sname, cname, P, C, skip_freqw=paired)
                if paired:                                                       # same records, merged per bubble
                    b = bmap.reindex(cols).to_numpy().astype(np.int64)
                    assert (b >= 0).all(), f"{ch} {cname}: records with no bubble code"
                    ub, li = np.unique(b, return_inverse=True); li = li.ravel()
                    M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
                    Pb = np.asarray((M @ P) > 0); Cb = np.asarray((M @ C) > 0)
                    bcls = cname.replace("allele", "event")
                    add(sname, bcls, Pb, Cb, skip_freqw=True)
                    add_freqw_paired(sname, cname, bcls, P, li, Pb, Cb)
        print(f"[{ch}] done", flush=True)

    np.savez_compressed(f"{OUT}/sv_sfs_event.npz", nproj=NPROJ, nproj_group=NPROJ_G, classes=np.array(CLASSES, dtype=str),
                        groups=np.array(GROUPS, dtype=str), sets=np.array(list(sets), dtype=str),
                        snp_keep=SNP_KEEP, **acc, **{f"n_{k}": v for k, v in nvar.items()},
                        **{f"freq_{k}": (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in freq.items()},
                        **{f"freqw_{k}": (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in freqw.items()},
                        **{f"freqw2_{k}": (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in freqw2.items()},
                        **{f"freqw12_{k}": (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in freqw12.items()})
    for s in sets:
        k = kproj(s, "all")
        print(f"\n{s}: variants (projected to {k} called founders), and % of them singletons in the projection")
        for c in CLASSES:
            v = acc[f"{s}_all_{c}"]
            print(f"  {c:<12}{nvar[f'{s}_all_{c}']:>9,}   singleton {100*v[1]/v[1:k].sum():5.1f}%   "
                  f"freq<=5% {100*v[1:4].sum()/v[1:k].sum():5.1f}%   mean freq {np.average(np.arange(k+1),weights=v)/k:.3f}")
    print(f"\n[wrote] {OUT}/sv_sfs_event.npz")


if __name__ == "__main__":
    main()
