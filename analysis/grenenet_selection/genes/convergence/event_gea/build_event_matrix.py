#!/usr/bin/env python
"""Inputs for an EVENT-level non-SNP GEA -- the multi-axis LFMM with SVs counted as events.

Why. The GEA tests every arch3 record. arch3 splits one SV into one record per assembly path,
so a single insertion appears as several "alleles" that differ by a few bp: its frequency is
split across them, alleles that are each rare never reach MAF 0.05 and are never tested, and
an SV whose alleles do reach it is tested several times over. The burden and purging work
collapses those records into truvari events (notebooks/sv_purging_fixed.ipynb). Counted that
way, 27,009 SV events clear MAF 0.05 -- against 21,805 events the record-level GEA touches --
and 5,369 of them are NEW: no single allele was common enough, their alleles together are.

What is identical to the GEA (r2_gea_nonsnp/phase1_replication/build_class_matrices.py +
build_lfmm_input.py): the pool frequency store (gen 9, 352 pools), the founding frequency p0,
Y = delta p = pool frequency - p0, MAF = min(p_bar, 1 - p_bar) over pools with MAF >= 0.05 (the only filter of the default
run; the 2026-09-30 re-run adds `--min-call-rate 0.50`, below), the z-scored climate axis, and the LFMM (ridge, K = 16) run downstream.

What changes is only the unit, reused from r1_sv_negative_selection/results/sv_adaptive/
nonsnp_event_units.npz (_build_event_units.py): short indels (|dlen| <= 50) and MNPs stay one
unit per record; SVs (|dlen| > 50, the GEA's own SV threshold) are one unit per truvari event.
A unit's pool frequency and p0 are the SUMS over its member records (members are carried by
disjoint founders), capped at 1.

Call rate (the 2026-09-30 re-run, `--min-call-rate 0.50`): a unit is CALLED at a founder if ANY member record is
called there (union), and must be called in >= ceil(0.5 * 231) = 116 founders. The record-level classes (snp,
smallindel) apply the same 116-founder floor per record, so record == unit for them; for SVs the event unit is
what the floor acts on, so the event and record arms keep different numbers of SV rows by design.

Two scans, as the GEA ran them: `nonsnp` (every unit, joint latent factors) and `sv` (SV events
only). The small-indel scan is unchanged at the unit level and is not repeated.

Check built in: a small-indel unit is a single record, so its delta p must equal the GEA's own
input column for that record exactly.

Writes event_gea/{cls}_gen9_units.csv and event_gea/{cls}_gen9_dp.npy (pools x units, float32).
env: kmate. Compute node (~20 GB RAM).
"""
import argparse, os, sys
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, GEA)
import lib                                                       # noqa: E402

UNITS = f"{lib.SV_ADAPTIVE}/nonsnp_event_units.npz"               # the panel version's event units (KMATE_PANEL)
CMDIR = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices_cr50{lib.VSUF}"   # for the CHECK only
MAF_MIN, SV_BP = 0.05, 50
N_FOUNDERS = 231


def unit_n_called(M: sp.csr_matrix, n_rec: int) -> np.ndarray:
    """Founders called at a unit = the UNION over its member records.

    The same convention as `_sv_event_catalog.py:49` (`(M @ var_called) > 0`), so the event
    GEA and the burden model agree on what "called" means for an event. Members of a truvari
    event sit in one graph bubble and almost always share a called set -- 42% of multi-allele
    events differ at all, and where they do the median spread is 2 founders -- so union, max
    and representative give nearly the same answer; union is chosen to match the catalogue.

    Built from the panel's var_called in af_store record order (the same positional
    correspondence `lib.founder_record_keep` asserts). Never joined on (chrom,pos,...).
    """
    blocks, off = [], 0
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    ch = idx["chrom"].astype(str)
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    for c in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        n = int((ch == c).sum())
        base = lib.var_pa_prefix(c)
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        panel = np.where(~((prl == 1) & (pal == 1)))[0]
        assert panel.size == n and np.array_equal(prl[panel], rl[off:off + n]), \
            f"{c}: af_store/panel non-SNP record order differs"
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()[:, panel]
        blocks.append((vc > 0).astype(np.int8))
        off += n
    VC = sp.hstack(blocks, format="csr")                    # 231 x n_rec
    assert VC.shape[1] == n_rec
    UQ = (M @ VC.T.tocsr()) > 0                             # units x 231, union of members
    return np.asarray(UQ.sum(1)).ravel().astype(int)


def main(out: str, min_call_rate: float):
    lib.check_version_dir(out, "--out")                 # never mix panel versions (KMATE_PANEL)
    os.makedirs(out, exist_ok=True)
    U = np.load(UNITS)
    ptr, ind, dlen, rep = U["indptr"], U["indices"], U["dlen"], U["rep"]
    nu = len(ptr) - 1; sizes = np.diff(ptr)
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    af = np.load(f"{lib.POOL_MATRICES}/pool_gen9_nonsnp_af.npy")   # pools x records (the panel version's)
    pm = pd.read_csv(f"{lib.POOL_MATRICES}/pool_gen9_nonsnp.meta.csv")
    assert len(pm) == af.shape[0] and not pm.site.astype(int).isin(lib.EXCLUDED_SITES).any(), "pool table / excluded site"
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    n_pools, n_rec = af.shape
    assert len(p0) == n_rec and ind.max() < n_rec
    print(f"{n_pools} pools x {n_rec:,} records -> {nu:,} units", flush=True)

    M = sp.csr_matrix((np.ones(len(ind)), (np.repeat(np.arange(nu), sizes), ind)), shape=(nu, n_rec))
    up0 = np.minimum(M @ p0.astype(np.float64), 1.0)
    # unit frequency per pool = sum of member frequencies; streamed over pools to bound memory
    pbar = np.zeros(nu)
    for i in range(n_pools):
        pbar += np.minimum(M @ af[i].astype(np.float64), 1.0)
    pbar /= n_pools
    maf = np.minimum(pbar, 1 - pbar)
    is_sv = np.abs(dlen) > SV_BP
    keep_all = maf >= MAF_MIN
    print(f"units with MAF >= {MAF_MIN}: {keep_all.sum():,} "
          f"({(keep_all & is_sv).sum():,} SV events, {(keep_all & ~is_sv).sum():,} short-indel/MNP units)", flush=True)

    # founder call-rate floor, the same filter the record-level GEA now applies
    # (r2_gea_nonsnp/phase1_replication/build_class_matrices.py). Collapsing records into
    # events does NOT rescue call rate -- members share a bubble, so they fail together:
    # of 50,280 SV records below 0.50, only 13 belong to an event with a better-called
    # member. The two filters are orthogonal and both are needed.
    ncall = unit_n_called(M, n_rec)
    if min_call_rate > 0:
        bar = int(np.ceil(min_call_rate * N_FOUNDERS))
        cr_keep = ncall >= bar
        print(f"call-rate floor {min_call_rate:.2f} (>= {bar}/{N_FOUNDERS} founders): "
              f"of the MAF-passing units it removes {int((keep_all & ~cr_keep).sum()):,} "
              f"({(keep_all & ~cr_keep & is_sv).sum():,} SV events)", flush=True)
        keep_all = keep_all & cr_keep
        print(f"units after both filters: {keep_all.sum():,} "
              f"({(keep_all & is_sv).sum():,} SV events, "
              f"{(keep_all & ~is_sv).sum():,} short-indel/MNP units)", flush=True)

    ch = idx["chrom"].astype(str); pos = idx["pos"]; rl = idx["ref_len"]; al = idx["alt_len"]
    for cls, keep in (("nonsnp", keep_all), ("sv", keep_all & is_sv)):
        k = np.where(keep)[0]
        Mk = M[k]
        dp = np.empty((n_pools, len(k)), dtype=np.float32)
        for i in range(n_pools):
            dp[i] = np.minimum(Mk @ af[i].astype(np.float64), 1.0) - up0[k]
        r = rep[k]
        T = pd.DataFrame(dict(unit=k, chrom=ch[r], pos=pos[r], ref_len=rl[r], alt_len=al[r],
                              dlen=dlen[k], n_alleles=sizes[k], is_sv=is_sv[k],
                              p_bar=pbar[k], MAF=maf[k], n_called=ncall[k]))
        T["block"] = lib.assign_ld_blocks(T.chrom.to_numpy(str), T.pos.to_numpy())
        np.save(f"{out}/{cls}_gen9_dp.npy", dp)
        T.to_csv(f"{out}/{cls}_gen9_units.csv", index=False)
        print(f"  {cls}: {len(k):,} units -> {cls}_gen9_dp.npy / _units.csv", flush=True)

        if cls == "nonsnp":
            # CHECK: single-record short-indel units must reproduce the GEA's own delta p
            R = pd.read_csv(f"{CMDIR}/nonsnp_gen9.records.csv", usecols=["col"])
            A = np.load(f"{CMDIR}/nonsnp_gen9_af.npy", mmap_mode="r")
            colpos = pd.Series(np.arange(len(R)), index=R["col"].to_numpy())
            single = np.where((sizes[k] == 1) & ~is_sv[k])[0]
            rng = np.random.default_rng(0)
            pick = rng.choice(single, min(2000, len(single)), replace=False)
            recs = ind[ptr[k[pick]]]
            have = colpos.reindex(recs).notna().to_numpy()
            gea = np.asarray(A[:, colpos.reindex(recs[have]).astype(int).to_numpy()], np.float64) - p0[recs[have]][None, :]
            mine = dp[:, pick[have]].astype(np.float64)
            print(f"  CHECK short-indel units vs GEA input: {have.sum()} compared, "
                  f"max |diff| = {np.abs(gea - mine).max():.2e}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--min-call-rate", type=float, default=0.0,
                    help="founder call-rate floor on the UNION of a unit's members; "
                         "0 = off (the original behaviour), 0.50 is the calibrated setting")
    a = ap.parse_args()
    main(a.out, a.min_call_rate)
