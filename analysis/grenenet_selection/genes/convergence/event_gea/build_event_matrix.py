#!/usr/bin/env python
"""Delta p per Truvari event: the input of the GEA's event arm (step 3 of _submit_gea.sh).

The panel splits one SV into one record per assembly path, so a single insertion appears as several alleles a few bp
apart: its frequency is split across them, and each may be too rare to test. Here non-SNP records are collapsed into
units (r1_sv_negative_selection/_build_event_units.py -> lib.SV_ADAPTIVE/nonsnp_event_units.npz): an SV
(|alt_len - ref_len| > 50 bp) is one unit per Truvari event; a short indel or MNP stays one unit per record.

  unit frequency per pool, and unit p0   sum over member records (members are carried by disjoint founders), capped at 1
  delta p                                unit pool frequency - unit p0
  MAF                                    min(p_bar, 1 - p_bar), p_bar = mean unit frequency over pools; kept if >= 0.05
  founder call rate                      a unit is called at a founder if any member is (the convention of
                                         _sv_event_catalog.py); kept if called in >= ceil(--min-call-rate * 231)

Everything else is the record arm's: the same generation-9 pools, p0, filters and LFMM (build_class_matrices.py,
build_lfmm_input.py). Two scans: `sv` (SV events only) and `nonsnp` (every unit).

Check: a single-record short-indel unit must reproduce the record arm's delta p for that record (--cmdir).

Writes in --out: {sv,nonsnp}_gen9_units.csv (one row per unit: unit,chrom,pos,ref_len,alt_len,dlen,n_alleles,is_sv,
p_bar,MAF,n_called,block) and {sv,nonsnp}_gen9_dp.npy (float32, pools x units, columns = rows of units.csv).
Environment: kmate. Compute node, ~120 GB.
"""
import argparse, os, sys
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))
import lib                                                       # noqa: E402

MAF_MIN, SV_BP, N_FOUNDERS = 0.05, 50, 231


def unit_n_called(M: sp.csr_matrix, n_rec: int) -> np.ndarray:
    """Founders called at each unit = the union over its member records, from the panel's var_called.

    The panel's non-SNP records are af_store's index_nonsnp records in the same order (asserted per chromosome);
    nothing is joined on (chrom, pos, ...).
    """
    blocks, off = [], 0
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    ch = idx["chrom"].astype(str)
    rl = idx["ref_len"].astype(np.int64)
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
    return np.asarray(((M @ VC.T.tocsr()) > 0).sum(1)).ravel().astype(int)


def check_short_indels(dp, k, ptr, ind, sizes, is_sv, p0, cmdir):
    """Single-record short-indel units against the record arm's smallindel matrix: delta p must agree."""
    R = pd.read_csv(f"{cmdir}/smallindel_gen9.records.csv", usecols=["col"])
    A = np.load(f"{cmdir}/smallindel_gen9_af.npy", mmap_mode="r")
    colpos = pd.Series(np.arange(len(R)), index=R["col"].to_numpy())
    single = np.where((sizes[k] == 1) & ~is_sv[k])[0]
    pick = np.random.default_rng(0).choice(single, min(2000, len(single)), replace=False)
    recs = ind[ptr[k[pick]]]
    have = colpos.reindex(recs).notna().to_numpy()          # all, unless --cmdir was built with another floor
    assert have.any(), f"no short-indel unit found in {cmdir}"
    rec_arm = np.asarray(A[:, colpos.reindex(recs[have]).astype(int).to_numpy()], np.float64) - p0[recs[have]][None, :]
    diff = np.abs(rec_arm - dp[:, pick[have]].astype(np.float64)).max()
    print(f"  CHECK short-indel units vs record arm: {have.sum()} compared, max |diff| = {diff:.2e}", flush=True)
    assert diff < 1e-6, "event-arm delta p disagrees with the record arm"


def main(out: str, min_call_rate: float, cmdir: str):
    lib.check_version_dir(out, "--out")                     # never mix panel versions (KMATE_PANEL)
    os.makedirs(out, exist_ok=True)
    U = np.load(f"{lib.SV_ADAPTIVE}/nonsnp_event_units.npz")
    ptr, ind, dlen, rep = U["indptr"], U["indices"], U["dlen"], U["rep"]
    nu = len(ptr) - 1; sizes = np.diff(ptr)
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    af = np.load(f"{lib.POOL_MATRICES}/pool_gen9_nonsnp_af.npy")   # pools x records
    pm = pd.read_csv(f"{lib.POOL_MATRICES}/pool_gen9_nonsnp.meta.csv")
    assert len(pm) == af.shape[0] and not pm.site.astype(int).isin(lib.EXCLUDED_SITES).any(), "pool table / excluded site"
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    n_pools, n_rec = af.shape
    assert len(p0) == n_rec and ind.max() < n_rec
    print(f"{n_pools} pools x {n_rec:,} records -> {nu:,} units", flush=True)

    M = sp.csr_matrix((np.ones(len(ind)), (np.repeat(np.arange(nu), sizes), ind)), shape=(nu, n_rec))
    up0 = np.minimum(M @ p0.astype(np.float64), 1.0)
    pbar = np.zeros(nu)
    for i in range(n_pools):                                # streamed over pools to bound memory
        pbar += np.minimum(M @ af[i].astype(np.float64), 1.0)
    pbar /= n_pools
    maf = np.minimum(pbar, 1 - pbar)
    is_sv = np.abs(dlen) > SV_BP
    ncall = unit_n_called(M, n_rec)
    floor = int(np.ceil(min_call_rate * N_FOUNDERS))
    keep_all = (maf >= MAF_MIN) & (ncall >= floor)
    print(f"units with MAF >= {MAF_MIN} and called in >= {floor}/{N_FOUNDERS} founders: {keep_all.sum():,} "
          f"({(keep_all & is_sv).sum():,} SV events, {(keep_all & ~is_sv).sum():,} short-indel/MNP units)", flush=True)

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
        if cls == "nonsnp":
            check_short_indels(dp, k, ptr, ind, sizes, is_sv, p0, cmdir)
        np.save(f"{out}/{cls}_gen9_dp.npy", dp)
        T.to_csv(f"{out}/{cls}_gen9_units.csv", index=False)
        print(f"  {cls}: {len(k):,} units -> {cls}_gen9_dp.npy / _units.csv", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-call-rate", type=float, default=0.50,
                    help="founder call-rate floor on the union of a unit's members (n_called / 231)")
    ap.add_argument("--cmdir", default=lib.GEA_CLASS_MATRICES,
                    help="the record arm's class matrices, for the short-indel check")
    a = ap.parse_args()
    main(a.out, a.min_call_rate, a.cmdir)
