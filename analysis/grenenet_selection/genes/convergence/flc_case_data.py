#!/usr/bin/env python
"""Assemble everything the FLC case study needs, into one small npz + CSVs.

The case: FLC was Bonferroni-significant in our pre-floor non-SNP GEA on bio6, and absent from the
phase-1 SNP GEA run with the same model on the same samples. The question the notebook has
to answer is whether that is a discovery or an artefact of who is in the panel at that locus.

Our scan is the production GEA, r2_gea_nonsnp/phase1_replication/results/multiaxis/tile_cr50
(MAF >= 0.05 AND founder call rate >= 0.50; raw LFMM p; snp/smallindel one row per record,
sv/nonsnp one row per Truvari EVENT `unit`). The hit this case examines was found in the
superseded pre-floor scan (_superseded_pre_cr50/wza_in_clq09_tile); FLC's records sit at
78/231 called (0.34) and fall below the 0.50 floor this case motivated, so under tile_cr50 the
FLC deletion is not tested at all -- region_pvals.csv says so (main() prints it).

Writes results/flc_case/:
  region_pvals.csv    our LFMM p across FLC +/-10 kb (tile_cr50), per class x axis, with the
                      founder-panel call rate of each row (event union for unit rows)
  region_pvals_pre_cr50.csv  the same window in the superseded pre-floor scan (records, no
                      call-rate floor): the hit the case examines, kept as its subject only
  phase1_pvals.csv    the phase-1 panel's own LFMM p (same model, K=16, 355 last-gen samples)
                      across the same window
  dp_three_ways.npz   per-pool Delta p at the 10 shared SNPs under three estimates:
                        kmate  = production, MAR over the 78 called founders
                        filled = kMate's OWN h, but the 153 genotypes filled from phase-1
                        hapfire= the phase-1 pipeline's own estimate
                      plus site and bio6 per pool
  panel_window.csv    every panel record in the window: position, size, carriers, called

env: kmate. Compute node.
"""
from __future__ import annotations
import os
import subprocess
import sys
import glob
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402

OUT = f"{HERE}/results/flc_case"
BCF = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools"
# production GEA (README_PRODUCTION.txt there); classes have different column sets, read by header
GEAIN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/tile_cr50"
# the pre-floor scan in which the FLC hit was found -- read ON PURPOSE, only to keep the case's
# own subject on record (region_pvals_pre_cr50.csv); never a current result
PRECR50 = (f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/"
           "_superseded_pre_cr50/wza_in_clq09_tile")
UNITS_NPZ = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/nonsnp_event_units.npz"
P1 = ("/global/scratch/users/tbellg/gea_grene-net/lfmm/lfmm_last_gen/lfmm_fullresults_all_k")
P1DP = "/global/scratch/users/tbellg/gea_grene-net/key_files/delta_p_maf05_mincount05_lastgensamples.csv"
Z = "/global/scratch/projects/fc_moilab/projects/grenenet-phase1/zenodo/greneNet_final_v1.1"

GENE = ("Chr5", 3173497, 3179448)            # FLC, minus strand; TSS at the high end
INTRON1 = (3175662, 3179154)
LO, HI = 3163497, 3189448                    # the +/-10 kb window
SNPS = [3170443, 3177330, 3178515, 3178749, 3178835, 3178838, 3178847, 3178855,
        3185704, 3185978]                    # the SNPs Bonferroni-significant in our scan
# their row indices in the phase-1 delta-p matrix (1-based data rows), same order
P1ROWS = [823635, 823692, 823706, 823709, 823712, 823713, 823714, 823716, 823782, 823784]
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1", "pc2", "pc3"]
FLCDEL_UNIT = 2055022                        # the 6 bp intron-1 deletion Chr5:3175803 (single-record unit)


def unit_info(units) -> pd.DataFrame:
    """unit -> n_members, union carriers and union called founders (nonsnp_event_units.npz)."""
    U = np.load(UNITS_NPZ, allow_pickle=True)
    u = np.asarray(units, np.int64)
    return pd.DataFrame(dict(unit=u, n_members=np.diff(U["indptr"])[u],
                             u_carriers=U["n"][u].astype(int), u_called=U["ncal"][u].astype(int)))


def row_units(rows) -> np.ndarray:
    """af_store index_nonsnp rows -> their unit id (each row belongs to exactly one unit)."""
    U = np.load(UNITS_NPZ, allow_pickle=True)
    r2u = np.empty(int(U["indices"].max()) + 1, np.int64)
    r2u[U["indices"]] = np.repeat(np.arange(len(U["indptr"]) - 1), np.diff(U["indptr"]))
    return r2u[np.asarray(rows, np.int64)]


def our_pvals(gea: str = GEAIN) -> pd.DataFrame:
    """LFMM p in the window from `gea` (default tile_cr50), every class x axis, with call rates.

    snp / smallindel rows are records; sv / nonsnp rows are Truvari events (`unit`), placed at
    the event's representative allele. Call rate: unit rows take the event's union over its
    member records (nonsnp_event_units.npz); SNP rows are matched to the panel record on
    (pos, 1, 1) -- they carry no unit id.
    """
    rows = []
    for f in sorted(glob.glob(f"{gea}/lfmm_*_gen9_*.csv")):
        b = os.path.basename(f)[len("lfmm_"):-len(".csv")]
        cls, ax = b.split("_gen9_")
        hdr = open(f).readline()
        out = subprocess.run(
            ["awk", "-F,", f'NR>1 && $1=="Chr5" && $2+0>={LO} && $2+0<={HI}',
             f], capture_output=True, text=True).stdout
        if not out.strip():
            continue
        d = pd.read_csv(pd.io.common.StringIO(hdr + out))
        d["cls"], d["axis"] = cls, ax
        # the production threshold is genome-wide for this class x axis
        n = sum(1 for _ in open(f)) - 1
        d["bonf"] = 0.05 / n
        rows.append(d)
    D = pd.concat(rows, ignore_index=True)
    # founder-panel call rate per record
    base = f"{lib.PROJ}/panel/arch3/chr5/var_pa_231_arch3_chr5"
    m = np.load(f"{base}.meta.npz", allow_pickle=True)
    vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
    sel = np.where((m["pos"] >= LO) & (m["pos"] <= HI))[0]
    cr = pd.DataFrame(dict(pos=m["pos"][sel], ref_len=m["ref_len"][sel].astype(int),
                           alt_len=m["alt_len"][sel].astype(int),
                           n_called=np.asarray((vc[:, sel] > 0).sum(0)).ravel(),
                           n_carriers=np.asarray(
                               (sp.load_npz(f"{base}.var_pa.npz").tocsc()[:, sel] > 0).sum(0)).ravel()))
    cr["call_rate"] = cr.n_called / 231
    cr.to_csv(f"{OUT}/panel_window.csv", index=False)
    # unit rows (smallindel / sv / nonsnp): the event's own union carriers + called founders
    D = D.drop(columns=["n_called"], errors="ignore")
    hasu = D["unit"].notna() if "unit" in D else pd.Series(False, index=D.index)
    D["n_called"] = np.nan; D["n_carriers"] = np.nan
    if hasu.any():
        UI = unit_info(D.loc[hasu, "unit"].astype(np.int64).to_numpy())
        D.loc[hasu, "n_members"] = UI.n_members.to_numpy()
        D.loc[hasu, "n_called"] = UI.u_called.to_numpy()
        D.loc[hasu, "n_carriers"] = UI.u_carriers.to_numpy()
    # rows without a unit id (SNPs; every row of the pre-floor scan): the panel record with the
    # same (pos, ref_len, alt_len), left NaN where that key is ambiguous (multiallelic trap)
    KEY = ["pos", "ref_len", "alt_len"]
    cs = cr.drop_duplicates(KEY, keep=False).set_index(KEY)
    k = pd.MultiIndex.from_frame(D.loc[~hasu, KEY].astype(int))
    D.loc[~hasu, "n_called"] = cs.n_called.reindex(k).to_numpy()
    D.loc[~hasu, "n_carriers"] = cs.n_carriers.reindex(k).to_numpy()
    D["call_rate"] = D.n_called / 231
    return D


def report_flc(D: pd.DataFrame) -> None:
    """Say plainly whether the FLC deletion and the 10 FLC SNPs are tested in tile_cr50."""
    b6 = D[D.axis == "bio6"]
    hit = b6[-np.log10(b6.pval) > -np.log10(b6.bonf)]
    print(f"  tile_cr50, bio6, FLC +/-10 kb: {b6.groupby('cls').size().to_dict()} rows; "
          f"Bonferroni hits {hit.groupby('cls').size().to_dict() or 0}")
    if "unit" in D and (D.unit == FLCDEL_UNIT).any():
        q = b6[b6.unit == FLCDEL_UNIT]
        print(f"  FLC 6 bp deletion (unit {FLCDEL_UNIT}): -log10 p bio6 = "
              f"{-np.log10(q.pval).round(2).to_dict()}")
    else:
        u = unit_info([FLCDEL_UNIT]).iloc[0]
        print(f"  FLC 6 bp deletion (unit {FLCDEL_UNIT}) is NOT in tile_cr50: called in "
              f"{u.u_called}/231 founders, below the 0.50 call-rate floor (116)")
    s = b6[(b6.cls == "snp") & b6.pos.isin(SNPS)]
    print(f"  the {len(SNPS)} FLC SNPs of the pre-floor scan: {len(s)} tested in tile_cr50")


def phase1_pvals() -> pd.DataFrame:
    """The phase-1 panel's own LFMM (same lfmm_ridge, K=16, 355 last-gen samples)."""
    ids = pd.read_csv(f"{P1}/lfmm_bio6_k16_results.csv", usecols=["snp_id"])
    ch = ids.snp_id.str.split("_", expand=True)
    ids["chr"], ids["pos"] = ch[0], ch[1].astype(int)
    m = (ids["chr"] == "5") & ids.pos.between(LO, HI)
    idx = np.where(m)[0]
    out = pd.DataFrame(dict(pos=ids.pos.values[idx]))
    n = len(ids)
    for ax in AXES:
        f = f"{P1}/wo_calibration_pvalue_full_genome_{ax}_k16.csv"
        if not os.path.exists(f):
            continue
        out[ax] = pd.read_csv(f, usecols=[ax]).iloc[idx, 0].values
    out["bonf"] = 0.05 / n
    out["n_tested"] = n
    return out


def dp_three_ways() -> None:
    """Per-pool Delta p at the 10 SNPs: production, genotype-filled, and hapFIRE."""
    base = f"{lib.PROJ}/panel/arch3/chr5/var_pa_231_arch3_chr5"
    m = np.load(f"{base}.meta.npz", allow_pickle=True)
    founders = m["founders"].astype(str)
    cols = np.array([np.where((m["pos"] == p) & (m["ref_len"] == 1) & (m["alt_len"] == 1))[0][0]
                     for p in SNPS])
    G78 = np.asarray((sp.load_npz(f"{base}.var_pa.npz").tocsc()[:, cols] > 0).todense()).astype(float)
    C78 = np.asarray((sp.load_npz(f"{base}.var_called.npz").tocsc()[:, cols] > 0).todense()).astype(float)
    alts = [str(m["alt"][c]) for c in cols]

    fam = pd.read_csv(Z + ".fam", sep=r"\s+", header=None)[0].astype(str).to_numpy()
    bim = pd.read_csv(Z + ".bim", sep="\t", header=None,
                      names=["chr", "id", "cm", "pos", "a1", "a2"])
    bim["row"] = np.arange(len(bim))
    b5 = bim[(bim["chr"] == 5) & bim.pos.isin(SNPS)].drop_duplicates("pos").set_index("pos")
    N = len(fam); nb = (N + 3) // 4
    LUT = np.array([[{0: 2, 1: np.nan, 2: 1, 3: 0}[(x >> (2 * k)) & 3] for k in range(4)]
                    for x in range(256)], float)
    order = [list(fam).index(f) for f in founders]
    G231 = np.full((231, len(SNPS)), np.nan)
    with open(Z + ".bed", "rb") as fh:
        fh.read(3)
        for j, p in enumerate(SNPS):
            r = b5.loc[p]
            fh.seek(3 + int(r.row) * nb)
            g = LUT[np.frombuffer(fh.read(nb), dtype=np.uint8)].ravel()[:N] / 2.0
            gg = g if r.a1 == alts[j] else (1 - g if r.a2 == alts[j] else np.full(N, np.nan))
            G231[:, j] = gg[order]

    meta = pd.read_csv(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp.meta.csv")
    PT = lib.pool_table(); PT = PT[PT.pool.isin(set(meta.pool))]
    def hmap(basedir):
        out = {}
        for f in glob.glob(f"{basedir}/*_Chr5.h_per_chrom.npz"):
            d = np.load(f, allow_pickle=True)
            out[os.path.basename(f).split("_Chr5")[0]] = d["Chr5"]
        return out
    HE, HS = hmap(lib.OUT), hmap(lib.SEEDMIX)
    pools = list(meta.pool)
    H = np.zeros((len(pools), 231))
    for i, pl in enumerate(pools):
        sub = PT[PT.pool == pl]
        w, hs = [], []
        for r in sub.itertuples():
            if r.sampleid in HE:
                w.append(max(r.flowerscollected, 1e-9)); hs.append(HE[r.sampleid])
        H[i] = np.average(np.array(hs), axis=0, weights=np.array(w))
    h0 = np.mean(np.array(list(HS.values())), axis=0)

    def mar(h):  return (h @ (G78 * C78)) / np.maximum(h @ C78, 1e-12)
    def full(h): return h @ np.nan_to_num(G231)
    dK = np.array([mar(H[i]) for i in range(len(pools))]) - mar(h0)
    dF = np.array([full(H[i]) for i in range(len(pools))]) - full(h0)

    hdr = open(f"{P1DP}", "r").readline().rstrip("\n").split(",")
    want = {r: i for i, r in enumerate(P1ROWS)}
    rows = {}
    with open(P1DP) as fh:
        fh.readline()
        for n, line in enumerate(fh, start=1):
            if n in want:
                rows[n] = np.array(line.rstrip("\n").split(","), dtype=float)
            if n > max(P1ROWS):
                break
    Hp = pd.DataFrame({SNPS[want[r]]: v for r, v in rows.items()}, index=hdr)
    shared = [p for p in pools if p in Hp.index]
    ix = [pools.index(p) for p in shared]
    np.savez(f"{OUT}/dp_three_ways.npz", snps=np.array(SNPS), pools=np.array(shared),
             kmate=dK[ix], filled=dF[ix], hapfire=Hp.loc[shared, SNPS].to_numpy(),
             site=meta.site.values[ix], bio6=meta.bio6.values[ix],
             h_mass_78=(H[ix] @ (C78[:, 0] > 0).astype(float)))
    print(f"  dp_three_ways: {len(shared)} shared pools x {len(SNPS)} SNPs")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    print("our per-record p in the window ...", flush=True)
    R = our_pvals()
    R.to_csv(f"{OUT}/region_pvals.csv", index=False)
    report_flc(R)
    print("the superseded pre-floor scan, same window (the case's subject) ...", flush=True)
    our_pvals(PRECR50).to_csv(f"{OUT}/region_pvals_pre_cr50.csv", index=False)
    print("phase-1 per-record p in the window ...", flush=True)
    phase1_pvals().to_csv(f"{OUT}/phase1_pvals.csv", index=False)
    print("three Delta-p estimates ...", flush=True)
    dp_three_ways()
    print(f"wrote {OUT}/")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------- appended: Delta p tracks
def dp_tracks() -> None:
    """Mean Delta p per RECORD across the window, for both panels -- the data under the p.

    A p-value says how well a record tracks climate; it hides whether the two panels even
    agree about how the allele moved. This writes, for every record in the window, the mean
    Delta p over the pools, from kMate (both classes) and from the phase-1 matrix.
    """
    rows = []
    for kind in ("snp", "nonsnp"):
        idx = np.load(f"{lib.AF_STORE}/index_{kind}.npz", allow_pickle=True)
        ch = idx["chrom"].astype(str)
        sel = np.where((ch == "Chr5") & (idx["pos"] >= LO) & (idx["pos"] <= HI))[0]
        af = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_{kind}_af.npy",
                     mmap_mode="r")
        p0 = np.load(f"{lib.AF_STORE}/p0_{kind}.npy")
        block = np.asarray(af[:, sel.min():sel.max() + 1], dtype=np.float64)
        sub = block[:, sel - sel.min()]
        rows.append(pd.DataFrame(dict(
            pos=idx["pos"][sel], ref_len=idx["ref_len"][sel].astype(int),
            alt_len=idx["alt_len"][sel].astype(int), kind=kind,
            mean_dp=sub.mean(0) - p0[sel], p0=p0[sel])))
        print(f"  kMate dp {kind}: {len(sel)} records", flush=True)
    pd.concat(rows, ignore_index=True).to_csv(f"{OUT}/dp_track_kmate.csv", index=False)

    ids = pd.read_csv(f"{P1}/lfmm_bio6_k16_results.csv", usecols=["snp_id"])
    c = ids.snp_id.str.split("_", expand=True)
    ids["chr"], ids["pos"] = c[0], c[1].astype(int)
    m = (ids["chr"] == "5") & ids.pos.between(LO, HI)
    want = {int(i) + 1: int(p) for i, p in zip(np.where(m)[0], ids.pos.values[m])}
    got = {}
    with open(P1DP) as fh:
        fh.readline()
        for n, line in enumerate(fh, start=1):
            if n in want:
                got[want[n]] = np.fromstring(line, sep=",").mean()
            if n > max(want):
                break
    pd.DataFrame(dict(pos=list(got), mean_dp=list(got.values()))).to_csv(
        f"{OUT}/dp_track_phase1.csv", index=False)
    print(f"  phase-1 dp: {len(got)} SNPs")


# ------------------------------------------------- appended: do the founder mixtures agree?
def h_compare() -> None:
    """kMate's inferred founder mixture vs the phase-1 pipeline's ecotype frequencies.

    The decisive test shows that filling the missing genotypes makes the two pipelines agree.
    That only isolates the PROJECTION if the thing being projected -- the founder mixture --
    already agrees. Here it is, per pool: kMate's Chr5 `h` (231 founders) against phase-1's
    `merged_ecotype_frequency.txt` for the same pool.

    Row order of that matrix is validated against `founder_ecotype_frequency.txt`
    (r = 0.263 for the row means vs the founding frequencies, against -0.007 after shuffling
    the founder list, so the rows are aligned; the modest value is three generations of
    selection moving the frequencies, not a misalignment).
    """
    ZD = "/global/scratch/projects/fc_moilab/projects/grenenet-phase1/zenodo"
    F = pd.read_csv(f"{ZD}/founder_ecotype_frequency.txt", sep="\t", header=None,
                    names=["eco", "p0"])
    M = pd.read_csv(f"{ZD}/merged_ecotype_frequency.txt", sep="\t")
    M.index = F.eco.astype(str).values
    meta = pd.read_csv(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp.meta.csv")
    PT = lib.pool_table(); PT = PT[PT.pool.isin(set(meta.pool))]
    founders = np.load(f"{lib.PROJ}/panel/arch3/chr5/var_pa_231_arch3_chr5.meta.npz",
                       allow_pickle=True)["founders"].astype(str)
    HE = {}
    for f in glob.glob(f"{lib.OUT}/*_Chr5.h_per_chrom.npz"):
        HE[os.path.basename(f).split("_Chr5")[0]] = np.load(f, allow_pickle=True)["Chr5"]
    rows, pairs = [], []
    for pl in meta.pool:
        if pl not in M.columns:
            continue
        sub = PT[PT.pool == pl]
        w, hs = [], []
        for r in sub.itertuples():
            if r.sampleid in HE:
                w.append(max(r.flowerscollected, 1e-9)); hs.append(HE[r.sampleid])
        if not hs:
            continue
        h = np.average(np.array(hs), axis=0, weights=np.array(w))
        p1 = M.loc[founders, pl].to_numpy(float)
        rows.append(dict(pool=pl, r=np.corrcoef(h, p1)[0, 1],
                         spearman=pd.Series(h).corr(pd.Series(p1), method="spearman")))
        pairs.append(np.c_[np.full(231, pl), founders, h, p1])
    pd.DataFrame(rows).to_csv(f"{OUT}/h_compare_pools.csv", index=False)
    P = pd.DataFrame(np.vstack(pairs), columns=["pool", "founder", "h_kmate", "f_phase1"])
    P[["h_kmate", "f_phase1"]] = P[["h_kmate", "f_phase1"]].astype(float)
    P.to_csv(f"{OUT}/h_compare_pairs.csv.gz", index=False)
    print(f"  h comparison: {len(rows)} pools, median r = {np.median([r['r'] for r in rows]):.3f}")


# --------------------------------------- appended: the regression coefficient on bio6
def coef_track() -> None:
    """Per-record slope of Delta p on bio6 -- the effect size behind each p-value.

    Mean Delta p (dp_tracks) is NOT what the GEA tests: a record can fall a lot on average
    without that fall tracking climate. The test is the association with the axis, so this
    is the quantity to put beside the p-values.

    Computed identically for both panels -- OLS slope of each record's per-pool Delta p on
    the pool's bio6, over the 352 pools both pipelines share -- so panels C and D are the
    same statistic. LFMM's own coefficient conditions on 16 latent factors, so it is not the
    same number; the stored LFMM B (ours for non-SNP from tile_cr50, joined on the event
    `unit`; phase-1's for SNPs) is carried alongside and the agreement reported.
    """
    meta = pd.read_csv(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp.meta.csv")
    bio6 = meta.bio6.to_numpy(float)
    pools = list(meta.pool)

    def slope(Y, x):                       # Y: pools x records
        xc = x - x.mean()
        return (Y - Y.mean(0)).T @ xc / (xc @ xc)

    rows = []
    for kind in ("snp", "nonsnp"):
        idx = np.load(f"{lib.AF_STORE}/index_{kind}.npz", allow_pickle=True)
        ch = idx["chrom"].astype(str)
        sel = np.where((ch == "Chr5") & (idx["pos"] >= LO) & (idx["pos"] <= HI))[0]
        af = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_{kind}_af.npy",
                     mmap_mode="r")
        p0 = np.load(f"{lib.AF_STORE}/p0_{kind}.npy")
        blk = np.asarray(af[:, sel.min():sel.max() + 1], dtype=np.float64)[:, sel - sel.min()]
        dp = blk - p0[sel]
        rows.append(pd.DataFrame(dict(
            pos=idx["pos"][sel], ref_len=idx["ref_len"][sel].astype(int),
            alt_len=idx["alt_len"][sel].astype(int), kind=kind,
            coef=slope(dp, bio6), mean_dp=dp.mean(0))))
        if kind == "nonsnp":
            rows_ns = sel
        print(f"  kMate coef {kind}: {len(sel)}", flush=True)
    K = pd.concat(rows, ignore_index=True)
    # the stored LFMM coefficient B, where we have it: tile_cr50's event-level nonsnp scan on
    # bio6, joined on `unit` (every member record of an event gets the event's B). The
    # comparison with the per-record slope is only like-for-like for single-record units.
    nk = K.kind == "nonsnp"
    K.loc[nk, "unit"] = row_units(rows_ns)
    UI = unit_info(K.loc[nk, "unit"].astype(np.int64).to_numpy())
    K.loc[nk, "n_members"] = UI.n_members.to_numpy()
    f = f"{GEAIN}/lfmm_nonsnp_gen9_bio6.csv"
    hdr = open(f).readline()
    out = subprocess.run(["awk", "-F,", f'NR>1 && $1=="Chr5" && $2+0>={LO - 50000} && $2+0<={HI + 50000}', f],
                         capture_output=True, text=True).stdout
    B = pd.read_csv(pd.io.common.StringIO(hdr + out), usecols=["unit", "B"])
    K = K.merge(B.assign(unit=B.unit.astype(float)), on="unit", how="left")
    ok = K.B.notna() & (K.n_members == 1)
    print(f"  stored LFMM B (tile_cr50) for {int(K.B.notna().sum())} non-SNP records "
          f"({K.loc[K.B.notna(), 'unit'].nunique()} units); slope vs B over the "
          f"{int(ok.sum())} single-record units: r = {K.loc[ok, 'coef'].corr(K.loc[ok, 'B']):.3f}")
    K.to_csv(f"{OUT}/coef_track_kmate.csv", index=False)

    ids = pd.read_csv(f"{P1}/lfmm_bio6_k16_results.csv", usecols=["snp_id"])
    c = ids.snp_id.str.split("_", expand=True)
    ids["chr"], ids["pos"] = c[0], c[1].astype(int)
    m = (ids["chr"] == "5") & ids.pos.between(LO, HI)
    want = {int(i) + 1: int(p) for i, p in zip(np.where(m)[0], ids.pos.values[m])}
    hdr = open(P1DP).readline().rstrip("\n").split(",")
    keep = [i for i, p in enumerate(hdr) if p in set(pools)]
    xp = meta.set_index("pool").loc[[hdr[i] for i in keep], "bio6"].to_numpy(float)
    got = {}
    with open(P1DP) as fh:
        fh.readline()
        for n, line in enumerate(fh, start=1):
            if n in want:
                v = np.fromstring(line, sep=",")[keep]
                xc = xp - xp.mean()
                got[want[n]] = float((v - v.mean()) @ xc / (xc @ xc))
            if n > max(want):
                break
    P = pd.DataFrame(dict(pos=list(got), coef=list(got.values())))
    Bp = pd.read_csv(f"{P1}/effect_sizes_simple_full_genome_bio6_k16.csv", usecols=["bio6"])
    P["B"] = Bp.iloc[[i for i in np.where(m)[0]], 0].values[:len(P)]
    P.to_csv(f"{OUT}/coef_track_phase1.csv", index=False)
    print(f"  phase-1 coef: {len(P)} SNPs; slope vs their LFMM B r = {P.coef.corr(P.B):.3f}")
