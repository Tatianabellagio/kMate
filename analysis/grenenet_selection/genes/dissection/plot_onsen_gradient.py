#!/usr/bin/env python
"""ONSEN locus (Col-0 copy AT1TE24850, Chr1:7,714,707-7,722,547): garden gradient + founder split.

The three SV-event GEA hits at this locus are deletions of (parts of) the Col-0 ONSEN copy, all
significant on pc3 only (the climate PC oriented to +bio10, warmest-quarter temperature):
  unit 536725  116 bp internal deletion   (af_store row 105658, single allele)  -- the lead
  unit 536720  444 bp 5'-end deletion     (105635, single allele)
  unit 536726  5,192 bp deletion          (105681 + 105682 + 105683: three alleles of one event)

Panel a-c: per garden, change in the deletion EVENT frequency (gen-9 pool frequency, i.e. each
plot's last sampled generation, minus the founding seed-mix frequency; event frequency = sum of its
member alleles' frequencies, as the GEA's event units define it), garden mean over plots, against
garden pc3, bio10 and bio4. Raw garden means -- no population-structure correction (the GEA's LFMM
uses K = 16 latent factors). Spearman rho over the 31 gardens is printed.

The founder split prints the home climate of founders that carry the lead deletion versus those
that carry the intact Col-0-like allele (founders not called at the lead are left out).

Data are the same kMate matrices the GEA and plot_variant_garden_grid.py read
(pool_matrices/pool_gen9_nonsnp_af.npy, af_store/p0_nonsnp.npy); pc1-3 are recomputed exactly as
r2_gea_nonsnp/phase1_replication/build_class_matrices.py does (SVD of the z-scored bio1-19 over
the gen-9 pools; PC1 oriented to +bio1, PC2+ so the dominant loading is positive).

Output: results/plots/ONSEN_AT1TE24850_gradient.{png,pdf}. No titles (project convention).
env: kmate or basic. Compute node.
"""
import os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib

OUT = f"{HERE}/results/plots"
PM = f"{lib.GEA}/common/results/pool_matrices"
R1 = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
UNITS = {536725: "116 bp internal deletion (lead)", 536720: "444 bp 5'-end deletion",
         536726: "5,192 bp deletion (3-allele event)"}
COL = {536725: "#B04A00", 536720: "#E08A3C", 536726: "#8C6A4A"}


def pick_cols(path, cols, block=8):
    """X[:, cols] of a big [pools x records] .npy, read in row blocks (mmap crawls on Lustre)."""
    with open(path, "rb") as f:
        v = np.lib.format.read_magic(f)
        rd = {(1, 0): np.lib.format.read_array_header_1_0, (2, 0): np.lib.format.read_array_header_2_0}[v]
        (nrow, ncol), fortran, dt = rd(f); assert not fortran
        out = np.empty((nrow, cols.size), np.float32); off = f.tell()
        for i in range(0, nrow, block):
            nb = min(block, nrow - i); f.seek(off + i * ncol * dt.itemsize)
            out[i:i + nb] = np.fromfile(f, dtype=dt, count=nb * ncol).reshape(nb, ncol)[:, cols]
    return out


def climate_pcs(meta):
    m = meta[[f"bio{i}" for i in range(1, 20)]].to_numpy(float)
    Z = (m - m.mean(0)) / m.std(0, ddof=0)
    U, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    for k in range(3):
        sc = U[:, k] * S[k]
        flip = (np.corrcoef(sc, m[:, 0])[0, 1] < 0) if k == 0 else (Vt[k][np.argmax(np.abs(Vt[k]))] < 0)
        sc = -sc if flip else sc
        meta[f"pc{k + 1}"] = (sc - sc.mean()) / sc.std(ddof=0)
    return meta


def main():
    EU = np.load(f"{R1}/nonsnp_event_units.npz")
    mem = {u: EU["indices"][EU["indptr"][u]:EU["indptr"][u + 1]] for u in UNITS}
    cols = np.unique(np.concatenate(list(mem.values())))
    X = pick_cols(f"{PM}/pool_gen9_nonsnp_af.npy", cols)
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")[cols]
    meta = climate_pcs(pd.read_csv(f"{PM}/pool_gen9_nonsnp.meta.csv"))
    ci = {c: i for i, c in enumerate(cols)}
    G = []
    for u in UNITS:
        j = [ci[c] for c in mem[u]]
        pe = X[:, j].sum(1); pe[np.isnan(X[:, j]).any(1)] = np.nan; pe = np.minimum(pe, 1)
        e0 = min(float(p0[j].sum()), 1.0)
        g = meta.assign(dp=pe - e0).groupby("site").agg(
            dp=("dp", "mean"), **{c: (c, "first") for c in ["bio1", "bio4", "bio7", "bio10", "pc3"]})
        G.append(g.reset_index().assign(unit=u, p0=e0))
        print(f"{UNITS[u]}: {len(j)} allele(s), founding frequency {e0:.3f}, "
              f"rose in {int((g.dp > 0).sum())}/{len(g)} gardens (mean change {g.dp.mean():+.3f})")
    G = pd.concat(G, ignore_index=True)

    print("\nSpearman rho (p), garden mean change in deletion frequency vs garden climate, 31 gardens")
    for u, lab in UNITS.items():
        g = G[G.unit == u]
        print(f"  {lab:36s}" + "  ".join(
            f"{v} {stats.spearmanr(g[v], g.dp)[0]:+.2f} ({stats.spearmanr(g[v], g.dp)[1]:.3f})"
            for v in ["pc3", "bio10", "bio4", "bio7", "bio1"]))

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), gridspec_kw=dict(wspace=0.32))
    XL = {"pc3": "garden pc3 (points to warmest-quarter temperature)",
          "bio10": "garden bio10, warmest-quarter temperature (°C)",
          "bio4": "garden bio4, temperature seasonality"}
    for a, (v, L) in zip(axs, zip(XL, "abc")):
        for u, lab in UNITS.items():
            h = G[G.unit == u]
            a.scatter(h[v], h.dp, s=22, color=COL[u], edgecolor="0.3", lw=0.3, label=lab, zorder=3)
            b = np.polyfit(h[v], h.dp, 1); xs = np.linspace(h[v].min(), h[v].max(), 10)
            a.plot(xs, np.polyval(b, xs), color=COL[u], lw=1.4)
        a.axhline(0, color="0.5", lw=0.8, ls=(0, (4, 2)))
        a.set_xlabel(XL[v], fontsize=8.5)
        a.set_ylabel("change in deletion frequency\n(gen 9 − founding mix), garden mean", fontsize=8)
        a.grid(color="0.9", lw=0.6); a.spines[["top", "right"]].set_visible(False)
        a.annotate(L, xy=(0.01, 1.03), xycoords="axes fraction", fontsize=10, fontweight="bold", color="0.3")
    axs[0].legend(fontsize=6.8, frameon=False, loc="lower left")
    for e in ("png", "pdf"):
        fig.savefig(f"{OUT}/ONSEN_AT1TE24850_gradient.{e}", dpi=170, bbox_inches="tight")
    print(f"\nwrote {OUT}/ONSEN_AT1TE24850_gradient.png/.pdf")

    # founders: carriers of the lead deletion vs carriers of the intact (Col-0-like) allele
    base = f"{lib.PROJ}/panel/arch3/chr1/var_pa_231_arch3_chr1"
    m = np.load(f"{base}.meta.npz", allow_pickle=True)
    rl, al = m["ref_len"].astype(int), m["alt_len"].astype(int)
    panel = np.where(~((rl == 1) & (al == 1)))[0]
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    off = int(np.flatnonzero(idx["chrom"].astype(str) == "Chr1")[0])
    rec = int(panel[105658 - off]); assert (rl[rec], al[rec], int(m["pos"][rec])) == (117, 1, 7716474)
    pa = sp.load_npz(f"{base}.var_pa.npz").tocsc()[:, rec].toarray().ravel() > 0
    called = sp.load_npz(f"{base}.var_called.npz").tocsc()[:, rec].toarray().ravel() > 0
    F = pd.DataFrame(dict(founder=[str(x) for x in m["founders"]], deletion=pa, called=called))
    FO = pd.read_csv(f"{R1}/sv_origin_partition_founders.csv", dtype={"founder": str})
    F = F.merge(FO[["founder", "bio1", "group"]], on="founder", how="left")
    c = F[F.called]
    print(f"\nFOUNDERS at the lead (116 bp deletion): {int(c.deletion.sum())} carry the deletion, "
          f"{int((~c.deletion).sum())} the intact allele, {int((~F.called).sum())} not called")
    for lab, d in (("deletion", c[c.deletion]), ("intact ONSEN copy", c[~c.deletion])):
        print(f"  {lab:18s} n={len(d):3d}  home bio1 median {d.bio1.median():.1f} °C  "
              f"origin groups {d.group.value_counts().to_dict()}")
    print(f"  Mann-Whitney home bio1, deletion vs intact: p = "
          f"{stats.mannwhitneyu(c[c.deletion].bio1.dropna(), c[~c.deletion].bio1.dropna()).pvalue:.3g}")


if __name__ == "__main__":
    main()
