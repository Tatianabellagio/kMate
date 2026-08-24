#!/usr/bin/env python
"""Build + execute notebooks/08_candidate_block.ipynb (run in the `basic` env).

Deep-dive on a single candidate LD block from the two-stage / WZA SV climate-GEA
(default 4_2781 = Chr4:15.47 Mb, inside AT4G31980 — the block that is a top hit in
BOTH the Δp and selection-coefficient statistics). Shows:
  - the block's SVs + their per-stat effects/p (table)
  - per-site allele-frequency TRAJECTORIES gen0(p0)->1->2->3, colored by climate
    (bio1) — the "rises-in-warm / falls-in-cold" readout
  - site-level Δp(gen3) vs bio1 scatter for the lead SV (the GEA at honest N=20)
  - all-block-SV Δp-vs-climate direction (heatmap)
  - gene-context map (AT4G31980 + neighbors, SV positions/sizes)
  - SNP-vs-SV independence: lead SV best r² with assembly-panel & short-read SNPs

Precomputes into a small cache so the notebook loads instantly.

  BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $BPY analysis/grenenet_selection/r2_gea_nonsnp/_build_candidate_nb.py            # default block 4_2781
"""
import os, glob, sys
import numpy as np
import pandas as pd
import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PROJ = "/global/scratch/users/tbellg/kmate"
GEA = lib.GEA; PM = f"{GEA}/common/results/pool_matrices"; STORE = lib.AF_STORE
NBDIR = f"{PROJ}/analysis/grenenet_selection/notebooks"
OUT = f"{NBDIR}/08_candidate_block.ipynb"
CACHE = f"{GEA}/r2_gea_nonsnp/results/gea/candidate_block_cache.npz"
CAND_BLOCK = "4_2781"
CHR_INT = int(CAND_BLOCK.split("_")[0])
CHROM = f"Chr{CHR_INT}"


def precompute():
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    size = np.abs(idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64))
    nc = np.asarray(np.load(sorted(glob.glob(f"{STORE}/nc_nonsnp/*.npy"))[0]))
    mask = (size > 50) & (nc >= 150); sv_idx = np.where(mask)[0]
    block = np.load(f"{GEA}/r2_gea_nonsnp/results/gea/twostage_blocks.npz", allow_pickle=True)["block"].astype(str)
    zd = np.load(f"{GEA}/r2_gea_nonsnp/results/gea/twostage_dp_bio1.npz", allow_pickle=True)
    zs = np.load(f"{GEA}/r2_gea_nonsnp/results/gea/twostage_scoef_bio1.npz", allow_pickle=True)
    sel = np.where(block == CAND_BLOCK)[0]
    nonsnp_idx = sv_idx[sel]
    p0 = zd["p0"][sel].astype(float)
    bm = dict(pos=zd["pos"][sel], ref_len=zd["ref_len"][sel], alt_len=zd["alt_len"][sel],
              sv_size=zd["sv_size"][sel], p0=p0,
              beta_dp=zd["beta"][sel], z_dp=zd["z_emp"][sel], p_dp=zd["p_perm"][sel],
              beta_sc=zs["beta"][sel], z_sc=zs["z_emp"][sel], p_sc=zs["p_perm"][sel])
    lead_local = int(np.nanargmax(np.abs(bm["z_sc"])))      # lead by |z_scoef|

    # per-gen per-site mean AF for the block SVs (gen0 = shared p0)
    site_means = {}; bio1 = {}; sites_by_gen = {}
    for g in (1, 2, 3):
        P = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy")[:, nonsnp_idx].astype(np.float32)
        m = pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv")
        sm = pd.DataFrame(P, columns=[f"sv{i}" for i in range(P.shape[1])])
        sm["site"] = m.site.to_numpy()
        gm = sm.groupby("site").mean()
        site_means[g] = gm
        sites_by_gen[g] = gm.index.to_numpy()
        bio1[g] = m.groupby("site").bio1.first()
    sites3 = sites_by_gen[3]
    bio1_3 = bio1[3].loc[sites3].to_numpy()

    nB = len(nonsnp_idx)
    # trajectory tensor [4 gens x nSites3 x nB]; gen0 = p0 broadcast
    traj = np.full((4, len(sites3), nB), np.nan, np.float32)
    traj[0] = p0[None, :]
    for gi, g in enumerate((1, 2, 3), start=1):
        gm = site_means[g]
        for si, s in enumerate(sites3):
            if s in gm.index:
                traj[gi, si] = gm.loc[s].to_numpy()
    dp3 = traj[3] - p0[None, :]                              # [nSites3 x nB] gen3 Δp

    # gene context: genes overlapping the locus window
    lo = int(bm["pos"].min()) - 12000; hi = int(bm["pos"].max()) + 12000
    genes = lib.load_genes()
    gc = genes[(genes.chrom == CHROM) & (genes.end >= lo) & (genes.start <= hi)]

    # SNP-vs-SV LD for the lead SV (assembly panel + short-read)
    lead_pos = int(bm["pos"][lead_local])
    ld = {}
    for tag, f in [("panel", f"{GEA}/r3_persite_gwas/results/sv_snp_ld/sv_snp_ld_panel_{CHROM}.npz"),
                   ("shortread", f"{GEA}/r3_persite_gwas/results/sv_snp_ld/sv_snp_ld_shortread_{CHROM}.npz")]:
        if os.path.exists(f):
            z = np.load(f, allow_pickle=True)
            j = np.argmin(np.abs(z["pos"] - lead_pos))
            if abs(int(z["pos"][j]) - lead_pos) <= 5:
                ld[tag] = (float(z["best_r2"][j]), int(z["best_snp_pos"][j]),
                           int(z["n_snp_window"][j]))

    np.savez(CACHE, block=CAND_BLOCK, chrom=CHROM, lead_local=lead_local,
             sites3=sites3, bio1_3=bio1_3, traj=traj, dp3=dp3, p0=p0,
             gene_start=gc.start.to_numpy(), gene_end=gc.end.to_numpy(),
             gene_id=gc.gene.to_numpy().astype("U16"),
             ld_panel_r2=ld.get("panel", (np.nan, -1, 0))[0],
             ld_panel_snp=ld.get("panel", (np.nan, -1, 0))[1],
             ld_short_r2=ld.get("shortread", (np.nan, -1, 0))[0],
             ld_short_snp=ld.get("shortread", (np.nan, -1, 0))[1],
             **{f"bm_{k}": v for k, v in bm.items()})
    print(f"precompute: block {CAND_BLOCK} ({CHROM}), {nB} SVs, lead@{lead_pos}, "
          f"{len(sites3)} gen3 sites, panel r2={ld.get('panel',('NA',))[0]} -> {CACHE}")


md_intro = fr"""# Candidate block {CAND_BLOCK} — Chr4:15.47 Mb (AT4G31980)

The one block that is a **top hit in BOTH** statistics of the two-stage SV climate-GEA
(endpoint Δp **and** the gen0→3 selection-coefficient), and survives the WZA block test
(scoef Z_pVal≈2e-6, up-in-warm). A tight ~3.3-kb cluster of SVs **inside the gene
AT4G31980**. Here we look at *what it actually does*: per-site frequency trajectories,
the climate gradient, gene context, and whether SNPs already capture it.
"""

code_load = r"""
import numpy as np, pandas as pd, sys
import matplotlib.pyplot as plt
from matplotlib import cm
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection"); import lib
C = np.load("/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/results/gea/candidate_block_cache.npz", allow_pickle=True)
bm = pd.DataFrame({k[3:]: C[k] for k in C.files if k.startswith("bm_")})
lead = int(C["lead_local"]); sites3 = C["sites3"]; bio1 = C["bio1_3"]
traj = C["traj"]; dp3 = C["dp3"]; p0 = C["p0"]
print(f"block {str(C['block'])} on {str(C['chrom'])}: {len(bm)} SVs, lead row {lead} "
      f"(pos {int(bm.pos[lead])}, {int(bm.sv_size[lead])}bp, p0={bm.p0[lead]:.3f})")
print(f"SNP-LD of lead SV: assembly-panel best r²={float(C['ld_panel_r2']):.2f}, "
      f"short-read best r²={float(C['ld_short_r2']):.2f}")
bm.assign(kind=np.where(bm.alt_len<bm.ref_len,"DEL",np.where(bm.alt_len>bm.ref_len,"INS","?")))\
  .sort_values("z_sc",key=np.abs,ascending=False)\
  [["pos","sv_size","kind","p0","beta_dp","z_dp","p_dp","beta_sc","z_sc","p_sc"]].head(12)
"""

code_traj = r"""
# Per-site allele-frequency trajectory of the LEAD SV: gen0(=shared p0) -> 1 -> 2 -> 3,
# one line per site, colored by site temperature (bio1). Up-in-warm = red lines rise.
gens=[0,1,2,3]; y=traj[:,:,lead]
norm=plt.Normalize(bio1.min(), bio1.max()); cmap=cm.coolwarm
fig, ax = plt.subplots(1,2, figsize=(13,4.6), gridspec_kw={"width_ratios":[1.4,1]})
for si in range(len(sites3)):
    ax[0].plot(gens, y[:,si], "-o", ms=3, lw=1.2, color=cmap(norm(bio1[si])), alpha=.85)
ax[0].axhline(p0[lead], color="k", ls=":", lw=1, label=f"founding p0={p0[lead]:.3f}")
ax[0].set_xlabel("generation"); ax[0].set_ylabel("site-mean alt-allele freq")
ax[0].set_xticks(gens); ax[0].set_title(f"Lead SV trajectory by site (n={len(sites3)})")
ax[0].legend(fontsize=8)
sm=cm.ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
fig.colorbar(sm, ax=ax[0], label="site bio1 (°C)")
# the headline GEA: site Δp(gen3) vs bio1
d=dp3[:,lead]; ok=np.isfinite(d)
ax[1].scatter(bio1[ok], d[ok], c=bio1[ok], cmap=cmap, s=70, edgecolor="k", lw=.4)
b,a=np.polyfit(bio1[ok], d[ok], 1); xs=np.array([bio1.min(),bio1.max()])
ax[1].plot(xs, a+b*xs, "k-", lw=1)
r=np.corrcoef(bio1[ok], d[ok])[0,1]
ax[1].axhline(0,color="grey",lw=.6); ax[1].set_xlabel("site bio1 (°C)")
ax[1].set_ylabel("Δp(gen3) = site-mean − p0"); ax[1].set_title(f"Lead SV Δp vs climate  (r={r:.2f}, slope={b:.4f}/°C)")
plt.tight_layout(); plt.show()
"""

code_block = r"""
# Whole block: gen3 Δp per site (rows) x SV (cols), sites ordered by bio1. Coherent
# vertical bands = the block's SVs move together with climate.
oo=np.argsort(bio1)
fig, ax = plt.subplots(figsize=(11,5))
im=ax.imshow(dp3[oo], aspect="auto", cmap="coolwarm",
             vmin=-np.nanmax(np.abs(dp3)), vmax=np.nanmax(np.abs(dp3)))
ax.set_yticks(range(len(sites3))); ax.set_yticklabels([f"s{int(s)} ({bio1[i]:.0f}°)" for i,s in zip(oo,sites3[oo])], fontsize=6)
ax.set_xticks(range(len(bm))); ax.set_xticklabels([f"{int(p/1e3)}k\n{int(sz)}bp" for p,sz in zip(bm.pos,bm.sv_size)], fontsize=5, rotation=90)
ax.set_xlabel("block SVs (pos / size)"); ax.set_ylabel("site (cold→warm)")
ax.set_title(f"block {str(C['block'])}: gen3 Δp (site × SV); sites sorted by bio1")
fig.colorbar(im, label="Δp"); plt.tight_layout(); plt.show()
"""

code_gene = r"""
# Gene context: AT4G31980 + neighbours, with the block's SVs drawn at their positions.
gs=C["gene_start"]; ge=C["gene_end"]; gid=C["gene_id"].astype(str)
fig, ax = plt.subplots(figsize=(12,2.8))
for s,e,g in zip(gs,ge,gid):
    ax.add_patch(plt.Rectangle((s,0.25),e-s,0.5, color="#9ab", ec="k", lw=.5))
    ax.text((s+e)/2, 0.8, g, fontsize=7, ha="center", rotation=0,
            color=("crimson" if g=="AT4G31980" else "k"))
for _,r in bm.iterrows():
    col="crimson" if r.alt_len<r.ref_len else "navy"
    ax.plot([r.pos],[0.05],"v",color=col,ms=6)
ax.set_ylim(0,1); ax.set_yticks([]); ax.set_xlabel(f"{str(C['chrom'])} position (bp)")
ax.set_title("Gene context (▼ SVs: red=DEL, blue=INS) — block sits inside AT4G31980")
plt.tight_layout(); plt.show()
print("AT4G31980 spans 15,464,905–15,469,204; the 22-SV block (15,466,240–15,469,500) is intragenic.")
"""

md_close = r"""## Reading this block

- The block sits **inside AT4G31980**; the lead signals are ~60–110 bp **deletions**
  (p0≈0.12) that **rise in warm gardens** across generations (B selection-coefficient
  z≈3.2, up-in-warm), plus a multi-allelic insertion site at the 3′ end.
- The trajectory + Δp-vs-bio1 panels are the honest readout at **N=20 sites** — judge the
  gradient by eye, not by a per-SV p-value.
- **SNP-independence:** compare the lead SV's best r² with assembly-panel vs short-read
  SNPs — if short-read r² is low, this is an SV a SNP-only GEA would (partly) miss.
- Caveats: still only ~20 sites of climate df; confirm with **cross-garden parallelism**
  and check AT4G31980's function before claiming adaptation.
"""

cells = [
    nbf.v4.new_markdown_cell(md_intro),
    nbf.v4.new_code_cell(code_load.strip()),
    nbf.v4.new_markdown_cell("## Per-site frequency trajectories + the climate gradient"),
    nbf.v4.new_code_cell(code_traj.strip()),
    nbf.v4.new_markdown_cell("## Whole-block Δp across sites (coherence)"),
    nbf.v4.new_code_cell(code_block.strip()),
    nbf.v4.new_markdown_cell("## Gene context"),
    nbf.v4.new_code_cell(code_gene.strip()),
    nbf.v4.new_markdown_cell(md_close),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec":
        {"name": "python3", "display_name": "Python 3"}})

if __name__ == "__main__":
    precompute()
    os.makedirs(NBDIR, exist_ok=True)
    ep = ExecutePreprocessor(timeout=900, kernel_name="python3", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": NBDIR}})
    with open(OUT, "w") as f:
        nbf.write(nb, f)
    print("wrote + executed", OUT)
