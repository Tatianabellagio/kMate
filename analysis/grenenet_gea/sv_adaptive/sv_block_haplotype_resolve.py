#!/usr/bin/env python
"""BLOCK vs HAPLOBLOCK: is each top-JOINT SV on a SELECTED haplotype (its own), or on a
NON-selected haplotype that merely shares the coarse clq0.9 block with the selected one?

My r^2-vs-lead audit conflated "not on the selected haplotype" with "passenger". The clq0.9
block is a coarse LD unit (can span kb, hold several hap-clusters). The founder-GWAS tests
EVERY hap-cluster in a block and the block's JOINT = the MAX cluster's chi2. So the right test:
for each SV, find the hap-cluster it actually TAGS (max founder-r^2), and read THAT cluster's
own GWAS chi2:
  - best-match cluster chi2 HIGH  -> the SV's own haplotype is itself selected  (NOT a passenger)
  - best-match cluster chi2 LOW   -> the SV sits on a non-selected sub-haplotype (passenger,
                                     the block scored on a DIFFERENT haplotype)
Uses only founder-r^2 + the GWAS per-cluster chi2 (raw['chi2_joint'] is per hap-cluster marker)
-> EMMAX-free, so it doesn't depend on my failed EMMAX self-check. Env: kmate.
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from founder_genotype import build_genotype

PANEL = "panel/arch3"; CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
MAC_MIN, MAC_GRM, S_DF = 3, 12, 30

raw = lib.multisite_gwas_raw("clq90_pc1")
chi2 = raw["chi2_joint"]; unit = raw["unit"]
# genome-wide chi2 reference points
null_med = stats.chi2.ppf(0.5, S_DF)                        # ~29.3 (a null cluster)
q = lib.bh(stats.chi2.sf(chi2, S_DF))
fdr_chi2 = chi2[q < 0.05].min() if (q < 0.05).any() else np.inf   # chi2 needed for genome-wide FDR<.05
top1_chi2 = np.quantile(chi2, 0.99)
print(f"reference chi2: null median {null_med:.1f} | genome-wide top-1% {top1_chi2:.1f} | FDR<.05 needs >= {fdr_chi2:.1f}")

os.environ["MEMB_TAG"] = "clq90"
G, founders, reg = build_genotype("clq90"); nF = len(founders)
cnt = G.sum(0); blk = reg.block.to_numpy()
order = np.lexsort((-cnt, blk)); sbk = blk[order]
is_ref = np.zeros(len(cnt), bool); is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
Gp = G[:, poly].astype(np.float64); regp = reg[poly].reset_index(drop=True)
regp_unit = (regp.chrom + ":" + regp.start.astype(str) + "-" + regp.end.astype(str)).to_numpy()
assert np.array_equal(regp_unit, unit), "marker order drift"
rows_by_unit = {u: np.where(regp_unit == u)[0] for u in np.unique(unit)}

up = pd.read_csv(f"{lib.GEA}/hapfreq/multisite_founder_gwas_clq90_pc1.csv").groupby("unit", as_index=False).p_joint.min()
L = pd.read_csv(f"{lib.GEA}/sv_adaptive/sv_landscape_clq0.9.csv").rename(columns={"block_id": "unit"})
up = up.merge(L[["unit", "chrom", "start_pos", "end_pos", "has_sv", "n_kept"]], on="unit")
up = up[up.n_kept >= 2]
ntop = int(round(0.005 * len(up)))
topsv = up.nsmallest(ntop, "p_joint"); topsv = topsv[topsv.has_sv == 1].reset_index(drop=True)
units = set(topsv.unit)
print(f"top-0.5% JOINT SV blocks: {len(units)}")

rows = []
for ch in CHROMS:
    cl = ch.lower(); meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
    pos = meta["pos"].astype(np.int64); dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
    vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
    vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
    na = np.asarray(vp.sum(0)).ravel(); ncl = np.asarray(vc.sum(0)).ravel()
    svmask = (dl > 50) & (na >= 12) & (na <= nF - 12) & (ncl / nF >= 0.9)
    blk = topsv[topsv.chrom == ch]
    for j in np.where(svmask)[0]:
        hit = blk[(blk.start_pos <= pos[j]) & (blk.end_pos >= pos[j])]
        if not len(hit):
            continue
        u = hit.iloc[0].unit; car = vp[:, j].toarray().ravel().astype(float)
        if car.sum() < 2 or car.sum() > nF - 2:
            continue
        cl_rows = rows_by_unit[u]
        r2 = np.array([np.corrcoef(car, Gp[:, r])[0, 1] ** 2 for r in cl_rows])
        cchi = chi2[cl_rows]
        best = int(np.argmax(r2)); lead = int(np.argmax(cchi))
        span = int(hit.iloc[0].end_pos - hit.iloc[0].start_pos)
        rows.append(dict(unit=u, sv_pos=int(pos[j]), block_span_bp=span, n_clusters=len(cl_rows),
                         lead_chi2=round(float(cchi[lead]), 1),
                         max_r2_anyclust=round(float(r2[best]), 2),
                         bestmatch_chi2=round(float(cchi[best]), 1),
                         bestmatch_is_lead=bool(best == lead),
                         r2_with_lead=round(float(r2[lead]), 2)))
df = pd.DataFrame(rows).sort_values("bestmatch_chi2", ascending=False)
pd.set_option("display.width", 190, "display.max_rows", 60)
print(df.to_string(index=False))
df.to_csv(f"{lib.GEA}/sv_adaptive/sv_block_haplotype_resolve.csv", index=False)

# ---- the decisive summary ----
sel = df.bestmatch_chi2 >= top1_chi2                        # SV's own best haplotype is itself top-1%-selected
tags = df.max_r2_anyclust >= 0.5                            # SV cleanly tags SOME cluster
print(f"\n=== does each SV's OWN best-matching haplotype carry selection? (block=coarse LD unit) ===")
print(f"  block spans: median {df.block_span_bp.median():.0f} bp (min {df.block_span_bp.min()}, max {df.block_span_bp.max()}); "
      f"clusters/block median {df.n_clusters.median():.0f}")
print(f"  SVs that cleanly tag a hap-cluster (max r2>=0.5): {int(tags.sum())}/{len(df)}")
print(f"  of those, best-matching cluster is itself top-1% selected (chi2>={top1_chi2:.0f}): "
      f"{int((sel & tags).sum())}/{int(tags.sum())}")
print(f"  best-matching cluster IS the block's lead (selected) cluster: {int(df.bestmatch_is_lead.sum())}/{len(df)}")
print(f"  SVs whose best haplotype is NON-selected (chi2 near null ~{null_med:.0f}, so the block scored on a "
      f"DIFFERENT haplotype): {int((~sel & tags).sum())}/{len(df)}")
print(f"\n  interpretation: 'best-matching cluster IS selected' => the SV's OWN haplotype is under selection "
      f"(NOT a mere passenger); 'best-matching cluster near null' => passenger on a co-block non-selected haplotype.")
