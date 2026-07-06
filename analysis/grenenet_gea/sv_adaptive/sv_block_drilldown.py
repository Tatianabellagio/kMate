#!/usr/bin/env python3
"""Drill into the SV-bearing blocks that drive the founder-GWAS JOINT enrichment (common SVs):
WHERE are they actually selected -- one garden, a consistent set, or up-here/down-there
heterogeneous (as the JOINT-not-GLOBAL result predicts)? And what genes/SVs are they?

For each top-JOINT SV-bearing block (common-SV landscape, MAF>=0.05) we take its lead-JOINT
hap-cluster marker's per-garden z vector (kinship-corrected founder-selection association at
each of the 30 gardens; z>0 = founders carrying the alt haplotype rose at that garden). Then:
  - how many gardens is it 'selected' at (|z|>2.5), and are the strong gardens same-sign
    (consistent) or mixed-sign (heterogeneous / local)?
  - which gardens, and their bio1 -- is there any climate pattern in the strong gardens?
  - what gene(s) does the block span, and what SV records (pos, size) does it contain?
Compare the SV top blocks to size-matched non-SV top blocks (are SV blocks more multi-garden /
more heterogeneous?). Run in `kmate` env (panel npz). Deterministic.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import lib

OUT = Path("results/grenenet_gea/sv_adaptive")
PANEL = Path("panel/arch3")
ZTHR = 2.5                       # per-garden |z| for "selected at that garden" (per-site lambda~1)
TOPFRAC = 0.01                   # define the enriched tail
MAC_MIN, CALLED_MIN, SV_BP = 12, 0.9, 50

raw = lib.multisite_gwas_raw("clq90_pc1")
Z, sites, bio1 = raw["Z"], raw["sites"], raw["bio1"]
chi2, unit = raw["chi2_joint"], raw["unit"]
# lead-JOINT marker per block -> its 30-garden z row
lead = pd.DataFrame({"unit": unit, "chi2": chi2, "row": np.arange(len(unit))})
lead = lead.loc[lead.groupby("unit")["chi2"].idxmax()]
u2row = dict(zip(lead.unit, lead.row))
u2chi = dict(zip(lead.unit, lead.chi2))

L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
L = L[L.block_id.isin(u2row)].copy()
L["chi2_joint"] = L.block_id.map(u2chi)
L["rank"] = L.chi2_joint.rank(ascending=False)
ntop = int(TOPFRAC * len(L))
L["top"] = L["rank"] <= ntop
genes = lib.load_genes()

def gene_span(chrom, s, e):
    g = genes[(genes.chrom == chrom) & (genes.start <= e) & (genes.end >= s)]
    return ";".join(dict.fromkeys(g.name.fillna(g.gene))) if len(g) else "-"

def profile(bid):
    z = Z[u2row[bid]]
    strong = np.abs(z) > ZTHR
    pos = [(int(sites[i]), float(z[i]), float(bio1[i])) for i in np.where(strong & (z > 0))[0]]
    neg = [(int(sites[i]), float(z[i]), float(bio1[i])) for i in np.where(strong & (z < 0))[0]]
    return z, pos, neg

# ---- summary: SV top blocks vs matched non-SV top blocks ----
def summarize(df, tag):
    ns, mixed, mono = [], 0, 0
    for bid in df.block_id:
        z, pos, neg = profile(bid)
        k = len(pos) + len(neg); ns.append(k)
        if pos and neg: mixed += 1
        if k == 1: mono += 1
    ns = np.array(ns)
    print(f"  {tag} (n={len(df)}): #gardens |z|>{ZTHR}: median {np.median(ns):.0f} mean {ns.mean():.2f} "
          f"(0:{(ns==0).sum()} 1:{(ns==1).sum()} 2:{(ns==2).sum()} 3+:{(ns>=3).sum()}) | "
          f"mixed-direction {mixed}/{len(df)} ({100*mixed/len(df):.0f}%) | mono-garden {mono}/{len(df)}")
    return ns

n_sv_top = int((L.top & (L.has_sv == 1)).sum())
print(f"{len(L):,} blocks | SV-bearing {int(L.has_sv.sum()):,} | top-{TOPFRAC:.0%} JOINT = {ntop} blocks "
      f"({n_sv_top} of them SV-bearing, vs ~{ntop*L.has_sv.mean():.0f} expected)")
sv_top = L[L.top & (L.has_sv == 1)].sort_values("chi2_joint", ascending=False)
snp_top = L[L.top & (L.has_sv == 0)]
print("\n=== per-garden selection breadth/direction (top-JOINT blocks) ===")
summarize(sv_top, "SV-bearing top")
# size-match the non-SV comparison to sv_top's n_kept distribution (coarse)
EDGES = [2,3,4,5,6,7,8,10,12,15,20,25,30,40,50,70,100,150,250,10**9]
sv_top["kb"] = pd.cut(sv_top.n_kept, EDGES, right=False, labels=False)
snp_top2 = snp_top.copy(); snp_top2["kb"] = pd.cut(snp_top2.n_kept, EDGES, right=False, labels=False)
rng = np.random.default_rng(0)
matched = pd.concat([snp_top2[snp_top2.kb == b].sample(min((sv_top.kb == b).sum() * 3, (snp_top2.kb == b).sum()),
              random_state=0) for b in sv_top.kb.unique() if (snp_top2.kb == b).any()])
summarize(matched, "non-SV top (size-matched)")

# ---- per-block detail for the top SV blocks ----
print(f"\n=== top {min(25,len(sv_top))} SV-bearing JOINT blocks: where selected + gene ===")
rows = []
for _, r in sv_top.head(25).iterrows():
    z, pos, neg = profile(r.block_id)
    gname = gene_span(r.chrom, r.start_pos, r.end_pos)
    ps = ",".join(f"{s}({zz:+.1f})" for s, zz, _ in sorted(pos, key=lambda x: -abs(x[1]))[:4])
    ng = ",".join(f"{s}({zz:+.1f})" for s, zz, _ in sorted(neg, key=lambda x: -abs(x[1]))[:4])
    rows.append(dict(block=r.block_id, gene=gname[:32], n_sv=int(r.n_sv),
                     chi2=round(r.chi2_joint,1), up_gardens=ps or "-", down_gardens=ng or "-"))
det = pd.DataFrame(rows)
pd.set_option("display.width", 200, "display.max_colwidth", 40, "display.max_rows", 60)
print(det.to_string(index=False))
det.to_csv(OUT / "sv_block_drilldown.csv", index=False)

# ---- actual SV records (pos,size,gene) inside the top SV blocks ----
print("\n=== SV records inside the top SV blocks (MAF>=0.05) ===")
chroms = sorted(sv_top.head(25).chrom.unique())
svrows = []
for ch in chroms:
    m = np.load(PANEL / ch.lower() / f"var_pa_231_arch3_{ch.lower()}.meta.npz", allow_pickle=True)
    pos = m["pos"].astype(np.int64); dl = m["alt_len"].astype(np.int64) - m["ref_len"].astype(np.int64)
    vp = sp.load_npz(PANEL / ch.lower() / f"var_pa_231_arch3_{ch.lower()}.var_pa.npz")
    vc = sp.load_npz(PANEL / ch.lower() / f"var_pa_231_arch3_{ch.lower()}.var_called.npz")
    na = np.asarray(vp.sum(0)).ravel().astype(float); nc = np.asarray(vc.sum(0)).ravel().astype(float)
    carr = np.divide(na, nc, out=np.zeros_like(na), where=nc > 0); maf = np.minimum(carr, 1 - carr)
    keep = (maf >= 0.05) & (nc / vp.shape[0] >= CALLED_MIN) & (np.abs(dl) > SV_BP)
    kp, kd, kmaf = pos[keep], dl[keep], maf[keep]
    for _, r in sv_top[sv_top.chrom == ch].head(25).iterrows():
        sel = (kp >= r.start_pos) & (kp <= r.end_pos)
        for p, d, mf in zip(kp[sel], kd[sel], kmaf[sel]):
            svrows.append(dict(block=r.block_id, chrom=ch, sv_pos=int(p), sv_dlen=int(d),
                               maf=round(float(mf),2), gene=gene_span(ch, int(p), int(p)+abs(int(d)))[:32]))
sv_recs = pd.DataFrame(svrows).sort_values(["block","sv_pos"])
print(sv_recs.to_string(index=False) if len(sv_recs) else "  (none mapped)")
sv_recs.to_csv(OUT / "sv_block_drilldown_svrecords.csv", index=False)
print(f"\n[done] -> {OUT}/sv_block_drilldown.csv + _svrecords.csv")
