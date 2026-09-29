#!/usr/bin/env python
"""Event-level GEA: SV hits compared with the allele-level GEA, and promoter-deletion direction.

Region per event = the SAME classifier (screen_sig_blocks.classify, 1 kb promoter) on the event's
representative record (the member with most carriers). Hit call: Bonferroni 0.05 / units per
class x axis, union over the two event scans (pooled nonsnp, sv-only), as the GEA defines it.
Direction: z > 0 = the deletion rises more in warmer gardens (temperature axes; pc1 -> +bio1,
pc3 -> +bio10). Contrast = promoter minus intergenic deletion events, mean z over all tested
events, stratified by size x MAF, block-bootstrap CI (gea_direction.contrast_mean_z).
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
GEA = os.path.abspath(os.path.join(CONV, "..", ".."))
sys.path.insert(0, CONV); sys.path.insert(0, os.path.join(GEA, "genes", "dissection")); sys.path.insert(0, GEA)
import screen_sig_blocks as S, gea_region_by_axis as GR, gea_direction as GD   # noqa: E402
TEMP = GD.TEMP

E = pd.read_csv(f"{HERE}/sv_gen9_units.csv")
C, _ = S.classify(E[["chrom", "pos", "ref_len", "alt_len"]].copy())
E["tier"] = C.tier.values
E["grp"] = E.tier.map({"3_promoter": "promoter", "7_proximal_intergenic": "intergenic",
                       "8_gene_desert": "intergenic"}).fillna("other")
E["kind"] = np.where(E.dlen < 0, "deletion", "insertion")
E["size"] = E.dlen.abs()
E["stratum"] = (pd.cut(E["size"], GR.SIZE_BINS, labels=False).astype(str) + "|"
                + pd.cut(E.MAF, GR.MAF_BINS, labels=False, include_lowest=True).astype(str))
# the allele-level GEA could test an event only if one of its alleles cleared MAF 0.05
U = np.load(f"{GEA}/r1_sv_negative_selection/results/sv_adaptive/nonsnp_event_units.npz")
L = pd.read_csv(f"{CONV}/results/gea_hit_landscape.csv.gz", usecols=["chrom", "pos", "ref_len", "alt_len", "hit", "cls"])
idx = np.load(f"{GEA}/../../data/af_store/index_nonsnp.npz", allow_pickle=True) if False else None

hits = set(); zs = {}
rng = np.random.default_rng(1); rows = []
for ax in GR.AXES:
    ax_hit = pd.Series(False, index=E.index)
    for cls in ("nonsnp", "sv"):
        T = pd.read_csv(f"{HERE}/lfmm_event_{cls}_gen9_{ax}.csv.gz", usecols=["unit", "pval", "z", "is_sv"])
        sig = set(T.unit[T.pval < 0.05 / len(T)])
        ax_hit |= E.unit.isin(sig)
        if cls == "sv":
            E[f"z_{ax}"] = E.unit.map(T.set_index("unit").z)
    E[f"hit_{ax}"] = ax_hit
    d = E[E.kind == "deletion"].assign(z=E[f"z_{ax}"], block=E.block)
    h = d[d[f"hit_{ax}"]]
    e, lo, hi, p = GD.contrast_mean_z(d, "promoter", "intergenic", rng)
    rows.append(dict(axis=ax, temperature=ax in TEMP, del_events=len(d), del_hits=len(h),
                     hits_promoter=int((h.grp == "promoter").sum()),
                     warm_hits_promoter=(h[h.grp == "promoter"].z > 0).mean() if (h.grp == "promoter").any() else np.nan,
                     hits_intergenic=int((h.grp == "intergenic").sum()),
                     warm_hits_intergenic=(h[h.grp == "intergenic"].z > 0).mean() if (h.grp == "intergenic").any() else np.nan,
                     dz=e, dz_lo=lo, dz_hi=hi, dz_p=p))
R = pd.DataFrame(rows)
E["hit"] = E[[f"hit_{a}" for a in GR.AXES]].any(axis=1)
E.drop(columns=[c for c in E.columns if c.startswith("z_") or c.startswith("hit_")]).to_csv(f"{HERE}/event_sv_hits.csv", index=False)
R.to_csv(f"{HERE}/event_direction.csv", index=False)
pd.set_option("display.width", 220)
print(f"SV events tested {len(E):,}; SV event hits (any axis, either scan) {int(E.hit.sum())} "
      f"({int((E.hit & (E.kind == 'deletion')).sum())} deletions, {int((E.hit & (E.kind == 'insertion')).sum())} insertions)")
print(f"  hits among events with >1 allele: {int((E.hit & (E.n_alleles > 1)).sum())}")
print("\nSV DELETION EVENTS, promoter vs intergenic (z > 0 = deletion rises in warmer gardens)")
print(R.round(3).to_string(index=False))
