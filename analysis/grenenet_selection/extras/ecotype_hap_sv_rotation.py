"""Adjudication: does my ~1.4x winning-haplotype SV-enrichment survive a SPATIAL (genome-
rotation) null? Another agent's audit showed the JOINT/TEMPORAL signal collapses under a
rotation null (fold~0.94) because BOTH the selection signal and the SV-tagging track are
spatially autocorrelated, and a dispersed covariate-matched null is anti-conservative
([[feedback-spatial-null-colocalization]]). My exact-k + pericentro/TE/gene null is exactly
that dispersed kind. Test it on MY hap table + MY kinship-corrected fitness axis.

Rotation null: order haplotypes by genomic position (chrom,start,cluster); keep the FITNESS
selection fixed (its spatial structure intact); circularly shift the SV-tag track within each
chromosome by a random offset (its spatial structure intact); recompute mean(tag) over the
fixed selected set. Both tracks keep their autocorrelation; only their relative phase breaks.

Reuses results/.../gwas/hap_sv_table.csv. Output -> hap_sv_rotation.csv. Env: kmate.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

G = f"{lib.GEA}/ecotype_fitness/gwas"
NROT = 5000
CELLS = [("zone_union", "fit_zone", 0.005), ("zone_union", "fit_zone", 0.01),
         ("global", "fit_glob", 0.01)]


def rotation_null(tag, chrom, sel_mask, nrot=NROT, seed=0):
    """mean(tag[sel]) vs circularly-shifted-tag null (shift within each chrom, sel fixed)."""
    rng = np.random.default_rng(seed)
    obs = tag[sel_mask].mean()
    # per-chrom index ranges (array already sorted by chrom,start,cluster)
    order = np.arange(len(tag))
    chrom_ids = np.unique(chrom)
    idx_by_chrom = {c: order[chrom == c] for c in chrom_ids}
    null = np.empty(nrot)
    sel_idx = np.where(sel_mask)[0]
    for r in range(nrot):
        rt = tag.copy()
        for c, ix in idx_by_chrom.items():
            sh = rng.integers(1, len(ix)) if len(ix) > 1 else 0
            rt[ix] = np.roll(tag[ix], sh)
        null[r] = rt[sel_idx].mean()
    fold = obs / null.mean() if null.mean() > 0 else np.nan
    p = (1 + (null >= obs).sum()) / (nrot + 1)
    return float(obs), float(null.mean()), float(fold), float(p)


def main():
    H = pd.read_csv(f"{G}/hap_sv_table.csv")
    H = H.sort_values(["chrom", "start", "k"]).reset_index(drop=True)
    chrom = H.chrom.to_numpy()
    rows = []
    for axis, fcol, q in CELLS:
        fc = H[fcol].to_numpy()
        sel = fc >= np.quantile(fc, 1 - q)
        for mac, tname in ((12, "tag12"), (24, "tag24"), (46, "tag46")):
            tag = H[tname].to_numpy()
            obs, nm, fold, p = rotation_null(tag, chrom, sel)
            rows.append(dict(axis=axis, top_q=q, sv_mac=mac, n_sel=int(sel.sum()),
                             obs_r2=round(obs, 4), rot_null_r2=round(nm, 4),
                             fold=round(fold, 3), p_rot=round(p, 4)))
    out = pd.DataFrame(rows)
    out.to_csv(f"{G}/hap_sv_rotation.csv", index=False)
    json.dump(rows, open(f"{G}/hap_sv_rotation.json", "w"), indent=2)
    print("=== ROTATION (spatial) null vs my earlier dispersed exact-k null ===")
    print("(earlier exact-k fold for zone0.5%: ~1.41/1.44/1.46; global1%: ~1.30/1.27/1.25)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
