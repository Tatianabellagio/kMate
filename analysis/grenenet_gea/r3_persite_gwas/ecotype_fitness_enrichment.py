"""Tier 2 — SV enrichment among ecotype-fitness-GWAS hits, at the founder unit.

Reads the fitness-GWAS z (ecotype_fitness_gwas.py: LOCO-kinship Z and naive-OLS Znaive) and
asks: are SVs over-represented among fitness-associated variants BEYOND SNPs, matched on
founder MAC? Three things make this the honest version:

  1) FREQ-MATCHED null: each SV is compared to MAC-binned SNP controls (not raw fractions) —
     the ×2.07 artifact came from unmatched size/frequency. [[feedback-test-at-signal-unit]]
  2) DIRECTION-AGNOSTIC 'selected in >=1 climate zone' hit = union of {cold,mid,hot} zone hits
     (matches CONDITIONAL NEUTRALITY; not the antagonistic slope). w_global = generalist axis.
  3) SNP-INDEPENDENCE split (best_r2 from sv_snp_ld panel): tagged SV (r2>=0.8, rides a SNP
     haplotype) vs independent SV (r2<0.2 — the ONLY place 'SVs reveal what SNPs miss' holds).

Also reports the LOCO-vs-NAIVE contrast (naive = no kinship): how much apparent SV signal is
just ecotype sorting / clade structure.

Output -> analysis/grenenet_gea/ecotype_fitness/gwas/enrichment.{csv,json}
Env: kmate (numpy/scipy/pandas). Light — run after the GWAS sbatch.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

G = f"{lib.GEA}/ecotype_fitness/gwas"
LD = f"{lib.GEA}/sv_snp_ld"
TOP_Q = [0.02, 0.01, 0.005]       # hit = top-q of |z| (per phenotype)
NPERM = 10000


def load_sv_r2():
    """(chrom,pos,sv_size) -> best_r2 with a local SNP, from the panel SV-SNP LD store."""
    d = {}
    for c in range(1, 6):
        p = f"{LD}/sv_snp_ld_panel_Chr{c}.npz"
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        for ch, po, sz, r2 in zip(z["chrom"], z["pos"], z["sv_size"], z["best_r2"]):
            d[(str(ch), int(po), int(sz))] = min(float(r2), 1.0)
    return d


def matched_fold(is_hit, mac, sv_idx, ctrl_mask, nperm=NPERM, seed=0):
    """P(hit|SV) vs P(hit|MAC-matched control). Controls drawn from ctrl_mask only.

    Returns (obs, null_mean, fold, p_one_sided). Each SV's null draws come from controls in
    its own MAC bin, so enrichment is frequency-matched.
    """
    rng = np.random.default_rng(seed)
    bins = np.digitize(mac, np.quantile(mac, np.linspace(0, 1, 11)[1:-1]))
    ctrl_by_bin = {b: np.where(ctrl_mask & (bins == b))[0] for b in np.unique(bins)}
    sv_idx = np.asarray(sv_idx)
    obs = is_hit[sv_idx].mean()
    null = np.empty((len(sv_idx), nperm))
    for k, i in enumerate(sv_idx):
        pool = ctrl_by_bin.get(bins[i])
        if pool is None or len(pool) == 0:
            null[k] = 0.0
        else:
            null[k] = is_hit[pool][rng.integers(0, len(pool), nperm)]
    nm = null.mean(0)
    fold = obs / nm.mean() if nm.mean() > 0 else np.nan
    p = (1 + (nm >= obs).sum()) / (nperm + 1)
    return float(obs), float(nm.mean()), float(fold), float(p)


def main():
    z = np.load(f"{G}/gwas_z.npz", allow_pickle=True)
    chrom = z["chrom"].astype(str); pos = z["pos"]; rl = z["ref_len"]; al = z["alt_len"]
    mac = z["mac"].astype(float); vclass = z["vclass"]
    phenos = list(z["pheno_names"].astype(str))
    ZZ = {"LOCO": z["Z"], "naive": z["Znaive"]}
    is_sv = vclass == 2; is_snp = vclass == 0
    sv_size = np.abs(al - rl)

    # SNP-independence per SV
    r2map = load_sv_r2()
    best_r2 = np.array([r2map.get((chrom[i], int(pos[i]), int(sv_size[i])), np.nan)
                        for i in np.where(is_sv)[0]])
    sv_all = np.where(is_sv)[0]
    have = np.isfinite(best_r2)
    sv_tagged = sv_all[have & (best_r2 >= 0.8)]
    sv_indep = sv_all[have & (best_r2 < 0.2)]
    print(f"variants={len(pos):,}  SV={is_sv.sum():,}  SNP={is_snp.sum():,}")
    print(f"SV with r2: {have.sum():,}  tagged(>=0.8)={len(sv_tagged):,}  "
          f"independent(<0.2)={len(sv_indep):,}")

    # zone-union (conditional-neutrality) and generalist axes, per flavour
    axes = {}
    for fl in ("rel", "cen"):
        axes[f"zone_union_{fl}"] = [f"w_{zt}_{fl}" for zt in ("cold", "mid", "hot")]
        axes[f"global_{fl}"] = [f"w_global_{fl}"]

    rows = []
    for model, Z in ZZ.items():
        col = {p: Z[:, phenos.index(p)] for p in phenos}
        for axis, members in axes.items():
            for q in TOP_Q:
                # union hit: top-q of |z| in ANY member phenotype (direction-agnostic)
                hit = np.zeros(len(pos), bool)
                for p in members:
                    zc = np.abs(col[p])
                    thr = np.quantile(zc, 1 - q)
                    hit |= zc >= thr
                for label, idx in (("SV_all", sv_all), ("SV_tagged", sv_tagged),
                                   ("SV_indep", sv_indep)):
                    if len(idx) < 20:
                        continue
                    obs, nm, fold, pv = matched_fold(hit, mac, idx, is_snp)
                    rows.append(dict(model=model, axis=axis, top_q=q, sv_set=label,
                                     n=len(idx), obs=round(obs, 4), null=round(nm, 4),
                                     fold=round(fold, 3), p=round(pv, 4)))
    df = pd.DataFrame(rows)
    df.to_csv(f"{G}/enrichment.csv", index=False)

    # headline: SV_indep on zone_union under LOCO (the decisive, novel number)
    print("\n=== decisive: SNP-independent SV enrichment, zone-union (conditional neutrality) ===")
    key = df[(df.sv_set == "SV_indep") & (df.axis.str.startswith("zone_union"))]
    print(key.to_string(index=False))
    print("\n=== LOCO vs naive (all SV, zone_union rel, how much is clade sorting) ===")
    con = df[(df.sv_set == "SV_all") & (df.axis == "zone_union_rel")]
    print(con.to_string(index=False))
    json.dump({"rows": rows}, open(f"{G}/enrichment.json", "w"), indent=2)
    print(f"\nwrote {G}/enrichment.csv")


if __name__ == "__main__":
    main()
