"""Tier 2 (the question the user actually wants): are the haploblocks that are FITNESS-RELATED
enriched for SVs?

NOT a causal-SV claim (SVs share a clq0.9 block with SNPs by construction, so per-SV
attribution is moot — user's point). The honest, useful question is block-level CO-OCCURRENCE /
ENRICHMENT: do the r2>=0.9 haploblocks whose SNPs associate with ecotype fitness carry SVs more
than size-matched non-fitness blocks?

Design (non-circular + size-controlled, the two things that made the retracted x2.07 an
artifact):
  - unit = clq0.9 block (blocks_mcf90/chrN_clq0.9_blocks_clq0.9.tsv), variants assigned by pos.
  - block FITNESS score = max |z| among the block's SNPs ONLY (>=3 SNPs), from the ecotype
    fitness-GWAS (ecotype_fitness_gwas.py). SNP-scored => the SV enrichment isn't circular.
    Axes: zone_union (max over cold/mid/hot = conditional neutrality) + global; LOCO and naive.
  - FITNESS-RELATED = top-frac genome-wide tail of that score.
  - enrichment = are fitness blocks more SV-bearing than SIZE-MATCHED (n_snp-binned) controls?
    two readouts: has_sv (>=1 SV in block) and sv_frac (SV share of block variants, size-robust).
    Size-matched permutation null (lib.matched_perm_test style).

Output -> analysis/grenenet_gea/ecotype_fitness/gwas/block_sv_enrichment.{csv,json}
Env: kmate. Light.
"""
from __future__ import annotations
import os, sys, glob, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

G = f"{lib.GEA}/ecotype_fitness/gwas"
BD = f"{lib.GEA}/blocks_mcf90"
TOP = [0.02, 0.01, 0.005]
MIN_SNP = 3
NPERM = 20000


def load_blocks():
    b = pd.concat([pd.read_csv(f, sep="\t")
                   for f in sorted(glob.glob(f"{BD}/chr*_clq0.9_blocks_clq0.9.tsv"))],
                  ignore_index=True)
    return b


def assign_block(chrom, pos, blocks):
    """Return per-variant block id ('Chrx:start'), '' if in no block. Disjoint blocks."""
    out = np.full(len(pos), "", dtype=object)
    for c, g in blocks.groupby("chrom"):
        m = np.where(chrom == c)[0]
        if len(m) == 0:
            continue
        st = g.start_pos.to_numpy(); en = g.end_pos.to_numpy()
        o = np.argsort(st); st, en = st[o], en[o]
        p = pos[m]
        j = np.searchsorted(st, p, "right") - 1
        ok = (j >= 0) & (j < len(st))
        inb = np.zeros(len(p), bool)
        inb[ok] = p[ok] <= en[j[ok]]
        ids = np.where(inb, [f"{c}:{s}" for s in np.where(ok, st[np.clip(j, 0, len(st)-1)], -1)], "")
        out[m] = ids
    return out


def matched(values, bins, idx, nperm=NPERM, seed=0):
    """mean(values[idx]) vs size(bins)-matched draws from ALL blocks. fold, p one-sided."""
    rng = np.random.default_rng(seed)
    members = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    obs = values[idx].mean()
    null = np.empty((len(idx), nperm))
    for k, i in enumerate(idx):
        pool = members[bins[i]]
        null[k] = values[pool][rng.integers(0, len(pool), nperm)]
    nm = null.mean(0)
    fold = obs / nm.mean() if nm.mean() > 0 else np.nan
    p = (1 + (nm >= obs).sum()) / (nperm + 1)
    return float(obs), float(nm.mean()), float(fold), float(p)


def main():
    z = np.load(f"{G}/gwas_z.npz", allow_pickle=True)
    chrom = np.array([c.capitalize() for c in z["chrom"].astype(str)])  # chr1 -> Chr1
    pos = z["pos"].astype(np.int64)
    vclass = z["vclass"]; phenos = list(z["pheno_names"].astype(str))
    Zt = {"LOCO": z["Z"], "naive": z["Znaive"]}
    is_snp = vclass == 0; is_sv = vclass == 2

    blocks = load_blocks()
    bid = assign_block(chrom, pos, blocks)
    inb = bid != ""
    print(f"variants={len(pos):,}  in a clq0.9 block={inb.sum():,} "
          f"({100*inb.mean():.0f}%)  SV in block={int((is_sv & inb).sum()):,}")

    df = pd.DataFrame({"bid": bid, "is_snp": is_snp, "is_sv": is_sv})
    df = df[inb].copy()
    # per-block composition
    grp = df.groupby("bid")
    comp = grp.agg(n=("is_snp", "size"), n_snp=("is_snp", "sum"), n_sv=("is_sv", "sum"))
    comp["has_sv"] = (comp.n_sv > 0).astype(float)
    comp["sv_frac"] = comp.n_sv / comp.n
    # abs-z per variant per axis, aggregated to block max over the block's SNPs
    absz = {}
    for model, Z in Zt.items():
        for ax, members in (("zone_union", ["w_cold_rel", "w_mid_rel", "w_hot_rel",
                                            "w_cold_cen", "w_mid_cen", "w_hot_cen"]),
                            ("global", ["w_global_rel", "w_global_cen"])):
            cols = [phenos.index(p) for p in members]
            absz[(model, ax)] = np.abs(Z[:, cols]).max(1)  # per-variant max|z| over axis members
    # block SNP-score = max |z| among block's SNPs
    snp_mask_inb = df.is_snp.to_numpy()
    idx_inb = np.where(inb)[0]
    rows = []
    n_snp = comp.n_snp.to_numpy()
    bins = np.digitize(np.log10(np.maximum(n_snp, 1)),
                       np.quantile(np.log10(np.maximum(n_snp, 1)), np.linspace(0, 1, 11)[1:-1]))
    order = comp.index.to_numpy()
    for (model, ax), az in absz.items():
        azb = az[idx_inb]                                 # aligned to df rows
        sc = pd.Series(np.where(snp_mask_inb, azb, np.nan), index=df.bid.to_numpy())
        block_score = sc.groupby(level=0).max().reindex(order).to_numpy()
        scored = np.isfinite(block_score) & (comp.n_snp.to_numpy() >= MIN_SNP)
        for q in TOP:
            thr = np.nanquantile(block_score[scored], 1 - q)
            fit_idx = np.where(scored & (block_score >= thr))[0]
            if len(fit_idx) < 10:
                continue
            for metric in ("has_sv", "sv_frac"):
                vals = comp[metric].to_numpy()
                obs, nm, fold, p = matched(vals, bins, fit_idx)
                rows.append(dict(model=model, axis=ax, top_q=q, metric=metric,
                                 n_fit_blocks=len(fit_idx), obs=round(obs, 4),
                                 null=round(nm, 4), fold=round(fold, 3), p=round(p, 4)))
    out = pd.DataFrame(rows)
    out.to_csv(f"{G}/block_sv_enrichment.csv", index=False)
    json.dump(rows, open(f"{G}/block_sv_enrichment.json", "w"), indent=2)
    print(f"\nblocks scored (>= {MIN_SNP} SNPs): {int((comp.n_snp>=MIN_SNP).sum()):,}  "
          f"of {len(comp):,};  %blocks w/ SV = {100*comp.has_sv.mean():.1f}")
    print("\n=== fitness-related haploblocks: SV enrichment (size-matched) ===")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
