#!/usr/bin/env python
"""Does cross-chromosome founder-frequency disagreement grow with GENERATION (recombination test)?

Coverage is ~uncorrelated with agreement, so the disagreement is not a depth artifact. Two remaining
hypotheses: (A) structural/EM -- each chromosome's panel k-mers give a different founder-collinearity
regime, present already at g0; (B) recombination -- selfing/outcrossing + chromosome-differential
selection make per-chromosome ancestry proportions diverge as generations proceed. The decisive test is
the g0 SEED MIX baseline (unrecombined whole genomes -> should agree tightly) vs the gen1/2/3 trend.

Light post-processing (the 2168-sample scan is cached): adds the 8 seed-mix (g0) samples, merges
generation, and summarizes agreement by generation. Writes cross_chrom_agreement_gen.{csv,png}. Env: kmate.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

WIN = "results/grenenet_kmate_window"; SEED = "results/grenenet_kmate_window_seedmix"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
H = "results/grenenet_gea/hapfreq"
iu = np.triu_indices(5, 1)


def agree(samp, base):
    """Full metric row for one sample (same columns as cross_chrom_agreement.csv)."""
    Hs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_blocks_per_chrom.npz"
        if not os.path.exists(f):
            return None
        Hs.append(np.load(f, allow_pickle=True)[f"{ch}_global_h"].astype(np.float64))
    Hs = np.vstack(Hs); R = np.corrcoef(Hs)
    tv = np.array([0.5 * np.abs(Hs[a] - Hs[b]).sum() for a, b in zip(*iu)])
    gw = Hs.mean(0); r_to_gw = [float(np.corrcoef(Hs[c], gw)[0, 1]) for c in range(5)]
    eff_n = 1.0 / (Hs ** 2).sum(1)
    return dict(sampleid=samp, mean_r=float(R[iu].mean()), min_r=float(R[iu].min()),
                mean_tv=float(tv.mean()), max_tv=float(tv.max()),
                eff_n_mean=float(eff_n.mean()), eff_n_cv=float(eff_n.std() / eff_n.mean()),
                **{f"r_{ch}": r_to_gw[c] for c, ch in enumerate(CHROMS)})


def main():
    df = pd.read_csv(f"{H}/cross_chrom_agreement.csv")     # cohort: full metrics + coverage
    pt = lib.pool_table().drop_duplicates("sampleid").set_index("sampleid")
    df["generation"] = df.sampleid.astype(str).map(pt["generation"].to_dict())

    # g0 seed-mix baseline (8 samples) with the SAME full metric set
    seeds = sorted(p.split("/")[-1].split("_Chr1")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz"))
    srows = [dict(agree(s, SEED), generation=0, coverage=np.nan) for s in seeds if agree(s, SEED)]
    sdf = pd.DataFrame(srows)
    allg = pd.concat([df, sdf], ignore_index=True)
    allg = allg[allg.generation.notna()].copy(); allg["generation"] = allg.generation.astype(int)
    # self-contained enriched table for the notebook (cohort + seed mix, all metrics + gen + coverage)
    allg.to_csv(f"{H}/cross_chrom_agreement_full.csv", index=False)
    allg[["sampleid", "mean_r", "min_r", "mean_tv", "generation"]].to_csv(f"{H}/cross_chrom_agreement_gen.csv", index=False)

    print("=== cross-chromosome agreement by generation ===")
    g = allg.groupby("generation").agg(n=("mean_r", "size"), mean_r_med=("mean_r", "median"),
                                       mean_tv_med=("mean_tv", "median"), min_r_med=("min_r", "median"))
    print(g.round(3).to_string())
    m = allg.generation >= 0
    rho_r = stats.spearmanr(allg.generation[m], allg.mean_r[m]).correlation
    rho_tv = stats.spearmanr(allg.generation[m], allg.mean_tv[m]).correlation
    print(f"\nSpearman(generation, mean_r)  = {rho_r:+.3f}   (negative => agreement falls with generation)")
    print(f"Spearman(generation, mean_tv) = {rho_tv:+.3f}   (positive => disagreement grows with generation)")
    print(f"\ng0 seed-mix (n={len(sdf)}): median mean_r={sdf.mean_r.median():.3f}, median TV={sdf.mean_tv.median():.3f}")
    print(f"gen3 evolved       : median mean_r={allg[allg.generation==3].mean_r.median():.3f}, "
          f"median TV={allg[allg.generation==3].mean_tv.median():.3f}")

    # figure: agreement + disagreement vs generation
    gens = sorted(allg.generation.unique())
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    for axi, col, lab, invert in [(ax[0], "mean_r", "mean pairwise r (agreement)", False),
                                  (ax[1], "mean_tv", "mean pairwise TV (disagreeing mass)", True)]:
        data = [allg[allg.generation == gg][col] for gg in gens]
        bp = axi.boxplot(data, labels=[f"g{gg}" + ("\n(seed mix)" if gg == 0 else "") for gg in gens],
                         showfliers=False, patch_artist=True)
        for b in bp["boxes"]:
            b.set(facecolor="#7aa0c4", alpha=.7)
        meds = [allg[allg.generation == gg][col].median() for gg in gens]
        axi.plot(range(1, len(gens) + 1), meds, "-o", color="firebrick", lw=1.6, ms=4)
        axi.set_ylabel(lab); axi.spines[["top", "right"]].set_visible(False)
    rho = stats.spearmanr(allg.generation, allg.mean_tv).correlation
    ax[0].set_title("(a) shape agreement falls with generation", loc="left", fontsize=10)
    ax[1].set_title(f"(b) disagreement grows with generation (Spearman \\u03c1={rho:+.2f})", loc="left", fontsize=10)
    fig.suptitle("Cross-chromosome founder-h disagreement vs generation — recombination signature", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f"{H}/cross_chrom_agreement_gen.png", dpi=150)
    print(f"\n[done] {H}/cross_chrom_agreement_gen.csv + .png")


if __name__ == "__main__":
    main()
