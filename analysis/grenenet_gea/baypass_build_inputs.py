#!/usr/bin/env python
"""Build BayPass TEMPORAL-contrast inputs on haploblocks (default site 4 smoke test).

⚠️ RETIRED INPUTS (2026-07-08): this reads `hapfreq/hapfreq_matrix.npy` +
`hapfreq_registry.csv`, which are STALE products of the retired hapfreq/Pipeline-B
chain (archived under archive/pipelineB_hapfreq_retired/; built on the pre-Kf_w window
`h`, not the corrected --unit chrom cohort). Do NOT run this as-is — the BayPass plan
(BAYPASS_TEMPORAL_PLAN.md) must first be repointed to a haploblock frequency table
regenerated under the corrected estimator, or retired.

Populations = 8 seedmix reps (founding, gen0) + site-SITE gen-3 plots (evolved).
Markers = testable haploblock one-vs-rest haplotypes, k-1 per block (drop the max-panel-freq
reference), same set logic as the founder GWAS. Per pop the allele counts are
  n1 = round(freq * Neff),  n2 = Neff - n1
with Neff = effective chromosomes (1/Neff = 1/(2*flowers) + 1/coverage; seedmix = fixed 200).

Writes (to OUTDIR):
  geno.txt        BayPass gfile: one row/locus, 2*npop space-sep counts
  poolsize.txt    one line, haploid pool size per pop
  contrast.txt    one line, -1 (founding) / +1 (evolved) per pop
  loci.csv        chrom,start,end per locus row (for the Manhattan later)
  pops.txt        population labels in column order
Env: kmate.  SITE / OUTDIR via env.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = "results/grenenet_gea/hapfreq"
SITE = int(os.environ.get("SITE", 4))
OUTDIR = os.environ.get("OUTDIR", f"{H}/baypass_temporal_site{SITE}")
# Neff = CENSUS (chromosomes collected = 2*flowers), NOT read depth: kMate freqs are EM
# estimates pooling genome-wide reads, so per-locus coverage (~5x) understates precision.
# Census is the real drift sampling. Founding seedmix = large well-mixed lot → high Neff.
SEED_NEFF = 500.0


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    mat = np.load(f"{H}/hapfreq_matrix.npy")                    # (2168, nhap) evolved freqs
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    reg = pd.read_csv(f"{H}/hapfreq_registry.csv")
    seedreps = np.load(f"{H}/hapfreq_p0_seedmix_reps.npy")      # (8, nhap) founding reps
    nhap = mat.shape[1]

    # --- evolved populations: site gen-3 plots (flower-weighted pool within plot) ---
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & (pt.generation == 3)].copy()
    s = s[s.sampleid.astype(str).isin(sidx)]
    evolved_F, evolved_Neff, evolved_lab = [], [], []
    for plot, g in s.groupby("plot"):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        fw = g.flowerscollected.to_numpy(float)
        fw = np.where(np.isfinite(fw) & (fw > 0), fw, 1.0)
        f = (mat[rows] * fw[:, None]).sum(0) / fw.sum()
        neff = 2.0 * fw.sum()                                  # census: chromosomes collected
        evolved_F.append(f); evolved_Neff.append(np.clip(neff, 10, 500)); evolved_lab.append(f"s{SITE}_g3_p{plot}")

    # --- founding population: ONE pop (8 seedmix reps are TECHNICAL reps of one founding
    #     frequency, not distinct populations; keeping them separate makes Omega singular) ---
    found_F = [seedreps.mean(0)]
    found_Neff = [SEED_NEFF]
    found_lab = ["seedmix"]

    F = np.vstack(found_F + evolved_F)                          # (npop, nhap), founding first
    Neff = np.array(found_Neff + evolved_Neff)
    labels = found_lab + evolved_lab
    contrast = np.array([-1] * len(found_F) + [+1] * len(evolved_F))
    npop = len(labels)

    # --- markers: testable, k-1 per block (drop max panel_freq reference) ---
    reg["unit"] = reg.chrom + ":" + reg.unit_start.astype(str) + "-" + reg.unit_end.astype(str)
    testable = (reg.covered & reg.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in reg[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]
    Fk = F[:, ki]                                               # (npop, M)
    M = len(ki)

    # --- counts -> gfile (loci rows, 2*npop cols) ---
    n1 = np.rint(Fk * Neff[:, None]).astype(int)               # (npop, M) allele=hap
    n2 = (np.rint(Neff)[:, None].astype(int) - n1)
    n2 = np.clip(n2, 0, None)
    # interleave per pop: [pop0_n1 pop0_n2 pop1_n1 pop1_n2 ...] per locus row
    G = np.empty((M, 2 * npop), int)
    G[:, 0::2] = n1.T
    G[:, 1::2] = n2.T

    np.savetxt(f"{OUTDIR}/geno.txt", G, fmt="%d")
    with open(f"{OUTDIR}/poolsize.txt", "w") as fh:
        fh.write(" ".join(str(int(x)) for x in np.rint(Neff)) + "\n")
    with open(f"{OUTDIR}/contrast.txt", "w") as fh:
        fh.write(" ".join(str(int(x)) for x in contrast) + "\n")
    with open(f"{OUTDIR}/pops.txt", "w") as fh:
        fh.write("\n".join(f"{l}\t{c}\t{int(round(n))}" for l, c, n in zip(labels, contrast, Neff)) + "\n")
    reg.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].to_csv(f"{OUTDIR}/loci.csv", index=False)

    print(f"[done] {OUTDIR}")
    print(f"  pops: {npop} ({len(found_F)} founding -1, {len(evolved_F)} evolved +1) | loci: {M:,}")
    print(f"  Neff founding={SEED_NEFF:.0f}, evolved median={np.median(evolved_Neff):.0f} "
          f"[{np.min(evolved_Neff):.0f}-{np.max(evolved_Neff):.0f}]")
    print(f"  geno.txt {G.shape} ; contrast {contrast.tolist()}")


if __name__ == "__main__":
    main()
