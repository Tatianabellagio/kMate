"""Every simulated pool and every kMate fit behind the paper's kMate figures, grouped by figure.

A pool is named by its parameters (p231_n50_g3_self97_s42_10x) and a fit by its pool and estimator settings
(p231_n50_g3_self97_s42_10x.p231.chrom). A pool or fit listed by two figures is simulated or fitted once.
`python runs.py` prints the number of pools and fits per figure.
"""
from dataclasses import dataclass

SEEDS5, SEEDS3 = (42, 43, 44, 45, 46), (42, 43, 44)


@dataclass(frozen=True)
class Pool:
    founders: str              # founder set (config.FOUNDER_SETS)
    n: int                     # individuals
    gen: int = 0               # generations of reproduction (0 = every individual is one whole founder)
    selfing: float = 0.0       # per-offspring probability of selfing instead of outcrossing
    seed: int = 42
    depth: int = 10            # total read depth of the pool
    dominant: float = 0.0      # pool share of individual 1 (0 = every individual the same share)
    founder_type: str = ""     # "long" / "short": 50 founders, 40 / 5 of them long-read (founder-type figure)

    @property
    def id(self):
        tags = [self.founders, f"n{self.n}", f"g{self.gen}"]
        tags += [f"self{round(self.selfing * 100)}"] if self.selfing else []
        tags += [f"dom{round(self.dominant * 100)}"] if self.dominant else []
        tags += [f"{self.founder_type}heavy"] if self.founder_type else []
        return "_".join(tags + [f"s{self.seed}", f"{self.depth}x"])


@dataclass(frozen=True)
class Fit:
    pool: Pool
    kmers: str = ""            # k-mer panel (config.KMER_PANELS); "" for hapFIRE
    variants: str = ""         # variant panel the allele frequencies are projected on ("" = founder frequencies only)
    unit: str = "chrom"        # chrom (one mixture per chromosome) | bp<N>k (N kb windows) | ld (r2 >= 0.1 LD blocks)
    normalize: str = "per_founder"   # EM M-step normalization: per_founder (production) | global
    tool: str = "kmate"        # kmate | hapfire (on its GrENE-Net panel; the other fields unused)
    timed: bool = False        # run alone on a whole node of one type (config.TIMED_RES), for the speed comparison

    @property
    def id(self):
        timed = ".timed" if self.timed else ""
        if self.tool == "hapfire":
            return f"{self.pool.id}.hapfire{timed}"
        norm = ".globalnorm" if self.normalize == "global" else ""
        return f"{self.pool.id}.{self.kmers}.{self.unit}{norm}{timed}"


def figures():
    """{figure: [Fit, ...]}"""
    F = {}
    # Main figure D, E and kMate vs hapFIRE: founder frequencies and run time by pool size and depth.
    # hapFIRE gets the same founders, its reads from its own SNP-only genomes.
    F["pool_size_depth"] = [f for n in (1, 2, 5, 20, 50, 150, 231) for d in (1, 5, 10) for s in SEEDS5 for f in
                            (Fit(Pool("p231", n, seed=s, depth=d), "p231", "p231"),
                             Fit(Pool("greneNet231", n, seed=s, depth=d), tool="hapfire"))]
    # Speed (kMate vs hapFIRE supplement C, D; caption): the 50-founder pools of pool_size_depth fitted again, each fit alone
    # on a whole node, so other users' jobs cannot slow either tool. Accuracy comes from pool_size_depth.
    F["timing"] = [f for d in (1, 5, 10) for s in SEEDS5 for f in
                   (Fit(Pool("p231", 50, seed=s, depth=d), "p231", "p231", timed=True),
                    Fit(Pool("greneNet231", 50, seed=s, depth=d), tool="hapfire", timed=True))]
    # Main figure B, C and the six-scenario supplement: allele-frequency accuracy (10x, seed 42).
    F["af_scenarios"] = [Fit(p, "p231", "p231") for p in AF_SCENARIOS.values()]
    # Founder type: per-founder frequencies with and without the private k-mer filter, with production (per-founder) and
    # legacy (global) normalization, in a uniform 231-founder pool and two 50-founder pools weighted to one founder type.
    F["founder_type"] = [Fit(p, k, normalize=nm)
                         for p in (Pool("p231", 231), Pool("p231", 50, founder_type="long"), Pool("p231", 50, founder_type="short"))
                         for k in ("p231_unfilt", "p231_filt") for nm in ("per_founder", "global")]
    # Estimation unit on the mixed panel: 50 founders, outcrossing or 97% selfing for 1, 3, 10 generations.
    F["units_p231"] = [Fit(Pool("p231", 50, g, sf, s), "p231", "p231", u)
                       for sf in (0.0, 0.97) for g in (1, 3, 10) for s in SEEDS3 for u in ("chrom", "bp10k", "bp50k", "ld")]
    # Long-read-only panel: pools of the 78 long-read founders, kMate on the 78-founder panel (with and without the
    # private k-mer filter) and on the production 231-founder panel; hapFIRE on its panel cut to the 78.
    F["longread_grid"] = [f for n in (1, 2, 5, 20, 50, 78) for d in (1, 5, 10) for s in SEEDS5 for f in
                          [Fit(Pool("p78", n, seed=s, depth=d), k, v) for k, v in (("p78_unfilt", "p78"), ("p78_filt", "p78"), ("p231", "p231"))]
                          + [Fit(Pool("greneNet78", n, seed=s, depth=d), tool="hapfire")]]
    # Estimation unit on the long-read panel, windows from 10 kb to 2 Mb.
    F["units_p78"] = [Fit(Pool("p78", 50, g, sf, s), "p78_unfilt", "p78", u)
                      for sf in (0.0, 0.97) for g in (1, 3, 10) for s in SEEDS3
                      for u in ("chrom", "bp10k", "bp25k", "bp50k", "bp100k", "bp250k", "bp500k", "bp1000k", "bp2000k", "ld")]
    return F


# Allele-frequency scenarios, by the names the figures use.
AF_SCENARIOS = {
    "n1_g0": Pool("p231", 1), "n2_g0": Pool("p231", 2), "n50_g0": Pool("p231", 50), "n231_g0": Pool("p231", 231),
    "n231_g1_self97": Pool("p231", 231, 1, 0.97), "n50_g1_self97": Pool("p231", 50, 1, 0.97),
    "n50_g3_self97": Pool("p231", 50, 3, 0.97), "n50_g3_dom50_self97": Pool("p231", 50, 3, 0.97, dominant=0.5),
}


def all_fits(names=None):
    """Unique fits of the named figures (all figures if None), in figure order."""
    out = {}
    for name, fits in figures().items():
        if names is None or name in names:
            out.update((f.id, f) for f in fits)
    return out


def all_pools(names=None):
    return {f.pool.id: f.pool for f in all_fits(names).values()}


if __name__ == "__main__":
    for name, fits in figures().items():
        print(f"{name:16s} {len({f.pool.id for f in fits}):4d} pools {len(fits):4d} fits")
    print(f"{'total (unique)':16s} {len(all_pools()):4d} pools {len(all_fits()):4d} fits")
