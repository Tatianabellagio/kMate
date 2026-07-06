#!/usr/bin/env python3
"""Build benchmarks/localonly_p231/FLOOR_DERIVATION.ipynb consolidating:
  - the local-only p231 benchmark recap (dynld vs 10kb),
  - the --min-kmers-per-block floor derivation (3 findings + figures).
Run, then execute with the `basic` env. House style: construct via nbformat."""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from pathlib import Path

REPO = "/global/scratch/users/tbellg/kmate"
OUT = Path(REPO) / "benchmarks/localonly_p231/FLOOR_DERIVATION.ipynb"
nb = new_notebook()
C = nb.cells


def md(s): C.append(new_markdown_cell(s.strip("\n")))
def co(s): C.append(new_code_cell(s.strip("\n")))


md(r"""
# Local-only window mode: the `--min-kmers-per-block` floor — derivation & benchmark

**Question.** kMate's window mode gates each genomic block on a minimum nonzero-k-mer
count (`--min-kmers-per-block`, default **200**, **50** in the bench). Below it a block
falls back to the chrom-wide *h* (anchored) or abstains → NaN (local-only). Two things
to settle:

1. Is the **variable** (a k-mer count) the right thing to gate on, or should it be a
   mathematical *identifiability* condition — "how many distinct k-mers to differentiate
   the founders/ecotypes"?
2. Is the **value 200** principled, or arbitrary?

This notebook consolidates (A) the local-only p231 benchmark that motivated the question
(dynld vs 10 kb units), and (B) a floor=1 per-block diagnostics experiment that answers
both. Production code is untouched; the experiment reuses kMate's own loaders, EM solver,
and Fisher/resolvability tools.
""")

co(r"""
import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from IPython.display import Image, display

REPO = Path("/global/scratch/users/tbellg/kmate")
RES  = REPO / "benchmarks/localonly_p231/results"
FD   = REPO / "benchmarks/localonly_p231/floor_diag"
SIMS = REPO / "benchmarks/p231/sims"
KEYS = ["chrom", "pos", "ref_len", "alt_len"]
POOLS = ["cov10_n50_g0_s42_hotspots_p231_chr1",
         "cov10_n231_g1_s42_self97_hotspots_p231_chr1",
         "cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1"]
SHORT = {"cov10_n50_g0_s42_hotspots_p231_chr1": "n50 g0 outcross",
         "cov10_n231_g1_s42_self97_hotspots_p231_chr1": "n231 g1 SELFING",
         "cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1": "n50 g3 dom-sel"}
pd.set_option("display.width", 160, "display.max_columns", 40)
print("repo:", REPO)
""")

md(r"""
## Part A — Local-only benchmark recap: dynld_K500 vs 10 kb windows

All 11 p231 sims, 5 modes, NaN-aware scoring. `call_rate` = fraction of records with a
finite estimate; accuracy metrics are computed only where the estimate is finite
("where it speaks, is it right?"). Below: the ALL-class headline by mode × unit.
""")

co(r"""
tab = pd.read_csv(RES / "localonly_p231_table.tsv", sep="\t")
allc = tab[tab.var_class == "ALL"].copy()
allc["mode_unit"] = allc["mode"] + " / " + allc["unit"].astype(str)
hd = (allc.groupby("mode_unit")
        .agg(mean_call=("call_rate", "mean"), mean_R2=("R2", "mean"),
             mean_RMSE=("RMSE", "mean"), n=("pool", "nunique"))
        .round(3).sort_values("mean_R2", ascending=False))
hd
""")

md(r"""
**dynld_K500 vs 10 kb (local-only), the contrast that prompted the floor question.**
10 kb wins on *both* axes — higher call rate AND higher R² — even though dynld abstains
more. So "abstain more → keep only easy records → lower error" is *false* here; both
axes are driven by the same hidden variable, per-block k-mer supply. Uniform 10 kb
windows spread the k-mer budget evenly and keep more blocks above the resolution floor;
LD-shaped dynld blocks have a heavy tail of k-mer-starved blocks that both abstain and,
where called, fit poorly.
""")

co(r"""
lo = allc[allc["mode"] == "localonly"].pivot_table(
        index="pool", columns="unit", values=["call_rate", "R2"]).round(3)
lo
""")

md(r"""
## Part B — Is the floor arbitrary? Per-block diagnostics experiment

**Design** (`scripts/block_floor_diag.py`). 3 sims spanning the founder-count
(50 vs 231 = identifiability difficulty) and regime (g0 / selfing / recombinant+selection)
axes, each fit local-only at **floor=1, no anchor, no HMM smoothing** — i.e. fit *every*
non-empty block on its own. One run gives the whole sweep in post: a record is "called at
floor T" iff its block has `nnz ≥ T`. Per block we logged the raw count plus principled
resolvability metrics derived from kMate's own `fisher_information_h` +
`_resolvability_from_J`:

| metric | meaning |
|---|---|
| `nnz` | nonzero-count k-mers in the block (the current proxy) |
| `effrank_design` | spectral-entropy rank of the founder-carriage matrix (coverage-free identifiability) |
| `effrank_fisher` / `cond` | resolvability at the data (Fisher info on the EM support) |
| `n_present` | # founders carrying ≥1 nonzero k-mer in the block |
""")

co(r"""
# helper: join per-record est to truth, attach the record's block diagnostics
def load_joined(pool, unit):
    diag = pd.read_csv(FD / f"{pool}_{unit}.blockdiag.tsv", sep="\t")
    est  = pd.read_csv(FD / f"{pool}_{unit}.recest.tsv", sep="\t")
    tr   = pd.read_csv(SIMS / pool / "recomb_truth_raw.tsv.gz", sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz", "n_present", "effrank_design", "effrank_fisher"]],
                on="block", how="left")
    m = m[np.isfinite(m.alt_freq.values)].copy()
    m["e2"] = (m.alt_freq - m.truth_af) ** 2
    return m

def per_block(m, min_rec=20):
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n":    m.groupby("block").alt_freq.size(),
                      "nnz":  m.groupby("block").nnz.first(),
                      "n_present": m.groupby("block").n_present.first(),
                      "effrank_design": m.groupby("block").effrank_design.first(),
                      "effrank_fisher": m.groupby("block").effrank_fisher.first()})
    return g[g.n >= min_rec]

joined = {(p, u): load_joined(p, u) for p in POOLS for u in ["w10kb", "dynldK500"]}
print("loaded", len(joined), "(pool,unit) joins")
""")

md(r"""
### Finding 1 — the *variable* is right: k-mer COUNT predicts AF error best

Spearman ρ of per-block AF RMSE vs each predictor (more negative = predicts lower error).
`nnz` is the strongest, most consistent predictor — **beating** the identifiability/rank
metrics.

**Why, against intuition:** the scored estimand is per-record **AF, not per-founder *h***.
AF cancels founder-collinearity (carriers in a block share their `var_pa` value, so the
non-identifiability bias cancels — the h-certainty result). So block AF error is governed
by **counting noise / total Fisher information ∝ nnz**, *not* by whether individual
founders are separable. Gating on a k-mer count is therefore the theoretically correct
choice of variable. (`effrank_design` even flips *positive*: conditional on supply, more
local-ancestry diversity = harder = more error — it measures mixture complexity, not
precision.)
""")

co(r"""
rows = []
for (p, u), m in joined.items():
    b = per_block(m)
    rows.append(dict(pool=SHORT[p], unit=u, n_blocks=len(b),
        **{k: round(spearmanr(b[k], b.rmse).correlation, 3)
           for k in ["nnz", "effrank_design", "effrank_fisher", "n_present"]}))
corr = pd.DataFrame(rows)
corr
""")

md(r"""
### Finding 2 — the *value* 200 is NOT where error plateaus

Pool the 3 sims; bin blocks by nonzero-k-mer count; report median + p90 per-block AF RMSE.
Median error keeps falling until **nnz ≈ 1000–2000**, where it hits the accuracy floor
(~0.048). At the current floors median error is ~2× the floor (≈0.10 at nnz=200) and the
p90 "worst-block" tail stays > 0.18 until nnz ≥ 1024. So **200 is a "minimum worth
attempting," not "enough to resolve the mixture."** It is the small sibling of the
~7000-k-mer *full-h* resolution figure — AF needs less than full *h*, but still ~1–2k,
far above 200.
""")

co(r"""
EDGES = [1, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 100000]
def binned(unit):
    b = pd.concat([per_block(joined[(p, unit)]) for p in POOLS], ignore_index=True)
    b["bin"] = pd.cut(b.nnz, EDGES, right=False)
    g = b.groupby("bin", observed=True).rmse.agg(
        nblk="size", median="median", p90=lambda x: x.quantile(.9)).round(4)
    return g
display(binned("w10kb"))
""")

md(r"""
#### Why the low-nnz bump? — it's the centromere (raw-data / swarm view)

The w10kb low-nnz bins are non-monotonic at the left (e.g. the 16–31 bin pops up). That's
**not noise — it's the Chr1 (peri)centromere.** A 10 kb window with only 16–31 nonzero
k-mers is k-mer-starved *despite* normal size, which happens in repeat-rich centromeric
sequence (most panel k-mers there are non-unique / repeat-guarded away). Those windows
also carry genuinely hard, structurally complex variation. Diagnostic for the low-nnz
w10kb bins (`%cen` = fraction in 12.5–17.5 Mb):
""")

co(r"""
def per_block_pos(pool, unit, min_rec=20):
    diag = pd.read_csv(FD / f"{pool}_{unit}.blockdiag.tsv", sep="\t")
    m = load_joined(pool, unit)
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n": m.groupby("block").alt_freq.size(),
                      "nnz": m.groupby("block").nnz.first()})
    g = g.merge(diag.set_index("block")[["start"]], left_index=True, right_index=True)
    g = g[g.n >= min_rec].copy()
    g["cen"] = (g.start >= 12_500_000) & (g.start <= 17_500_000)
    return g

w = pd.concat([per_block_pos(p, "w10kb") for p in POOLS], ignore_index=True)
w["bin"] = pd.cut(w.nnz, EDGES, right=False)
diag_tbl = w.groupby("bin", observed=True).agg(
    n=("rmse", "size"), med_rmse=("rmse", "median"),
    p90_rmse=("rmse", lambda x: x.quantile(.9)),
    pct_centromere=("cen", lambda x: round(100 * x.mean(), 0))).head(6).round(3)
diag_tbl
""")

co(r"""
# raw-data swarm/strip view; centromere blocks in red
display(Image(filename=str(RES / "block_floor_strip.png")))
""")

md(r"""
The low-nnz bins are 50–75% centromeric (vs ~46% by 128–255), peaking exactly at the bumpy
16–31 bin, and the worst windows sit at ~14.4–15.3 Mb (dead centre of CEN1) with high true
AF. So the bump is a real structural signal, and it *strengthens* the floor argument: the
blocks below the floor are not random — they are pathological centromere/repeat windows
where abstaining (local-only) or falling back to global (anchored) is exactly the right
move. (Ties to the repeat-contamination guard, which is *why* those windows are k-mer-poor.)

### Finding 2b — it's k-mer COUNT, not DENSITY (why dynld ≈ 10 kb)

The curves above collapse dynld onto w10kb when plotted vs *raw* nnz. Is that because
total supply is the sufficient statistic (so block size / variant content are irrelevant
once you fix the count), or would a *density* (k-mers per variant, per kb) explain more?
Per-block `h` is fit **once** from all the block's k-mers, and every variant projects
through that same `h` — so accuracy should depend on the *total* k-mers feeding the EM,
not on how densely they're packed. Test it: Spearman ρ(RMSE, metric) for **dynld**
(variable block size, so the comparison is meaningful):
""")

co(r"""
def per_block_dense(pool, unit, min_rec=20):
    diag = pd.read_csv(FD / f"{pool}_{unit}.blockdiag.tsv", sep="\t")
    m = load_joined(pool, unit)
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n_var": m.groupby("block").alt_freq.size()})
    g = g.merge(diag.set_index("block")[["nnz", "start", "end"]],
                left_index=True, right_index=True)
    g = g[g.n_var >= min_rec].copy()
    g["kmers_per_var"] = g.nnz / g.n_var
    g["kmers_per_kb"]  = g.nnz / ((g.end - g.start + 1) / 1000.0)
    return g

rows = []
for p in POOLS:
    for u in ["dynldK500", "w10kb"]:
        b = per_block_dense(p, u)
        rows.append(dict(pool=SHORT[p], unit=u,
            **{k: round(spearmanr(b[k], b.rmse).correlation, 3)
               for k in ["nnz", "kmers_per_var", "kmers_per_kb", "n_var"]}))
pd.DataFrame(rows)
""")

md(r"""
Raw `nnz` (ρ≈−0.6 on dynld) beats k-mers/variant (−0.27) and k-mers/kb (−0.13). The
within-band check is even sharper: holding nnz fixed, density does **not** reduce error
(ρ≈0 or wrong-signed). So **at equal total k-mers, a wide variant-rich block and a tiny
one have the same AF error** — total supply is the sufficient statistic, which is why the
two units lie on one curve vs nnz. The figure below shows only the *raw-count* panel
collapses; the density panels do not.
""")

co(r"""
display(Image(filename=str(RES / "block_floor_normalized.png")))
""")

md(r"""
### Finding 3 — no clean elbow: it's a coverage/accuracy dial

The full sweep (call-rate & R²-on-called vs the floor) is smooth and concave: each
doubling of the floor buys ~+0.01–0.02 R² at a rising coverage cost. Marginal
R²-per-coverage is roughly flat up to ~200 then falls off above ~400 — so **200 sits at a
defensible soft knee, and 50 is a valid point** (≈6 pts more coverage for a small accuracy
cost). Neither is a discoverable optimum. To actually reach the accuracy floor you'd need
a floor (~1000) so high it abstains on ~70% of records — windows don't hold that many
k-mers at 10× / 5–10 kb — which is exactly why production uses **global/anchored**, not a
higher local floor.
""")

co(r"""
sw = pd.read_csv(RES / "block_floor_summary.tsv", sep="\t")
sw["pool"] = sw.pool.map(SHORT)
g0 = sw[(sw.gate == "nnz") & (sw.unit == "w10kb") & (sw.pool == "n50 g0 outcross")]
g0[["floor", "call_rate", "R2", "RMSE"]].round(3).reset_index(drop=True)
""")

co(r"""
# detailed multi-panel figures (pre-rendered by analyze_block_floor.py)
for png in ["block_floor_sweep.png", "block_floor_derivation.png"]:
    display(Image(filename=str(RES / png)))
""")

md(r"""
## Part C — Panel & the private-k-mer filter (p80 vs p231)

The 231 production panel applies **filt2inv**: drop `ac==1` private singletons (error/
repeat-prone in a *mixed* long+short-read panel) and `ac==F` invariants. Two questions:
does an all-long-read 80-founder panel resolve better, and — the clean within-panel test —
does *keeping* private k-mers help when the founders are all long-read? Same 3 scenarios,
w10kb, floor=1. Three series: p231 filt2inv / p80 filt2 (private dropped) / p80 unfiltered
(private kept). (`run_block_floor_p80.sbatch`, job 35291298.)
""")

co(r"""
def pb_panel(fd, pool, variant, sims, truthfn, min_rec=20):
    base = REPO / "benchmarks/localonly_p231" / fd
    tag = f"{pool}_w10kb" + (f"_{variant}" if variant else "")
    if not (base / f"{tag}.blockdiag.tsv").exists():
        return None
    diag = pd.read_csv(base / f"{tag}.blockdiag.tsv", sep="\t")
    est = pd.read_csv(base / f"{tag}.recest.tsv", sep="\t")
    tr = pd.read_csv(REPO / sims / pool / truthfn, sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz", "start"]], on="block", how="left")
    m = m[np.isfinite(m.alt_freq.values)].copy(); m["e2"] = (m.alt_freq - m.truth_af) ** 2
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n": m.groupby("block").alt_freq.size(),
                      "nnz": m.groupby("block").nnz.first(),
                      "start": m.groupby("block").start.first()})
    g = g[g.n >= min_rec].copy()
    g["cen"] = (g.start >= 12_500_000) & (g.start <= 17_500_000)
    return g

SCEN = ["cov10_n50_g0_s42_hotspots", "cov10_n231_g1_s42_self97_hotspots",
        "cov10_n50_g3_s42_hotspots_dom500nr"]
SER = [("p231 filt2inv (231 mixed)", "floor_diag", "_p231_chr1", "", "benchmarks/p231/sims", "recomb_truth_raw.tsv.gz"),
       ("p80 filt2 (80 long-read)", "floor_diag_p80", "_p80_chr1", "filt2", "benchmarks/p80/sims", "recomb_truth.tsv.gz"),
       ("p80 unfiltered (keep private)", "floor_diag_p80", "_p80_chr1", "unfiltered", "benchmarks/p80/sims", "recomb_truth.tsv.gz")]
rows = []
for lab, fd, suf, var, sims, tfn in SER:
    parts = [pb_panel(fd, s + suf, var, sims, tfn) for s in SCEN]
    df = pd.concat([p for p in parts if p is not None], ignore_index=True)
    hi = df[df.nnz >= 512]
    rows.append(dict(series=lab, arm_RMSE=round(hi[~hi.cen].rmse.median(), 3),
                     cen_RMSE=round(hi[hi.cen].rmse.median(), 3),
                     cen_penalty=round(hi[hi.cen].rmse.median() - hi[~hi.cen].rmse.median(), 3)))
pd.DataFrame(rows)
""")

md(r"""
**Distribution view (strip).** Same three series as full per-block RMSE distributions per
nnz bin (one panel each; centromere blocks red; black bar = median). Note in *p80
unfiltered* the low-nnz bins hold far fewer blocks (keeping private k-mers shifts blocks
rightward to higher supply) and the distributions are visibly tighter.
""")

co(r"""
display(Image(filename=str(RES / "panel_compare_strip.png")))
""")

md(r"""
**Joint-density view.** The strip/median plots show *where the error is* but hide *how many
blocks sit at each supply level*. Here one jointplot per series: the main panel is a hexbin
of nnz (log-x) vs per-block RMSE coloured by **block count** (log scale, **reversed** — the
sparse tail gets the strong colour so it doesn't get lost; the dense floor blob goes white);
the **top marginal** is the nnz supply histogram and the **right marginal** the RMSE
histogram (series colour = all blocks, **grey** = centromere). What this makes legible:
in *p80 unfiltered* the top marginal is **shifted right** with almost nothing at low nnz —
keeping private k-mers both adds supply and moves blocks up the supply axis — and the hexbin
bulk is **tightest against the floor**. The centromere (grey) is the elevated-RMSE mass that
survives in every panel.

The grey dotted lines are candidate `--min-kmers-per-block` floors (100 / 200 / 500). Each
**vertical** line is labelled with the **% of scored blocks it would drop to NaN** *and* the
**mean per-block AF RMSE of the blocks that survive** the floor; the matching **horizontal**
dotted line draws that kept-mean level rightward over the retained region. Even 200 drops only ~7–9%
(500 ~13–16%), so the floor is cheap in coverage — and *p80 unfiltered* drops the fewest at
every threshold (the rightward supply shift again). The payoff is **modest and diminishing**:
raising the floor 0→500 pulls the mean only 0.073→0.062 (p231) — and the *median* barely
moves (0.058→0.054), because the floor removes exactly the high-error low-nnz tail that the
mean feels and the median doesn't.
""")

co(r"""
display(Image(filename=str(RES / "panel_compare_joint.png")))
""")

md(r"""
**Read-out.** (1) **Keeping private k-mers helps** (green lowest, and lower *at matched
nnz* — private ac=1 k-mers each pin one founder, so they raise quality-per-k-mer, not just
supply; arm error drops to 0.044, below the filtered ~0.05 floor). (2) The long-read panel
helps **specifically in the centromere** (penalty +0.016 vs +0.027) but not the arm
(coverage-limited either way) — note the **founder-count confound** (p80 also has fewer
founders). (3) The centromere penalty **persists even unfiltered** (+0.014) → its residual
is bias (coverage / missingness / single-h), not uniqueness.

> ⚠️ **Closed-loop caveat.** Sim reads are generated *from* the panel, so every panel
> private k-mer is real — there are **no short-read assembly artifacts**, the exact failure
> mode `filt2inv` exists to suppress. So this is the **best case** for keeping private
> k-mers. Honest reading: private k-mers carry real signal *when the panel's private k-mers
> are trustworthy* (plausible for all-long-read p80; not for the mixed 231 panel — hence the
> filter). The repeat-contamination guard was also off (consistent across all series).

## Takeaway for the paper

- **Gate on a k-mer count — and say why.** A k-mer count is the sufficient statistic for
  AF precision *because AF is collinearity-robust* (ties to the h-uncertainty framework).
  This is the principled justification the floor was missing; an identifiability/rank gate
  would be correct only if the estimand were per-founder *h*.
- **200 is an operating point, not a resolution requirement.** The AF resolution
  requirement is ~1–2k k-mers/block (≪ the ~7k for full *h*), and is unreachable at 10× in
  5–10 kb windows. 200 / 50 are points on a smooth coverage/accuracy dial near a soft knee.
- **Reinforces the standing decision:** global/anchored for production; local-only only
  when a global prior would mask the signal (e.g. recombinant MAGIC), at a real cost in
  coverage and precision.

*Sources:* `FLOOR_DERIVATION.md`; scripts `block_floor_diag.py`, `analyze_block_floor.py`,
`plot_floor_curve.py`, `plot_floor_box.py`, `analyze_block_density.py`; sbatch job 35286184.
Figures in `results/`.
""")

OUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, str(OUT))
print("wrote", OUT, "with", len(C), "cells")
