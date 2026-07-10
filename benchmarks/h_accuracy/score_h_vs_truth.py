#!/usr/bin/env python
"""
score_h_vs_truth.py  --  reusable PER-FOUNDER h benchmark for kMate.

The production AF benchmark (benchmarks/build_benchmark_table.py) scores allele
frequency only (R2 ~ 0.99). kMate's *per-founder h* is a separate, weaker
estimand: on the seed mix / n231_g0 it absorbs many founders to ~0 and, in the
raw arm, over-credits long-read (cactus) founders.  This scorer scores h the
same way for any candidate estimate against a simulation truth (pool_weights.tsv)
so production can be baselined and prototype fixes can be scored with one call.

Run it (compute + plot) in the `basic` mamba env (has numpy/pandas/scipy/matplotlib;
matplotlib is fine here -- it is the `plotting` env that hangs).

    conda activate basic
    python benchmarks/h_accuracy/score_h_vs_truth.py \
        --est   <h_per_chrom.npz | h.tsv> \
        --truth benchmarks/p231/sims/<sim>/pool_weights.tsv \
        --label my_run \
        --out   benchmarks/h_accuracy

Outputs (under --out, prefixed by --label):
    <label>_metrics.json      scalar metrics
    <label>_per_founder.csv   aligned truth/est per founder + flags
    <label>_diag.png          ordered (expected=truth) + log-log scatter diagnostics

METRICS reported:
  * per-founder h RMSE, MAE, Pearson/Spearman r (linear AND log space)
  * ABSORBED founders: est < absorb_thr while truth >= truth_present_thr (count + ids)
  * worst-under-called and worst-over-called founders
  * LONG-READ (cactus) vs SHORT-READ (PG) mass balance: sum(est) vs sum(truth)
    over each class -- catches any fix that reintroduces long-read over-credit.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# loaders
# --------------------------------------------------------------------------- #
def load_h_estimate(path, chrom=None):
    """Return a dict founder_id(str) -> h(float), normalized to sum 1.

    Accepts:
      * .npz  with 'founders' array + one or more per-chrom float arrays
              (keys like 'Chr1'..'Chr5'); h averaged over the chrom keys
              present (or the single --chrom if given).
      * .tsv/.csv with a founder-id column and an h/weight column
              (auto-detected).
    """
    if path.endswith(".npz"):
        d = np.load(path, allow_pickle=True)
        founders = [str(x) for x in d["founders"]]
        chrom_keys = [k for k in d.files if k != "founders"]
        if chrom is not None:
            if chrom not in chrom_keys:
                sys.exit(f"[err] chrom {chrom} not in {path} (have {chrom_keys})")
            chrom_keys = [chrom]
        mats = np.vstack([np.asarray(d[k], float) for k in chrom_keys])
        h = mats.mean(axis=0)   # per-founder h averaged across chroms present
        est = dict(zip(founders, h))
    else:
        sep = "\t" if path.endswith((".tsv", ".txt")) else ","
        df = pd.read_csv(path, sep=sep)
        fcol = _pick_col(df, ["founder", "founder_id", "ecotype", "id", "sample"])
        hcol = _pick_col(df, ["h", "weight", "freq", "frequency", "est", "value"])
        est = {str(f): float(v) for f, v in zip(df[fcol], df[hcol])}
    s = sum(est.values())
    if s > 0:
        est = {k: v / s for k, v in est.items()}
    return est


def load_truth(path):
    """pool_weights.tsv -> dict founder_id(str) -> weight(float), sum-1 normalized."""
    df = pd.read_csv(path, sep="\t")
    fcol = _pick_col(df, ["founder", "founder_id", "ecotype", "id"])
    wcol = _pick_col(df, ["weight", "freq", "count", "h"])
    truth = {str(f): float(w) for f, w in zip(df[fcol], df[wcol])}
    s = sum(truth.values())
    if s > 0:
        truth = {k: v / s for k, v in truth.items()}
    return truth


def _pick_col(df, candidates):
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in low:
            return low[cand]
    # fall back: first column for id, last numeric for value
    raise SystemExit(f"[err] none of {candidates} in columns {list(df.columns)}")


def load_split(path):
    d = json.load(open(path))
    return {"cactus": set(map(str, d["cactus"])), "PG": set(map(str, d["PG"]))}


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #
def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).correlation)


def pearson(a, b):
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def score(est, truth, split, absorb_thr=1e-3, present_thr=1e-3, eps=1e-6):
    founders = sorted(set(est) | set(truth), key=lambda x: (len(x), x))
    e = np.array([est.get(f, 0.0) for f in founders])
    t = np.array([truth.get(f, 0.0) for f in founders])
    n = len(founders)

    diff = e - t
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
    le, lt = np.log10(e + eps), np.log10(t + eps)

    # absorbed: present in truth but read as ~0 by the estimate
    present = t >= present_thr
    absorbed_mask = present & (e < absorb_thr)
    absorbed = [(founders[i], float(t[i]), float(e[i]))
                for i in np.where(absorbed_mask)[0]]
    absorbed.sort(key=lambda x: -x[1])

    order_under = np.argsort(-(t - e))   # truth >> est
    order_over = np.argsort(-(e - t))    # est  >> truth
    worst_under = [(founders[i], float(t[i]), float(e[i]),
                    "cactus" if founders[i] in split["cactus"] else "PG")
                   for i in order_under[:15]]
    worst_over = [(founders[i], float(t[i]), float(e[i]),
                   "cactus" if founders[i] in split["cactus"] else "PG")
                  for i in order_over[:15]]

    # long-read (cactus) vs short-read (PG) mass balance
    def _mass(fset, vec):
        return float(sum(vec[i] for i, f in enumerate(founders) if f in fset))
    lr_est, lr_tru = _mass(split["cactus"], e), _mass(split["cactus"], t)
    sr_est, sr_tru = _mass(split["PG"], e), _mass(split["PG"], t)

    metrics = {
        "n_founders": n,
        "n_present_in_truth": int(present.sum()),
        "rmse": rmse,
        "mae": mae,
        "pearson_linear": pearson(e, t),
        "spearman_linear": spearman(e, t),
        "pearson_log": pearson(le, lt),
        "spearman_log": spearman(le, lt),
        "absorb_thr": absorb_thr,
        "present_thr": present_thr,
        "n_absorbed": len(absorbed),
        "absorbed_founders": [{"founder": f, "truth": tr, "est": es}
                              for f, tr, es in absorbed],
        "worst_undercalled": [{"founder": f, "truth": tr, "est": es, "class": c}
                              for f, tr, es, c in worst_under],
        "worst_overcalled": [{"founder": f, "truth": tr, "est": es, "class": c}
                             for f, tr, es, c in worst_over],
        "mass_balance": {
            "cactus_longread": {"truth": lr_tru, "est": lr_est,
                                "est_minus_truth": lr_est - lr_tru},
            "PG_shortread": {"truth": sr_tru, "est": sr_est,
                             "est_minus_truth": sr_est - sr_tru},
        },
    }
    per_founder = pd.DataFrame({
        "founder": founders,
        "class": ["cactus" if f in split["cactus"] else
                  ("PG" if f in split["PG"] else "unknown") for f in founders],
        "truth": t,
        "est": e,
        "diff_est_minus_truth": diff,
        "absorbed": absorbed_mask,
    })
    return metrics, per_founder


# --------------------------------------------------------------------------- #
# figure
# --------------------------------------------------------------------------- #
def make_figure(per_founder, metrics, label, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = per_founder.copy()
    cmap = {"cactus": "#1b7837", "PG": "#7fbf7b", "unknown": "#999999"}
    colors = df["class"].map(cmap)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (1) ordered-by-kMate expected=truth diagnostic
    ax = axes[0]
    dsort = df.sort_values("est", ascending=False).reset_index(drop=True)
    x = np.arange(len(dsort))
    ax.scatter(x, dsort["est"], s=14,
               c=dsort["class"].map(cmap), label="_kmate", zorder=3)
    ax.scatter(x, dsort["truth"], s=10, facecolors="none",
               edgecolors="#d6604d", linewidths=0.8, label="truth", zorder=2)
    ax.axhline(0, color="#cccccc", lw=0.5)
    ax.set_xlabel("founder rank (ordered by kMate ĥ, descending)")
    ax.set_ylabel("founder frequency")
    ax.set_title(f"ordered by kMate, expected = truth\n"
                 f"RMSE={metrics['rmse']:.4g}  absorbed={metrics['n_absorbed']}")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], marker="o", ls="", mfc=cmap["cactus"],
                  mec=cmap["cactus"], label="ĥ cactus (long-read)"),
           Line2D([0], [0], marker="o", ls="", mfc=cmap["PG"],
                  mec=cmap["PG"], label="ĥ PG (short-read)"),
           Line2D([0], [0], marker="o", ls="", mfc="none",
                  mec="#d6604d", label="truth")]
    ax.legend(handles=leg, frameon=False, fontsize=8)

    # (2) log-log est vs truth
    ax = axes[1]
    eps = 1e-6
    ax.scatter(df["truth"] + eps, df["est"] + eps, s=16, c=colors, alpha=0.8)
    lim_lo, lim_hi = eps, max(df["truth"].max(), df["est"].max()) * 1.5
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "--", color="#999999", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("truth h (+1e-6)")
    ax.set_ylabel("kMate ĥ (+1e-6)")
    mb = metrics["mass_balance"]
    ax.set_title("kMate ĥ vs truth (log-log)\n"
                 f"LR mass est {mb['cactus_longread']['est']:.3f} / "
                 f"truth {mb['cactus_longread']['truth']:.3f}   |   "
                 f"SR est {mb['PG_shortread']['est']:.3f} / "
                 f"truth {mb['PG_shortread']['truth']:.3f}")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.suptitle(f"per-founder h benchmark  —  {label}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--est", required=True,
                    help="kMate per-founder h: .npz (founders+chrom) or .tsv/.csv")
    ap.add_argument("--truth", required=True, help="sim pool_weights.tsv")
    ap.add_argument("--split",
                    default="data/founder_split_cactus_pg.json",
                    help="cactus/PG founder split json")
    ap.add_argument("--chrom", default=None,
                    help="restrict npz to one chrom key (default: mean of all)")
    ap.add_argument("--label", default="run", help="run label / output prefix")
    ap.add_argument("--out", default="benchmarks/h_accuracy", help="output dir")
    ap.add_argument("--absorb-thr", type=float, default=1e-3)
    ap.add_argument("--present-thr", type=float, default=1e-3)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    est = load_h_estimate(args.est, chrom=args.chrom)
    truth = load_truth(args.truth)
    split = load_split(args.split)
    metrics, per_founder = score(est, truth, split,
                                 absorb_thr=args.absorb_thr,
                                 present_thr=args.present_thr)
    metrics["label"] = args.label
    metrics["est_path"] = os.path.abspath(args.est)
    metrics["truth_path"] = os.path.abspath(args.truth)

    os.makedirs(args.out, exist_ok=True)
    pfx = os.path.join(args.out, args.label)
    json.dump(metrics, open(pfx + "_metrics.json", "w"), indent=2)
    per_founder.to_csv(pfx + "_per_founder.csv", index=False)
    if not args.no_plot:
        make_figure(per_founder, metrics, args.label, pfx + "_diag.png")

    # ------- console report -------
    m = metrics
    mb = m["mass_balance"]
    print(f"\n=== per-founder h benchmark : {args.label} ===")
    print(f"  founders={m['n_founders']}  present_in_truth={m['n_present_in_truth']}")
    print(f"  RMSE={m['rmse']:.5g}  MAE={m['mae']:.5g}")
    print(f"  Pearson r  lin={m['pearson_linear']:.4f}  log={m['pearson_log']:.4f}")
    print(f"  Spearman r lin={m['spearman_linear']:.4f}  log={m['spearman_log']:.4f}")
    print(f"  ABSORBED (truth>={m['present_thr']:g}, est<{m['absorb_thr']:g}): "
          f"{m['n_absorbed']}")
    if m["absorbed_founders"]:
        ids = ", ".join(a["founder"] for a in m["absorbed_founders"][:20])
        print(f"    ids: {ids}" + (" ..." if m["n_absorbed"] > 20 else ""))
    print("  mass balance (est vs truth):")
    print(f"    cactus/long-read : est {mb['cactus_longread']['est']:.4f}  "
          f"truth {mb['cactus_longread']['truth']:.4f}  "
          f"Δ {mb['cactus_longread']['est_minus_truth']:+.4f}")
    print(f"    PG/short-read    : est {mb['PG_shortread']['est']:.4f}  "
          f"truth {mb['PG_shortread']['truth']:.4f}  "
          f"Δ {mb['PG_shortread']['est_minus_truth']:+.4f}")
    print("  worst over-called:")
    for w in m["worst_overcalled"][:5]:
        print(f"    {w['founder']:>7} [{w['class']}] truth {w['truth']:.4f} "
              f"est {w['est']:.4f}  (+{w['est']-w['truth']:.4f})")
    print("  worst under-called:")
    for w in m["worst_undercalled"][:5]:
        print(f"    {w['founder']:>7} [{w['class']}] truth {w['truth']:.4f} "
              f"est {w['est']:.4f}  ({w['est']-w['truth']:+.4f})")
    print(f"  wrote: {pfx}_metrics.json  {pfx}_per_founder.csv  "
          f"{'(no plot)' if args.no_plot else pfx+'_diag.png'}")


if __name__ == "__main__":
    main()
