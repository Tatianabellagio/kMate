#!/usr/bin/env python3
"""kMate accuracy scenario grids — canonical aesthetic, 3x2, raw arm.

Two SCENARIOS, each a self-contained grid with all its regimes (NOT a
side-by-side comparison):
  outcross : [[n231_g0, n50_g0],[n231_g1, n50_g1],[n50_g3, n50_g3_dom500]]
             (forced outcrossing; dom500 winner is a recombinant)
  selfing  : [[n231_g0, n50_g0],[n231_g1_self97, n50_g1_self97],
              [n50_g3_self97, n50_g3_dom500_self97]]
             (97% selfing; g0 == outcross g0 (no recomb); dom500 winner ~pure)

Two MODES: global (kmate_global_*) and window (kmate_window_*, the recombinant
estimator). Missing filter via env KMATE_MISS_THR=0.5 (-> _miss50 suffix).

Usage: plot_scenario_grids.py            # all scenarios x modes x panels
       KMATE_MISS_THR=0.5 plot_scenario_grids.py
"""
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _accuracy_panel import grid_panel, load_cell

ROOT = Path("/global/scratch/users/tbellg/kmate/benchmarks")
_THR = float(os.environ["KMATE_MISS_THR"]) if os.environ.get("KMATE_MISS_THR") else None
SUFFIX = f"_miss{int(_THR*100)}" if _THR is not None else ""
NOTE = f"missing-filtered n_called>={int((1-_THR)*100)}%" if _THR is not None else "all records, no filter"

SCENARIOS = {
    "outcross": [["n231_g0", "n50_g0"], ["n231_g1", "n50_g1"], ["n50_g3", "n50_g3_dom500"]],
    "selfing":  [["n231_g0", "n50_g0"], ["n231_g1_self97", "n50_g1_self97"],
                 ["n50_g3_self97", "n50_g3_dom500_self97"]],
}
# per panel: F, truth filename, est dir templates {mode}, est prefix
# pfx is mode-dependent: the p80 window driver tags files 'filt2mbW' (note W).
PANEL = {
    "p231": dict(F=231, truth="recomb_truth_raw.tsv.gz",
                 pfx={"global": "p231_filt2mb_raw", "window": "p231_filt2mb_raw"},
                 edir={"global": "kmate_global_filt2mb_raw", "window": "kmate_window_filt2mb_raw"},
                 title="arch3 raw: SNP+indel+SV"),
    "p80":  dict(F=80,  truth="recomb_truth.tsv.gz",
                 pfx={"global": "p80_filt2mb", "window": "p80_filt2mbW"},
                 edir={"global": "cactus_em_global_filt2_mb", "window": "cactus_em_window_filt2_mb"},
                 title="var_pa_p80: SNP+indel+SV"),
}


def subdir(reg, panel):
    """regime name -> sim dir basename for this panel."""
    self_tag = "_self97" if reg.endswith("_self97") else ""
    core = reg[:-7] if self_tag else reg          # strip _self97
    if core == "n50_g3_dom500":
        return f"cov10_n50_g3_s42{self_tag}_hotspots_dom500_{panel}_chr1"
    n, g = core.split("_g")
    return f"cov10_{n}_g{g}_s42{self_tag}_hotspots_{panel}_chr1"


def make(panel, scenario, mode):
    p = PANEL[panel]
    grid = SCENARIOS[scenario]
    cells = []
    for row in grid:
        for reg in row:
            tp = ROOT / f"{panel}/sims" / subdir(reg, panel) / p["truth"]
            ep = ROOT / f"{panel}/results" / p["edir"][mode] / reg / f"{p['pfx'][mode]}_{reg}_cov10_s42.tsv"
            if not (tp.exists() and ep.exists()):
                print(f"  [skip cell] {panel} {scenario} {mode} {reg}: missing truth or est")
                cells.append((f"{reg}\n(missing)", [0.0], [0.0])); continue
            cells.append((reg, *load_cell(tp, ep, miss_thr=_THR, F=p["F"])))
    out = ROOT / f"{panel}/results/plots/{panel}_{scenario}_{mode}{SUFFIX}.png"
    grid_panel(cells, f"{panel} kMate filt2 + 1/m_b {mode.upper()} ({p['title']}) — "
               f"FULLY {scenario.upper()} — {NOTE}", out)


if __name__ == "__main__":
    for panel in ("p231", "p80"):
        for scenario in ("outcross", "selfing"):
            for mode in ("global", "window"):
                make(panel, scenario, mode)
