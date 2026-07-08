#!/usr/bin/env python3
"""kMate accuracy scenario grids for the NEW --unit ld (r2=0.1) production estimator.

Variant of plot_scenario_grids.py pointed at the new kmate_ldr01_* output dirs
(GLOBAL / --unit ld only; window is out of scope for the 2026-07-07 refresh).
Two scenarios (outcross, selfing) x two panels (p231 raw, p80 filt2inv in-house).

Usage: plot_scenario_grids_ldr01.py            # all records
       KMATE_MISS_THR=0.5 plot_scenario_grids_ldr01.py    # missing-filtered
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
# NEW estimator: --unit ld --ld-r2 0.1 + per_founder + uniform, production in-house panel.
PANEL = {
    "p231": dict(F=231, truth="recomb_truth_raw.tsv.gz",
                 edir="kmate_ldr01_raw", pfx="p231_ldr01_raw",
                 title="prod in-house filt2inv, --unit ld r2=0.1: SNP+indel+SV"),
    "p80":  dict(F=80,  truth="recomb_truth.tsv.gz",
                 edir="kmate_ldr01_p80_filt2inv", pfx="p80_ldr01_filt2inv",
                 title="p80 in-house filt2inv, --unit ld r2=0.1: SNP+indel+SV"),
}


def subdir(reg, panel):
    self_tag = "_self97" if reg.endswith("_self97") else ""
    core = reg[:-7] if self_tag else reg
    if core == "n50_g3_dom500":
        return f"cov10_n50_g3_s42{self_tag}_hotspots_dom500_{panel}_chr1"
    n, g = core.split("_g")
    return f"cov10_{n}_g{g}_s42{self_tag}_hotspots_{panel}_chr1"


def make(panel, scenario):
    p = PANEL[panel]
    cells = []
    for row in SCENARIOS[scenario]:
        for reg in row:
            tp = ROOT / f"{panel}/sims" / subdir(reg, panel) / p["truth"]
            ep = ROOT / f"{panel}/results" / p["edir"] / reg / f"{p['pfx']}_{reg}_cov10_s42.tsv"
            if not (tp.exists() and ep.exists()):
                print(f"  [skip cell] {panel} {scenario} {reg}: missing truth or est ({ep.name})")
                cells.append((f"{reg}\n(missing)", [0.0], [0.0])); continue
            cells.append((reg, *load_cell(tp, ep, miss_thr=_THR, F=p["F"])))
    out = ROOT / f"{panel}/results/plots/{panel}_{scenario}_ldr01{SUFFIX}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    grid_panel(cells, f"{panel} kMate --unit ld r2=0.1 ({p['title']}) — "
               f"FULLY {scenario.upper()} — {NOTE}", out)
    print(f"  wrote {out}")


if __name__ == "__main__":
    for panel in ("p231", "p80"):
        for scenario in ("outcross", "selfing"):
            make(panel, scenario)
