#!/usr/bin/env python3
"""Classic truth-vs-estimate density panel: global vs dynld-window, p231 + p80
(outcross g3). 2x2, canonical kMate aesthetic."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # benchmarks/
from _accuracy_panel import grid_panel

D = "benchmarks/ldblock_window_test"
def load(est_tsv, truth_gz):
    o = pd.read_csv(est_tsv, sep="\t"); t = pd.read_csv(truth_gz, sep="\t")
    assert len(o) == len(t), (est_tsv, len(o), len(t))
    snp = ((o.ref_len == 1) & (o.alt_len == 1)).values
    e = o.alt_freq.values.astype(np.float32); tr = t.truth_af.values.astype(np.float32)
    m = snp & np.isfinite(e) & np.isfinite(tr)
    return tr[m], e[m]

TR231 = "benchmarks/p231/sims/cov10_n50_g3_s42_hotspots_p231_chr1/recomb_truth_raw.tsv.gz"
TR80  = "benchmarks/p80/sims/cov10_n50_g3_s42_hotspots_p80_chr1/recomb_truth.tsv.gz"
cells = [
    ("p231 — GLOBAL (1 h / chrom)", *load("benchmarks/p231/results/kmate_global_filt2mb_raw/n50_g3/p231_filt2mb_raw_n50_g3_cov10_s42.tsv", TR231)),
    ("p231 — dynld-WINDOW (per-block h)", *load(f"{D}/cov10_n50_g3_s42_hotspots_p231_chr1_p231_chr1_units_dynld_K500_mkb50.tsv", TR231)),
    ("p80 — GLOBAL (1 h / chrom)", *load("benchmarks/p80/results/cactus_em_global_filt2_mb/n50_g3/p80_filt2mb_n50_g3_cov10_s42.tsv", TR80)),
    ("p80 — dynld-WINDOW (per-block h)", *load(f"{D}/cov10_n50_g3_s42_hotspots_p80_chr1_p80_chr1_units_dynld_K500_mkb50.tsv", TR80)),
]
grid_panel(cells, "", f"{D}/blockpanel_dynld_outcross_g3.png", ncols=2, cell=4.2)
