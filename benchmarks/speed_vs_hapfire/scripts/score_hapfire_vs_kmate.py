#!/usr/bin/env python3
"""hapFIRE vs kMate, each on its OWN native panel (the corrected, fair design --
see ../PANEL_MISMATCH_BUG.md). kMate: arch3 panel, arch3-simulated reads
(existing benchmarks/p231 results, no rerun). hapFIRE: greneNet_final_v1.1
panel, reads simulated from greneNet-derived founder FASTAs
(sims_greneNet/, results/greneNet_fair/). Same founders-meta + seed for both ->
identical pool composition (founder draw), so h-accuracy is a fair, matched-pool
comparison. Same grid as kMate's poolsize_depth sweep: N in {2,5,20,50,150,231} x
depth in {1,10}, plus N=50 x depth in {30,50}, x seeds {42-46}.

Scope: h (founder-mixture) accuracy + speed/compute only. AF (per-variant)
accuracy is NOT computed here -- kMate's arch3 truth and hapFIRE's greneNet
truth come from two different variant-calling lineages (arch3 = full-assembly
SNP+indel+SV calls, greneNet = independent short-read SNP-only calls), so a
"true AF" at a shared position isn't a single unambiguous number across panels.
The old (pre-fix) hapfire_vs_kmate_table.tsv's af_R2/af_RMSE columns depended on
hapFIRE and kMate sharing the SAME reads against a common truth, which no
longer holds once each tool gets its own native-panel reads. Re-deriving a
greneNet-lineage AF truth (a var_pa-style matrix from the greneNet VCF + a
compute_recomb_truth.py-style projection) is possible but out of scope here.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/score_hapfire_vs_kmate.py
"""
import os, re
import numpy as np
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
BASE = f"{ROOT}/benchmarks/speed_vs_hapfire"
P231 = f"{ROOT}/benchmarks/p231"
HF_SIMS = f"{BASE}/sims_greneNet"
HF_RES = f"{BASE}/results/greneNet_fair"
OUT = f"{BASE}/results"
os.makedirs(OUT, exist_ok=True)

POOL_SIZES = [2, 5, 20, 50, 150, 231]
# 30x/50x only exist at N=50 (coverage-gradient diagnostic); other N x {30,50}
# combos will just [skip] as missing output, which is fine.
DEPTHS = [1, 10, 30, 50]
SEEDS = [42, 43, 44, 45, 46]
REL_EPS = 1e-6


def r2_rmse(est, truth):
    d = est - truth
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    degenerate = ss_tot < (REL_EPS * truth.mean()) ** 2 * len(truth)
    r2 = np.nan if degenerate else 1 - np.sum(d ** 2) / ss_tot
    rmse = float(np.sqrt(np.mean(d ** 2)))
    sd = truth.std()
    rmse_norm = np.nan if degenerate else rmse / sd
    return float(r2), rmse, float(rmse_norm)


TIME_RE = {
    "elapsed": re.compile(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)"),
    "user": re.compile(r"User time \(seconds\):\s*([\d.]+)"),
    "sys": re.compile(r"System time \(seconds\):\s*([\d.]+)"),
    "maxrss_kb": re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)"),
}


def parse_time_v(path):
    text = open(path).read()
    out = {}
    em = TIME_RE["elapsed"].search(text)
    if em:
        parts = [float(p) for p in em.group(1).split(":")]
        while len(parts) < 3:
            parts.insert(0, 0)
        h, m, s = parts[-3:]
        out["elapsed_s"] = h * 3600 + m * 60 + s
    um, sm, rm = (TIME_RE[k].search(text) for k in ("user", "sys", "maxrss_kb"))
    out["cpu_s"] = (float(um.group(1)) if um else 0) + (float(sm.group(1)) if sm else 0)
    out["max_rss_mb"] = float(rm.group(1)) / 1024 if rm else None
    return out


kmate_speed = pd.read_csv(f"{ROOT}/benchmarks/poolsize_depth/results/kmate_speed_table.tsv", sep="\t")
kmate_speed = kmate_speed[kmate_speed.panel == "p231"]

rows = []
for n in POOL_SIZES:
    for cov in DEPTHS:
        for seed in SEEDS:
            hf_dir = f"{HF_RES}/n{n}_cov{cov}_s{seed}"
            hf_sample = f"greneNet_hapfire_n{n}_cov{cov}_s{seed}"
            kmate_dir = f"{P231}/results/kmate_chrom_poolsize_depth/n{n}_cov{cov}_s{seed}"
            kmate_sample = f"p231_chrom_psd_n{n}_cov{cov}_s{seed}"
            truth_path = f"{HF_SIMS}/cov{cov}_n{n}_g0_s{seed}_greneNet_chr1/pool_weights.tsv"

            ef_path = f"{hf_dir}/{hf_sample}_ecotype_frequency.txt"
            kmate_npz = f"{kmate_dir}/{kmate_sample}.h_per_chrom.npz"
            if not all(os.path.exists(p) for p in (ef_path, kmate_npz, truth_path)):
                print(f"  [skip] n={n} cov={cov} s={seed}: missing output")
                continue

            truth_w = pd.read_csv(truth_path, sep="\t")
            truth_w["founder"] = truth_w["founder"].astype(str)
            tw = dict(zip(truth_w["founder"], truth_w["weight"]))
            true_set = set(truth_w[truth_w.weight > 0]["founder"])
            n_truth_founders = len(true_set)

            # "founders recovered" = overlap between the true pool members and
            # the top-K estimated founders (K = n_truth_founders). Neither
            # tool's raw estimate is a clean zero/nonzero mask -- both spread a
            # small residual weight across off-pool founders (LD/relatedness),
            # so top-K-vs-truth (precision at the true cardinality) is the
            # threshold-free way to count recovery + false positives.

            # --- h accuracy: hapFIRE (greneNet panel, greneNet-sim truth) ---
            ef = pd.read_csv(ef_path, sep="\t", header=None, names=["founder", "freq"])
            ef["founder"] = ef["founder"].astype(str)
            ef["truth"] = ef["founder"].map(lambda f: tw.get(f, 0.0))
            hf_h_r2, hf_h_rmse, hf_h_rmse_norm = r2_rmse(ef["freq"].values, ef["truth"].values)
            hf_top = set(ef.sort_values("freq", ascending=False).head(n_truth_founders)["founder"])
            hf_n_found = len(hf_top & true_set)

            # --- h accuracy: kMate (arch3 panel, arch3-sim truth) ---
            hd = np.load(kmate_npz, allow_pickle=True)
            fo, h = hd["founders"].astype(str), hd["Chr1"].astype(float)
            truth_k = np.array([tw.get(f, 0.0) for f in fo])
            km_h_r2, km_h_rmse, km_h_rmse_norm = r2_rmse(h, truth_k)
            km_top = set(fo[np.argsort(-h)[:n_truth_founders]])
            km_n_found = len(km_top & true_set)

            # --- speed/compute ---
            align_t = parse_time_v(f"{hf_dir}/align_time.txt")
            hapfire_t = parse_time_v(f"{hf_dir}/hapfire_time.txt")
            hf_total_elapsed = align_t["elapsed_s"] + hapfire_t["elapsed_s"]
            hf_total_cpu = align_t["cpu_s"] + hapfire_t["cpu_s"]
            hf_peak_rss = max(align_t["max_rss_mb"] or 0, hapfire_t["max_rss_mb"] or 0)

            ks = kmate_speed[(kmate_speed.N == n) & (kmate_speed.depth == cov) & (kmate_speed.seed == seed)]
            km_elapsed = ks["elapsed_s"].iloc[0] if len(ks) else None
            km_cpu = ks["cpu_s"].iloc[0] if len(ks) else None
            km_rss = ks["max_rss_mb"].iloc[0] if len(ks) else None

            rows.append(dict(
                N=n, depth=cov, seed=seed, n_truth_founders=n_truth_founders,
                kmate_h_R2=km_h_r2, kmate_h_RMSE=km_h_rmse, kmate_h_RMSE_norm=km_h_rmse_norm,
                kmate_n_found=km_n_found,
                hapfire_h_R2=hf_h_r2, hapfire_h_RMSE=hf_h_rmse, hapfire_h_RMSE_norm=hf_h_rmse_norm,
                hapfire_n_found=hf_n_found,
                kmate_elapsed_s=km_elapsed, kmate_cpu_s=km_cpu, kmate_max_rss_mb=km_rss,
                hapfire_align_elapsed_s=align_t["elapsed_s"], hapfire_align_cpu_s=align_t["cpu_s"],
                hapfire_proper_elapsed_s=hapfire_t["elapsed_s"], hapfire_proper_cpu_s=hapfire_t["cpu_s"],
                hapfire_total_elapsed_s=hf_total_elapsed, hapfire_total_cpu_s=hf_total_cpu,
                hapfire_max_rss_mb=hf_peak_rss,
            ))
            print(f"n={n:>3} cov={cov:>2}x s={seed}: "
                  f"h R2 kmate={km_h_r2:.3f} hapfire={hf_h_r2:.3f} | "
                  f"found kmate={km_n_found}/{n_truth_founders} hapfire={hf_n_found}/{n_truth_founders} | "
                  f"wall kmate={km_elapsed:.0f}s hapfire={hf_total_elapsed:.0f}s "
                  f"(align={align_t['elapsed_s']:.0f}+hf={hapfire_t['elapsed_s']:.0f})")

df = pd.DataFrame(rows)
table_path = f"{OUT}/hapfire_vs_kmate_table.tsv"
df.to_csv(table_path, sep="\t", index=False)
print(f"\nwrote {table_path} ({len(df)} rows)")
