#!/usr/bin/env python3
"""hapFIRE vs kMate on the SAME p231 pool-size x depth grid (N in
{2,5,20,50,150}, depth in {1,10}, seeds 42-46) -- accuracy (h, AF) and
speed/compute, side by side. p231 ONLY (hapFIRE's panel is the greneNet
231-founder SNP set; there is no equivalent panel for p80).

Accuracy:
  - h (founder-mixture) accuracy: both tools scored against the same
    pool_weights.tsv truth kMate's poolsize_depth sweep uses.
  - AF accuracy: hapFIRE has no indel/SV calls (SNP-only VCF panel) and its
    snp_frequency.txt carries no REF/ALT identity, so multi-allelic truth
    positions (ambiguous without allele identity) are dropped. Both tools are
    scored on the SAME bi-allelic-SNP-only, hapFIRE-covered position set --
    this makes kMate's AF numbers here narrower than (and not comparable to)
    poolsize_depth_table.tsv's whole-panel af_R2/af_RMSE.

Speed/compute:
  - kMate: benchmarks/poolsize_depth/results/kmate_speed_table.tsv (sacct
    Elapsed/MaxRSS, no rerun needed).
  - hapFIRE: align_time.txt + hapfire_time.txt (/usr/bin/time -v) per
    condition, reported both separately (align vs hapFIRE-proper) and summed
    (fastq-in to AF-out, the fair comparison against kMate's single fastq-to-AF
    step, since kMate is alignment-free and hapFIRE is not).

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/score_hapfire_vs_kmate.py
"""
import os, re, glob
import numpy as np
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
P231 = f"{ROOT}/benchmarks/p231"
HF_RES = f"{ROOT}/benchmarks/speed_vs_hapfire/results/poolsize_depth"
OUT = f"{ROOT}/benchmarks/speed_vs_hapfire/results"
os.makedirs(OUT, exist_ok=True)

POOL_SIZES = [2, 5, 20, 50, 150]
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
    # RMSE/sd(truth) -- matches the original paper's "Normalized RMSE" convention
    # (Fig S13). Moves in the SAME direction as R2 (both relative to truth's own
    # spread), unlike plain RMSE which shrinks in absolute terms as the target
    # itself shrinks toward 0 with growing N -- same degenerate guard as R2.
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


kmate_tbl = pd.read_csv(f"{ROOT}/benchmarks/poolsize_depth/results/poolsize_depth_table.tsv", sep="\t")
kmate_tbl = kmate_tbl[kmate_tbl.panel == "p231"]
kmate_speed = pd.read_csv(f"{ROOT}/benchmarks/poolsize_depth/results/kmate_speed_table.tsv", sep="\t")
kmate_speed = kmate_speed[kmate_speed.panel == "p231"]

rows = []
for n in POOL_SIZES:
    for cov in DEPTHS:
        for seed in SEEDS:
            sim = f"{P231}/sims/cov{cov}_n{n}_g0_s{seed}_hotspots_p231_chr1"
            hf_dir = f"{HF_RES}/n{n}_cov{cov}_s{seed}"
            hf_sample = f"p231_hapfire_psd_n{n}_cov{cov}_s{seed}"
            kmate_dir = f"{P231}/results/kmate_chrom_poolsize_depth/n{n}_cov{cov}_s{seed}"
            kmate_sample = f"p231_chrom_psd_n{n}_cov{cov}_s{seed}"

            ef_path = f"{hf_dir}/{hf_sample}_ecotype_frequency.txt"
            sf_path = f"{hf_dir}/{hf_sample}_snp_frequency.txt"
            kmate_tsv = f"{kmate_dir}/{kmate_sample}.tsv"
            kmate_npz = f"{kmate_dir}/{kmate_sample}.h_per_chrom.npz"
            if not all(os.path.exists(p) for p in (ef_path, sf_path, kmate_tsv, kmate_npz)):
                print(f"  [skip] n={n} cov={cov} s={seed}: missing output")
                continue

            truth_w = pd.read_csv(f"{sim}/pool_weights.tsv", sep="\t")
            truth_w["founder"] = truth_w["founder"].astype(str)
            tw = dict(zip(truth_w["founder"], truth_w["weight"]))

            # --- h accuracy: hapFIRE ---
            ef = pd.read_csv(ef_path, sep="\t", header=None, names=["founder", "freq"])
            ef["founder"] = ef["founder"].astype(str)
            ef["truth"] = ef["founder"].map(lambda f: tw.get(f, 0.0))
            hf_h_r2, hf_h_rmse, hf_h_rmse_norm = r2_rmse(ef["freq"].values, ef["truth"].values)

            # --- h accuracy: kMate (same truth) ---
            hd = np.load(kmate_npz, allow_pickle=True)
            fo, h = hd["founders"].astype(str), hd["Chr1"].astype(float)
            truth_k = np.array([tw.get(f, 0.0) for f in fo])
            km_h_r2, km_h_rmse, km_h_rmse_norm = r2_rmse(h, truth_k)

            # --- AF accuracy: bi-allelic-SNP-only, common to both tools ---
            truth_af = pd.read_csv(f"{sim}/recomb_truth_raw.tsv.gz", sep="\t").dropna(subset=["truth_af"])
            is_snp = (truth_af.ref_len == 1) & (truth_af.alt_len == 1)
            truth_snp = truth_af[is_snp]
            truth_snp_uniq = truth_snp[~truth_snp.pos.duplicated(keep=False)][["pos", "truth_af"]]

            hf_snp = pd.read_csv(sf_path, sep="\t", header=None, names=["chrom", "pos", "freq"])
            common = truth_snp_uniq.merge(hf_snp[["pos", "freq"]], on="pos", how="inner")

            kmate_est = pd.read_csv(kmate_tsv, sep="\t")
            kmate_snp = kmate_est[(kmate_est.ref_len == 1) & (kmate_est.alt_len == 1)][["pos", "alt_freq"]]
            common = common.merge(kmate_snp, on="pos", how="inner")

            hf_af_r2, hf_af_rmse, _ = r2_rmse(common["freq"].values, common["truth_af"].values)
            km_af_r2, km_af_rmse, _ = r2_rmse(common["alt_freq"].values, common["truth_af"].values)
            n_common = len(common)

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
                N=n, depth=cov, seed=seed,
                kmate_h_R2=km_h_r2, kmate_h_RMSE=km_h_rmse, kmate_h_RMSE_norm=km_h_rmse_norm,
                hapfire_h_R2=hf_h_r2, hapfire_h_RMSE=hf_h_rmse, hapfire_h_RMSE_norm=hf_h_rmse_norm,
                kmate_af_R2=km_af_r2, kmate_af_RMSE=km_af_rmse,
                hapfire_af_R2=hf_af_r2, hapfire_af_RMSE=hf_af_rmse,
                af_n_common=n_common,
                kmate_elapsed_s=km_elapsed, kmate_cpu_s=km_cpu, kmate_max_rss_mb=km_rss,
                hapfire_align_elapsed_s=align_t["elapsed_s"], hapfire_align_cpu_s=align_t["cpu_s"],
                hapfire_proper_elapsed_s=hapfire_t["elapsed_s"], hapfire_proper_cpu_s=hapfire_t["cpu_s"],
                hapfire_total_elapsed_s=hf_total_elapsed, hapfire_total_cpu_s=hf_total_cpu,
                hapfire_max_rss_mb=hf_peak_rss,
            ))
            print(f"n={n:>3} cov={cov:>2}x s={seed}: "
                  f"h R2 kmate={km_h_r2:.3f} hapfire={hf_h_r2:.3f} | "
                  f"AF(SNP,n={n_common:,}) R2 kmate={km_af_r2:.4f} hapfire={hf_af_r2:.4f} | "
                  f"wall kmate={km_elapsed:.0f}s hapfire={hf_total_elapsed:.0f}s "
                  f"(align={align_t['elapsed_s']:.0f}+hf={hapfire_t['elapsed_s']:.0f})")

df = pd.DataFrame(rows)
table_path = f"{OUT}/hapfire_vs_kmate_table.tsv"
df.to_csv(table_path, sep="\t", index=False)
print(f"\nwrote {table_path} ({len(df)} rows)")
