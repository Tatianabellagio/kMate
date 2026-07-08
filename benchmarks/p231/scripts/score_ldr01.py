"""Score the new --unit ld (r2=0.1) p231 AF outputs vs sim truth, by regime x
var-class. Outputs results/kmate_ldr01_summary.tsv. Est TSV and truth share the
arch3 var_pa record order (same length) — aligned by row after a length + pos check."""
import glob, os
import numpy as np, pandas as pd

CTRL = "/global/scratch/users/tbellg/kmate/benchmarks/p231"
SUB = {"n50_g0": "cov10_n50_g0_s42_hotspots_p231_chr1",
       "n231_g0": "cov10_n231_g0_s42_hotspots_p231_chr1",
       "n50_g1": "cov10_n50_g1_s42_hotspots_p231_chr1",
       "n231_g1": "cov10_n231_g1_s42_hotspots_p231_chr1",
       "n50_g3": "cov10_n50_g3_s42_hotspots_p231_chr1",
       "n50_g3_dom500": "cov10_n50_g3_s42_hotspots_dom500_p231_chr1"}

rows = []
for cnvar in ["raw", "atomized"]:
    truth_name = f"recomb_truth_{cnvar}.tsv.gz"
    for reg, sub in SUB.items():
        est_glob = glob.glob(f"{CTRL}/results/kmate_ldr01_{cnvar}/{reg}/*.tsv")
        if not est_glob:
            continue
        est = pd.read_csv(est_glob[0], sep="\t", usecols=["pos", "ref_len", "alt_len", "alt_freq"])
        tr = pd.read_csv(f"{CTRL}/sims/{sub}/{truth_name}", sep="\t",
                         usecols=["pos", "ref_len", "alt_len", "truth_af"])
        assert len(est) == len(tr), f"{reg}/{cnvar}: len {len(est)} != {len(tr)}"
        assert np.array_equal(est["pos"].values, tr["pos"].values), f"{reg}/{cnvar}: pos order mismatch"
        e = est["alt_freq"].values.astype(float); t = tr["truth_af"].values.astype(float)
        rl = tr["ref_len"].values.astype(int); al = tr["alt_len"].values.astype(int)
        cls = np.where(np.maximum(rl, al) >= 50, "SV", np.where((rl == 1) & (al == 1), "SNP", "indel"))
        fin = np.isfinite(e) & np.isfinite(t)
        for c in ["ALL", "SNP", "indel", "SV"]:
            m = fin if c == "ALL" else (fin & (cls == c))
            if m.sum() == 0:
                continue
            d = e[m] - t[m]
            rows.append(dict(var_pa=cnvar, regime=reg, cls=c, n=int(m.sum()),
                             finite_pct=round(100 * fin.mean(), 2),
                             MAE=round(float(np.mean(np.abs(d))), 5),
                             RMSE=round(float(np.sqrt(np.mean(d ** 2))), 5),
                             R2=round(float(np.corrcoef(e[m], t[m])[0, 1] ** 2), 5),
                             outlier_gt10pct=round(float(np.mean(np.abs(d) > 0.10)), 5)))
        del est, tr, e, t

df = pd.DataFrame(rows)
out = f"{CTRL}/results/kmate_ldr01_summary.tsv"
df.to_csv(out, sep="\t", index=False)
print(f"wrote {out}  ({len(df)} rows)\n")
# headline: raw, ALL + classes, per regime
print("=== RAW var_pa, per regime ===")
print(df[(df.var_pa == "raw")].pivot_table(index="regime", columns="cls", values="MAE")
      .reindex(columns=["ALL", "SNP", "indel", "SV"]).round(4).to_string())
print("\n=== RAW var_pa R2, ALL ===")
print(df[(df.var_pa == "raw") & (df.cls == "ALL")][["regime", "n", "finite_pct", "MAE", "R2", "outlier_gt10pct"]].to_string(index=False))
