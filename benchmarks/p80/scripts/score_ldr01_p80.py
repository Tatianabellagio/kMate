"""Score p80 --unit ld outputs (raw vs filt2inv arms) vs sim truth, by regime x
class. The neutrality test: on the homogeneous 80-cactus panel the private-k-mer
(ac=1) drop should be ~neutral, so raw ≈ filt2inv. Writes kmate_ldr01_p80_summary.tsv."""
import glob
import numpy as np, pandas as pd

CTRL = "/global/scratch/users/tbellg/kmate/benchmarks/p80"
SUB = {"n50_g0": "cov10_n50_g0_s42_hotspots_p80_chr1",
       "n231_g0": "cov10_n231_g0_s42_hotspots_p80_chr1",
       "n50_g1": "cov10_n50_g1_s42_hotspots_p80_chr1",
       "n231_g1": "cov10_n231_g1_s42_hotspots_p80_chr1",
       "n50_g3": "cov10_n50_g3_s42_hotspots_p80_chr1",
       "n50_g3_dom500": "cov10_n50_g3_s42_hotspots_dom500_p80_chr1"}

rows = []
for arm in ["raw", "filt2inv"]:
    for reg, sub in SUB.items():
        g = glob.glob(f"{CTRL}/results/kmate_ldr01_p80_{arm}/{reg}/*.tsv")
        if not g:
            continue
        est = pd.read_csv(g[0], sep="\t", usecols=["pos", "ref_len", "alt_len", "alt_freq"])
        tr = pd.read_csv(f"{CTRL}/sims/{sub}/recomb_truth.tsv.gz", sep="\t",
                         usecols=["pos", "ref_len", "alt_len", "truth_af"])
        assert len(est) == len(tr) and np.array_equal(est["pos"].values, tr["pos"].values), \
            f"{arm}/{reg}: order/len mismatch"
        e = est["alt_freq"].values.astype(float); t = tr["truth_af"].values.astype(float)
        rl = tr["ref_len"].values.astype(int); al = tr["alt_len"].values.astype(int)
        cls = np.where(np.maximum(rl, al) >= 50, "SV", np.where((rl == 1) & (al == 1), "SNP", "indel"))
        fin = np.isfinite(e) & np.isfinite(t)
        for c in ["ALL", "SNP", "indel", "SV"]:
            m = fin if c == "ALL" else (fin & (cls == c))
            if m.sum() == 0:
                continue
            d = e[m] - t[m]
            rows.append(dict(arm=arm, regime=reg, cls=c, n=int(m.sum()),
                             finite_pct=round(100 * fin.mean(), 2),
                             MAE=round(float(np.mean(np.abs(d))), 5),
                             RMSE=round(float(np.sqrt(np.mean(d ** 2))), 5),
                             R2=round(float(np.corrcoef(e[m], t[m])[0, 1] ** 2), 5)))
        del est, tr

df = pd.DataFrame(rows)
out = f"{CTRL}/results/kmate_ldr01_p80_summary.tsv"
df.to_csv(out, sep="\t", index=False)
print(f"wrote {out}\n")
piv = df[df.cls == "ALL"].pivot_table(index="regime", columns="arm", values="MAE")
piv["Δ(filt2inv-raw)"] = (piv.get("filt2inv") - piv.get("raw")).round(5)
print("=== ALL-class MAE: raw vs filt2inv (neutrality test) ===")
print(piv.round(5).to_string())
print("\n=== per-arm ALL finite% + R2 ===")
print(df[df.cls == "ALL"][["arm", "regime", "finite_pct", "MAE", "R2"]].to_string(index=False))
