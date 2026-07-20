"""Combine persample_SEEDMIX_S*.npz -> seedmix_allsamples_h.npz + summary.json
(the format plot_seedmix.py expects). Aggregate h = mean over the 5 chromosomes."""
import json, glob, os
import numpy as np

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/seedmix_validation/fix_seedmix"
SPLIT = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CACTUS = set(map(str, SPLIT["cactus"]))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]

files = sorted(glob.glob(f"{OUT}/persample_SEEDMIX_S*.npz"))
if not files:
    raise SystemExit("no persample npz found")
fo0 = None; save = {}; summ = {}
samples = []
for f in files:
    sid = os.path.basename(f).replace("persample_", "").replace(".npz", "")
    d = np.load(f, allow_pickle=True)
    fo = d["founders"].astype(str)
    if fo0 is None:
        fo0 = fo; F = len(fo); u = 1.0 / F
        is_cac = np.array([x in CACTUS for x in fo])
    assert np.array_equal(fo, fo0)
    Hp = d["perchrom_prod"]; Hf = d["perchrom_fix"]     # 5 x F
    save[f"perchrom_prod__{sid}"] = Hp
    save[f"perchrom_fix__{sid}"] = Hf
    for tag, Hm in [("prod", Hp), ("fix", Hf)]:
        hbar = Hm.mean(0); hbar = hbar / hbar.sum()
        save[f"agg_{tag}__{sid}"] = hbar
        summ[f"{sid}_{tag}"] = dict(n_abs=int((hbar < 1e-3).sum()),
            eff_n=float(1/np.sum(hbar**2)),
            rmse_vs_uniform=float(np.sqrt(np.mean((hbar-u)**2))),
            cac_ratio=float(hbar[is_cac].sum()/is_cac.mean()))
    samples.append(sid)

save["founders"] = fo0; save["chroms"] = np.array(CHROMS)
save["expected"] = np.full(len(fo0), 1.0/len(fo0))
np.savez_compressed(f"{OUT}/seedmix_allsamples_h.npz", **save)
json.dump(summ, open(f"{OUT}/seedmix_allsamples_summary.json", "w"), indent=2)
print(f"aggregated {len(samples)} samples: {samples}")
print("\n=== aggregate (chrom-mean) n_absorbed per sample: prod -> fix ===")
for sid in samples:
    print(f"  {sid}: prod {summ[f'{sid}_prod']['n_abs']:3d} -> fix {summ[f'{sid}_fix']['n_abs']:3d}  "
          f"| rmse_vs_1/231 {summ[f'{sid}_prod']['rmse_vs_uniform']:.2e} -> {summ[f'{sid}_fix']['rmse_vs_uniform']:.2e}")
