"""Combine persample_SEEDMIX_S*.npz -> allsamples_h.npz + summary.json.
Aggregate h = mean over the 5 chromosomes (per sample), old vs new Kf_w."""
import json, glob, os
import numpy as np

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/seedmix_validation/fix_kfw_fullpanel"
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
    Ho = d["perchrom_old"]; Hn = d["perchrom_new"]     # 5 x F
    save[f"perchrom_old__{sid}"] = Ho
    save[f"perchrom_new__{sid}"] = Hn
    for tag, Hm in [("old", Ho), ("new", Hn)]:
        hbar = Hm.mean(0); hbar = hbar / hbar.sum()
        save[f"agg_{tag}__{sid}"] = hbar
        summ[f"{sid}_{tag}"] = dict(n_abs=int((hbar < 1e-3).sum()),
            eff_n=float(1/np.sum(hbar**2)),
            rmse_vs_uniform=float(np.sqrt(np.mean((hbar-u)**2))),
            cac_ratio=float(hbar[is_cac].sum()/is_cac.mean()))
    samples.append(sid)

save["founders"] = fo0; save["chroms"] = np.array(CHROMS)
save["expected"] = np.full(len(fo0), 1.0/len(fo0))
np.savez_compressed(f"{OUT}/allsamples_h.npz", **save)
json.dump(summ, open(f"{OUT}/allsamples_summary.json", "w"), indent=2)
print(f"aggregated {len(samples)} samples: {samples}")
print("\n=== aggregate (chrom-mean) n_absorbed per sample: old -> new ===")
for sid in samples:
    print(f"  {sid}: old {summ[f'{sid}_old']['n_abs']:3d} -> new {summ[f'{sid}_new']['n_abs']:3d}  "
          f"| rmse_vs_1/231 {summ[f'{sid}_old']['rmse_vs_uniform']:.2e} -> {summ[f'{sid}_new']['rmse_vs_uniform']:.2e}")
