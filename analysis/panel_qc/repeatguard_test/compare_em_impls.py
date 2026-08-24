import numpy as np, scipy.sparse as sp, json, pandas as pd, sys, time
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/src")
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/qc/seedmix_validation/fix_norm")
from kmate.em_solver import solve_em as prod_solve_em
from em_variants import solve_em as proto_solve_em

ROOT = "/global/scratch/users/tbellg/kmate"
RAW = f"{ROOT}/benchmarks/p231/data/kmer_pa_p231/kmer_pa_Chr1"
SIM = f"{ROOT}/benchmarks/p231/sims/cov10_n231_g0_s42_hotspots_p231_chr1"
COUNTS_CACHE = f"{ROOT}/analysis/grenenet_selection/qc/seedmix_validation/fix_ablation/raw_cov10_n231_g0_counts.npy"

t0 = time.time()
meta = np.load(f"{RAW}.meta.npz", allow_pickle=True)
FO = np.asarray(meta["founders"]).astype(str)
counts = np.load(COUNTS_CACHE)
K = sp.load_npz(f"{RAW}.kmer_pa.npz").astype(np.float32).tocsr()
F, Kn = K.shape
ac = np.asarray(K.sum(axis=0)).ravel().astype(np.int32)
keep = (ac >= 2) & (ac <= F - 1)   # FILT2INV
Kf = K.tocsc()[:, keep].tocsr()
cf = counts[keep].astype(np.float32)
print(f"F={F} Kn_filt2inv={keep.sum():,} nonzero_counts={(cf>0).sum():,}/{len(cf):,}  [{time.time()-t0:.0f}s load]", flush=True)

pw = pd.read_csv(f"{SIM}/pool_weights.tsv", sep="\t")
wmap = dict(zip(pw["founder"].astype(str), pw["weight"]))
h_true = np.array([wmap.get(f, 0.0) for f in FO]); h_true /= h_true.sum()

t1 = time.time()
h_proto, info_proto = proto_solve_em(cf, Kf, mode="poisson", omega=None, max_iter=300, tol=1e-7)
n_abs_proto = int((h_proto < 1e-3).sum())
print(f"PROTOTYPE (full-panel Kf_w, unfiltered-to-nz K): n_absorbed={n_abs_proto}  iters={info_proto['iters']}  "
      f"h_rmse={np.sqrt(np.mean((h_proto-h_true)**2)):.5f}  [{time.time()-t1:.0f}s]", flush=True)

t2 = time.time()
h_prod_full, info_pf = prod_solve_em(cf, Kf, 1.0, omega=None, normalize="per_founder", max_iter=300, tol=1e-7)
n_abs_pf = int((h_prod_full < 1e-3).sum())
print(f"PRODUCTION code, NOT pre-filtered to nz: n_absorbed={n_abs_pf}  h_rmse={np.sqrt(np.mean((h_prod_full-h_true)**2)):.5f}  [{time.time()-t2:.0f}s]", flush=True)

t3 = time.time()
nz = cf > 0
h_prod_real, info_pr = prod_solve_em(cf[nz], Kf[:, nz], 1.0, omega=None, normalize="per_founder", max_iter=300, tol=1e-7)
n_abs_pr = int((h_prod_real < 1e-3).sum())
print(f"PRODUCTION code, pre-filtered to nz (REAL pipeline behavior): n_absorbed={n_abs_pr}  h_rmse={np.sqrt(np.mean((h_prod_real-h_true)**2)):.5f}  [{time.time()-t3:.0f}s]", flush=True)

print(f"\nmax|h_proto - h_prod_full| = {np.max(np.abs(h_proto-h_prod_full)):.2e}", flush=True)
print(f"max|h_prod_full - h_prod_real| = {np.max(np.abs(h_prod_full-h_prod_real)):.2e}", flush=True)
print(f"TOTAL {time.time()-t0:.0f}s", flush=True)
