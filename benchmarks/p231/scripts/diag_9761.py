#!/usr/bin/env python3
"""Test why founder 9761 is the EM identifiability sink, on the raw p231 kmer_pa.

TEST 1  most-common-kmers: per founder, the average 'commonness' (# founders
        sharing each of its k-mers) = (K @ ac)/Kf. Rank 9761.
TEST 2  NNLS sink: reconstruct each founder's k-mer profile as a non-negative
        blend of the other 230 (via the 231x231 overlap Gram). Rank by R2.
"""
import numpy as np, json, time
from scipy.sparse import load_npz
from scipy.optimize import nnls

ROOT = "/global/scratch/users/tbellg/kmate"
pre  = f"{ROOT}/benchmarks/p231/data/kmer_pa_p231/kmer_pa"
fo   = np.asarray(np.load(f"{ROOT}/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz",
                          allow_pickle=True)["founders"]).astype(str)
LONG = set(map(str, json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
i = list(fo).index("9761")

t = time.time(); print("loading K (csc, float32)...", flush=True)
K = load_npz(f"{pre}_Chr1.kmer_pa.npz").astype(np.float32).tocsc()
F, Kn = K.shape
ac = np.asarray(K.sum(axis=0)).ravel()        # carriers per k-mer
Kf = np.asarray(K.sum(axis=1)).ravel()        # k-mers per founder
print(f"  K={K.shape}  nnz={K.nnz:,}  ({time.time()-t:.0f}s)", flush=True)

# ---- TEST 1: most-common-kmers (no Gram needed) ----
Kcsr = K.tocsr()
mean_ac = (Kcsr @ ac) / Kf
order = np.argsort(mean_ac)[::-1]
rk = int(np.where(order == i)[0][0])
print("\n[TEST 1] mean k-mer commonness = avg #founders sharing each of a founder's k-mers")
print(f"  9761: mean_ac={mean_ac[i]:.2f}   rank {rk+1}/{F} (1 = most common)")
print(f"  median={np.median(mean_ac):.2f}  max={mean_ac.max():.2f}")
print(f"  TOP 5 most-common: {[(fo[j], round(float(mean_ac[j]),1), 'LR' if fo[j] in LONG else 'SR') for j in order[:5]]}")
print(f"  9761 percentile: {100*(mean_ac < mean_ac[i]).mean():.0f}th")

# ---- Gram G = K K^T (231x231), block-wise over columns ----
t = time.time(); print("\ncomputing overlap Gram block-wise...", flush=True)
G = np.zeros((F, F), dtype=np.float64)
B = 2_000_000
for s in range(0, Kn, B):
    blk = K[:, s:s+B].toarray()            # 231 x b  float32
    G += blk @ blk.T
print(f"  Gram done ({time.time()-t:.0f}s)", flush=True)

# ---- TEST 2: NNLS reconstructability (sink) ----
def recon_R2(f):
    oth = [j for j in range(F) if j != f]
    Goo = G[np.ix_(oth, oth)] + 1e-6*np.eye(F-1)
    g = G[oth, f]; G99 = G[f, f]
    U = np.linalg.cholesky(Goo).T          # Goo = U^T U
    b = np.linalg.solve(U.T, g)
    w, _ = nnls(U, b)
    resid = G99 - 2*g@w + w@Goo@w
    return 1 - resid/G99

R2 = np.array([recon_R2(f) for f in range(F)])
order2 = np.argsort(R2)[::-1]
rk2 = int(np.where(order2 == i)[0][0])
print("\n[TEST 2] NNLS reconstructability by the other 230 founders (sink = high R2)")
print(f"  9761: R2={R2[i]:.4f}   rank {rk2+1}/{F} (1 = most reconstructable = worst sink)")
print(f"  median R2={np.median(R2):.4f}  max={R2.max():.4f}")
print(f"  TOP 5 sinks: {[(fo[j], round(float(R2[j]),4), 'LR' if fo[j] in LONG else 'SR') for j in order2[:5]]}")
print(f"  9761 percentile: {100*(R2 < R2[i]).mean():.0f}th")

# cross-check: correlation of the two rankings with the observed spurious over-credit isn't
# computed here, but print both metrics for 9761 vs the panel for the writeup.
print("\nsummary: 9761  mean_ac rank", rk+1, " | NNLS-sink rank", rk2+1, "of", F)
