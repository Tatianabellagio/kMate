import numpy as np, scipy.sparse as sp, sys
sys.path.insert(0, "src")
import pandas as pd
from kmate.em_solver import haploblock_collapse_indices
from kmate.block_em import BlockSpec, assign_kmers_to_blocks

pfx = "data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
K = sp.load_npz(pfx + ".kmer_pa.npz").tocsc()
m = np.load(pfx + ".meta.npz", allow_pickle=True)
bid = np.asarray(m["bubble_id"]).astype(np.int64)
bchrom = np.asarray(m["bubble_chrom"]).astype(str)
bstart = np.asarray(m["bubble_start"]).astype(np.int64)
bend = np.asarray(m["bubble_end"]).astype(np.int64)
bt = pd.read_csv("benchmarks/p231/data/ld_blocks_r2_0.10_Chr1.tsv", sep="\t")
blocks = [BlockSpec("Chr1", int(s), int(e)) for s, e in zip(bt.start_pos, bt.end_pos)]
kmer_block = assign_kmers_to_blocks(bid, bchrom, bstart, bend, blocks)
d = np.load("benchmarks/p231/results/kmate_ldr01_raw/n231_g0/"
            "p231_ldr01_raw_n231_g0_cov10_s42.h_blocks_per_chrom.npz", allow_pickle=True)
hb = d["Chr1_h_blocks"]
print("blk region(Mb)     n_kmers   K_b  K_b/231  informative(ac2..229)  fit_effn  verdict")
for b in range(len(blocks)):
    cols = np.flatnonzero(kmer_block == b)
    sub = np.asarray(K[:, cols].todense(), dtype=np.float32)
    lab, reps, csize, Kb = haploblock_collapse_indices(sub, eps=0.0)
    ac = sub.sum(0)
    informative = int(((ac >= 2) & (ac <= 229)).sum())
    effn = 1.0 / np.sum(hb[b] ** 2)
    s, e = blocks[b].start, blocks[b].end
    verdict = "NON-IDENT (K_b<200)" if Kb < 200 else ("EM-DRIFT (K_b hi,effn lo)" if effn < 180 else "ok")
    print(f"{b:3d} {s/1e6:5.1f}-{e/1e6:5.1f} {len(cols):>9,} {Kb:>4d}  {Kb/231:.3f}  {informative:>10,}  {effn:6.1f}   {verdict}")
