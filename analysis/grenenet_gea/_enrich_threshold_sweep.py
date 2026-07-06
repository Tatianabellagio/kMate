import numpy as np, pandas as pd
ST = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/site_temporal"
tab = pd.read_parquet(f"{ST}/site4_clq90_blocks.parquet")
scored = tab[tab.n_snp >= 3].copy()
rng = np.random.default_rng(0)
non_all = scored  # recompute selected per threshold
print(f"{'top_frac':>9} {'thr':>6} {'n_sel':>6} {'hasSV_sel':>9} {'null_med':>8} {'fold':>5} {'p_emp':>6}")
for tf in (0.02, 0.01, 0.005, 0.001):
    thr = scored.block_score.quantile(1 - tf)
    sel = scored[scored.block_score >= thr]
    non = scored[scored.block_score < thr]
    f_sel = sel.has_sv.mean()
    non_by_size = {nv: sub.index.to_numpy() for nv, sub in non.groupby("n_var")}
    sizes = np.array(sorted(non_by_size))
    hs = non.has_sv
    nperm = 3000
    null = np.empty(nperm)
    sel_sizes = sel.n_var.to_numpy()
    for b in range(nperm):
        picks = [rng.choice(non_by_size[sizes[np.argmin(np.abs(sizes - nv))]]) for nv in sel_sizes]
        null[b] = hs.loc[picks].mean()
    p = (np.sum(null >= f_sel) + 1) / (nperm + 1)
    fold = f_sel / np.median(null)
    print(f"{tf:>9} {thr:>6.3f} {len(sel):>6} {f_sel:>9.3f} {np.median(null):>8.3f} {fold:>5.2f} {p:>6.3f}")
