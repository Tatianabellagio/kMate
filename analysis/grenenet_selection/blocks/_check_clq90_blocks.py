import numpy as np, pandas as pd, glob, os
BD = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/blocks_mcf90"
blocks = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{BD}/chr*_clq0.9_blocks_clq0.9.tsv"))])
print("total clq0.9 blocks:", len(blocks), "| chroms:", sorted(blocks.chrom.unique()))
# disjoint per chrom?
ov = 0
for c, g in blocks.groupby("chrom"):
    g = g.sort_values("start_pos")
    ov += int((g.start_pos.values[1:] < g.end_pos.values[:-1]).sum())
print("overlapping adjacent block pairs:", ov)
print("span bp: median", int((blocks.end_pos - blocks.start_pos).median()),
      "| n_variants/block median", int(blocks.n_variants.median()),
      "| total covered Mb", round((blocks.end_pos - blocks.start_pos).sum()/1e6, 1))

# assign per-variant to block, per class, count SNP coverage
ST = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/site_temporal"
def load(c):
    z = np.load(f"{ST}/site4_scoef_{c}.npz", allow_pickle=True); k = z["keep"]
    return pd.DataFrame(dict(chrom=z["chrom"][k].astype(str), pos=z["pos"][k].astype(np.int64), cls=c))
V = pd.concat([load(c) for c in ("snp", "smallindel", "sv")], ignore_index=True)

def assign(V, blocks):
    bid = np.full(len(V), -1, np.int64)
    blocks = blocks.reset_index(drop=True)
    for c, g in blocks.groupby("chrom"):
        g = g.sort_values("start_pos")
        st = g.start_pos.to_numpy(); en = g.end_pos.to_numpy(); idx = g.index.to_numpy()
        vi = np.where(V.chrom.values == c)[0]
        p = V.pos.values[vi]
        j = np.searchsorted(st, p, "right") - 1
        ok = (j >= 0) & (j < len(st))
        inb = ok & (p <= en[np.clip(j, 0, len(en)-1)])
        bid[vi[inb]] = idx[j[inb]]
    return bid
V["bid"] = assign(V, blocks)
print("\nvariants assigned to a clq0.9 block:", (V.bid >= 0).mean().round(3),
      "| unassigned (inter-block gaps):", int((V.bid < 0).sum()))
comp = V[V.bid >= 0].groupby(["bid", "cls"]).size().unstack(fill_value=0)
for c in ("snp", "smallindel", "sv"):
    if c not in comp: comp[c] = 0
nb = len(comp)
print(f"blocks containing reachable variants: {nb:,}")
print(f"  with >=1 SNP: {(comp.snp>0).mean():.3f} | with >=1 SV: {(comp.sv>0).mean():.3f} "
      f"| SV-only (no snp,no indel): {((comp.sv>0)&(comp.snp==0)&(comp.smallindel==0)).sum()}")
print(f"  median n_snp/block: {int(comp.snp.median())} | blocks with 0 SNP: {(comp.snp==0).sum()}")
