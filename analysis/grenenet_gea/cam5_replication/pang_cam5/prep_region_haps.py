#!/usr/bin/env python
"""Collapse the 231 founders to UNIQUE HAPLOTYPES over a focused region (by drawn
ALT pattern), keeping the whole-gene coordinate origin so the gene model aligns.
Emits gg_rhapseqs.tsv / gg_rhapfeats.tsv for plot_gggenomes_region.R.

Usage: prep_region_haps.py [REGION_LO REGION_HI]   (default = 3' end of CAM5)
"""
import sys, numpy as np, pandas as pd, pysam

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/cam5_replication/pang_cam5"
WIN_LO = 11531800                                  # same origin as the whole-gene fig
LO = int(sys.argv[1]) if len(sys.argv) > 2 else 11533810   # starts after the near-fixed cluster
HI = int(sys.argv[2]) if len(sys.argv) > 2 else 11534333

bio1 = pd.read_csv("/tmp/ecotype_bio1.tsv", sep="\t").set_index("ecotype")["bio1"].to_dict()
vcf = pysam.VariantFile(f"{HERE}/cam5_region.vcf.gz")
samp = [int(s) for s in vcf.header.samples]

wsites, Gw = [], []
for r in vcf.fetch("Chr2", LO, HI):
    rl, al = len(r.ref), len(r.alts[0]); p0 = r.pos - WIN_LO
    if rl == 1 and al == 1:                       # SNP
        vtype, size, s0, e0 = "snp", 0, p0, p0 + 1
    elif rl == al:                                # equal-length multi-base substitution
        vtype, size, s0, e0 = "mnp", 0, p0, p0 + rl
    elif rl > al:                                 # deletion: bar spans the deleted bases
        size = rl - al; vtype, s0, e0 = "del", p0 + al, p0 + al + size
    else:                                         # insertion: point, sized by inserted bp
        size = al - rl; vtype, s0, e0 = "ins", p0 + rl, p0 + rl
    wsites.append((s0, e0, vtype, size))
    Gw.append([1 if r.samples[s].get("GT", (None,))[0] == 1 else 0 for s in vcf.header.samples])
Gw = np.array(Gw).T
pat, inv = np.unique(Gw, axis=0, return_inverse=True)
b1 = np.array([bio1.get(e, np.nan) for e in samp])

hseq, hfeat = [], []
for k in range(pat.shape[0]):
    m = np.where(inv == k)[0]
    mb = np.nanmean(b1[m]) if np.any(~np.isnan(b1[m])) else np.nan
    hid = f"rhap{k:03d}"
    hseq.append((hid, hid, 2533, "hap", mb, len(m)))
    for j in np.where(pat[k] == 1)[0]:
        s0, e0, vtype, size = wsites[j]
        hfeat.append((hid, s0, e0, vtype, size))
hs = pd.DataFrame(hseq, columns=["seq_id", "bin_id", "length", "kind", "mean_bio1", "n_eco"])
hs = hs.sort_values("mean_bio1", kind="mergesort", na_position="last")
iso = pd.DataFrame([(i, i, 2533, "isoform", np.nan, 0) for i in
                    ["AT2G27030.1", "AT2G27030.2", "AT2G27030.3"]], columns=hs.columns)
pd.concat([iso, hs]).to_csv(f"{HERE}/gg_rhapseqs.tsv", sep="\t", index=False)
pd.DataFrame(hfeat, columns=["seq_id", "start", "end", "vtype", "size"]).to_csv(
    f"{HERE}/gg_rhapfeats.tsv", sep="\t", index=False)

# write the region bounds (gene coords) + overall founder mean bio1 (diverging midpoint)
pd.DataFrame({"name": ["region_lo", "region_hi", "bio1_mean"],
              "x": [LO - WIN_LO, HI - WIN_LO, round(float(np.nanmean(b1)), 3)]}
             ).to_csv(f"{HERE}/gg_region.tsv", sep="\t", index=False)
print(f"region {LO}-{HI} ({HI-LO}bp, {len(wsites)} sites): {pat.shape[0]} unique haplotypes "
      f"from {len(samp)} founders; max group {int(np.bincount(inv).max())}")
print("wrote gg_rhapseqs.tsv, gg_rhapfeats.tsv, gg_region.tsv")
