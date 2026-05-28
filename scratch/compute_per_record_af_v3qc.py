"""Compute per-record carrier-AF for v3qc vs GrENE-Net (Chr1 SNPs, matched on
(pos, ref, alt)). Saves npz with gn_af, v3qc_af, diff. Clone of
compute_per_record_af_panels.py task=1 (v3 all 231) but pointing at v3qc artifacts.
"""
import time, gzip
import numpy as np, pysam
from scipy.sparse import load_npz
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')
CN_VAR = ROOT / 'data/cn_var_231_v3qc.cn_var.npz'
CN_VAR_META = ROOT / 'data/cn_var_231_v3qc.meta.npz'
REFALT = ROOT / 'data/cn_var_231_v3qc.ref_alt.tsv.gz'
OUT_NPZ = ROOT / 'scratch/v3qc_vs_gn_per_record_af.npz'
GN_VCF = '/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf'


def log(*a, **kw): print(*a, **kw, flush=True)
t0 = time.time()

log('Loading GN Chr1 SNP carriers...')
gn = pysam.VariantFile(GN_VCF)
gn_samples = list(gn.header.samples)
F_gn = len(gn_samples)
keys, mats = [], []
for rec in gn:
    if rec.chrom != '1': continue
    if rec.alts is None or len(rec.alts) != 1: continue
    if len(rec.ref) != 1 or len(rec.alts[0]) != 1: continue
    carriers = np.zeros(F_gn, dtype=bool)
    for i, s in enumerate(gn_samples):
        gt = rec.samples[s]['GT']
        if gt is None: continue
        if any(a is not None and a > 0 for a in gt): carriers[i] = True
    keys.append((rec.pos, rec.ref, rec.alts[0]))
    mats.append(carriers)
gn_matrix = np.array(mats, dtype=bool)
gn_key_to_idx = {k: i for i, k in enumerate(keys)}
log(f'  GN: {gn_matrix.shape}, {time.time()-t0:.0f}s')

log('Loading v3qc cn_var + meta + ref_alt...')
cn = load_npz(CN_VAR)
meta = np.load(CN_VAR_META, allow_pickle=True)
founders = list(np.asarray(meta['founders']).astype(str))
chrom = np.asarray(meta['chrom']).astype(str)
pos = np.asarray(meta['pos']).astype(np.int64)
rl = np.asarray(meta['ref_len']).astype(np.int32)
al = np.asarray(meta['alt_len']).astype(np.int32)

chr1_snp_idx = np.where((chrom == 'Chr1') & (rl == 1) & (al == 1))[0]
log(f'  v3qc Chr1 SNP records: {len(chr1_snp_idx):,}')

refs = [None] * len(chrom); alts = [None] * len(chrom)
with gzip.open(REFALT, 'rt') as f:
    for i, line in enumerate(f):
        parts = line.rstrip().split('\t')
        refs[i] = parts[2]; alts[i] = parts[3]

panel_local, gn_idx = [], []
for li, g in enumerate(chr1_snp_idx):
    idx = gn_key_to_idx.get((int(pos[g]), refs[g], alts[g]))
    if idx is None: continue
    panel_local.append(li); gn_idx.append(idx)
panel_local = np.array(panel_local, dtype=np.int64)
gn_idx = np.array(gn_idx, dtype=np.int64)
log(f'  matched on (pos, ref, alt): {len(panel_local):,}')

f_to_gn = np.array([gn_samples.index(f) if f in gn_samples else -1 for f in founders])
overlap = f_to_gn >= 0
F_overlap = int(overlap.sum())
log(f'  founders overlap with GN: {F_overlap}/{len(founders)}')

panel_global = chr1_snp_idx[panel_local]
panel_carriers = (cn[:, panel_global].toarray() > 0)            # (F, n_matched)
panel_carriers_overlap = panel_carriers[overlap]                # (F_overlap, n_matched)
gn_carriers = gn_matrix[gn_idx]                                 # (n_matched, F_gn)
gn_indices_for_panel = f_to_gn[overlap]
gn_carriers_overlap = gn_carriers[:, gn_indices_for_panel].T    # (F_overlap, n_matched)

v3qc_af = panel_carriers_overlap.sum(axis=0).astype(np.float64) / F_overlap
gn_af   = gn_carriers_overlap.sum(axis=0).astype(np.float64)    / F_overlap
diff    = v3qc_af - gn_af

log(f'\n=== v3qc-vs-GN ===')
log(f'  n_records: {len(diff):,}  F_overlap: {F_overlap}')
log(f'  mean diff: {diff.mean():+.5f}  median: {np.median(diff):+.5f}')
log(f'  >+0.1 (over-call): {(diff>0.1).sum():,}  <-0.1 (under-call): {(diff<-0.1).sum():,}')

np.savez_compressed(OUT_NPZ, gn_af=gn_af, v3qc_af=v3qc_af, diff=diff,
                    label='v3qc (231)', F_overlap=F_overlap)
log(f'  saved {OUT_NPZ}')
log(f'  total wall: {time.time()-t0:.0f}s')
