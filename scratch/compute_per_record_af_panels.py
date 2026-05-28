"""Compute per-record carrier-AF for one panel vs GrENE-Net (Chr1 SNPs, matched on
(pos, ref, alt)). Saves npz with gn_af, panel_af, diff. Mirrors one_panel_compare.py
geometry but writes per-record arrays instead of summary counts.

Args: <task_id 1..3>
"""
import sys, time, json, gzip
import numpy as np, pysam
from scipy.sparse import load_npz
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')
task = int(sys.argv[1])

PANELS = {
    1: dict(label='v3 (all 231)',
            cn=ROOT/'data/cn_var_231_v3.cn_var.npz',
            meta=ROOT/'data/cn_var_231_v3.meta.npz',
            refalt=ROOT/'data/cn_var_231_v3.ref_alt.tsv.gz',
            subset=None,
            out_npz=ROOT/'scratch/per_record_af_v3.npz'),
    2: dict(label='v3 cactus-only (80)',
            cn=ROOT/'data/cn_var_231_v3.cn_var.npz',
            meta=ROOT/'data/cn_var_231_v3.meta.npz',
            refalt=ROOT/'data/cn_var_231_v3.ref_alt.tsv.gz',
            subset='cactus',
            out_npz=ROOT/'scratch/per_record_af_v3_cactus_only.npz'),
    3: dict(label='p82 (82 cactus, no PG)',
            cn=ROOT/'control_p82/data/cn_var_p82.cn_var.npz',
            meta=ROOT/'control_p82/data/cn_var_p82.meta.npz',
            refalt=ROOT/'control_p82/data/cn_var_p82.ref_alt.tsv',
            subset=None,
            out_npz=ROOT/'scratch/per_record_af_p82.npz'),
}
cfg = PANELS[task]

def log(*a, **kw): print(*a, **kw, flush=True)
t0 = time.time()

log('Loading GN Chr1 SNP carriers...')
gn = pysam.VariantFile('/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf')
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
    keys.append((rec.pos, rec.ref, rec.alts[0])); mats.append(carriers)
gn_matrix = np.array(mats, dtype=bool)
gn_key_to_idx = {k: i for i, k in enumerate(keys)}
log(f'  GN: {gn_matrix.shape}, {time.time()-t0:.0f}s')

log(f'\nLoading {cfg["label"]}...')
cn = load_npz(cfg['cn'])
meta = np.load(cfg['meta'], allow_pickle=True)
founders = list(np.asarray(meta['founders']).astype(str))
chrom = np.asarray(meta['chrom']).astype(str)
pos = np.asarray(meta['pos']).astype(np.int64)
rl = np.asarray(meta['ref_len']).astype(np.int32)
al = np.asarray(meta['alt_len']).astype(np.int32)

if cfg['subset'] == 'cactus':
    split = json.load(open(ROOT/'data/founder_split_cactus_pg.json'))
    cactus_set = set(split['cactus'])
    keep = np.array([f in cactus_set for f in founders])
    cn = cn[np.where(keep)[0], :]
    founders = [f for f, k in zip(founders, keep) if k]
    log(f'  subset to cactus: {cn.shape}, {len(founders)} founders')

chr1_snp_idx = np.where((chrom=='Chr1') & (rl==1) & (al==1))[0]
log(f'  Chr1 SNP records: {len(chr1_snp_idx):,}')
refs = [None]*len(chrom); alts = [None]*len(chrom)
open_fn = (lambda p: gzip.open(p, 'rt')) if str(cfg['refalt']).endswith('.gz') else open
i = 0
with open_fn(cfg['refalt']) as f:
    for line in f:
        parts = line.rstrip().split('\t')
        refs[i] = parts[2]; alts[i] = parts[3]; i += 1

panel_local, gn_idx = [], []
for li, g in enumerate(chr1_snp_idx):
    idx = gn_key_to_idx.get((int(pos[g]), refs[g], alts[g]))
    if idx is None: continue
    panel_local.append(li); gn_idx.append(idx)
panel_local = np.array(panel_local, dtype=np.int64)
gn_idx = np.array(gn_idx, dtype=np.int64)
log(f'  matched: {len(panel_local):,}')

f_to_gn = np.array([gn_samples.index(f) if f in gn_samples else -1 for f in founders])
overlap = f_to_gn >= 0
F_overlap = int(overlap.sum())
log(f'  overlap with GN: {F_overlap}/{len(founders)}')

panel_global = chr1_snp_idx[panel_local]
panel_carriers = (cn[:, panel_global].toarray() > 0)            # (F_panel, n_matched)
panel_carriers_overlap = panel_carriers[overlap]                # (F_overlap, n_matched)
gn_carriers = gn_matrix[gn_idx]                                 # (n_matched, F_gn)
gn_indices_for_panel = f_to_gn[overlap]
gn_carriers_overlap = gn_carriers[:, gn_indices_for_panel].T    # (F_overlap, n_matched)

panel_af = panel_carriers_overlap.sum(axis=0).astype(np.float64) / F_overlap
gn_af    = gn_carriers_overlap.sum(axis=0).astype(np.float64)    / F_overlap
diff     = panel_af - gn_af

log(f'\n=== {cfg["label"]} ===')
log(f'  n_records: {len(diff):,}  F_overlap: {F_overlap}')
log(f'  mean diff: {diff.mean():+.5f}  median: {np.median(diff):+.5f}')
log(f'  >+0.1: {(diff>0.1).sum():,}  <-0.1: {(diff<-0.1).sum():,}')

np.savez_compressed(cfg['out_npz'], gn_af=gn_af, panel_af=panel_af, diff=diff,
                    label=cfg['label'], F_overlap=F_overlap)
log(f'  saved {cfg["out_npz"]}')
log(f'  total wall: {time.time()-t0:.0f}s')
