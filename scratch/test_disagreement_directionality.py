"""Characterize where v3 vs GrENE-Net carrier disagreements come from:
direction (v3 over- vs under-assigns), founder class, INDEL distance, AF tier.

Inputs: cn_var_v3 (atomic Chr1 records + carrier matrix), GrENE-Net VCF (Chr1
biallelic SNPs), founder split.

Output: tables answering — is v3 systematically over-calling ALT (discovery
signal), under-calling ALT (panel-build loss), or symmetric (noise)?

Run as SLURM (32G mem, ~10 min wall).
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import load_npz

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
SRC_VCF = ROOT / 'pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz'
GN_VCF  = Path('/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf')
CN_VAR  = ROOT / 'poolfreq/data/cn_var_231_v3.cn_var.npz'
CN_VAR_META = ROOT / 'poolfreq/data/cn_var_231_v3.meta.npz'
REFALT  = ROOT / 'poolfreq/data/cn_var_231_v3.ref_alt.tsv.gz'
SPLIT_FILE = ROOT / 'data/founder_split_cactus_pg.json'
OUT_TSV = ROOT / 'scratch/disagreement_directionality.tsv'


def log(*a, **kw): print(*a, **kw, flush=True)


def main():
    t0 = time.time()
    log('Loading cn_var_v3 + meta...')
    cn_var = load_npz(CN_VAR)
    meta = np.load(CN_VAR_META, allow_pickle=True)
    rec_chrom = np.asarray(meta['chrom']).astype(str)
    rec_pos = np.asarray(meta['pos']).astype(np.int64)
    rec_ref_len = np.asarray(meta['ref_len']).astype(np.int32)
    rec_alt_len = np.asarray(meta['alt_len']).astype(np.int32)
    founders = list(np.asarray(meta['founders']).astype(str))
    F = len(founders)
    chr1 = np.where(rec_chrom == 'Chr1')[0]
    log(f'  cn_var Chr1 records: {len(chr1):,}')

    log('Loading REF/ALT sidecar...')
    rec_ref = [None] * len(rec_chrom)
    rec_alt = [None] * len(rec_chrom)
    with pd.read_csv(REFALT, sep='\t', header=None, names=['chrom','pos','ref','alt'],
                      chunksize=500_000) as reader:
        idx = 0
        for ch in reader:
            for r, a in zip(ch['ref'].values, ch['alt'].values):
                rec_ref[idx] = r; rec_alt[idx] = a; idx += 1

    # Distance from each Chr1 SNP to nearest INDEL (any v3 record where ref_len != alt_len)
    log('Computing distance to nearest INDEL on Chr1...')
    chr1_pos = rec_pos[chr1]
    chr1_is_snp = (rec_ref_len[chr1] == 1) & (rec_alt_len[chr1] == 1)
    chr1_is_indel = ~chr1_is_snp
    indel_positions = np.sort(chr1_pos[chr1_is_indel])
    log(f'  Chr1 INDELs in v3: {len(indel_positions):,}')

    # Founder class split
    split = json.load(open(SPLIT_FILE))
    cactus_set = set(split['cactus']); pg_set = set(split['PG'])
    is_cactus = np.array([f in cactus_set for f in founders])
    is_pg = np.array([f in pg_set for f in founders])

    # SNP indices in chr1
    snp_idx_in_chr1 = np.where(chr1_is_snp)[0]
    log(f'  Chr1 SNPs in v3: {len(snp_idx_in_chr1):,}')

    # 2. Load GrENE-Net carriers per SNP key (pos, REF, ALT) — only those that match
    log('Loading GrENE-Net Chr1 biallelic SNP carriers...')
    t = time.time()
    gn_vcf = pysam.VariantFile(str(GN_VCF))
    gn_samples = list(gn_vcf.header.samples)
    f_to_gn = np.array([gn_samples.index(f) if f in gn_samples else -1 for f in founders])

    gn_carriers_by_key = {}
    n = 0
    for rec in gn_vcf:
        if rec.chrom != '1':
            continue
        if rec.alts is None or len(rec.alts) != 1: continue
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1: continue
        key = (rec.pos, rec.ref, rec.alts[0])
        carriers = np.zeros(F, dtype=bool)
        for v3_i, gn_i in enumerate(f_to_gn):
            if gn_i < 0: continue
            gt = rec.samples[gn_samples[gn_i]]['GT']
            if gt is None: continue
            if any(a is not None and a > 0 for a in gt):
                carriers[v3_i] = True
        gn_carriers_by_key[key] = carriers
        n += 1
        if n % 200000 == 0:
            log(f'    parsed {n:,} GN biallelic SNPs')
    log(f'  GN biallelic SNPs total: {n:,}, {time.time()-t:.0f}s')

    # 3. For each v3 Chr1 SNP, look up GN carriers, compute 2x2 confusion
    log('Computing 2x2 confusion per (founder, matched SNP)...')
    t = time.time()
    cn_var_chr1 = cn_var[:, chr1].tocsc()

    # Accumulators (vectorized would be ideal, but doing per-SNP loop for clarity)
    n_matched = 0
    n_v3_only_total = 0          # v3 says 1, GN says 0
    n_gn_only_total = 0          # GN says 1, v3 says 0
    n_both_total   = 0          # both 1
    n_neither_total = 0

    # By founder class
    cls_counts = {
        'cactus': dict(both=0, v3_only=0, gn_only=0, neither=0),
        'pg':     dict(both=0, v3_only=0, gn_only=0, neither=0),
    }

    # By distance-to-INDEL bucket
    dist_buckets = [(0,10),(10,50),(50,100),(100,500),(500,5000),(5000,10_000_000)]
    dist_counts = {bk: dict(both=0, v3_only=0, gn_only=0, neither=0) for bk in dist_buckets}

    # By GN AF tier (rare / mid / common)
    af_buckets = [(0.0, 0.05),(0.05,0.2),(0.2,0.5),(0.5,1.0)]
    af_counts = {bk: dict(both=0, v3_only=0, gn_only=0, neither=0) for bk in af_buckets}

    is_cactus_idx = np.where(is_cactus)[0]
    is_pg_idx = np.where(is_pg)[0]

    n_processed = 0
    for s_local in snp_idx_in_chr1:
        s_global = chr1[s_local]
        key = (int(rec_pos[s_global]), rec_ref[s_global], rec_alt[s_global])
        gn = gn_carriers_by_key.get(key)
        if gn is None:
            continue
        v3 = (cn_var_chr1[:, s_local].toarray().flatten() > 0)
        n_matched += 1

        # Per-founder confusion at this SNP
        both = (v3 & gn)
        v3_only = (v3 & ~gn)
        gn_only = (~v3 & gn)
        neither = (~v3 & ~gn)

        n_both_total    += int(both.sum())
        n_v3_only_total += int(v3_only.sum())
        n_gn_only_total += int(gn_only.sum())
        n_neither_total += int(neither.sum())

        # By class
        cls_counts['cactus']['both']    += int(both[is_cactus].sum())
        cls_counts['cactus']['v3_only'] += int(v3_only[is_cactus].sum())
        cls_counts['cactus']['gn_only'] += int(gn_only[is_cactus].sum())
        cls_counts['cactus']['neither'] += int(neither[is_cactus].sum())
        cls_counts['pg']['both']    += int(both[is_pg].sum())
        cls_counts['pg']['v3_only'] += int(v3_only[is_pg].sum())
        cls_counts['pg']['gn_only'] += int(gn_only[is_pg].sum())
        cls_counts['pg']['neither'] += int(neither[is_pg].sum())

        # By distance to INDEL
        pos = int(rec_pos[s_global])
        i_left = np.searchsorted(indel_positions, pos) - 1
        i_right = i_left + 1
        d_left = abs(pos - indel_positions[max(i_left,0)]) if len(indel_positions) else 1e9
        d_right = abs(indel_positions[min(i_right, len(indel_positions)-1)] - pos) if len(indel_positions) else 1e9
        d = min(d_left, d_right)
        for bk in dist_buckets:
            if bk[0] <= d < bk[1]:
                dist_counts[bk]['both']    += int(both.sum())
                dist_counts[bk]['v3_only'] += int(v3_only.sum())
                dist_counts[bk]['gn_only'] += int(gn_only.sum())
                dist_counts[bk]['neither'] += int(neither.sum())
                break

        # By GN AF
        af = gn.sum() / max((f_to_gn >= 0).sum(), 1)
        for bk in af_buckets:
            if bk[0] <= af < bk[1]:
                af_counts[bk]['both']    += int(both.sum())
                af_counts[bk]['v3_only'] += int(v3_only.sum())
                af_counts[bk]['gn_only'] += int(gn_only.sum())
                af_counts[bk]['neither'] += int(neither.sum())
                break

        n_processed += 1
        if n_processed % 50000 == 0:
            log(f'    matched {n_processed:,}; '
                f'v3_over={n_v3_only_total:,} v3_under={n_gn_only_total:,} '
                f'ratio={n_v3_only_total/max(n_gn_only_total,1):.2f}')

    log(f'\n=== Done. Matched SNPs: {n_matched:,} ({time.time()-t:.0f}s) ===')
    log(f'\n=== Overall direction (per-founder, per-SNP cell counts across {n_matched:,} matched SNPs × 231 founders) ===')
    total_cells = n_both_total + n_v3_only_total + n_gn_only_total + n_neither_total
    log(f'  total cells:            {total_cells:>14,}')
    log(f'  both agree carrier:     {n_both_total:>14,}  ({100*n_both_total/total_cells:5.2f}%)')
    log(f'  both agree non-carrier: {n_neither_total:>14,}  ({100*n_neither_total/total_cells:5.2f}%)')
    log(f'  v3 OVER (v3=1, gn=0):   {n_v3_only_total:>14,}  ({100*n_v3_only_total/total_cells:5.2f}%)')
    log(f'  v3 UNDER (v3=0, gn=1):  {n_gn_only_total:>14,}  ({100*n_gn_only_total/total_cells:5.2f}%)')
    if n_gn_only_total > 0:
        log(f'\n  *** RATIO v3_over / v3_under = {n_v3_only_total/n_gn_only_total:.3f} ***')
        log(f'  Reading:')
        log(f'    >> 1: v3 systematically calls MORE carriers than GN (long-read discovery? or over-call?)')
        log(f'    ≈ 1: symmetric disagreement (noise without directionality)')
        log(f'    << 1: v3 systematically calls FEWER carriers than GN (atomization loses signal)')

    log(f'\n=== By founder class ===')
    log(f'{"class":>8s}  {"both":>10s}  {"v3_OVER":>10s}  {"v3_UNDER":>10s}  '
        f'{"ratio_O/U":>10s}  {"%_disagree":>10s}')
    for cls in ['cactus', 'pg']:
        d = cls_counts[cls]
        tot = sum(d.values())
        disagree = d['v3_only'] + d['gn_only']
        ratio = d['v3_only'] / max(d['gn_only'], 1)
        log(f'{cls:>8s}  {d["both"]:>10,}  {d["v3_only"]:>10,}  {d["gn_only"]:>10,}  '
            f'{ratio:>10.3f}  {100*disagree/max(tot,1):>9.2f}%')

    log(f'\n=== By distance to nearest INDEL ===')
    log(f'{"dist (bp)":>15s}  {"both":>12s}  {"v3_OVER":>12s}  {"v3_UNDER":>12s}  '
        f'{"ratio_O/U":>10s}  {"%_disagree":>10s}')
    for bk in dist_buckets:
        d = dist_counts[bk]
        tot = sum(d.values())
        if tot == 0: continue
        disagree = d['v3_only'] + d['gn_only']
        ratio = d['v3_only'] / max(d['gn_only'], 1)
        lab = f'{bk[0]}–{bk[1]}'
        log(f'{lab:>15s}  {d["both"]:>12,}  {d["v3_only"]:>12,}  {d["gn_only"]:>12,}  '
            f'{ratio:>10.3f}  {100*disagree/max(tot,1):>9.2f}%')

    log(f'\n=== By GN MAF ===')
    log(f'{"AF":>15s}  {"both":>12s}  {"v3_OVER":>12s}  {"v3_UNDER":>12s}  '
        f'{"ratio_O/U":>10s}  {"%_disagree":>10s}')
    for bk in af_buckets:
        d = af_counts[bk]
        tot = sum(d.values())
        if tot == 0: continue
        disagree = d['v3_only'] + d['gn_only']
        ratio = d['v3_only'] / max(d['gn_only'], 1)
        lab = f'{bk[0]:.2f}–{bk[1]:.2f}'
        log(f'{lab:>15s}  {d["both"]:>12,}  {d["v3_only"]:>12,}  {d["gn_only"]:>12,}  '
            f'{ratio:>10.3f}  {100*disagree/max(tot,1):>9.2f}%')

    # Save flat TSV
    rows = []
    rows.append(dict(stratum='OVERALL', class_or_bucket='all',
                     both=n_both_total, v3_over=n_v3_only_total,
                     v3_under=n_gn_only_total, neither=n_neither_total))
    for cls in ['cactus', 'pg']:
        d = cls_counts[cls]
        rows.append(dict(stratum='founder_class', class_or_bucket=cls,
                         both=d['both'], v3_over=d['v3_only'],
                         v3_under=d['gn_only'], neither=d['neither']))
    for bk, d in dist_counts.items():
        rows.append(dict(stratum='dist_to_indel', class_or_bucket=f'{bk[0]}-{bk[1]}',
                         both=d['both'], v3_over=d['v3_only'],
                         v3_under=d['gn_only'], neither=d['neither']))
    for bk, d in af_counts.items():
        rows.append(dict(stratum='gn_af', class_or_bucket=f'{bk[0]}-{bk[1]}',
                         both=d['both'], v3_over=d['v3_only'],
                         v3_under=d['gn_only'], neither=d['neither']))
    pd.DataFrame(rows).to_csv(OUT_TSV, sep='\t', index=False)
    log(f'\nSaved {OUT_TSV}')
    log(f'\nTotal wall: {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
