"""Per-LD-block SNP and SV density audit.

For each hapFIRE LD block (3.2M-SNP-panel, density-0.5 partition):
  - n_snp_xwu   : SNPs in xwu's GrENE-Net VCF (the OLD less-dense panel; already in block index as block_n_snps)
  - n_snp_dense : SNPs in xwu ∪ cactus pangenome SNP positions (the DENSE panel we'd hand hapFIRE)
  - n_sv        : SVs from the production 231-founder panel (cactus 80 + Beagle-imputed 151), size_change >= 50

Produces summary tables for SV-rich / SNP-poor / SNP-empty blocks — the cells
where direct SV evidence in cactus_em could be load-bearing.
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd

_ROOT = Path(__file__).resolve().parents[2]
BLOCK_INDEX = str(_ROOT / 'poolfreq/data/hapfire_block_index.npz')
XWU_POS     = '/global/home/users/tbellg/scratch/freqk_gr/panel_overlap_test/data/xwu231_positions.tsv'   # chrom='1','2',...
CACTUS_POS  = '/global/home/users/tbellg/scratch/freqk_gr/panel_overlap_test/data/our82_snp_positions.tsv' # chrom='Chr1',...
CN_VAR_META = str(_ROOT / 'poolfreq/data/cn_var_231_v2.meta.npz')  # production 231-founder panel
OUT_TSV     = str(_ROOT / 'results/blocks_snp_sv_density.tsv.gz')
SV_MIN_SIZE_CHANGE = 50  # standard SV definition

# Chromosome name normalization: block index uses '1','2',... — map to 'Chr1',...
def to_chr(c):  # '1' -> 'Chr1'; 'Chr1' -> 'Chr1'
    s = str(c)
    return s if s.startswith('Chr') else f'Chr{s}'

def load_positions(path, has_header=False):
    """Returns dict chrom -> sorted np.int64 array of positions."""
    df = pd.read_csv(path, sep='\t', header=None if not has_header else 0,
                     names=['chrom','pos'], dtype={'chrom':str,'pos':np.int64})
    df['chrom'] = df['chrom'].map(to_chr)
    out = {}
    for c, g in df.groupby('chrom'):
        out[c] = np.sort(g['pos'].to_numpy())
    return out

def main():
    t0 = time.time()
    print(f'Loading block index: {BLOCK_INDEX}', flush=True)
    bi = np.load(BLOCK_INDEX, allow_pickle=True)
    block_chrom = np.asarray([to_chr(c) for c in bi['block_chrom']])
    block_pos_start = bi['block_pos_start'].astype(np.int64)
    block_pos_end   = bi['block_pos_end'].astype(np.int64)
    block_n_snp_xwu = bi['block_n_snps'].astype(np.int32)   # already counted from xwu VCF
    block_n_uniq    = bi['block_n_uniq'].astype(np.int32)
    n_blocks = len(block_chrom)
    print(f'  n_blocks={n_blocks:,}', flush=True)

    print(f'Loading cactus SNP positions: {CACTUS_POS}', flush=True)
    cact = load_positions(CACTUS_POS)
    print(f'  total cactus SNP positions: {sum(len(v) for v in cact.values()):,}', flush=True)
    for c, v in sorted(cact.items()):
        print(f'    {c}: {len(v):,}')

    print(f'Loading production cn_var meta (231-founder panel): {CN_VAR_META}', flush=True)
    m = np.load(CN_VAR_META, allow_pickle=True)
    sv_all = pd.DataFrame({
        'chrom': m['chrom'],
        'pos': m['pos'].astype(np.int64),
        'ref_len': m['ref_len'],
        'alt_len': m['alt_len'],
    })
    sv_all['chrom'] = sv_all['chrom'].map(to_chr)
    sv_all['size_change'] = (sv_all.alt_len.astype(int) - sv_all.ref_len.astype(int)).abs()
    sv = sv_all[sv_all.size_change >= SV_MIN_SIZE_CHANGE].copy()
    sv_pos = {c: np.sort(g['pos'].to_numpy()) for c, g in sv.groupby('chrom')}
    print(f'  total records in cn_var: {len(sv_all):,}')
    print(f'  SVs (size_change >= {SV_MIN_SIZE_CHANGE}): {len(sv):,}')
    for c, v in sorted(sv_pos.items()):
        print(f'    {c}: {len(v):,}')

    print(f'Loading xwu SNP positions (for sanity check): {XWU_POS}', flush=True)
    xwu = load_positions(XWU_POS)
    print(f'  total xwu SNP positions: {sum(len(v) for v in xwu.values()):,}', flush=True)

    # Per-block counts via searchsorted (positions arrays are sorted within chrom)
    n_snp_dense = np.zeros(n_blocks, dtype=np.int32)   # xwu ∪ cactus per block
    n_snp_cactus_only = np.zeros(n_blocks, dtype=np.int32)  # cactus minus xwu (extras)
    n_snp_xwu_recount = np.zeros(n_blocks, dtype=np.int32)  # sanity
    n_sv = np.zeros(n_blocks, dtype=np.int32)
    block_bp_len = (block_pos_end - block_pos_start + 1).astype(np.int64)

    print('Counting per-block intersections...', flush=True)
    for b in range(n_blocks):
        c = block_chrom[b]
        s, e = int(block_pos_start[b]), int(block_pos_end[b])
        if s == 0 and e == 0:
            continue
        if c in cact:
            arr = cact[c]
            n_cact = int(np.searchsorted(arr, e, side='right') - np.searchsorted(arr, s, side='left'))
        else:
            n_cact = 0
        if c in xwu:
            arr = xwu[c]
            n_xwu = int(np.searchsorted(arr, e, side='right') - np.searchsorted(arr, s, side='left'))
        else:
            n_xwu = 0
        if c in sv_pos:
            arr = sv_pos[c]
            n_sv[b] = int(np.searchsorted(arr, e, side='right') - np.searchsorted(arr, s, side='left'))
        # for "dense", count UNION via merge — fast via concatenate+unique
        if c in cact and c in xwu:
            cact_in = cact[c][(cact[c] >= s) & (cact[c] <= e)]
            xwu_in = xwu[c][(xwu[c] >= s) & (xwu[c] <= e)]
            n_snp_dense[b] = int(np.unique(np.concatenate([cact_in, xwu_in])).size)
            n_snp_cactus_only[b] = int(np.setdiff1d(cact_in, xwu_in, assume_unique=False).size)
        else:
            n_snp_dense[b] = max(n_cact, n_xwu)
        n_snp_xwu_recount[b] = n_xwu
        if (b + 1) % 2000 == 0:
            print(f'  {b+1}/{n_blocks} [{time.time()-t0:.0f}s]', flush=True)

    # Sanity check: xwu recount should match block_n_snps from index
    diff = (n_snp_xwu_recount - block_n_snp_xwu).astype(np.int64)
    print(f'\nSanity: xwu recount vs index block_n_snps — '
          f'mean Δ {diff.mean():.2f}, max |Δ| {abs(diff).max()}, '
          f'agreement {(diff==0).mean()*100:.1f}%', flush=True)

    df = pd.DataFrame({
        'chrom': block_chrom,
        'pos_start': block_pos_start,
        'pos_end': block_pos_end,
        'bp_len': block_bp_len,
        'n_snp_xwu': block_n_snp_xwu,
        'n_snp_xwu_recount': n_snp_xwu_recount,
        'n_snp_dense': n_snp_dense,
        'n_snp_cactus_only': n_snp_cactus_only,
        'n_sv': n_sv,
        'n_uniq_haps_xwu': block_n_uniq,
    })

    # Headline tables
    print('\n=== Per-block SNP / SV density ===\n', flush=True)
    print(f'Total blocks: {len(df):,}')
    print(f'Total bp covered: {df.bp_len.sum():,}')
    print(f'Total xwu SNPs in blocks: {df.n_snp_xwu.sum():,}')
    print(f'Total dense SNPs in blocks: {df.n_snp_dense.sum():,}')
    print(f'Total cactus-only extras: {df.n_snp_cactus_only.sum():,} '
          f'(+{df.n_snp_cactus_only.sum()/max(1,df.n_snp_xwu.sum())*100:.1f}% over xwu)')
    print(f'Total SVs in blocks: {df.n_sv.sum():,}')

    print('\n--- Distribution of n_snp_xwu per block ---')
    for q in [0,5,25,50,75,95,99,100]:
        print(f'  p{q:3d}: {np.percentile(df.n_snp_xwu, q):.0f}')
    print(f'  mean: {df.n_snp_xwu.mean():.1f}')

    print('\n--- Distribution of n_snp_dense per block ---')
    for q in [0,5,25,50,75,95,99,100]:
        print(f'  p{q:3d}: {np.percentile(df.n_snp_dense, q):.0f}')

    print('\n--- Distribution of n_sv per block ---')
    for q in [0,50,75,90,95,99,100]:
        print(f'  p{q:3d}: {np.percentile(df.n_sv, q):.0f}')
    print(f'  mean: {df.n_sv.mean():.2f}')
    print(f'  blocks with ≥1 SV: {(df.n_sv >= 1).sum():,} ({(df.n_sv >= 1).mean()*100:.1f}%)')

    print('\n--- Joint counts (xwu panel, the one hapFIRE actually uses) ---')
    bins_snp = [-1, 0, 9, 99, 999, 1e9]
    bins_lab = ['0', '1-10', '11-100', '101-1000', '>1000']
    df['snp_xwu_bin'] = pd.cut(df.n_snp_xwu, bins=bins_snp, labels=bins_lab)
    bins_sv = [-1, 0, 4, 19, 1e9]
    bins_sv_lab = ['0', '1-5', '6-20', '>20']
    df['sv_bin'] = pd.cut(df.n_sv, bins=bins_sv, labels=bins_sv_lab)
    ct = pd.crosstab(df.snp_xwu_bin, df.sv_bin, margins=True)
    print(ct.to_string())

    print('\n--- Joint counts (dense panel) ---')
    df['snp_dense_bin'] = pd.cut(df.n_snp_dense, bins=bins_snp, labels=bins_lab)
    ct2 = pd.crosstab(df.snp_dense_bin, df.sv_bin, margins=True)
    print(ct2.to_string())

    # The cells we care about
    print('\n=== CASES WHERE SV EVIDENCE WOULD BE LOAD-BEARING ===')

    # SV-only blocks (SNP-empty in xwu)
    sv_only_xwu = df[(df.n_snp_xwu == 0) & (df.n_sv >= 1)]
    print(f'\n[xwu] Blocks with 0 SNPs but ≥1 SV: {len(sv_only_xwu):,} '
          f'({len(sv_only_xwu)/len(df)*100:.2f}%)')
    if len(sv_only_xwu):
        print(f'  carrying {sv_only_xwu.n_sv.sum():,} SVs total')
        print(f'  median bp_len: {sv_only_xwu.bp_len.median():.0f}')

    sv_only_dense = df[(df.n_snp_dense == 0) & (df.n_sv >= 1)]
    print(f'\n[dense] Blocks with 0 SNPs but ≥1 SV: {len(sv_only_dense):,} '
          f'({len(sv_only_dense)/len(df)*100:.2f}%)')
    if len(sv_only_dense):
        print(f'  carrying {sv_only_dense.n_sv.sum():,} SVs total')

    # SV-rich / SNP-poor (SVs ≥ SNPs in count, i.e. SV evidence dominates)
    snp_poor_xwu = df[(df.n_snp_xwu < 50) & (df.n_sv >= 1)]
    print(f'\n[xwu] Blocks with <50 SNPs and ≥1 SV: {len(snp_poor_xwu):,} '
          f'({len(snp_poor_xwu)/len(df)*100:.2f}%) carrying {snp_poor_xwu.n_sv.sum():,} SVs')

    snp_poor_dense = df[(df.n_snp_dense < 50) & (df.n_sv >= 1)]
    print(f'[dense] Blocks with <50 SNPs and ≥1 SV: {len(snp_poor_dense):,} '
          f'({len(snp_poor_dense)/len(df)*100:.2f}%) carrying {snp_poor_dense.n_sv.sum():,} SVs')

    # The "SV evidence proportionally significant" regime
    # — k-mer share of SV ≈ 32 per SV, k-mer share per SNP ≈ 1
    df['sv_kmer_share'] = (32 * df.n_sv) / (df.n_snp_xwu + 32 * df.n_sv).clip(lower=1)
    df['sv_kmer_share_dense'] = (32 * df.n_sv) / (df.n_snp_dense + 32 * df.n_sv).clip(lower=1)
    print('\n--- SV-k-mer share of total block evidence (assuming 1 k-mer/SNP, 32 k-mers/SV) ---')
    print(f'[xwu]   blocks where SV evidence ≥ 30% of block evidence: '
          f'{(df.sv_kmer_share >= 0.30).sum():,} ({(df.sv_kmer_share >= 0.30).mean()*100:.2f}%)')
    print(f'[xwu]   blocks where SV evidence ≥ 50%: '
          f'{(df.sv_kmer_share >= 0.50).sum():,} ({(df.sv_kmer_share >= 0.50).mean()*100:.2f}%)')
    print(f'[dense] blocks where SV evidence ≥ 30%: '
          f'{(df.sv_kmer_share_dense >= 0.30).sum():,} ({(df.sv_kmer_share_dense >= 0.30).mean()*100:.2f}%)')
    print(f'[dense] blocks where SV evidence ≥ 50%: '
          f'{(df.sv_kmer_share_dense >= 0.50).sum():,} ({(df.sv_kmer_share_dense >= 0.50).mean()*100:.2f}%)')

    # Where do SV-dominated blocks live?
    sv_dom = df[df.sv_kmer_share >= 0.30]
    if len(sv_dom):
        print('\n[xwu] SV-evidence-dominated blocks (≥30% SV share) by chrom:')
        print(sv_dom.groupby('chrom').agg(
            n_blocks=('n_sv','size'),
            n_sv_total=('n_sv','sum'),
            mean_bp_len=('bp_len','mean'),
        ).to_string())

    print(f'\nWriting {OUT_TSV}', flush=True)
    from pathlib import Path
    Path(OUT_TSV).parent.mkdir(parents=True, exist_ok=True)
    df.drop(columns=['snp_xwu_bin','snp_dense_bin','sv_bin']).to_csv(
        OUT_TSV, sep='\t', index=False, compression='gzip')
    print(f'Done [{time.time()-t0:.0f}s]', flush=True)


if __name__ == '__main__':
    main()
