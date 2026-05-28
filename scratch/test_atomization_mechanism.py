"""Test the hypothesis: are missing-carrier founders (GN says yes, v3 says no
at the same SNP) enriched for INDEL-co-occurrence within the same cactus snarl?

If yes → consistent with `bcftools norm` decomposition dropping founders whose
chosen bubble allele was indel-bearing from the atomic SNP record.

If no → mechanism is something else; rethink.

Built-in control: founders correctly listed in both panels at the same SNPs.
Compute P(indel-co-occurrence within snarl | missing) vs same | correctly-listed.

Chr1 only (tractable).
"""
from __future__ import annotations
import sys, time
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


def log(*a, **kw): print(*a, **kw, flush=True)


def main():
    t0 = time.time()
    # 1. Load cn_var_v3 + meta + ref/alt sidecar (Chr1 slice)
    log('Loading cn_var_v3...')
    cn_var = load_npz(CN_VAR)
    meta = np.load(CN_VAR_META, allow_pickle=True)
    rec_chrom = np.asarray(meta['chrom']).astype(str)
    rec_pos = np.asarray(meta['pos']).astype(np.int64)
    rec_ref_len = np.asarray(meta['ref_len']).astype(np.int32)
    rec_alt_len = np.asarray(meta['alt_len']).astype(np.int32)
    founders = list(np.asarray(meta['founders']).astype(str))
    F = len(founders)
    log(f'  founders: {F}')

    chr1 = np.where(rec_chrom == 'Chr1')[0]
    log(f'  Chr1 records: {len(chr1):,}')

    # Load ref/alt for these records — keep as Python lists; one SV ALT is ~100K chars,
    # which would force a fixed-width numpy array to allocate TB. Lists are fine.
    log('Loading REF/ALT sidecar...')
    rec_ref = [None] * len(rec_chrom)
    rec_alt = [None] * len(rec_chrom)
    with pd.read_csv(REFALT, sep='\t', header=None, names=['chrom','pos','ref','alt'],
                      chunksize=500_000) as reader:
        idx = 0
        for ch in reader:
            for r, a in zip(ch['ref'].values, ch['alt'].values):
                rec_ref[idx] = r; rec_alt[idx] = a; idx += 1

    # 2. Pull snarl IDs from source VCF for Chr1 records
    log('Pulling snarl IDs from source VCF (Chr1)...')
    t = time.time()
    snarl_id_by_pos_ref_alt = {}
    vcf = pysam.VariantFile(str(SRC_VCF))
    for rec in vcf.fetch('Chr1'):
        if rec.alts is None or len(rec.alts) != 1:
            continue
        snarl_id_by_pos_ref_alt[(rec.pos, rec.ref, rec.alts[0])] = rec.id
    log(f'  snarl IDs loaded: {len(snarl_id_by_pos_ref_alt):,} in {time.time()-t:.0f}s')

    # Map each Chr1 cn_var record to its snarl ID
    log('Mapping cn_var records → snarls...')
    snarl_per_chr1_rec = np.array([
        snarl_id_by_pos_ref_alt.get((int(rec_pos[i]), rec_ref[i], rec_alt[i]), None)
        for i in chr1
    ])
    n_mapped = (snarl_per_chr1_rec != None).sum()
    log(f'  mapped {n_mapped:,} / {len(chr1):,} Chr1 records to snarls')

    # 3. Identify SNP records on Chr1
    chr1_snp_mask = (rec_ref_len[chr1] == 1) & (rec_alt_len[chr1] == 1)
    snp_idx_in_chr1 = np.where(chr1_snp_mask)[0]
    log(f'  Chr1 SNPs in v3: {len(snp_idx_in_chr1):,}')

    # 4. Build position-keyed Chr1 SNP carrier sets from v3
    cn_var_chr1 = cn_var[:, chr1]  # F × |chr1|
    log(f'  cn_var Chr1 shape: {cn_var_chr1.shape}')

    # 5. Build GN VCF Chr1 SNP carrier sets, matched by (pos, REF, ALT)
    log('Loading GN Chr1 SNP carriers (matched by pos/REF/ALT)...')
    t = time.time()
    gn_vcf = pysam.VariantFile(str(GN_VCF))
    gn_samples = list(gn_vcf.header.samples)
    # GrENE-Net uses chrom '1', not 'Chr1'
    gn_samples_set = set(gn_samples)
    # Map founder index in v3 → GN sample index (or -1 if missing)
    f_to_gn = []
    for f in founders:
        if f in gn_samples_set:
            f_to_gn.append(gn_samples.index(f))
        else:
            f_to_gn.append(-1)
    f_to_gn = np.array(f_to_gn)
    log(f'  v3 founders matched in GN: {(f_to_gn >= 0).sum()} / {F}')

    # Build a dict: (pos, ref, alt) → set of v3 founder IDs (carrier set in GN)
    gn_chrom = '1'
    gn_carrier_by_key = {}
    n_gn_rec = 0
    n_gn_match = 0
    # GN VCF has no index; just stream the whole file and filter on chrom
    for rec in gn_vcf:
        if rec.chrom != gn_chrom:
            continue
        n_gn_rec += 1
        if rec.alts is None or len(rec.alts) != 1:
            continue
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        key = (rec.pos, rec.ref, rec.alts[0])
        # only build dict for v3 SNP keys we care about
        carriers = []
        for v3_f_idx, gn_idx in enumerate(f_to_gn):
            if gn_idx < 0:
                continue
            gt = rec.samples[gn_samples[gn_idx]]['GT']
            if gt is None: continue
            if any(a is not None and a > 0 for a in gt):
                carriers.append(v3_f_idx)
        gn_carrier_by_key[key] = set(carriers)
        n_gn_match += 1
        if n_gn_match % 100000 == 0:
            log(f'    GN parsed {n_gn_rec:,} records, {n_gn_match:,} biallelic SNPs')
    log(f'  GN total Chr1 records: {n_gn_rec:,}, biallelic SNPs: {n_gn_match:,}, '
        f'{time.time()-t:.0f}s')

    # 6. Match v3 Chr1 SNPs to GN Chr1 SNPs by (pos, REF, ALT)
    log('Matching v3 SNP records to GN by (pos, REF, ALT)...')
    snp_global_indices = chr1[snp_idx_in_chr1]
    v3_snp_keys = [
        (int(rec_pos[g]), rec_ref[g], rec_alt[g])
        for g in snp_global_indices
    ]
    matched_v3_snp_local_idx = []  # indices into snp_idx_in_chr1
    matched_gn_carriers = []        # corresponding GN carrier set
    for local_i, key in enumerate(v3_snp_keys):
        if key in gn_carrier_by_key:
            matched_v3_snp_local_idx.append(local_i)
            matched_gn_carriers.append(gn_carrier_by_key[key])
    log(f'  matched v3↔GN SNPs: {len(matched_v3_snp_local_idx):,}')

    # 7. For each matched SNP, compute v3 carrier set and the asymmetry vs GN
    log('Computing per-SNP carrier disagreements and indel-co-occurrence...')
    # For each matched SNP:
    #   missing = GN_carriers - V3_carriers (founders GN says yes, v3 says no)
    #   correct_listed = GN_carriers ∩ V3_carriers (founders both say yes)
    # For each missing carrier F at SNP s:
    #   snarl = snarl_per_chr1_rec[s_idx_in_chr1]
    #   indel_in_snarl = any chr1 record sharing this snarl with ref_len != alt_len
    #   F-indel-co-occurrence = does F have v3 carrier=1 at any INDEL in the snarl?

    # Pre-group records by snarl: snarl → list of chr1 record local indices
    snarl_to_local_idxs = {}
    for local_j, sn in enumerate(snarl_per_chr1_rec):
        if sn is None:
            continue
        snarl_to_local_idxs.setdefault(sn, []).append(local_j)
    log(f'  unique snarls on Chr1: {len(snarl_to_local_idxs):,}')

    # Per-snarl: indel record indices
    indel_records_in_snarl = {}
    for sn, idxs in snarl_to_local_idxs.items():
        indel_idxs = [j for j in idxs
                       if rec_ref_len[chr1[j]] != rec_alt_len[chr1[j]]]
        if indel_idxs:
            indel_records_in_snarl[sn] = indel_idxs
    log(f'  snarls containing ≥1 INDEL atomic record: {len(indel_records_in_snarl):,}')

    # We need fast access to v3 carrier status — convert to CSR with chr1 cols
    cn_var_chr1_csc = cn_var_chr1.tocsc()  # F × |chr1|

    # Aggregate counts
    missing_total = 0
    missing_with_indel_co = 0
    correct_total = 0
    correct_with_indel_co = 0
    # Stratify by founder class (cactus/PG)
    import json
    split = json.load(open(ROOT / 'data/founder_split_cactus_pg.json'))
    cactus_set = set(split['cactus']); pg_set = set(split['PG'])
    is_cactus = np.array([f in cactus_set for f in founders])
    is_pg = np.array([f in pg_set for f in founders])

    miss_cac = miss_pg = 0
    miss_co_cac = miss_co_pg = 0
    corr_cac = corr_pg = 0
    corr_co_cac = corr_co_pg = 0

    # Also record per-SNP disagreement counts for distance bucketing later
    snp_disagree_records = []

    for matched_i, local_in_chr1_snps in enumerate(matched_v3_snp_local_idx):
        gn_carriers = matched_gn_carriers[matched_i]
        s_local = snp_idx_in_chr1[local_in_chr1_snps]  # local idx into chr1
        s_global = chr1[s_local]
        # v3 carriers at this SNP
        col = cn_var_chr1_csc[:, s_local]
        v3_carriers_arr = col.toarray().flatten()
        v3_carriers = set(np.where(v3_carriers_arr > 0)[0].tolist())

        missing = gn_carriers - v3_carriers
        correct = gn_carriers & v3_carriers
        extra = v3_carriers - gn_carriers

        # Snarl ID for this SNP record
        sn = snarl_per_chr1_rec[s_local]
        indel_locals = indel_records_in_snarl.get(sn, []) if sn is not None else []
        # For each founder in missing, check carriership at any INDEL in snarl
        if indel_locals:
            indel_carriers_in_snarl = np.zeros(F, dtype=bool)
            for il in indel_locals:
                col_i = cn_var_chr1_csc[:, il].toarray().flatten()
                indel_carriers_in_snarl |= (col_i > 0)
        else:
            indel_carriers_in_snarl = np.zeros(F, dtype=bool)

        snp_disagree_records.append(dict(
            pos=int(rec_pos[s_global]),
            n_missing=len(missing),
            n_correct=len(correct),
            n_extra=len(extra),
            has_indel_in_snarl=bool(indel_locals),
        ))

        for F_idx in missing:
            missing_total += 1
            if is_cactus[F_idx]: miss_cac += 1
            elif is_pg[F_idx]:   miss_pg += 1
            if indel_carriers_in_snarl[F_idx]:
                missing_with_indel_co += 1
                if is_cactus[F_idx]: miss_co_cac += 1
                elif is_pg[F_idx]:   miss_co_pg += 1
        for F_idx in correct:
            correct_total += 1
            if is_cactus[F_idx]: corr_cac += 1
            elif is_pg[F_idx]:   corr_pg += 1
            if indel_carriers_in_snarl[F_idx]:
                correct_with_indel_co += 1
                if is_cactus[F_idx]: corr_co_cac += 1
                elif is_pg[F_idx]:   corr_co_pg += 1

        if matched_i % 100000 == 0 and matched_i:
            log(f'    processed {matched_i:,} / {len(matched_v3_snp_local_idx):,} matched SNPs '
                f'({time.time()-t0:.0f}s); '
                f'missing_total so far: {missing_total:,}')

    # Save raw per-SNP stats
    df = pd.DataFrame(snp_disagree_records)
    df.to_csv(ROOT / 'scratch/atomization_per_snp.tsv.gz', sep='\t', index=False,
              compression='gzip')

    # === REPORT ===
    log(f'\n=== Mechanism test result ===')
    log(f'Matched v3↔GN Chr1 SNPs:     {len(matched_v3_snp_local_idx):,}')
    log(f'Missing-carrier (founder, SNP) cases:   {missing_total:,}')
    log(f'Correct-listed (founder, SNP) cases:    {correct_total:,}')
    log(f'')
    log(f'P(carrier at INDEL in same snarl | missing)        = '
        f'{missing_with_indel_co:>10,} / {missing_total:>10,} = '
        f'{100*missing_with_indel_co/max(missing_total,1):6.2f}%')
    log(f'P(carrier at INDEL in same snarl | correctly-listed) = '
        f'{correct_with_indel_co:>10,} / {correct_total:>10,} = '
        f'{100*correct_with_indel_co/max(correct_total,1):6.2f}%')

    or_ = (missing_with_indel_co * (correct_total - correct_with_indel_co)) / max(
        (correct_with_indel_co * (missing_total - missing_with_indel_co)), 1)
    log(f'  Odds ratio (missing vs correct): {or_:.2f}×')
    log(f'  If hypothesis is right, this should be MUCH greater than 1.')
    log(f'  If close to 1, the indel-bearing-allele mechanism is NOT it.')

    log(f'\nStratified by founder class:')
    log(f'  cactus founders:')
    log(f'    P(indel-co | missing)  = {miss_co_cac:>8,} / {miss_cac:>8,} = '
        f'{100*miss_co_cac/max(miss_cac,1):6.2f}%')
    log(f'    P(indel-co | correct)  = {corr_co_cac:>8,} / {corr_cac:>8,} = '
        f'{100*corr_co_cac/max(corr_cac,1):6.2f}%')
    log(f'  PG founders:')
    log(f'    P(indel-co | missing)  = {miss_co_pg:>8,} / {miss_pg:>8,} = '
        f'{100*miss_co_pg/max(miss_pg,1):6.2f}%')
    log(f'    P(indel-co | correct)  = {corr_co_pg:>8,} / {corr_pg:>8,} = '
        f'{100*corr_co_pg/max(corr_pg,1):6.2f}%')

    log(f'\nTotal wall: {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
