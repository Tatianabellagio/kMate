# Seedmix contamination test (1141-panel hapFIRE re-run)

## Question

In the GrENE-Net seedmix, was there contamination from ecotypes outside the
expected 231 founders?

## Method

Re-run hapFIRE on the 8 SEEDMIX BAMs against an **expanded 1141-ecotype
founder VCF** (Xing's pre-subset 1001G + 80pilot + israel + regmap panel
which contains all 231 GrENE-Net founders + 910 other 1001G ecotypes).

If the seedmix is pure, hapFIRE's CVXPY ecotype projection should put nearly
all the mass on the 231 expected founders. Mass on the 910 non-231 ecotypes
is the contamination signal — modulo CVXPY identifiability slack at F=1141
(see "Known caveat" below).

## Inputs

| | Path |
|---|---|
| Panel VCF (1141 samples, all 5 chroms) | `/carnegie/nobackup/scratch/xwu/GrENE_net/vcf/1001G_80pilot_israel_regmap_overlapping.vcf.gz` |
| Reference (TAIR10, chrom names `1..5,Mt,Pt`) | `/home/tbellagio/scratch/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa` |
| Seedmix BAMs (S1..S8, already aligned + filtered) | `/carnegie/nobackup/scratch/xwu/GrENE_net/seed_mix/filtered_bam/filtered_seeds-{1..8}.bam` |
| 231 GrENE-Net founder list | `/home/tbellagio/scratch/hapfire_sv/data/vcf_samples_231.txt` |

Per-chrom panel splits (built from the 1141 VCF, biallelic SNPs only,
fully-phased non-missing only): `vcf/panel_1141_chr{1..5}.vcf`.

Sample list (1141 IDs): `vcf/samples_1141.txt`.

## Pipeline

1. **Split panel VCF** by chrom → 5 plain-text VCFs that hapFIRE.py can `open()`.
   `scripts/split_panel_per_chrom.sh <chrom>`, submitted as
   `scripts/sbatch_split.sh` (array 1-5).

2. **Run hapFIRE per (sample, chrom)** = 40 jobs total.
   `scripts/run_hapfire_perchrom.sh <SEEDMIX_S{1..8}> <chrom 1..5>`.
   - Resources per job: 4 CPUs, 128 GB, 10 h, memex partition.
   - Output: `results/<SAMPLE>_chr<N>/<SAMPLE>_chr<N>_ecotype_frequency_selected.txt`
     (1141 rows, sums to 1).

3. **Job stats**: `scripts/collect_stats.sh` pulls SLURM accounting (Elapsed, MaxRSS, NodeList, Start/End) for every `hapfire_1141` job and writes `results/job_stats.tsv`. Re-run anytime to refresh; running it after fan-out completes produces final MaxRSS for all 40 jobs. Includes a per-chrom summary line at the bottom.

4. **Aggregate + diagnose**:
   `scripts/aggregate_and_diagnose.py` — averages per-chrom ecotype vectors
   per sample (simple mean across chroms, mirrors hapFIRE.py's internal
   averaging), splits mass into `in231 / non231`, ranks the 910 non-231
   ecotypes by mean mass.

   Outputs:
   - `results/per_sample_ecotype_freq_1141.tsv` — full 1141 × 8 matrix
   - `results/contamination_summary.tsv` — per-sample headline numbers
   - `results/non231_top.tsv` — ranked non-231 ecotype contributions

## Submission

```bash
cd /home/tbellagio/scratch/hapfire_sv/contamination_test
# 1. split panel (5 jobs, ~10 min each)
sbatch scripts/sbatch_split.sh

# 2a. smoke (S1 chr5 only) before fan-out
bash scripts/submit_all.sh smoke

# 2b. fan-out to all (8 samples x 5 chroms = 40 jobs)
bash scripts/submit_all.sh all

# 3. refresh job stats (anytime)
bash scripts/collect_stats.sh

# 4. aggregate ecotype frequencies
/home/tbellagio/miniforge3/envs/hapfm/bin/python scripts/aggregate_and_diagnose.py
```

## Known caveat: identifiability slack at F=1141

CVXPY's ecotype projection at F=1141 has more rank-deficiency slack than at
F=231. Many 1001G ecotypes are genetically very similar to GrENE-Net
founders; CVXPY will distribute some mass to those "twin" ecotypes even
when no real contamination is present.

A **null calibration** (run the same pipeline on a synthetic uniform-231
pool simulated by VISOR) would establish the spillover floor. Decision:
deferred — first see what the real numbers look like before deciding if the
null sim is needed for interpretation.
