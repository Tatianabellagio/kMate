# Pipeline FASTQ preprocessing — runbook for future big runs

**Status:** Established 2026-05-20 after audit of FASTQ-trim provenance and methodological-consistency review.
**Scope:** Applies to all pool-seq runs (SEEDMIX, GrENE-Net evolution samples) projected through cactus_em.

---

## TL;DR — two data classes, two different Trimmomatic configs

**This is intentional, not a bug.** xwu's GrENE-Net pipeline uses two different Trimmomatic configs depending on data class:

| Data class | What it is | Trimmomatic config | Why |
|---|---|---|---|
| **Pool-seq** | SEEDMIX (S1..S8); GrENE evolution pools | `ILLUMINACLIP + SLIDINGWINDOW:4:20 + LEADING:5 + TRAILING:5 + MINLEN:36` | PCR-free Lucigen libraries are high-Q; stricter trimming is affordable and removes more error-prone tail bases |
| **Single-individual short reads** | 1001G ecotype FASTQs (for PG-genotyping founder substitutes; the 151 PG_153 panel) | `ILLUMINACLIP + LEADING:5 + TRAILING:5 + MINLEN:36`  **(no SLIDINGWINDOW)** | 1001G data has lower-quality 3' tails; SLIDINGWINDOW drops ~4.7% of reads below MINLEN, losing too much coverage |

**Then both classes** get the same downstream step:

| Step | Tool | Params |
|---|---|---|
| Clumpify dedup | `clumpify.sh` (BBTools) | `dedupe=t dupesubs=0 optical=f` |

**Critical hygiene point**: Clumpify dedup was **missing** from prior SEEDMIX cactus_em runs. Adding it removes ~35% PCR-duplicate k-mer inflation. This is the **principal fix**, not the Trimmomatic params (which were already correct for SEEDMIX via xwu's re-trim).

The 151 founder side already gets clumpify via `panel/pangenie_genotyping/scripts/preprocess_one.sh`.

---

## 1. Trimmomatic — exact config to use

**Pool-seq data (SEEDMIX, evolution pools)**: use xwu's GrENE pool-seq config.
Source: `/carnegie/nobackup/scratch/xwu/GrENE_net/seed_mix/re-trimmed/commands.sh`

```
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True
SLIDINGWINDOW:4:20
LEADING:5
TRAILING:5
MINLEN:36
```

**Single-individual ecotype data (1001G short reads, founder substitution for PG_153 panel)**:
use xwu's 1001G/SRA config — **NO SLIDINGWINDOW**.
Source: `/carnegie/nobackup/scratch/xwu/GrENE_net/vcf/sra/commands.sh`

```
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True
LEADING:5
TRAILING:5
MINLEN:36
```

This matches `panel/pangenie_genotyping/scripts/preprocess_one.sh`.

**Why different params for pool-seq vs ecotype**: xwu's deliberate choice; the 1001G data has lower-quality 3' tails, and SLIDINGWINDOW drops too many reads (4.7% empirically) → loses coverage at edges. Pool-seq data (PCRfreeLucigen) is higher quality, can afford stricter trimming.

**Why this doesn't break apples-to-apples** vs kmer_pa: kmer_pa's founder k-mers are **FASTA-derived** (consensus sequences from VCFs), not FASTQ-derived. The trim choice on founder FASTQs only propagates indirectly through PanGenie GT calls for the 151 PG founders. See `memory/project_kmer_pipeline_provenance.md`.

**Tool path**:
- jar: `/home/tbellagio/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/trimmomatic.jar`
- adapter PE: `/home/tbellagio/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa`
- env: `sequencing_pipeline`

**Resources**: 4 CPU, 4 GB RAM, ~30 min per sample at typical SEEDMIX coverage.

---

## 2. Clumpify — exact config to use

Match `panel/pangenie_genotyping/scripts/preprocess_one.sh:118-125`:

```bash
clumpify.sh in=$TRIM_R1 in2=$TRIM_R2 \
    out=$DEDUP_R1 out2=$DEDUP_R2 \
    dedupe=t dupesubs=0 optical=f
```

**Flag rationale**:
- `dedupe=t` — enable deduplication
- `dupesubs=0` — exact-match dedup (no allowance for sequencing errors creating "near-duplicate" pairs; conservative)
- `optical=f` — disable optical-duplicate detection; we want **all PCR duplicates** (typically ~30-40%), not just optical (~1-2%)

**Tool path**: `/home/tbellagio/miniforge3/envs/pang/bin/clumpify.sh` (BBTools install in `pang` conda env)

**Resources**: 4 CPU, 32 GB RAM, ~30-60 min per sample.

**Expected dedup rate**:
- SEEDMIX: ~35% (measured: 2.9 GB → 1.9 GB on S1)
- 1001G ecotypes: variable, typically 10-25% depending on library prep

---

## 3. kmer_pa rebuild — when it's needed (and when it's not)

**kmer_pa does NOT need to be rebuilt when**:
- SEEDMIX FASTQs are reprocessed (clumpify dedup added). kmer_pa is founder-side only.
- New pool-seq samples are added. kmer_pa is sample-agnostic.
- Trimmomatic params change for pool-seq. Same reason.

**kmer_pa DOES need to be rebuilt when**:
- The founder panel changes (different ecotypes added/removed).
- The cactus pangenome graph is rebuilt (different bubble structure → different consensus).
- The TAIR10 reference changes.
- The 151 PG founder VCF changes (PanGenie re-genotyped with different reads or different graph).
- The k-mer size (default 31) changes.

**If kmer_pa DOES need rebuilding** — example: new founder substituted into PG_153:
1. Re-process that founder's raw FASTQs through `preprocess_one.sh` (ecotype Trimmomatic + clumpify).
2. Re-run PanGenie genotyping (`panel/pangenie_genotyping/scripts/pangenie_one.sh`) → new per-sample VCF.
3. Re-merge into the 153-sample PG-panel VCF, then re-merge with cactus_78 (see `scratch/arch3_chr1/jobA1-A4` for current Arch 3 merge pipeline).
4. Re-run `build_kmer_pa.py` (founder k-mer matrix from VCF + TAIR10 ref) → new kmer_pa.
5. Re-run `build_var_pa.py` (founder × variant carrier matrix from biallelic VCF) → new var_pa.

---

## 4. Re-running cactus_em SEEDMIX on dedup'd FASTQs

**Goal**: refit `h` vector for each SEEDMIX_S{1..8} using clumpify-dedup'd FASTQs, then re-project against existing var_pa (or new Arch 3 var_pa).

**Prerequisite**: B1 array (trim+clumpify) complete. Outputs at
`/home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S{N}_{1,2}.dedup.fq.gz`.

**Steps**:

1. Update SEEDMIX SLURM script(s) to point at the dedup'd FASTQ paths:
   ```bash
   READS=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup
   R1=$READS/SEEDMIX_S${S}_1.dedup.fq.gz
   R2=$READS/SEEDMIX_S${S}_2.dedup.fq.gz
   ```
   Already done for `scripts/archive/run_seedmix_v3qc_v3_mixedloose_chr1.sh` on 2026-05-20.
   For genome-wide runs, update the corresponding `_genomewide.sh` scripts.

2. Submit cactus_em per-sample, per-chrom (one SLURM array task per (sample, chrom)):
   - 8 samples × 5 chroms = 40 tasks
   - Per task: 8 CPU, 80 GB, ~2 h
   - Outputs: `<OUT_DIR>/SEEDMIX_S{N}.tsv` (per-record AF) and
     `<OUT_DIR>/SEEDMIX_S{N}.h_per_chrom.npz` (per-chrom h vectors, dict keyed by 'Chr1'..'Chr5' + 'founders')

3. Compare new h vectors against existing (pre-dedup) h to quantify the dedup effect.
   Compare new AF against xwu hapFIRE (canonical 231-panel reference at
   `xwu/GrENE_net/hapFIRE_frequencies/seed_mix/s{N}_snp_frequency.txt`) on the
   4-tuple `(chrom, pos, ref, alt)` — see `scratch/arch3_chr1/jobA7_compare_vs_hapfire.sh`
   for the join template.

**What to expect after dedup**:
- Total k-mer counts drop ~35% uniformly (PCR duplicates removed).
- Per-record AF should change MOST at SNPs in PCR-duplicate hotspots (typically
  high-GC, near repeats, near rRNA-like loci). MAE shift vs hapFIRE is hard to
  predict — could be 0.001-0.01 globally.
- Effective coverage drops from ~50× to ~33× per sample. Still well above the
  ~20× threshold where cactus_em becomes simplex-limited.

---

## 5. Order of operations for a clean future run

For a full reprocessing of SEEDMIX from scratch:

```
1. Trimmomatic (pool-seq config)        ← per_sample SLURM array, ~30 min × 8 = 4 h wall (parallel)
2. Clumpify dedup                       ← per_sample SLURM array, ~30 min × 8 = 4 h wall (parallel)
3. cactus_em SEEDMIX per-chrom          ← per_(sample,chrom) SLURM array, ~2 h × 40 = 80 h CPU
4. Compare vs hapFIRE on 4-tuple join   ← single job, ~30 min
```

Steps 1-2 can be combined into one SLURM script per sample (as in
`scratch/arch3_chr1/jobB1_seedmix_trim_clumpify.sh`, currently set up to
**skip** step 1 because xwu's re-trimmed FASTQs are byte-identical to what
xwu's Trimmomatic pool-seq config would produce — so we just clumpify
those directly. For TRULY fresh data (no xwu-trimmed intermediate available),
restore the full trim step from the earlier version of jobB1 / from
`preprocess_one.sh`'s trim block).

Total wall time for one round on 8 samples × 5 chroms: ~1 day with full SLURM parallelism.

---

## 6. Provenance — exactly where FASTQs come from RIGHT NOW

### SEEDMIX pool-seq (S1..S8)

| Tier | Path on disk | Trimmomatic? | Clumpify? | Compute-node accessible? | Status |
|---|---|---|---|---|---|
| **Raw** | `/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/ena_seeds/illumina_ST-J00101_flowcellA_SampleIdDNAfromseeds{N}_RunId0102_LaneId3/PCRfreeLucigen_S{N}_L003_R{1,2}_001.fastq.gz` | NO | NO | **NO — Carnegie-DPB share not mounted on memex compute nodes** | Cannot use directly |
| **Trim-only** (xwu's "re-trimmed" pool-seq output) | `/home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S{N}-1.{1,2}_P.fq.gz` | YES (xwu's pool-seq config, with SLIDINGWINDOW:4:20) | NO | YES | **What every prior SEEDMIX cactus_em run actually used (v3, v3qc, v3qc_v3, etc.)** — verified byte-identical (md5) to `/carnegie/nobackup/scratch/xwu/GrENE_net/seed_mix/re-trimmed/S{N}-1.{1,2}_P.fq.gz` |
| **Trim + clumpify-dedup** (NEW, 2026-05-20) | `/home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S{N}_{1,2}.dedup.fq.gz` | YES (xwu's pool-seq config, inherited) | YES (added by jobB1) | YES | What going forward all SEEDMIX cactus_em runs should use |

Since the Carnegie-DPB raw share is invisible to compute nodes, our jobB1 (`scratch/arch3_chr1/jobB1_seedmix_trim_clumpify.sh`) reads from the **trim-only** tier and only adds clumpify on top. The Trimmomatic step is effectively skipped because xwu already did it with the exact pool-seq config we'd run.

### 151 PG founder ecotypes (single-individual 1001G short reads)

| Tier | Path on disk | Trimmomatic? | Clumpify? | Compute-node accessible? | Status |
|---|---|---|---|---|---|
| **Raw** (ENA-downloaded or xwu_BAM-derived) | `panel/pangenie_genotyping/data/raw_fastqs/{ECOTYPE}/{RUN}_{1,2}.fastq.gz` (ENA), or `panel/pangenie_genotyping/data/raw_fastqs/{ECOTYPE}/{ECOTYPE}_{1,2}.fastq.gz` (xwu_BAM, lane-concat'd) | NO | NO | YES | Source for `preprocess_one.sh` |
| **Trim + clumpify-dedup** | `panel/pangenie_genotyping/data/preprocessed/{ECOTYPE}_{1,2}.dedup.fq.gz` | YES (xwu's 1001G config, NO SLIDINGWINDOW) | YES (`dedupe=t dupesubs=0 optical=f`) | YES | What PanGenie consumes via `pangenie_one.sh`; what built the 151-PG-founder VCFs that feed kmer_pa |

The founder pipeline already does trim + clumpify correctly — no fix needed there.

### Where the two streams meet

- **Founders' k-mers in kmer_pa**: built FASTA-side from VCFs, NOT directly from these FASTQs. The founder FASTQ → PanGenie → VCF chain produces clean consensus sequences which then get k-merized by `build_kmer_pa.py`. So founder FASTQ trim choices propagate through PanGenie GT calls only (modest effect; PanGenie has internal noise tolerance).
- **SEEDMIX k-mers**: counted directly from FASTQs via jellyfish inside `per_sample_per_chrom.py`. SEEDMIX trim/dedup choices directly affect which k-mers enter the EM. This is where clumpify matters most.

---

## ✅ 2026-05-20 04:10 — RESOLVED: SEEDMIX libraries are PCR-FREE, no dedup needed

### What I initially thought

The 35% file-size reduction in `seed_mix_trimmed/dedup/seeds-1_1_dedup.fq.gz` (1.9 GB) vs `seed_mix_trimmed/seeds-1.1.fastq.gz` (2.89 GB) looked like a 35% PCR-duplicate rate that B1's clumpify dedup should recover.

### What's actually going on

Counted reads in BOTH files directly:

```
seeds-1.1.fastq.gz        :  39,126,655 reads  (2.7 GB)
seeds-1_1_dedup.fq.gz     :  38,458,036 reads  (1.8 GB)
  → real read-count dedup :   1.71%
  → file-size shrink      :  34.24%
```

**The 34% size reduction is entirely from clumpify reordering reads for better gzip compression, NOT actual deduplication.**

B1's 1.33-1.70% dedup rates across S1-S8 are the **true PCR-duplicate rates** for this data, not artifacts of SLIDINGWINDOW-driven length variation. My earlier "SLIDINGWINDOW breaks exact-match dedup" hypothesis was wrong; clumpify is finding all the duplicates that exist.

### Why so few duplicates

The raw FASTQ filenames contain **`PCRfreeLucigen_S{N}_*`** — these are **PCR-free Lucigen library prep**. PCR-free libraries lack the PCR amplification step that creates duplicates. The residual ~1.5% is optical duplicates (tile-edge) + chance sequence collisions, exactly what's expected for a PCR-free Illumina library.

### Implication for the pipeline

**Clumpify dedup on SEEDMIX is essentially a no-op.** The trim-only FASTQs at `seed_mix/S{N}-1.{1,2}_P.fq.gz` and the trim+clumpify FASTQs at `seed_mix_trimdedup/SEEDMIX_S{N}_{1,2}.dedup.fq.gz` differ by <2% in actual content. Any cactus_em h-estimate difference between them will be at noise level.

**Recommended action**:
- For SEEDMIX (PCR-free Lucigen libs): no need to repopulate the production script with clumpify-deduped paths. The current trim-only path was fine all along. The 35% PCR-duplicate concern that motivated this whole exercise didn't actually exist.
- The B1-generated dedup files at `seed_mix_trimdedup/` can stay on disk as a 1.7%-deduped variant — useful as a control but not load-bearing.
- **Revert** the `run_seedmix_v3qc_v3_mixedloose_chr1.sh` path change if you want to stay on the version that has actual prior-run h-vector compatibility. Or leave it — barely matters.

### What this means for evolution pool-seq samples (NOT PCR-free)

The GrENE-Net evolution pool-seq samples may use PCR-based library prep, in which case PCR duplicates DO exist (typically 20-40%) and clumpify dedup IS necessary. For those samples:
- Confirm library prep type from sample metadata before running.
- If PCR-based: dedup matters, and the SLIDINGWINDOW-breaks-exact-match-dedup concern (which I raised then ruled out for SEEDMIX) may matter. Test dedup rate on a single sample first to verify clumpify is finding duplicates.
- If PCR-free (like Lucigen): skip clumpify; it's a no-op.

### What I got right vs wrong

- **Right**: The kmer_pa FASTA-derived k-mer provenance, the xwu trim-by-data-class split, the lack of SLIDINGWINDOW being intentional for ecotype data.
- **Wrong**: Assumed the 35% size delta in `seed_mix_trimmed/dedup/` represented 35% read-count dedup. Should have counted reads first. Sloppy.
- **The whole "PCR duplicates are inflating SEEDMIX k-mer counts" worry**: false premise; PCR-free libraries don't have meaningful PCR duplicates.

---

## 7. Open empirical questions for future investigation

1. **PG founder GT-call sensitivity to SLIDINGWINDOW** — does PanGenie's k-mer-based genotyping miscall founders when reads are not SLIDINGWINDOW-trimmed? Empirical Q-analysis shows TRAILING:5 already eats most low-Q tail, so net effect likely small, but unmeasured.
2. **Spatial structure of duplicates** — is PCR duplication uniform across the genome, or hotspot-concentrated? Per-locus dedup rate would tell us whether the AF shift is uniform (then h-invariant) or biased.
3. **Founder-side clumpify benefit** — `preprocess_one.sh` already clumpifies founders; verify the dedup rate is comparable to what we see on SEEDMIX (~35%) or whether 1001G libraries have different duplicate profiles.
