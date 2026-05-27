# Imputation framework — extending G_SV from 80 → 231 founders

**Goal**: impute SV genotypes for the 151 GrENE-Net founders that don't have long-read assemblies, using SV-SNP LD trained on the 80 we have.

**Why**: closes the 65% missing-mass gap that bounds SEEDMIX r² to ~0.81. After imputation, our pool-seq frequency model estimates over all 231 founders (matching hapFIRE's setup, where they reach r²=0.99 against the same recipe-based truth).

---

## Why imputation is the right answer

The pool-seq frequency model fits founder freqs `h` on the simplex spanned by the panel founders. Mass on founders NOT in the panel has nowhere to go and gets falsely attributed to in-panel founders. With 65% of seedmix mass on the 151 unobserved GrENE-Net founders, the in-panel `h` estimates are inflated by ~3× — this is the structural cause of the r=−0.15 we saw on real SEEDMIX data.

hapFIRE side-steps this because they estimate `h` over all 231 GrENE-Net founders (their SNP panel covers all 231). To match their setup, we need our cactus SV panel to also cover 231. We don't have long-read assemblies for the 151 missing founders, but we DO have their SNP genotypes from the GrENE-Net VCF. **Imputation lets us infer their SV genotypes from SNP-haplotype linkage.**

The biological assumption: **SV genotypes are correlated with surrounding SNP haplotypes within recombination scale**. If two founders share a haplotype in a region, they should also share the SV genotype within that region. Beagle's HMM exploits this exactly.

---

## Inputs

| File | Description | Path |
|---|---|---|
| Cactus pangenome VCF | 82 founders × 4.45M variants (SNPs + SVs), TAIR10-based | `/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz` |
| GrENE-Net SNP VCF | 231 founders × ~3.2M SNPs, TAIR10-based, chroms "1/2/3/4/5" | `/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_updatedVCF/greneNet_final_v1.1.recode.vcf` |
| Assembly→1001G mapping | `Assembly_ID → 1001G Accession_ID` | `data/sv_panel_to_accession_id.tsv` |
| GrENE-Net 231 sample list | 1001G IDs for the 231 founders | `data/vcf_samples_231.txt` |

**Sample overlap**: 80 of 82 cactus assemblies map to 1001G IDs in GrENE-Net 231. The 2 non-overlap (cactus accessions not in GrENE-Net) are dropped — they have no GrENE-Net SNP haplotype to anchor the imputation.

---

## Reference-genome verification

Both VCFs are TAIR10-based, just with different chromosome naming conventions:

| | Cactus | GrENE-Net |
|---|---|---|
| Chrom name | `Chr1`, ..., `Chr5` | `1`, ..., `5` |
| Chr1 contig length declared | 30,427,671 ✓ (matches TAIR10) | not declared, but positions are TAIR10-consistent (verified by SNP overlap) |
| Reference field | not present | not present |

We rename cactus chroms `Chr1→1` etc. to align with GrENE-Net convention.

**SNP overlap test (Chr1 first 1Mb)** confirmed coordinate consistency:
- Cactus: 28,748 biallelic SNPs in 1Mb
- GrENE-Net: 21,488 biallelic SNPs in 1Mb
- Position overlap: 14,953 (~70% of GrENE)
- Position + REF + ALT exact match: 14,135 (95% of position-matched)

The 5% allele mismatch at shared positions reflects calling differences (long-read vs short-read) — small enough to ignore. The overall agreement confirms both pipelines use the same reference coordinates.

---

## Architecture — Beagle 5.5 with two parallel data streams

```
Reference panel (80 founders, all phased, no missing GTs):
    grene_80 SNPs (3.2M)   +   cactus_80 SVs (240K, 50bp–50kb)
                                 ↓ concat + sort + filter F_MISSING
                                ref_80.vcf.gz

Target panel (151 GrENE-only founders):
    grene_151 SNPs (3.2M)  +   target_151 SVs MASKED (./.)
                                 ↓ concat + sort
                                target_151.vcf.gz

                                 ↓ Beagle imputation per chromosome
                                imputed_151.vcf.gz
```

**Why this works architecturally**:

1. **Both panels share the SAME GrENE-Net SNPs at the SAME positions.** SNPs come from one source VCF (just different sample subsets), so haplotype patterns are perfectly comparable between ref and target.

2. **Beagle traces SNP haplotypes from the ref to the target.** For each target sample, Beagle finds which ref founders' SNP haplotypes match best in each window, then assigns the matching ref founder's SV genotype to the target sample. This is the Li–Stephens HMM applied per locus.

3. **Cactus SVs are at distinct positions from GrENE SNPs** — they don't compete or conflict. The merged panel just has more rows.

4. **We use ONLY cactus SVs (50bp–50kb), not cactus SNPs.** Cactus SNPs would be redundant with GrENE SNPs (already cover all 231) and could create allele-conflict at the ~5% of position-matched-but-allele-mismatched sites.

---

## Pipeline scripts (in execution order)

| Script | Purpose | Output |
|---|---|---|
| `01_prepare_inputs.sh` | rename samples (Assembly→1001G), rename chroms (Chr→1), subset to 80 GrENE-overlap × SVs (50bp–50kb) | `cactus_svs_renamed.vcf.gz` |
| `02_build_merged_panel.sh` | build `ref_80` and `target_151` VCFs | `ref_80.vcf.gz`, `target_151.vcf.gz` |
| `03_run_imputation.sh` | run Beagle per chrom, concat results | `imputed_151.vcf.gz` |
| `04_validate_loo.sh` | leave-one-out validation on 10 random founders | concordance report |
| `05_build_imputed_cn.py` | merge cactus_svs (80) + imputed_151 → 231 × N_SVs cn matrix (variant-level) | `cn_imputed_231.cn.npz` |
| `06_build_231_cn_matrices.sh` | merge VCFs, rename chroms back, build full 231-founder `cn_kmer` and `cn_var` for the production pool-seq pipeline | `cn_full_231/cn_Chr*.cn.npz`, `cn_var_231.cn_var.npz` |

---

## Two non-obvious bugs we found and fixed

These would have killed silent runs at scale; smoke testing on Chr1:1–2Mb caught them:

### 1. `bcftools +setGT -i` requires `-t q`

What we tried (wrong):
```bash
bcftools +setGT input.vcf.gz -- -t a -i 'TYPE="indel"' -n .
```
Fails because `-i` (include filter) requires `-t q` (target = filtered subset). For our use case (mask all genotypes in a SV-only VCF), the right form is:
```bash
bcftools +setGT input.vcf.gz -- -t a -n .
```

### 2. Beagle reference must be FULLY non-missing

Beagle 5.x phased-reference panels reject any missing genotype (`./.` or `.|.`) anywhere in the reference. Cactus output has missing GTs at sites with multi-path-conflict (CONFLICT INFO field). Fix:
```bash
bcftools view -e 'F_MISSING > 0' merged.vcf.gz
```
Drops records where ANY sample has missing GT. Dropped a small fraction; everything else was fine.

---

## Validation strategy

**Smoke test (passed 2026-04-27)**: Chr1:1–2Mb, 9 ref founders + 1 leave-out → **96.8% concordance** (2,083/2,152 SVs match cactus truth, 0 missing). Beagle runtime: 3 seconds.

**Full leave-one-out (`04_validate_loo.sh`)**: 10 random founders × full genome. Each founder's SVs masked and re-imputed using the other 79 + GrENE SNPs. Aggregate per-SV concordance and per-founder dosage R².

Expected at full scale:
- Common SVs (AC ≥ 10 in panel): >99% concordance (dense SNP linkage)
- Rare SVs (AC ≤ 5 in panel): 70–90% concordance (sparser linkage)
- Overall mean: comparable to the 96.8% smoke test result, possibly slightly higher because larger ref panel = denser haplotype graph

---

## Integration with the pool-seq pipeline

After imputation, the production cn matrices are rebuilt over 231 founders (script `06`):

| Old (82-founder) | New (231-founder) |
|---|---|
| `cn_full_<chrom>.cn.npz` | `cn_full_231/cn_<chrom>.cn.npz` |
| `cn_var_82.cn_var.npz` | `cn_var_231.cn_var.npz` |

Both have shape `(231, K)` instead of `(82, K)`. The same set of K is used (PanGenie unique-k-mer set is panel-independent — the cn just has 151 more rows for the imputed founders).

The pool-seq pipeline (block solver, EM, projection) runs unchanged. The simplex constraint just operates over a 231-dim founder space instead of 82-dim.

**Expected outcome on real SEEDMIX**: r vs recipe should jump from −0.15 (82-founder, missing-mass distorted) to ~0.99 (231-founder, no missing-mass) — matching hapFIRE.

---

## Outputs (final)

- `cn_imputed_231.cn.npz` — 231 × n_SVs sparse matrix of imputed SV dosages (0/1, integer)
- `cn_full_231/cn_<chrom>.cn.npz` — per-chrom 231 × K_kmer matrices for the EM solver
- `cn_var_231.cn_var.npz` — 231 × N_records matrix for projection to per-record alt-allele frequencies
- LOO concordance report from `04_validate_loo.sh`
- Per-sample alt-allele freq TSVs (from `per_sample_driver.py` running over the 231-founder cn)
