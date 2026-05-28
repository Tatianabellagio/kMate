#!/bin/bash
# =============================================================================
# Submit the full p80 control experiment from scratch, chained via SLURM
# dependencies. Idempotent: each step skips if its outputs already exist.
#
# Pipeline (arch decomposition):
#   Shared prerequisite (NOT in this submit chain): arch3/chr1/jobA1_*.sh
#     -> chr1_135_annotated.sorted.vcf.gz + chr1_135_annotated_biallelic.sorted.vcf.gz
#     Built once on the full 135-asm pangenome + GFA. Already present.
#
#   01 -- build_biallelic_p80: subset 135-asm annotated to 80 Asm_IDs, run
#         convert-to-biallelic against arch3 biallelic catalog, fill-tags,
#         drop AC=0, reheader Asm_ID -> Acc_ID
#         -> canonical pangenome_p80_chr1.vcf.gz (biallelic, 80 samples)
#   03 -- cn_full_p80 (uses pang_135 PG-index from pangenie_genotyping/data/)
#   04 -- cn_var_p80
#   05 -- 80 founder consensus FASTAs
#   06 -- recomb sims for {n50_g1, n50_g3}
#   07 -- cactus_em {global, star2} per regime
#
# Usage:
#   bash scripts/submit_all_p80.sh
# =============================================================================
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/control_p80

# Verify the shared arch3 A1 prerequisite exists
A1_ANNOT=/global/scratch/users/tbellg/kmate/arch3/chr1/chr1_135_annotated.sorted.vcf.gz
A1_BIAL=/global/scratch/users/tbellg/kmate/arch3/chr1/chr1_135_annotated_biallelic.sorted.vcf.gz
for f in $A1_ANNOT $A1_BIAL; do
    [ -s "$f" ] || { echo "ERROR: missing arch3 A1 output $f" >&2; exit 1; }
done

A1=$(sbatch --parsable scripts/01_build_biallelic_p80.sh)
A3=$(sbatch --dependency=afterok:$A1 --parsable scripts/03_build_cn_full_p80.sh)
A4=$(sbatch --dependency=afterok:$A1 --parsable scripts/04_build_cn_var_p80.sh)
A5=$(sbatch --dependency=afterok:$A1 --parsable scripts/05_build_fastas_p80.sh)

B1=$(sbatch --dependency=afterok:$A3:$A4:$A5 --parsable scripts/06_run_sim_p80.sh 50 1)
B3=$(sbatch --dependency=afterok:$A3:$A4:$A5 --parsable scripts/06_run_sim_p80.sh 50 3)

C1G=$(sbatch --dependency=afterok:$B1 --parsable scripts/07_run_cactus_em_p80.sh n50_g1 global)
C1S=$(sbatch --dependency=afterok:$B1 --parsable scripts/07_run_cactus_em_p80.sh n50_g1 star2)
C3G=$(sbatch --dependency=afterok:$B3 --parsable scripts/07_run_cactus_em_p80.sh n50_g3 global)
C3S=$(sbatch --dependency=afterok:$B3 --parsable scripts/07_run_cactus_em_p80.sh n50_g3 star2)

# Optional: v3 panel on p80 reads (for apples-to-apples comparison vs v3).
# Uncomment to include in the chain.
# V3G1=$(sbatch --dependency=afterok:$B1 --parsable scripts/08_run_v3panel_on_p80reads.sh n50_g1 global)
# V3S1=$(sbatch --dependency=afterok:$B1 --parsable scripts/08_run_v3panel_on_p80reads.sh n50_g1 star2)
# V3G3=$(sbatch --dependency=afterok:$B3 --parsable scripts/08_run_v3panel_on_p80reads.sh n50_g3 global)
# V3S3=$(sbatch --dependency=afterok:$B3 --parsable scripts/08_run_v3panel_on_p80reads.sh n50_g3 star2)

cat <<EOF
Submitted job chain (arch decomposition; A1 reused from arch3/chr1/):
  A1   = $A1     -- subset to 80, convert-to-biallelic, fill-tags, drop AC=0, reheader
  A3   = $A3     -- cn_full_p80 (uses pang_135 PG-index)
  A4   = $A4     -- cn_var_p80
  A5   = $A5     -- 80 founder FASTAs
  B1   = $B1     -- sim n50_g1
  B3   = $B3     -- sim n50_g3
  C1G  = $C1G    -- cactus_em global on n50_g1
  C1S  = $C1S    -- cactus_em star2  on n50_g1
  C3G  = $C3G    -- cactus_em global on n50_g3
  C3S  = $C3S    -- cactus_em star2  on n50_g3

Monitor:
  squeue -u \$USER -t PD,R -o "%.10i %.18j %.10T %.8M %R"
EOF
