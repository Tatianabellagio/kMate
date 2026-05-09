#!/bin/bash
#SBATCH --job-name=panel_split
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/contamination_test/logs/split_chr%a.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/contamination_test/logs/split_chr%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-5

set -eo pipefail
CH="${SLURM_ARRAY_TASK_ID}"
/home/tbellagio/scratch/hapfire_sv/contamination_test/scripts/split_panel_per_chrom.sh "${CH}"
