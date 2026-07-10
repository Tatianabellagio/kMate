#!/bin/bash
#SBATCH --job-name=vg_packcov
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=4:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -eo pipefail
# Raw-pack-coverage SV AF: ref/alt READ SUPPORT straight from `vg pack`, NO
# `vg call` genotype-conditioning (the open vg#4443 problem). Allele node sets
# come from `vg paths -a` (the _alt_ variant paths), in the pack's node space.
#   1. vg paths -a -A/-F  -> _alt_ path node walks + sequences
#   2. vg paths -p Chr1 -A -> reference walk (for ref-allele node spans)
#   3. stage1_gaf -> map _alt_ -> panel svidx by ALT seq; emit node list + JSON
#   4. vg pack -d -N       -> per-node coverage (cov10 AND cov50 packs)
#   5. compute-af          -> AF = alt_support/(alt_support+ref_support)
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg   # vg env has NO python
PY=~/miniforge3/envs/basic/bin/python
export PYTHONPATH=/global/scratch/users/tbellg/kmate/src:$PYTHONPATH
T=${SLURM_CPUS_PER_TASK:-16}
W=benchmarks/accuracy_vs_competitors/work
SC=benchmarks/accuracy_vs_competitors/scripts
XGL=$W/c_cov10_n80_g0_s42_hotspots_p80_chr1.alts.xg
SVVCF=$W/panel_sv_p80_Chr1.vcf.gz

TMP=${TMPDIR:-/tmp}/packcov.$SLURM_JOB_ID; mkdir -p "$TMP"
echo "[$(date)] host=$(hostname) staging xg -> $TMP"
cp "$XGL" "$TMP"/g.xg
LXG=$TMP/g.xg

ALTGAF=$W/altpaths.gaf; ALTFA=$W/altpaths.fa; CHR1GAF=$W/chr1_ref.gaf
NODES=$W/packcov_nodes.txt; JSON=$W/packcov_rec.json

echo "[$(date)] STEP1 vg paths -a -A (node walks)"
[ -s "$ALTGAF" ] || vg paths -x "$LXG" -a -A -t "$T" > "$ALTGAF"
echo "[$(date)] STEP1b vg paths -a -F (sequences)"
[ -s "$ALTFA" ]  || vg paths -x "$LXG" -a -F > "$ALTFA"
echo "[$(date)] STEP1c vg paths (reference walk) — select Chr1 via names-file"
echo Chr1 > "$TMP/chr1name.txt"
[ -s "$CHR1GAF" ] || vg paths -x "$LXG" -p "$TMP/chr1name.txt" -A > "$CHR1GAF"
echo "[$(date)] paths done: gaf=$(wc -l <$ALTGAF) fa=$(grep -c '^>' $ALTFA)"

echo "[$(date)] STEP2 stage1 (map _alt_ -> svidx, build node sets)"
$PY "$SC/vg_packcov_stage1_gaf.py" --svlen 50 \
  --alt-gaf "$ALTGAF" --alt-fa "$ALTFA" --chr1-gaf "$CHR1GAF" \
  --panel-sv-vcf "$SVVCF" --out-nodelist "$NODES" --out-json "$JSON"

for COV in cov10 cov50; do
  if [ "$COV" = "cov10" ]; then PACK=$W/cov10_n80_g0_s42_hotspots_p80_chr1_construct.pack
  else PACK=$W/cov50_n80_g0_construct.pack; fi
  TAB=$W/packcov_${COV}_table.tsv
  OUT=$W/vg_packcov_sv_${COV}_n80_g0.tsv
  echo "[$(date)] STEP3 vg pack -d ($COV) over $(wc -l < $NODES) nodes"
  vg pack -d -x "$LXG" -i "$PACK" -N "$NODES" -t "$T" > "$TAB"
  echo "[$(date)] STEP4 compute-af ($COV)"
  $PY "$SC/vg_cov_sv.py" --svlen 50 compute-af \
    --pack-table "$TAB" --in-json "$JSON" --out-tsv "$OUT"
done
rm -rf "$TMP"
echo "[$(date)] DONE"
