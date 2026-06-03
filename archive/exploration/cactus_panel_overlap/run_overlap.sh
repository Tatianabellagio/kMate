#!/bin/bash -l
#SBATCH --job-name=cactus_overlap
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cactus_overlap_%j.out
#SBATCH --error=logs/cactus_overlap_%j.err

# =============================================================================
# Cactus pangenome (82 founders) ↔ xwu GrENE-Net (231 founders) SNP overlap.
#
# Mirrors freqk_gr/panel_overlap_test architecture:
#   1. Filter xwu 231-panel to 80 and 82 cactus-overlap subsets, drop monomorphic
#   2. Extract SNP position lists from cactus VCF and from xwu subsets
#   3. Compute overlap (# shared sites, % cactus in xwu, % xwu in cactus)
#
# Question: how much do cactus's natural SNPs add over xwu's GrENE-Net SNPs?
# This decides whether to include cactus SNPs in the imputation reference panel.
# =============================================================================
mkdir -p logs
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pang

BASE=/global/scratch/users/tbellg/hapfire_sv/cactus_panel_overlap
DATA=$BASE/data
KEEP_80=$DATA/ecotype_ids_to_keep.txt        # 80 raw-ID matches
KEEP_82=$DATA/ecotype_ids_to_keep_82.txt     # 82 with Ct-1/No-0 overrides

XWU_VCF=/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf
XWU_BGZ=/global/scratch/users/tbellg/hapfire_sv/imputation/work/grene_all.vcf.gz

CACTUS_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz

# 1) Reuse already-bgzipped xwu VCF from imputation work dir
if [[ ! -s $XWU_BGZ ]]; then
  echo "[$(date)] bgzipping xwu VCF (would take 10 min)"
  bgzip -c -@ 4 "$XWU_VCF" > "$XWU_BGZ"
  tabix -p vcf "$XWU_BGZ"
fi
echo "[$(date)] xwu source: $XWU_BGZ"

# 2) Subset xwu to 80 cactus-overlap, drop monomorphic
SUB80=$DATA/xwu_subset80.vcf.gz
if [[ ! -s $SUB80 ]]; then
  echo "[$(date)] subset xwu to 80 cactus-overlap accessions"
  bcftools view -S "$KEEP_80" --force-samples "$XWU_BGZ" -Ou \
    | bcftools view -e 'AC=0 || AC=AN' -Oz -o "$SUB80"
  tabix -p vcf "$SUB80"
fi

# 3) Subset xwu to 82 cactus-overlap (with Ct-1/No-0 overrides)
SUB82=$DATA/xwu_subset82.vcf.gz
if [[ ! -s $SUB82 ]]; then
  echo "[$(date)] subset xwu to 82 cactus-overlap accessions"
  bcftools view -S "$KEEP_82" --force-samples "$XWU_BGZ" -Ou \
    | bcftools view -e 'AC=0 || AC=AN' -Oz -o "$SUB82"
  tabix -p vcf "$SUB82"
fi

# 4) Position lists. xwu chroms are "1..5"; cactus chroms are "Chr1..Chr5".
echo "[$(date)] dump position lists"
bcftools query -f '%CHROM\t%POS\n' "$XWU_BGZ" > $DATA/xwu231_positions.tsv
bcftools query -f '%CHROM\t%POS\n' "$SUB80" > $DATA/xwu80_positions.tsv
bcftools query -f '%CHROM\t%POS\n' "$SUB82" > $DATA/xwu82_positions.tsv
echo "  xwu 231: $(wc -l < $DATA/xwu231_positions.tsv) sites"
echo "  xwu 80-subset: $(wc -l < $DATA/xwu80_positions.tsv) sites"
echo "  xwu 82-subset: $(wc -l < $DATA/xwu82_positions.tsv) sites"

# 5) Cactus SNP positions only (filter to SNPs: REF and ALT both length 1).
#    cactus chroms are Chr1..Chr5 — strip "Chr" so they match xwu's 1..5
echo "[$(date)] dump cactus SNP positions"
bcftools view -v snps "$CACTUS_VCF" -Ou \
  | bcftools query -f '%CHROM\t%POS\n' \
  | sed 's/^Chr//' \
  | sort -u > $DATA/cactus82_snp_positions.tsv
echo "  cactus 82 SNP sites (unique chrom,pos): $(wc -l < $DATA/cactus82_snp_positions.tsv)"

# 6) Cactus indel/SV position list too (for completeness)
echo "[$(date)] dump cactus indel/SV positions"
bcftools view -V snps "$CACTUS_VCF" -Ou \
  | bcftools query -f '%CHROM\t%POS\n' \
  | sed 's/^Chr//' \
  | sort -u > $DATA/cactus82_indel_sv_positions.tsv
echo "  cactus 82 non-SNP sites: $(wc -l < $DATA/cactus82_indel_sv_positions.tsv)"

# 7) Compute overlap inline
echo ""
echo "[$(date)] === OVERLAP COMPUTATION ==="
python3 <<'PYEOF'
import os
DATA = "/global/scratch/users/tbellg/hapfire_sv/cactus_panel_overlap/data"

def load_pos(path):
    s = set()
    with open(path) as f:
        for line in f:
            c, p = line.rstrip("\n").split("\t")
            s.add((c, p))
    return s

print("  loading positions ...")
xwu231 = load_pos(f"{DATA}/xwu231_positions.tsv")
xwu80  = load_pos(f"{DATA}/xwu80_positions.tsv")
xwu82  = load_pos(f"{DATA}/xwu82_positions.tsv")
cactus = load_pos(f"{DATA}/cactus82_snp_positions.tsv")
print(f"  xwu231: {len(xwu231):,}")
print(f"  xwu80:  {len(xwu80):,}")
print(f"  xwu82:  {len(xwu82):,}")
print(f"  cactus SNPs: {len(cactus):,}")

print()
print(f"=== Cactus 82 SNPs vs xwu 231 ===")
shared = xwu231 & cactus
print(f"  shared positions: {len(shared):,}")
print(f"  % of xwu in cactus: {100*len(shared)/len(xwu231):.2f}%")
print(f"  % of cactus in xwu: {100*len(shared)/len(cactus):.2f}%")
print(f"  cactus-only (NOT in xwu): {len(cactus - xwu231):,}  ← these are EXTRA SNPs we get from cactus")
print(f"  xwu-only (NOT in cactus): {len(xwu231 - cactus):,}  ← what we'd lose if we dropped xwu and used only cactus")

print()
print(f"=== Cactus vs xwu 82-subset (apples-to-apples on the 82 cactus founders) ===")
shared82 = xwu82 & cactus
print(f"  shared positions: {len(shared82):,}")
print(f"  % of xwu82 in cactus: {100*len(shared82)/len(xwu82):.2f}%")
print(f"  % of cactus in xwu82: {100*len(shared82)/len(cactus):.2f}%")
print(f"  cactus-only (NOT in xwu82): {len(cactus - xwu82):,}")
print(f"  xwu82-only (NOT in cactus): {len(xwu82 - cactus):,}")

print()
print(f"=== Implications for ref_80 (imputation reference panel) ===")
print(f"  Current ref_80 has only xwu's SNPs: {len(xwu80):,} SNPs (drops {len(xwu231 - xwu80):,} that became monomorphic in 80-subset)")
print(f"  Adding cactus's SNPs (union) would give: {len(xwu80 | cactus):,} unique positions")
print(f"  → +{100*len(xwu80 | cactus)/len(xwu80) - 100:.1f}% more SNPs in the ref panel")
print(f"  Cactus-only sites add: {len(cactus - xwu80):,} new positions to the imputation backbone")
PYEOF

echo ""
echo "[$(date)] DONE"
