#!/bin/bash
#SBATCH --job-name=block_pc1
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/wza/investigation/block_pc1_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/wza/investigation/block_pc1_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_selection/wza/investigation                 # PC1 scripts now live here

for CLS in snp smallindel sv; do
  echo "===== $CLS ====="
  $PY -u $PR/run_block_pc1.py --class $CLS --gen 9 --climate bio1
done

echo "== COMPARE: block-PC1 vs Kendall->WZA (deg7), CAM5 & 4_2519 =="
$PY - <<'EOF'
import pandas as pd, numpy as np
PC="analysis/grenenet_selection/wza/investigation/results/pc1"            # PC1 outputs (relocated)
WZ="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/wza"           # WZA inputs from the replication
def wza(cls):
    w=pd.read_csv(f"{WZ}/wza_kendall_{cls}_gen9_bio1_deg7nocap.csv").rename(columns={"index":"block"})
    w["block"]=w["block"].astype(str); w=w[w.Z_pVal.notna()].sort_values("Z_pVal").reset_index(drop=True)
    return w
for cls in ["snp","smallindel","sv"]:
    pc=pd.read_csv(f"{PC}/pc1_{cls}_gen9_bio1.csv"); pc["block"]=pc["block"].astype(str)
    pc=pc.sort_values("pval").reset_index(drop=True)
    w=wza(cls)
    print(f"\n--- {cls}: PC1 test = {len(pc)} blocks | WZA = {len(w)} blocks ---")
    print(f"    PC1 BH q<0.05: {(pc.q<0.05).sum()} | median PC1 var-explained {pc.pc1_var_explained.median():.2f}")
    for blk,name in [("2_1265","CAM5"),("4_2519","Chr4-12Mb"),("4_2781","Chr4-15.5Mb")]:
        rp=pc.index[pc.block==blk]; rw=w.index[w.block==blk]
        ps=(f"PC1 rank {int(rp[0])+1}/{len(pc)} p={pc.loc[rp[0],'pval']:.2e} q={pc.loc[rp[0],'q']:.3f} "
            f"VE={pc.loc[rp[0],'pc1_var_explained']:.2f} nSNP={int(pc.loc[rp[0],'n_snps'])}" if len(rp) else "PC1 absent")
        ws=(f"WZA rank {int(rw[0])+1}/{len(w)}" if len(rw) else "WZA absent")
        print(f"    {name:11s} {blk}: {ps} | {ws}")
    # rank concordance on shared blocks
    m=pc[["block","pval"]].merge(w[["block","Z_pVal"]],on="block")
    from scipy.stats import spearmanr
    print(f"    Spearman(PC1 -log10p, WZA -log10p) = {spearmanr(-np.log10(m.pval.clip(1e-300)), -np.log10(m.Z_pVal.clip(1e-300))).statistic:.3f}")
    # top PC1 blocks
    print("    top 5 PC1 blocks:")
    print(pc.head(5)[["block","chrom","pos","n_snps","pc1_var_explained","tau","pval","q"]].to_string(index=False))
EOF
echo "== ALLDONE_PC1 =="
