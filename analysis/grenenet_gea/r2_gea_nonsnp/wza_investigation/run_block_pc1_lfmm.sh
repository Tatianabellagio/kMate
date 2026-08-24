#!/bin/bash
#SBATCH --job-name=pc1_lfmm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/pc1_lfmm_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/pc1_lfmm_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_gea/wza_investigation                  # PC1 scripts now live here
P1R=analysis/grenenet_gea/phase1_replication                # run_lfmm_lastgen.R stays in the replication
PCDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/results/pc1
K=16

for CLS in snp smallindel sv; do
  echo "===== $CLS ====="
  echo "== build [pools x blocks] PC1 matrix =="
  $PY -u $PR/build_block_pc1_matrix.py --class $CLS --gen 9 --climate bio1
  echo "== LFMM K=$K on PC1 matrix (structure correction) =="
  $R $P1R/run_lfmm_lastgen.R "$PCDIR/pc1lfmm_${CLS}_gen9" $K "$PCDIR/pc1lfmm_${CLS}_gen9_calibp.csv"
  echo "== join -> per-block result =="
  $PY - "$CLS" <<'EOF'
import sys, pandas as pd, numpy as np
cls=sys.argv[1]; PCDIR="analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/results/pc1"
b=pd.read_csv(f"{PCDIR}/pc1lfmm_{cls}_gen9_blocks.csv")
pv=pd.read_csv(f"{PCDIR}/pc1lfmm_{cls}_gen9_calibp.csv")["pval"].to_numpy()
assert len(pv)==len(b),(len(pv),len(b))
b["pval"]=pv
p=b["pval"].to_numpy(); n=len(p); o=np.argsort(p); q=np.empty(n)
q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; b["q"]=np.clip(q,0,1)
b.sort_values("pval").to_csv(f"{PCDIR}/pc1lfmm_{cls}_gen9_bio1.csv",index=False)
print("  wrote pc1lfmm result | BH q<0.05:",int((b.q<0.05).sum()),"| min p %.2e"%b.pval.min())
EOF
done

echo "== COMPARE: PC1-LFMM (structure-corrected) vs PC1-Kendall vs WZA =="
$PY - <<'EOF'
import pandas as pd, numpy as np
from scipy.stats import chi2
PCDIR="analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/results/pc1"; WZ="analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/wza"
def lam(p): p=np.clip(np.asarray(p,float),1e-300,1); return np.median(chi2.isf(p,1))/chi2.isf(0.5,1)
for cls in ["snp","smallindel","sv"]:
    pl=pd.read_csv(f"{PCDIR}/pc1lfmm_{cls}_gen9_bio1.csv"); pl["block"]=pl["block"].astype(str)
    pk=pd.read_csv(f"{PCDIR}/pc1_{cls}_gen9_bio1.csv"); pk["block"]=pk["block"].astype(str)
    pl=pl.sort_values("pval").reset_index(drop=True)
    print(f"\n--- {cls}: n={len(pl)} blocks ---")
    print(f"    lambda: PC1-Kendall={lam(pk.pval):.2f} -> PC1-LFMM(K=16)={lam(pl.pval):.2f}  | PC1-LFMM BH q<0.05: {(pl.q<0.05).sum()}")
    for blk,name in [("2_1265","CAM5"),("4_2519","Chr4-12Mb"),("4_2781","Chr4-15.5Mb")]:
        r=pl.index[pl.block==blk]
        if len(r): print(f"    {name:11s} {blk}: PC1-LFMM rank {int(r[0])+1}/{len(pl)} p={pl.loc[r[0],'pval']:.2e} q={pl.loc[r[0],'q']:.3f}")
        else: print(f"    {name:11s} {blk}: absent")
    print("    top 6 PC1-LFMM blocks:")
    print(pl.head(6)[["block","chrom","pos","n_snps","pc1_var_explained","pval","q"]].to_string(index=False))
EOF
echo "== ALLDONE_PC1LFMM =="
