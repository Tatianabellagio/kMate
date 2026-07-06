#!/bin/bash
#SBATCH --job-name=lfmm_snp_v_nonsnp_k16
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/snp_v_nonsnp_k16_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/snp_v_nonsnp_k16_%j.out
# STAGE B: raw LFMM (K=16, GIF-calibrated) for the pooled non-SNP scan, then a
# first "where do the hits land" comparison vs the existing SNP K=16 scan.
# SNP per-record calibrated p already exists (lfmm_snp_gen9_bio1_calibp.csv); we
# only need to run nonsnp, then summarize both at record + LD-block level.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_gea/phase1_replication
LDIR=results/grenenet_gea/phase1_replication/lfmm
CMDIR=results/grenenet_gea/phase1_replication/class_matrices
OUT=results/grenenet_gea/gea_newpanel/snp_vs_nonsnp
mkdir -p "$OUT"
K=16

echo "== nonsnp LFMM ridge K=$K + gif =="
$R $PR/run_lfmm_lastgen.R "$LDIR/lfmm_nonsnp_gen9" $K "$LDIR/lfmm_nonsnp_gen9_bio1_calibp.csv"

echo "== join calibrated p -> per-record tables (snp reuse; nonsnp new) =="
$PY - <<'EOF'
import pandas as pd, numpy as np
LDIR="results/grenenet_gea/phase1_replication/lfmm"
CMDIR="results/grenenet_gea/phase1_replication/class_matrices"
for cls in ["snp","nonsnp"]:
    rec=pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
    pv=pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_bio1_calibp.csv")["pval"].to_numpy()
    assert len(pv)==len(rec),(cls,len(pv),len(rec))
    rec=rec.assign(pval=pv).rename(columns={"maf":"MAF"})
    rec[["chrom","pos","ref_len","alt_len","MAF","block","pval"]].to_csv(
        f"{LDIR}/lfmm_{cls}_gen9_bio1.csv",index=False)
    print(f"  {cls}: {len(rec):,} records, min p={np.nanmin(pv):.2e}")
EOF

echo "== WHERE DO THE HITS LAND: record + LD-block comparison =="
$PY - <<'EOF'
import pandas as pd, numpy as np, scipy.stats as st
LDIR="results/grenenet_gea/phase1_replication/lfmm"
OUT="results/grenenet_gea/gea_newpanel/snp_vs_nonsnp"
def gif(p):
    p=np.clip(np.asarray(p,float),1e-300,1); return np.median(st.chi2.isf(p,1))/st.chi2.isf(0.5,1)

D={}
for cls in ["snp","nonsnp"]:
    d=pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_bio1.csv")
    d=d[np.isfinite(d.pval)].copy()
    d["sv_size"]=(d.alt_len-d.ref_len).abs()
    D[cls]=d
    print(f"\n### {cls}: {len(d):,} records | GIF(calib)={gif(d.pval):.2f} | "
          f"min p={d.pval.min():.2e} | p<1e-5: {(d.pval<1e-5).sum():,}")
    top=d.nsmallest(20,"pval")[["chrom","pos","sv_size","MAF","block","pval"]]
    print(top.to_string(index=False))
    top.to_csv(f"{OUT}/top20_records_{cls}.csv",index=False)

# LD-block aggregation: min p (== -log10 peak) and record count per block, per scan
def blockmin(d):
    g=d.groupby("block").agg(min_p=("pval","min"), n=("pval","size"),
                             chrom=("chrom","first"), pos=("pos","median")).reset_index()
    g["nlp"]=-np.log10(g.min_p.clip(1e-300))
    return g
bs=blockmin(D["snp"]); bn=blockmin(D["nonsnp"])
m=bs.merge(bn,on="block",suffixes=("_snp","_nonsnp"))
print(f"\n### LD blocks: snp={len(bs):,} nonsnp={len(bn):,} shared={len(m):,}")
r=st.spearmanr(m.nlp_snp,m.nlp_nonsnp)
print(f"block peak -log10p Spearman (shared blocks): rho={r.correlation:.3f} p={r.pvalue:.1e}")
# top-100 block overlap by peak
def topset(g,n=100): return set(g.nsmallest(n,"min_p").block)
for n in (25,50,100,200):
    a,b=topset(bs,n),topset(bn,n); j=len(a&b)/len(a|b)
    print(f"  top-{n} block overlap: shared={len(a&b)}  Jaccard={j:.3f}")
m.to_csv(f"{OUT}/block_peak_snp_vs_nonsnp.csv",index=False)
# blocks that are top in nonsnp but NOT top in snp (candidate SV/indel-only signal)
mtop=m.assign(rank_snp=m.min_p_snp.rank(),rank_nonsnp=m.min_p_nonsnp.rank())
only_nonsnp=mtop[(mtop.rank_nonsnp<=50)&(mtop.rank_snp>500)].sort_values("min_p_nonsnp")
print(f"\n### top-50 nonsnp blocks that are weak (>rank500) in snp: {len(only_nonsnp)}")
print(only_nonsnp[["block","chrom_snp","pos_snp","min_p_snp","min_p_nonsnp","n_nonsnp"]].head(20).to_string(index=False))
only_nonsnp.to_csv(f"{OUT}/nonsnp_only_top_blocks.csv",index=False)
print("\nwrote:",OUT)
EOF
echo DONE_SNP_VS_NONSNP
