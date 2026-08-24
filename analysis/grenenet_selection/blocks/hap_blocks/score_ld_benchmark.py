"""Score the block-size benchmark: per-config per-record AF vs recomb_truth.
Answers 'error vs block size' + 'coverage lost when blocks too small'.
Configs: global baseline + r2={0.1,0.2,0.3,0.4} local-only per-block EM.
Metrics per config: coverage (finite-AF frac), AF RMSE/MAE/Pearson over finite
records, block count + median block width, block resolve/NaN status.
"""
import numpy as np, pandas as pd, json, os
HB="/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/hap_blocks"
OUT=f"{HB}/bench_g0_231"
SIM="/global/scratch/users/tbellg/kmate/benchmarks/p231/sims/cov10_n231_g0_s42_hotspots_p231_chr1"

KEY=["chrom","pos","ref_len","alt_len"]
truth=pd.read_csv(f"{SIM}/recomb_truth_raw.tsv.gz",sep="\t")
truth["chrom"]=truth["chrom"].astype(str).str.replace("Chr","",regex=False)
truth=truth[["chrom","pos","ref_len","alt_len","truth_af"]]

part=json.load(open(f"{HB}/ld_partitions_summary.json"))
psize={f"{p['r2']:.2f}":p for p in part}

def score(tsv):
    d=pd.read_csv(tsv,sep="\t")
    d["chrom"]=d["chrom"].astype(str).str.replace("Chr","",regex=False)
    m=d.merge(truth,on=KEY,how="inner")
    fin=np.isfinite(m["alt_freq"].values)
    cov=fin.mean()
    a=m["alt_freq"].values[fin]; t=m["truth_af"].values[fin]
    rmse=float(np.sqrt(np.mean((a-t)**2))); mae=float(np.mean(np.abs(a-t)))
    r=float(np.corrcoef(a,t)[0,1]) if fin.sum()>2 else np.nan
    return dict(n_records=len(m),coverage=float(cov),rmse=rmse,mae=mae,pearson=r)

def block_status(cfg):
    h=f"{OUT}/{cfg}.tsv".replace(".tsv",".h_blocks_per_chrom.npz")
    if not os.path.exists(h): return None
    z=np.load(h,allow_pickle=True)
    st=z["Chr1_status"] if "Chr1_status" in z.files else None
    if st is None: return None
    return dict(n_blocks=int(len(st)),resolved=int((st==0).sum()),
                nan_lowk=int((st==1).sum()),empty=int((st==2).sum()))

rows=[]
# global baseline
if os.path.exists(f"{OUT}/global.tsv"):
    s=score(f"{OUT}/global.tsv"); s.update(config="global",r2=np.nan,block_kb=np.nan,n_blocks=1)
    rows.append(s)
for r2 in ["0.10","0.20","0.30","0.40"]:
    tsv=f"{OUT}/alltype_r2_{r2}.tsv"
    if not os.path.exists(tsv): continue
    s=score(tsv); s.update(config=f"r2_{r2}",r2=float(r2),
        block_kb=psize.get(r2,{}).get("width_med_kb",np.nan),
        n_blocks=psize.get(r2,{}).get("n_blocks",np.nan))
    bs=block_status(f"r2_{r2}")
    if bs: s.update(resolved=bs["resolved"],nan_lowk=bs["nan_lowk"],empty=bs["empty"])
    rows.append(s)

df=pd.DataFrame(rows)
cols=["config","r2","block_kb","n_blocks","coverage","rmse","mae","pearson"]
extra=[c for c in ["resolved","nan_lowk","empty","n_records"] if c in df.columns]
print(df[cols+extra].to_string(index=False,float_format=lambda x:f"{x:.4f}"))
df.to_csv(f"{OUT}/scores.csv",index=False)
print(f"\nwrote {OUT}/scores.csv")
print("\nRead: global = all-231 baseline. As blocks shrink (r2 up), watch RMSE and")
print("coverage: the block size where RMSE stops improving / coverage drops = the floor.")
