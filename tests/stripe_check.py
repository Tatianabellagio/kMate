import numpy as np
from pathlib import Path
from scipy.sparse import load_npz
ROOT=str(Path(__file__).resolve().parents[2]); CV=f"{ROOT}/panel/arch3/chr1/cn_var_231_arch3_chr1"
cn_var=load_npz(f"{CV}.cn_var.npz").tocsc(); cn_called=load_npz(f"{CV}.cn_var_called.npz").tocsc()
fo=np.asarray(np.load(f"{CV}.meta.npz",allow_pickle=True)["founders"]).astype(str); F=len(fo)
fidx={f:i for i,f in enumerate(fo)}
sim="g0_n50_rep101_pg"
present=np.zeros(F,bool)
for ln in open(f"{ROOT}/sims/visor_freqk/g0_reps/{sim}/h_truth.tsv"):
    p=ln.split()
    if p and p[0]!="founder" and float(p[1])>0 and p[0] in fidx: present[fidx[p[0]]]=True
print(f"present founders in pool: {present.sum()}")
# per-record: # PRESENT founders called, and truth_AF = present_carriers/present_called
called_present = np.asarray(cn_called[present].sum(axis=0)).flatten()   # of the 50 present, how many called
carr_present   = np.asarray(cn_var[present].sum(axis=0)).flatten()
ncalled_all    = np.asarray(cn_called.sum(axis=0)).flatten()            # of all 231
valid=called_present>0
truth=np.where(valid,carr_present/np.maximum(called_present,1),np.nan)
print(f"\ndistribution of #present-founders-called (of {present.sum()}):")
for lo,hi in [(1,5),(6,10),(11,20),(21,35),(36,50)]:
    m=valid&(called_present>=lo)&(called_present<=hi)
    print(f"  {lo:>2}-{hi:<2} present-called: {m.sum():>9,} records  ({100*m.mean():.1f}%)")
print("\nDo STRIPE records (truth_AF at coarse fractions) have FEW present-called?")
for v,lab in [(0.20,"~0.20"),(0.25,"~0.25"),(0.333,"~0.333"),(0.50,"~0.50")]:
    near=valid&(np.abs(truth-v)<0.004)
    if near.sum()>0:
        print(f"  truth_AF {lab}: n={near.sum():>8,}  median present-called={np.median(called_present[near]):.0f}  median %missing(all231)={100*(1-ncalled_all[near].mean()/F):.0f}%")
print("\nHIGH present-called (>=40) records — truth should be SMOOTH (no stripes):")
hi=valid&(called_present>=40); print(f"  n={hi.sum():,}  unique truth values in [0.18,0.22]: {len(np.unique(np.round(truth[hi&(truth>0.18)&(truth<0.22)],3)))} (many=smooth)")
lo=valid&(called_present<=5);  print(f"  LOW present-called (<=5): n={lo.sum():,}  unique truth in [0.18,0.22]: {len(np.unique(np.round(truth[lo&(truth>0.18)&(truth<0.22)],3)))} (few=stripes)")
