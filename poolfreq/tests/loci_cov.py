import numpy as np, json
from pathlib import Path
from scipy.sparse import load_npz
ROOT=str(Path(__file__).resolve().parents[2])
CAC=set(map(str,json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
cn=load_npz(f"{ROOT}/poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1.cn.npz").tocsr()
meta=np.load(f"{ROOT}/poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1.meta.npz",allow_pickle=True)
fo=np.asarray(meta["founders"]).astype(str); bid=np.asarray(meta["bubble_id"]).astype(np.int64)
nbub=int(bid.max())+1
import numpy as np
ncov=np.zeros(len(fo)); nkm=np.zeros(len(fo))
for f in range(len(fo)):
    s,e=cn.indptr[f],cn.indptr[f+1]; cols=cn.indices[s:e]
    nkm[f]=len(cols); ncov[f]=len(np.unique(bid[cols]))
cac=np.array([x in CAC for x in fo])
print(f"total bubbles in panel: {nbub:,}")
print(f"{'class':<8}{'n':>4}{'med loci covered':>18}{'med kmers':>12}{'med kmers/locus':>16}")
for lab,m in [("cactus",cac),("PG",~cac)]:
    lc=ncov[m]; km=nkm[m]
    print(f"{lab:<8}{m.sum():>4}{np.median(lc):>18.0f}{np.median(km):>12.0f}{np.median(km/np.maximum(lc,1)):>16.2f}")
print(f"\ncactus/PG median loci-covered ratio = {np.median(ncov[cac])/np.median(ncov[~cac]):.2f}")
print(f"cactus/PG median kmers/locus ratio  = {np.median((nkm/np.maximum(ncov,1))[cac])/np.median((nkm/np.maximum(ncov,1))[~cac]):.2f}")
