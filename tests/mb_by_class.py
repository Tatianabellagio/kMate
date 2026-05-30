import numpy as np, json
from pathlib import Path
from scipy.sparse import load_npz
ROOT=str(Path(__file__).resolve().parents[2])
CAC=set(map(str,json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
kmer_pa=load_npz(f"{ROOT}/data/kmer_pa_231_v3qc_v3_filt2/kmer_pa_Chr1.kmer_pa.npz").tocsr()
meta=np.load(f"{ROOT}/data/kmer_pa_231_v3qc_v3_filt2/kmer_pa_Chr1.meta.npz",allow_pickle=True)
fo=np.asarray(meta["founders"]).astype(str); bid=np.asarray(meta["bubble_id"]).astype(np.int64)
ac=np.asarray(kmer_pa.sum(0)).flatten()             # AC per k-mer
m_b=np.bincount(bid)[bid].astype(float)        # bubble size per k-mer
cac=np.array([x in CAC for x in fo])
# mean m_b and mean AC of the k-mers each class CARRIES (averaged over carried entries)
for lab,m in [("cactus",cac),("PG",~cac)]:
    rows=np.where(m)[0]; cols=np.concatenate([kmer_pa.indices[kmer_pa.indptr[f]:kmer_pa.indptr[f+1]] for f in rows])
    print(f"{lab}: carried-kmer mean m_b={m_b[cols].mean():.2f} median m_b={np.median(m_b[cols]):.0f} | mean AC={ac[cols].mean():.2f} frac ac=2:{(ac[cols]==2).mean():.3f}")
