import numpy as np, json
from scipy.sparse import load_npz
ROOT="/carnegie/nobackup/scratch/tbellagio/hapfire_sv"
CAC=set(map(str,json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
cn=load_npz(f"{ROOT}/poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1.cn.npz").tocsr()
meta=np.load(f"{ROOT}/poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1.meta.npz",allow_pickle=True)
fo=np.asarray(meta["founders"]).astype(str); bid=np.asarray(meta["bubble_id"]).astype(np.int64)
ac=np.asarray(cn.sum(0)).flatten(); m_b=np.bincount(bid)[bid].astype(float)
is_cac=np.array([x in CAC for x in fo])
print("LINK 1: are LOW-AC (rare) k-mers in HIGH-m_b (multi-allelic) bubbles?")
for lo,hi,lab in [(2,2,"ac=2"),(3,4,"ac3-4"),(5,10,"ac5-10"),(11,50,"ac11-50"),(51,999,"ac>50")]:
    m=(ac>=lo)&(ac<=hi); print(f"  {lab:<8}: mean m_b={m_b[m].mean():.1f}  median={np.median(m_b[m]):.0f}  (n={m.sum():,})")
print("\nLINK 2: among ac=2 (discriminating) k-mers, does cactus carry MORE / higher-m_b ones?")
a2=ac==2
cac_car=np.asarray(cn[is_cac].sum(0)).flatten()>0
pg_car =np.asarray(cn[~is_cac].sum(0)).flatten()>0
print(f"  ac=2 carried by cactus: n={ (a2&cac_car).sum():,}  mean m_b={m_b[a2&cac_car].mean():.1f}")
print(f"  ac=2 carried by PG    : n={ (a2&pg_car).sum():,}  mean m_b={m_b[a2&pg_car].mean():.1f}")
