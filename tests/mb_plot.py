import numpy as np, json
from pathlib import Path
from scipy.sparse import load_npz
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=str(Path(__file__).resolve().parents[2])
CAC=set(map(str,json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
cn=load_npz(f"{ROOT}/data/cn_full_231_v3qc_v3_filt2/cn_Chr1.cn.npz").tocsr()
meta=np.load(f"{ROOT}/data/cn_full_231_v3qc_v3_filt2/cn_Chr1.meta.npz",allow_pickle=True)
fo=np.asarray(meta["founders"]).astype(str); bid=np.asarray(meta["bubble_id"]).astype(np.int64)
is_cac=np.array([x in CAC for x in fo])
m_b=np.bincount(bid)[bid].astype(np.int64)                 # per-kmer bubble size
# per-bubble size (each bubble once)
mb_per_bubble=np.bincount(bid)
# which class carries each k-mer (any founder of that class)
cac_car=np.asarray(cn[is_cac].sum(0)).flatten()>0
pg_car =np.asarray(cn[~is_cac].sum(0)).flatten()>0
fig,ax=plt.subplots(1,2,figsize=(14,5))
b=np.linspace(0,300,60)
ax[0].hist(np.clip(mb_per_bubble,0,300),bins=b,color="#555")
ax[0].set_title(f"k-mers per bubble (all {len(mb_per_bubble):,} bubbles)\nmedian={np.median(mb_per_bubble):.0f} mean={mb_per_bubble.mean():.1f}")
ax[0].set_xlabel("m_b (k-mers in bubble, clipped 300)"); ax[0].set_ylabel("# bubbles")
ax[1].hist(np.clip(m_b[cac_car],0,300),bins=b,histtype="step",lw=2,density=True,color="#1f77b4",
           label=f"cactus-carried (med m_b={np.median(m_b[cac_car]):.0f})")
ax[1].hist(np.clip(m_b[pg_car],0,300),bins=b,histtype="step",lw=2,density=True,color="#ff7f0e",
           label=f"PG-carried (med m_b={np.median(m_b[pg_car]):.0f})")
ax[1].set_title("bubble size m_b of k-mers carried by each class\n(density-normalized)")
ax[1].set_xlabel("m_b (clipped 300)"); ax[1].set_ylabel("density"); ax[1].legend()
fig.tight_layout(); out=f"{ROOT}/notebooks/plots/kmers_per_bubble_class.png"
fig.savefig(out,dpi=110); print("wrote",out)
print(f"per-bubble m_b: median={np.median(mb_per_bubble):.0f} mean={mb_per_bubble.mean():.1f} max={mb_per_bubble.max()}")
