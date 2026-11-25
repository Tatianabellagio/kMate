import numpy as np, pandas as pd
G="results/sv_adaptive"
D=pd.read_csv(f"{G}/founder_burden_event.csv",dtype={"founder":str})
TERMS=["snp_CDS","snp_TE","snp_promoter","snp_intergenic",
       "ins_CDS","ins_TE","ins_promoter","ins_intergenic",
       "del_CDS","del_TE","del_promoter","del_intergenic"]
rows=[]
for t in TERMS:
    k,x=t.split("_")
    num=D[f"{k}_shared_{x}_all"].to_numpy(float); den=D[f"called_{k}_shared_{x}_all"].to_numpy(float)
    r=num/den
    rows.append(dict(term=t,sd_rate=r.std(),called_median=np.median(den),
                     variants_per_sd=r.std()*np.median(den),carried_median=np.median(num)))
V=pd.DataFrame(rows)
old=pd.read_csv(f"{G}/sv_burden_full_model_event_CDS+promoter.csv").set_index("term")
V["variants_per_sd_published"]=[old.variants_per_sd.get(t,np.nan) for t in V.term]
V["ratio_check"]=V.variants_per_sd/V.variants_per_sd_published
print("recomputed variants_per_sd = sd(rate) x median(called), checked against the published table:")
print(V.round(4).to_string(index=False))
V[["term","sd_rate","called_median","carried_median","variants_per_sd"]].to_csv(f"{G}/variants_per_sd_event_SV6.csv",index=False)
print(f"\n[wrote] {G}/variants_per_sd_event_SV6.csv")
