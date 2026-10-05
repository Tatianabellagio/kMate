import sys, os
import pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); CONV=os.path.dirname(HERE); sys.path.insert(0,HERE); sys.path.insert(0,CONV)
import units as UN, tfbs_turnover as TF
R=f"{CONV}/results"
N=pd.read_csv(f"{R}/candidate_list_v2_events_cr50/v2_layer_tfbs.csv"); N=N[N.skipped_reason.isna()]
O=pd.read_csv(f"{R}/tfbs_pool_summary.csv"); E=pd.read_csv(f"{R}/evidence_matrix.csv",usecols=["store_row","chrom","pos","ref_len","alt_len"])
O=O.merge(E,on="store_row")
A=UN.rep_alleles(N.unit.values)[["unit","chrom","pos","ref_len","alt_len","ref","alt"]]
M=N.merge(A,on="unit").merge(O,on=["chrom","pos","ref_len","alt_len"],suffixes=("","_old"))
print("key matches (unit x old row):",len(M))
cmp=["tfbs_n_ref","tfbs_lost","tfbs_gained","tfbs_repeat_hits","tfbs_lost_motifs","tfbs_gained_motifs","tfbs_n_families","tfbs_tandem","tfbs_sites_per_motif"]
same=0;diff=[];noseq=0
for r in M.itertuples():
    ro,ao=TF.fetch_allele(r.chrom,int(r.pos),int(r.ref_len),int(r.alt_len))
    if (ro,ao)!=(r.ref.upper(),r.alt.upper()): noseq+=1; continue
    ok=True
    for c in cmp:
        a,b=getattr(r,c),getattr(r,c+"_old")
        if not ((pd.isna(a) and pd.isna(b)) or a==b): ok=False
    if (str(r.tfbs_lost_families) if pd.notna(r.tfbs_lost_families) else "")!=(str(r.tfbs_lost_families_old) if pd.notna(r.tfbs_lost_families_old) else ""): ok=False
    if ok: same+=1
    else: diff.append(r.unit)
print("sequence-identical pairs:",same+len(diff),"| key match but different seq:",noseq,"| exact agree:",same,"| disagree units:",diff)
