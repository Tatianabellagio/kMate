import numpy as np, pandas as pd
L=pd.read_csv("results/sv_adaptive/onestage/long.csv",dtype={"founder":str})
F=L.drop_duplicates("founder").set_index("founder")
P=["snp_CDS","ins_CDS","ins_TE","ins_promoter","ins_intergenic","del_CDS","del_TE","del_promoter","del_intergenic","cactus"]
X=F[P].astype(float)
print(f"founders: {len(X)}   predictors: {len(P)}\n")
print("correlation matrix of the fixed-effect burden predictors (231 founders):")
print(X.corr().round(2).to_string())
print("\nvariance inflation factors (VIF = 1/(1-R2) regressing each predictor on all others):")
Z=(X-X.mean())/X.std()
for j,c in enumerate(P):
    y=Z[c].to_numpy(); A=Z.drop(columns=c).to_numpy()
    A1=np.c_[np.ones(len(A)),A]
    b,*_=np.linalg.lstsq(A1,y,rcond=None); r2=1-((y-A1@b)**2).sum()/((y-y.mean())**2).sum()
    print(f"   {c:<16} R2 = {r2:.3f}   VIF = {1/max(1-r2,1e-12):7.2f}")
C=np.corrcoef(Z.T.to_numpy()); ev=np.linalg.eigvalsh(C)
print(f"\ncondition number of the predictor correlation matrix = {ev.max()/ev.min():.1f}"
      f"   (eigenvalues {np.round(ev,3).tolist()})")
print("\nsame, for the OLD SV4 model without the two TE terms:")
P4=[p for p in P if not p.endswith("_TE")]
Z4=Z[P4]
for c in P4:
    y=Z4[c].to_numpy(); A=Z4.drop(columns=c).to_numpy(); A1=np.c_[np.ones(len(A)),A]
    b,*_=np.linalg.lstsq(A1,y,rcond=None); r2=1-((y-A1@b)**2).sum()/((y-y.mean())**2).sum()
    print(f"   {c:<16} VIF = {1/max(1-r2,1e-12):7.2f}")
C4=np.corrcoef(Z4.T.to_numpy()); e4=np.linalg.eigvalsh(C4)
print(f"   condition number = {e4.max()/e4.min():.1f}")
