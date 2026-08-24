"""Self-contained r2 LD-block partitions from kMate's OWN 231-founder panel.
NO external SNP panel: markers = var_pa records (founder x record presence/absence,
maf>=0.05 common set). Uses hapFIRE's own CompleteLDPartition algorithm (r2 cutoff
+/-100 window) so the ONLY thing swept is the r2 threshold. Blocks are tiled to
genomic midpoints so every downstream k-mer lands in exactly one block.
Writes one blocks-TSV per cutoff (chrom start_pos end_pos n_variants) for
kmate --block-mode window --blocks-tsv, plus a stats summary.
"""
import sys, time, json
import numpy as np, scipy.sparse as sp
from sklearn import preprocessing

# ---------------------------------------------------------------------------
# find_ld + CompleteLDPartition copied VERBATIM from
#   external/HapFIRE/haplotype_generation.py  (lines 30-109)
# to avoid that module's top-level `import cvxpy` (unneeded here; these two fns
# are pure-numpy). Faithful to hapFIRE's own r2 independent-LD-block algorithm.
# ---------------------------------------------------------------------------
def find_ld(i,snps,cutoff,window_size):
	n_inds,n_snps = snps.shape
	left = max(i - window_size, 0)
	right = min(i + window_size,n_snps)
	left_snps_window = snps[:,left:(i+1)]
	right_snps_window = snps[:,i:(right+1)]
	left_cor = np.matmul(np.transpose(left_snps_window),snps[:,i])/n_inds
	left_cor_rev = np.flip(left_cor)
	right_cor = np.matmul(np.transpose(right_snps_window),snps[:,i])/n_inds
	left_list_ = np.where(left_cor_rev**2 > cutoff)[0]
	for j in range(len(left_list_)-1):
		if left_list_[j+1] - left_list_[j] > 25:
			left_list_ = left_list_[:j+1]
			break
	left_list_ = np.flip(left_list_) * -1
	right_list_ = np.where(right_cor**2 > cutoff)[0]
	for j in range(len(right_list_)-1):
		if right_list_[j+1] - right_list_[j] > 25:
			right_list_ = right_list_[:j+1]
			break
	SNPinLD_index = np.unique(np.concatenate((left_list_, right_list_))) + i
	return(SNPinLD_index)

def CompleteLDPartition(standardized_genotype_matrix,cutoff,window_size):
	n_inds,n_snps = standardized_genotype_matrix.shape
	snp_list = {}
	cummax_list = []
	max_list = []
	boundary = []
	alone_SNPs_index = []
	for i in range(n_snps):
		snp_list[i] = find_ld(i,snps=standardized_genotype_matrix,cutoff=0.1,window_size=50)
		if len(snp_list[i]) == 1:
			alone_SNPs_index.append(i)
	for i in range(n_snps):
		snp_list[i] = find_ld(i,snps=standardized_genotype_matrix,cutoff=cutoff,window_size=window_size)
	for i in range(len(snp_list)):
		if len(snp_list[i]) > 0:
			max_list.append(np.max(snp_list[i]))
		else:
			max_list.append(i)
	cummax_list.append(max_list[0])
	for i in range(1,len(max_list)):
		if max_list[i] > cummax_list[i-1]:
			cummax_list.append(max_list[i])
		else:
			cummax_list.append(cummax_list[i-1])
	idx = np.where( cummax_list - np.array(range(n_snps)) == 0)
	idx = idx[0]
	boundary_ = np.concatenate(([-1],np.array(idx)))
	for i in range(len(boundary_)-1):
		left = boundary_[i]+1
		right = boundary_[i+1]
		if right - left > 0:
			boundary.append([left,right])
		else:
			boundary.append([left])
	j = 0
	while j < len(boundary):
		if len(boundary[j]) == 1:
			if j ==0:
				boundary[j+1][0] = boundary[j][0]
				del boundary[j]
			else:
				boundary[j-1][1] = boundary[j][0]
				del boundary[j]
		else:
			j += 1
	return(boundary,alone_SNPs_index)

ROOT="/global/scratch/users/tbellg/kmate"
CHR="Chr1"; CHRNUM="Chr1"
VP=f"{ROOT}/panel/arch3/chr1/var_pa_231_arch3_chr1"
OUT=f"{ROOT}/analysis/grenenet_gea/hap_blocks"
CUTOFFS=[0.1,0.2,0.3,0.4]; WINDOW=100; MAF=0.05; CALLRATE_MIN=0.9
CHR1_LEN=30427671   # TAIR10 Chr1, so blocks tile to the true chromosome end

t0=time.time()
vp=sp.load_npz(VP+".var_pa.npz").tocsr(); F=vp.shape[0]
vc=sp.load_npz(VP+".var_called.npz").tocsr()                    # 1 = genotype called
meta=np.load(VP+".meta.npz",allow_pickle=True); pos=meta["pos"].astype(np.int64)
# --- MISSINGNESS-AWARE marker filtering ---
# var_pa=0 conflates "carries reference" with "uncalled". Use var_called to
# (a) require a call-rate floor, (b) compute maf on CALLED genotypes only.
called=np.asarray(vc.sum(axis=0)).ravel().astype(np.float64)
alt=np.asarray(vp.sum(axis=0)).ravel().astype(np.float64)
callrate=called/F
with np.errstate(invalid="ignore", divide="ignore"):
    freqc=np.where(called>0, alt/called, 0.0)                   # called-based allele freq
common=np.flatnonzero((callrate>=CALLRATE_MIN)&(freqc>MAF)&(freqc<1-MAF))
common=common[np.argsort(pos[common], kind="stable")]          # position-sorted (LD window == genomic)
cpos=pos[common]
print(f"[{time.time()-t0:.0f}s] common (callrate>={CALLRATE_MIN}, called-maf>={MAF}): "
      f"{len(common):,} markers  pos {cpos.min()}-{cpos.max()}",flush=True)

# genotype matrix with MISSING mean-imputed per marker (to called allele freq),
# so uncalled entries center to ~0 and contribute nothing to r (standard LD handling).
G=np.asarray(vp[:,common].todense(), dtype=np.float64)          # 231 x n_common (uncalled currently 0)
Cs=np.asarray(vc[:,common].todense(), dtype=bool)               # called mask (1 byte/elt, cheap)
fc=freqc[common]
miss=~Cs; rows,cols=np.where(miss)                              # in-place impute (avoid full copy)
G[rows,cols]=fc[cols]                                           # uncalled -> marker called-freq
del Cs,miss,rows,cols
Gs=preprocessing.scale(G); del G                                # column-standardized (imputed -> ~0)
chrom_end=CHR1_LEN
print(f"[{time.time()-t0:.0f}s] standardized {Gs.shape} (missing mean-imputed)",flush=True)

summary=[]
for r2 in CUTOFFS:
    tb=time.time()
    boundary,alone=CompleteLDPartition(standardized_genotype_matrix=Gs,cutoff=r2,window_size=WINDOW)
    # boundary: list of [left_idx,right_idx] in common-marker index space
    edges=[]
    for blk in boundary:
        l=blk[0]; r=blk[-1] if len(blk)>1 else blk[0]
        edges.append((int(l),int(r)))
    edges.sort()
    # tile to genomic midpoints so blocks are contiguous (no unassigned k-mers)
    starts=[1]
    for k in range(len(edges)-1):
        r_pos=cpos[edges[k][1]]; nxt_pos=cpos[edges[k+1][0]]
        mid=(int(r_pos)+int(nxt_pos))//2
        starts.append(mid+1)
    rows=[]
    for k,(l,r) in enumerate(edges):
        s=starts[k]; e=(starts[k+1]-1) if k+1<len(starts) else chrom_end
        rows.append((CHRNUM,s,e,r-l+1))
    path=f"{OUT}/ld_blocks_r2_{r2:.2f}.tsv"
    with open(path,"w") as fh:
        fh.write("chrom\tstart_pos\tend_pos\tn_variants\n")
        for c,s,e,n in rows: fh.write(f"{c}\t{s}\t{e}\t{n}\n")
    widths=np.array([e-s+1 for _,s,e,_ in rows],float)
    nv=np.array([n for *_,n in rows])
    rec=dict(r2=r2,n_blocks=len(rows),alone=len(alone),
             width_med_kb=float(np.median(widths)/1e3),width_p10_kb=float(np.percentile(widths,10)/1e3),
             width_p90_kb=float(np.percentile(widths,90)/1e3),width_max_mb=float(widths.max()/1e6),
             nvar_med=float(np.median(nv)),nvar_p10=float(np.percentile(nv,10)))
    summary.append(rec)
    print(f"[{time.time()-t0:.0f}s] r2={r2}: {len(rows)} blocks  width(kb) med={rec['width_med_kb']:.0f} "
          f"p10={rec['width_p10_kb']:.1f} p90={rec['width_p90_kb']:.0f} max={rec['width_max_mb']:.1f}Mb  "
          f"nvar/block med={rec['nvar_med']:.0f}  [+{time.time()-tb:.0f}s]",flush=True)

json.dump(summary,open(f"{OUT}/ld_partitions_summary.json","w"),indent=1)
print(f"\nwrote {len(CUTOFFS)} partitions + summary to {OUT}  ({time.time()-t0:.0f}s)")
