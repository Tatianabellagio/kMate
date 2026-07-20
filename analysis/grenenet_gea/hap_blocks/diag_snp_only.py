"""Diagnostic: does the var_pa r2=0.1 fragmentation (25 blocks vs hapFIRE's 2)
come from SV/pangenome markers? Re-run the SAME CompleteLDPartition at r2=0.1 on
SNP-ONLY common markers (ref_len==1 & alt_len==1) and compare block count.
"""
import numpy as np, scipy.sparse as sp, time
from sklearn import preprocessing
import importlib.util
# reuse vendored CompleteLDPartition from gen_ld_partitions
spec=importlib.util.spec_from_file_location("gen","/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/hap_blocks/gen_ld_partitions.py")
# can't exec (it runs the pipeline). Re-vendor the two fns instead:
def find_ld(i,snps,cutoff,window_size):
	n_inds,n_snps=snps.shape
	left=max(i-window_size,0); right=min(i+window_size,n_snps)
	lcor=np.matmul(np.transpose(snps[:,left:(i+1)]),snps[:,i])/n_inds
	rcor=np.matmul(np.transpose(snps[:,i:(right+1)]),snps[:,i])/n_inds
	lrev=np.flip(lcor)
	ll=np.where(lrev**2>cutoff)[0]
	for j in range(len(ll)-1):
		if ll[j+1]-ll[j]>25: ll=ll[:j+1]; break
	ll=np.flip(ll)*-1
	rl=np.where(rcor**2>cutoff)[0]
	for j in range(len(rl)-1):
		if rl[j+1]-rl[j]>25: rl=rl[:j+1]; break
	return np.unique(np.concatenate((ll,rl)))+i
def CompleteLDPartition(G,cutoff,window_size):
	n,m=G.shape; sl={}; cum=[]; mx=[]; bd=[]; alone=[]
	for i in range(m):
		sl[i]=find_ld(i,G,0.1,50)
		if len(sl[i])==1: alone.append(i)
	for i in range(m): sl[i]=find_ld(i,G,cutoff,window_size)
	for i in range(m): mx.append(np.max(sl[i]) if len(sl[i])>0 else i)
	cum.append(mx[0])
	for i in range(1,len(mx)): cum.append(mx[i] if mx[i]>cum[i-1] else cum[i-1])
	idx=np.where(cum-np.array(range(m))==0)[0]
	b_=np.concatenate(([-1],idx))
	for i in range(len(b_)-1):
		l=b_[i]+1; r=b_[i+1]
		bd.append([l,r] if r-l>0 else [l])
	j=0
	while j<len(bd):
		if len(bd[j])==1:
			if j==0: bd[j+1][0]=bd[j][0]; del bd[j]
			else: bd[j-1][1]=bd[j][0]; del bd[j]
		else: j+=1
	return bd,alone

t0=time.time()
VP="/global/scratch/users/tbellg/kmate/panel/arch3/chr1/var_pa_231_arch3_chr1"
vp=sp.load_npz(VP+".var_pa.npz").tocsr(); F=vp.shape[0]
m=np.load(VP+".meta.npz",allow_pickle=True)
pos=m["pos"].astype(np.int64); rl=m["ref_len"].astype(np.int64); al=m["alt_len"].astype(np.int64)
freq=np.asarray(vp.sum(0)).ravel()/F
is_snp=(rl==1)&(al==1)
common=(freq>0.05)&(freq<0.95)
sel=np.flatnonzero(common & is_snp)
sel=sel[np.argsort(pos[sel],kind="stable")]
print(f"[{time.time()-t0:.0f}s] SNP-only common markers: {len(sel):,} "
      f"(vs 645,882 all-type). pos {pos[sel].min()}-{pos[sel].max()}",flush=True)
G=preprocessing.scale(np.asarray(vp[:,sel].todense(),dtype=np.float64))
bd,alone=CompleteLDPartition(G,0.1,100)
cpos=pos[sel]
ends=[round(cpos[b[-1] if len(b)>1 else b[0]]/1e6,2) for b in sorted([(x[0],x[-1] if len(x)>1 else x[0]) for x in bd])]
print(f"[{time.time()-t0:.0f}s] SNP-only r2=0.1: {len(bd)} blocks; ends(Mb)={ends if len(ends)<=15 else ends[:8]+['...']+ends[-3:]}")
print("hapFIRE (SNP panel) r2=0.1: 2 blocks [24.02, 30.42]")
print("all-type var_pa r2=0.1: 25 blocks (from gen_ld_partitions)")
