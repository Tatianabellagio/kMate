"""Correctness check: our var_pa-based r2 partition vs hapFIRE's SNP-based blocks.
Same CompleteLDPartition algorithm, different (self-contained) marker set. At
r2=0.1 the arm-scale structure should be comparable to hapFIRE's independent LD
blocks (chr1 = 24.0Mb + 6.4Mb, breakpoint ~24Mb). We compare block COUNT and the
main breakpoint position(s), not exact identity (different panels).
"""
import numpy as np
ROOT="/global/scratch/users/tbellg/kmate"
OUT=f"{ROOT}/analysis/grenenet_selection/blocks/hap_blocks"
HAPF="/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_independent_genomewide_partition.txt"

def load_ours(r2, chrom="1"):
    rows=[]
    for i,ln in enumerate(open(f"{OUT}/ld_blocks_r2_{r2:.2f}.tsv")):
        if i==0: continue
        c,s,e,n=ln.split()
        if c==chrom: rows.append((int(s),int(e),int(n)))
    return rows

def load_hapfire(chrom="1"):
    rows=[]
    for ln in open(HAPF):
        p=ln.split()
        if p[0]==chrom: rows.append((int(p[1]),int(p[2])))
    return rows

hf=load_hapfire()
print("=== hapFIRE s1 independent LD blocks (SNP panel, r2=0.1), chr1 ===")
print(f"  {len(hf)} blocks; breakpoints(Mb): {[round(e/1e6,2) for s,e in hf]}")
print()
for r2 in [0.1,0.2,0.3,0.4]:
    try: ours=load_ours(r2)
    except FileNotFoundError:
        print(f"  r2={r2}: (not generated yet)"); continue
    bounds=[round(e/1e6,2) for s,e,n in ours]
    tag=""
    if r2==0.1:
        # main breakpoint proximity to hapFIRE's ~24Mb arm boundary
        hf_bp=hf[0][1]/1e6
        near=min([abs(e/1e6-hf_bp) for s,e,n in ours]) if ours else 999
        tag=f"   [our nearest breakpoint to hapFIRE {hf_bp:.1f}Mb: {near:.2f}Mb away]"
    show = bounds if len(bounds)<=12 else bounds[:6]+["...",]+bounds[-3:]
    print(f"  our var_pa r2={r2}: {len(ours):3d} blocks  ends(Mb)={show}{tag}")
print()
print("Interpretation: at r2=0.1 we expect FEW arm-scale blocks with a breakpoint")
print("near ~24Mb; more/smaller blocks are fine at higher r2. Hundreds of blocks at")
print("r2=0.1 would indicate the var_pa LD is distorted (check SVs/standardization).")
