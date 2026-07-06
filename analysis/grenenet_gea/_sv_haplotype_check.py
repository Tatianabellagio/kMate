#!/usr/bin/env python
"""Reframe at the HAPLOBLOCK (hap-cluster) level -- the founder-GWAS's actual unit.

Resolves two things:
 (A) Why a "clq0.9 block" shows low block-WIDE variant r² (median ~0.38) even though 0.9 was
     the LD threshold: BigLD builds LD REGIONS that accrete MULTIPLE founder haplotypes
     (hap-clusters). 0.9 is the edge/graph threshold used to GROW a block, NOT a guarantee
     that all pairs in it are r²>=0.9. Large blocks (where SVs live -- has_sv~block-size
     corr 0.45) span 2-3 haplotypes, so block-wide r² is low BY DESIGN. Report n hap-clusters
     and within-block variant r² vs block size.
 (B) The right question for an SV: does it TAG one of the block's haplotypes (hap-cluster),
     and is THAT haplotype the selected one? Compute r²(SV, each hap-cluster) -> best-tagged
     cluster; compare its founder-fitness gap to the block's most-fit cluster.
Env: kmate.
"""
import os, sys, glob
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype

SEED = "results/grenenet_kmate_window_seedmix"; CHROMS = ["Chr1","Chr2","Chr3","Chr4","Chr5"]
PANEL = "panel/arch3"; FG = f"{lib.GEA}/hapfreq/multisite_founder_gwas_clq90_pc1"
SVL = f"{lib.GEA}/sv_adaptive/sv_landscape_clq0.9.csv"


def genome_h(s, base):
    gs = []
    for c in CHROMS:
        f = f"{base}/{s}_{c}.h_blocks_per_chrom.npz"
        if not os.path.exists(f): return None
        gs.append(np.load(f, allow_pickle=True)[f"{c}_global_h"].astype(np.float64))
    return np.mean(gs, 0)


def per_founder_fitness():
    cache = np.load(f"{lib.GEA}/fitness/sample_genome_h.npz", allow_pickle=True)
    H = cache["H"]; samps = cache["samples"].astype(str); hmap = {s:i for i,s in enumerate(samps)}
    seeds = sorted({p.split("/")[-1].split("_Chr")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    h0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(hmap)]
    Sg = []
    for site, sd in pt.groupby("site"):
        cell = {}
        for gen, g in sd.groupby("generation"):
            if int(gen) not in (1,2,3): continue
            hs, ws = [], []
            for _, r in g.iterrows():
                hs.append(H[hmap[str(r.sampleid)]])
                wv = r.flowerscollected if np.isfinite(r.flowerscollected) and r.flowerscollected>0 else 1.0
                ws.append(wv)
            ws=np.asarray(ws); cell[int(gen)]=(np.vstack(hs)*ws[:,None]).sum(0)/ws.sum()
        pres=sorted(cell)
        if 1 not in pres: continue
        t=np.array([0.0]+[float(g) for g in pres]); tc=t-t.mean()
        Y=np.vstack([h0]+[np.clip(cell[g],0,1) for g in pres])
        Sg.append((tc[:,None]*Y).sum(0)/(tc@tc))
    S=np.vstack(Sg); return S, S.mean(0)


def r2vec(a, B):
    a=a-a.mean(); Bc=B-B.mean(0)
    num=(a[:,None]*Bc).sum(0)**2; den=(a@a)*(Bc**2).sum(0)
    return np.divide(num,den,out=np.zeros_like(num),where=den>0)


def main():
    S, fit = per_founder_fitness()
    G, gf, reg = build_genotype()                    # MEMB_TAG=clq90 must be set
    reg["unit"] = reg.chrom+":"+reg.start.astype(str)+"-"+reg.end.astype(str)
    col_by_unit = {u:g.index.to_numpy() for u,g in reg.groupby("unit")}
    csv = pd.read_csv(f"{FG}.csv"); up = csv.groupby("unit",as_index=False).p_joint.min()
    L = pd.read_csv(SVL).rename(columns={"block_id":"unit"})
    up = up.merge(L[["unit","chrom","start_pos","end_pos","has_sv","n_kept"]],on="unit")
    up = up[up.n_kept>=2]; ntop=int(round(0.005*len(up)))
    topsv = up.nsmallest(ntop,"p_joint"); topsv=topsv[topsv.has_sv==1]

    rows=[]
    for ch in CHROMS:
        blk=topsv[topsv.chrom==ch]
        if not len(blk): continue
        cl=ch.lower()
        meta=np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz",allow_pickle=True)
        pos=meta["pos"].astype(np.int64); dl=np.abs(meta["alt_len"].astype(np.int64)-meta["ref_len"].astype(np.int64))
        vp=sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc=sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na=np.asarray(vp.sum(0)).ravel(); nc=np.asarray(vc.sum(0)).ravel()
        sv=(dl>50)&(na>=12)&(na<=219)&(nc/231>=0.9)
        for _,b in blk.iterrows():
            cols=col_by_unit.get(b.unit)
            if cols is None: continue
            Gc=G[:,cols].astype(float)                # 231 x n_clusters (hap membership)
            # fitness gap per cluster (per-site absmax) + which cluster is most fit
            capgap=[np.abs(S[:,Gc[:,k]==1].mean(1)-S[:,Gc[:,k]==0].mean(1)).max() if 1<=Gc[:,k].sum()<=229 else 0
                    for k in range(Gc.shape[1])]
            best_clu_gap=max(capgap) if capgap else np.nan
            inb=(pos>=b.start_pos)&(pos<=b.end_pos)
            for j in np.where(inb&sv)[0]:
                a=vp[:,j].toarray().ravel().astype(float)
                rc=r2vec(a,Gc)                        # r2 of SV vs each hap-cluster
                kbest=int(rc.argmax())
                sv_gap=np.abs(S[:,a==1].mean(1)-S[:,a==0].mean(1)).max()
                rows.append(dict(unit=b.unit,pos=int(pos[j]),n_clusters=Gc.shape[1],
                                 sv_r2_bestcluster=float(rc.max()),
                                 bestcluster_is_mostfit=bool(capgap[kbest]==best_clu_gap),
                                 sv_persite_absmax=float(sv_gap),
                                 bestcluster_gap=float(capgap[kbest]),
                                 block_mostfit_gap=float(best_clu_gap)))
    R=pd.DataFrame(rows)
    print(f"top-JOINT SVs analysed at HAPLOTYPE level: {len(R)}  (blocks median {R.n_clusters.median():.0f} hap-clusters)")
    print("\n(A) does the SV cleanly TAG a haplotype (hap-cluster)?")
    print(f"    r2(SV, best-matching hap-cluster): median {R.sv_r2_bestcluster.median():.3f} | "
          f">=0.9 in {int((R.sv_r2_bestcluster>=0.9).sum())}/{len(R)} | >=0.5 in {int((R.sv_r2_bestcluster>=0.5).sum())}/{len(R)}")
    print("\n(B) is the SV's haplotype the SELECTED (most-fit) one in its block?")
    print(f"    SV's best-tagged cluster == block's most-fit cluster: {int(R.bestcluster_is_mostfit.sum())}/{len(R)}")
    print(f"    SV's tagged-cluster fitness gap {R.bestcluster_gap.median():.4f}  vs  "
          f"block's most-fit cluster {R.block_mostfit_gap.median():.4f}")
    R.to_csv(f"{lib.GEA}/sv_adaptive/sv_haplotype_check.csv",index=False)

    # (C) block-wide variant r2 vs block size, chr1 -- explains the "0.9 block but 0.38 median"
    cl="chr1"; meta=np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz",allow_pickle=True)
    pos=meta["pos"].astype(np.int64); vp=sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
    vc=sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
    na=np.asarray(vp.sum(0)).ravel(); nc=np.asarray(vc.sum(0)).ravel(); common=(na>=12)&(na<=219)&(nc/231>=0.9)
    bl=pd.read_csv(f"{lib.GEA}/blocks_mcf90/chr1_clq0.9_blocks_clq0.9.tsv",sep="\t")
    import numpy.random as npr; rng=npr.default_rng(0)
    print("\n(C) within-block variant founder-r² vs block size (chr1 sample):")
    for lo,hi,lab in [(2,8,"small (2-8 var)"),(8,25,"med (8-25)"),(25,80,"large (25-80)"),(80,10**9,"huge (80+)")]:
        sel=bl[(bl.n_variants>=lo)&(bl.n_variants<hi)]
        if not len(sel): continue
        sel=sel.sample(min(40,len(sel)),random_state=0)
        meds=[]; fr=[]
        for _,b in sel.iterrows():
            jj=np.where((pos>=b.start_pos)&(pos<=b.end_pos)&common)[0]
            if len(jj)<2: continue
            M=vp[:,jj].toarray().astype(float); C=np.corrcoef(M.T)**2
            iu=np.triu_indices(len(jj),1)
            if len(iu[0]): meds.append(np.median(C[iu])); fr.append((C[iu]>=0.9).mean())
        if meds:
            print(f"    {lab:16s}: median pairwise r² {np.median(meds):.2f} | frac pairs r²>=0.9 {np.median(fr):.2f}  (n={len(meds)} blocks)")


if __name__=="__main__":
    main()
