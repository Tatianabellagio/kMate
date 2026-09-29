#!/usr/bin/env python
"""Is the origin-climate x inserted-sequence relationship a PANGENOME CONSTRUCTION artifact?

The earlier divergence control used `n_snp` = ALT count at *common* panel SNPs (MAC>=12).
That is a weak instrument for this question, for three reasons:
  - it is MAC-filtered, so the rare and private variation that actually marks divergent
    lineages is excluded (median 284k vs 524k unfiltered);
  - it is a SNP count, while the worry is about how much of a genome the graph can express;
  - it does not measure the real construction asymmetry at all. The graph was built from 82
    assemblies. What should matter for a founder -- above all a PanGenie one -- is how close
    it is to THOSE 82. A founder far from every assembly has insertions the graph cannot
    contain, regardless of its SNP count.

Four instruments, in increasing relevance:

  D1  UNFILTERED DIVERGENCE   n_snp_all / n_alt_all -- ALT over every segregating record, no
      MAC floor. A genuine divergence proxy rather than a common-allele-load statistic.

  D2  GRAPH REPRESENTATION    per founder, kinship (K_snp) to the 80 cactus/assembly
      founders, self-excluded: max (nearest assembly) and mean (overall representation).
      This is the variable the construction worry is actually about.

  D3  WITHIN THE PANGENIE HALF ONLY. Cactus founders are trivially at distance zero from the
      assembly set -- they ARE it -- so pooling both halves blurs D2. The 151 PG founders,
      whose insertions can only be genotyped from what other accessions contributed, are
      where a representation artifact would show.

  D4  ASSEMBLY SIZE (cactus only, n<=80). The most direct genome-content control available:
      total chromosome-level assembly length from the .fai indexes. Tests both "do
      cold-origin accessions simply have bigger genomes?" and "is inserted sequence just
      assembly size?". Only possible for founders that HAVE an assembly.

Env: kmate. Reads founder_sv_content.npz, founder_climate_confound.npz, class_grms.npz,
the assembly .fai indexes, and data/request_assemblies_for_Moi.csv.
Writes results/sv_adaptive/founder_graph_representation.{npz,csv} + printed report.
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PROJ = lib.PROJ
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
VAREXP = f"{lib.GEA}/r3_persite_gwas/results/varexp"
ASM_DIR = "/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only"
CHROMS = {"Chr1", "Chr2", "Chr3", "Chr4", "Chr5"}


def partial_spearman(x, y, Z):
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    RZ = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    A = np.column_stack([np.ones(len(rx)), RZ])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None); by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex = rx - A @ bx; ey = ry - A @ by
    r = float(np.corrcoef(ex, ey)[0, 1]); dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof))


def rp(x, y):
    s = stats.spearmanr(x, y)
    return s.correlation, s.pvalue


def main():
    C = np.load(f"{OUT}/founder_climate_confound.npz", allow_pickle=True)
    S = np.load(f"{OUT}/founder_sv_content.npz", allow_pickle=True)
    GR = np.load(f"{VAREXP}/class_grms.npz", allow_pickle=True)
    keep = C["keep"].astype(bool)
    fid = C["founders"].astype("U6")[keep]
    home1 = C["home_bio1"][keep]; home18 = C["home_bio18"][keep]
    g1 = C["gamma_bio1"][keep]
    cac = C["is_cactus"][keep].astype(bool)
    M = {k: S[k][keep].astype(float) for k in S.files if k != "founders"}
    KB = M["bp_ins"] / 1000.0
    HOME = np.column_stack([home1, home18])
    rows = []

    # ------------------------------------------------------------------ D1
    print("=" * 92)
    print("D1  Unfiltered divergence (rare + private variants included)")
    for nm in ("n_snp", "n_snp_all", "n_alt_all"):
        v = M[nm]
        r_h, p_h = rp(v, home1); r_g, p_g = rp(v, g1); r_k, p_k = rp(v, KB)
        tag = "  (MAC>=12, the OLD proxy)" if nm == "n_snp" else ""
        print(f"  {nm:<12} median {np.median(v):>10,.0f} | vs origin bio1 {r_h:+.3f} (p={p_h:.3f})"
              f" | vs gamma {r_g:+.3f} (p={p_g:.3f}) | vs kb inserted {r_k:+.3f}{tag}")
        rows.append(dict(test="D1", var=nm, rho_home=r_h, p_home=p_h, rho_gamma=r_g,
                         p_gamma=p_g, rho_kb=r_k))
    for nm in ("n_snp_all", "n_alt_all"):
        r, p = partial_spearman(KB, g1, M[nm][:, None])
        r2, p2 = partial_spearman(KB, home1, M[nm][:, None])
        print(f"    kb~gamma  | {nm} held fixed: {r:+.3f} (p={p:.1e})   "
              f"kb~origin | {nm}: {r2:+.3f} (p={p2:.1e})")
        rows.append(dict(test="D1_partial", var=nm, rho_gamma=r, p_gamma=p,
                         rho_home=r2, p_home=p2))

    # ------------------------------------------------------------------ D2
    print("\n" + "=" * 92)
    print("D2  Graph representation: kinship to the assembly (cactus) founders")
    K = GR["K_snp"].astype(float); K = K / np.mean(np.diag(K))
    Ks = K[np.ix_(keep, keep)]
    cidx = np.where(cac)[0]
    max_kin = np.empty(len(fid)); mean_kin = np.empty(len(fid))
    for i in range(len(fid)):
        others = cidx[cidx != i]
        max_kin[i] = Ks[i, others].max()
        mean_kin[i] = Ks[i, others].mean()
    print(f"  max kinship to an assembly founder : cactus median {np.median(max_kin[cac]):.3f} | "
          f"PG median {np.median(max_kin[~cac]):.3f}")
    for nm, v in (("max_kin_to_assembly", max_kin), ("mean_kin_to_assembly", mean_kin)):
        r_h, p_h = rp(v, home1); r_k, p_k = rp(v, KB); r_g, p_g = rp(v, g1)
        print(f"  {nm:<22} vs origin bio1 {r_h:+.3f} (p={p_h:.3f}) | vs kb inserted "
              f"{r_k:+.3f} (p={p_k:.3f}) | vs gamma {r_g:+.3f} (p={p_g:.3f})")
        rows.append(dict(test="D2", var=nm, rho_home=r_h, p_home=p_h, rho_kb=r_k,
                         p_kb=p_k, rho_gamma=r_g, p_gamma=p_g))
    r, p = partial_spearman(KB, home1, np.column_stack([max_kin, mean_kin]))
    r2, p2 = partial_spearman(KB, g1, np.column_stack([max_kin, mean_kin]))
    print(f"    kb~origin | graph representation held fixed: {r:+.3f} (p={p:.1e})  "
          f"(was {rp(KB, home1)[0]:+.3f})")
    print(f"    kb~gamma  | graph representation held fixed: {r2:+.3f} (p={p2:.1e})  "
          f"(was {rp(KB, g1)[0]:+.3f})")
    rows.append(dict(test="D2_partial", var="graph_repr", rho_home=r, p_home=p,
                     rho_gamma=r2, p_gamma=p2))

    # ------------------------------------------------------------------ D3
    print("\n" + "=" * 92)
    print(f"D3  Within the PanGenie half only (n={int((~cac).sum())}, no assemblies)")
    pg = ~cac
    for nm, v in (("max_kin_to_assembly", max_kin), ("mean_kin_to_assembly", mean_kin),
                  ("n_snp_all", M["n_snp_all"])):
        r_h, p_h = rp(v[pg], home1[pg]); r_k, p_k = rp(v[pg], KB[pg])
        print(f"  {nm:<22} vs origin bio1 {r_h:+.3f} (p={p_h:.3f}) | vs kb inserted {r_k:+.3f} (p={p_k:.3f})")
        rows.append(dict(test="D3", var=nm, rho_home=r_h, p_home=p_h, rho_kb=r_k, p_kb=p_k))
    r_pg, p_pg = rp(KB[pg], home1[pg])
    rc, pc = partial_spearman(KB[pg], home1[pg],
                              np.column_stack([max_kin[pg], mean_kin[pg], M["n_snp_all"][pg]]))
    print(f"  kb~origin within PG: raw {r_pg:+.3f} (p={p_pg:.1e})  ->  "
          f"| graph repr + unfiltered divergence: {rc:+.3f} (p={pc:.1e})")
    rows.append(dict(test="D3_partial", var="kb_origin_pg", rho_home=r_pg, p_home=p_pg,
                     rho_kb=rc, p_kb=pc))

    # ------------------------------------------------------------------ D4
    print("\n" + "=" * 92)
    print("D4  Assembly size (cactus founders only)")
    req = pd.read_csv(f"{PROJ}/data/request_assemblies_for_Moi.csv")
    req["Assembly_ID"] = req["Assembly_ID"].astype(str)
    req["Accession_ID"] = req["Accession_ID"].astype(str)
    size_by_asm = {}
    for fai in glob.glob(f"{ASM_DIR}/*.chr.fa.fai"):
        asm = os.path.basename(fai).split(".")[0]
        tot = 0
        for line in open(fai):
            f0 = line.split("\t")
            if f0[0] in CHROMS:
                tot += int(f0[1])
        if tot > 0:
            size_by_asm[asm] = tot
    print(f"  {len(size_by_asm)} assemblies indexed in {os.path.basename(ASM_DIR)}")
    # accession -> assembly size, preferring the flagged-best assembly
    best = req[req.Best_for_Moi == "yes"]
    acc2size = {}
    for src in (best, req):
        for _, r0 in src.iterrows():
            a = r0["Accession_ID"]
            if a not in acc2size and r0["Assembly_ID"] in size_by_asm:
                acc2size[a] = size_by_asm[r0["Assembly_ID"]]
    asm_mb = np.array([acc2size.get(f, np.nan) for f in fid]) / 1e6
    have = np.isfinite(asm_mb)
    print(f"  matched {int(have.sum())} founders to an assembly "
          f"({int((have & cac).sum())} of the {int(cac.sum())} cactus founders)")
    if have.sum() >= 20:
        print(f"  assembly size (Mb, 5 chroms): median {np.nanmedian(asm_mb[have]):.2f} "
              f"range {np.nanmin(asm_mb[have]):.2f}-{np.nanmax(asm_mb[have]):.2f}  "
              f"(TAIR10 = 119.15)")
        for nm, v in (("ecotype origin bio1", home1), ("kb inserted", KB),
                      ("gamma_bio1", g1), ("n_snp_all", M["n_snp_all"])):
            r, p = rp(asm_mb[have], v[have])
            print(f"    assembly size vs {nm:<22} rho={r:+.3f} p={p:.4f}  [n={int(have.sum())}]")
            rows.append(dict(test="D4", var=nm, rho=r, p=p, n=int(have.sum())))
        r, p = partial_spearman(KB[have], home1[have], asm_mb[have][:, None])
        r2, p2 = partial_spearman(KB[have], g1[have], asm_mb[have][:, None])
        print(f"    kb~origin | assembly size held fixed: {r:+.3f} (p={p:.1e})  "
              f"(raw within this subset {rp(KB[have], home1[have])[0]:+.3f})")
        print(f"    kb~gamma  | assembly size held fixed: {r2:+.3f} (p={p2:.1e})  "
              f"(raw within this subset {rp(KB[have], g1[have])[0]:+.3f})")
        rows.append(dict(test="D4_partial", var="kb_vs_origin", rho=r, p=p))
        rows.append(dict(test="D4_partial", var="kb_vs_gamma", rho=r2, p=p2))
    else:
        print("  too few matched assemblies to test")

    np.savez_compressed(f"{OUT}/founder_graph_representation.npz",
                        founders=fid, max_kin=max_kin, mean_kin=mean_kin,
                        asm_mb=asm_mb, is_cactus=cac, kb_ins=KB,
                        n_snp_all=M["n_snp_all"], n_alt_all=M["n_alt_all"])
    pd.DataFrame(rows).to_csv(f"{OUT}/founder_graph_representation.csv", index=False)
    print(f"\n[wrote] {OUT}/founder_graph_representation.npz + .csv")


if __name__ == "__main__":
    main()
