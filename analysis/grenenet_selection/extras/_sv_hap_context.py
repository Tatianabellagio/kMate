#!/usr/bin/env python
"""Is the common-SV / garden-fitness haplotype enrichment a REGIONAL confound?

exact-mac matching (sv_hap_freqrobust) ruled out haplotype frequency. This tests the three
remaining regional confounds on the JOINT high-MAF cell (fitness-selected hap-clusters tag a
common SV ~1.7-2x): pericentromere/genomic context, repetitiveness (TE density), and gene
content. Each is added to the matched null (mac + covariate strata); the signal must survive
all three. Also an arms-only rerun and a gene-content readout of the SVs that drive it.

Env: kmate (MEMB_TAG=clq90). Writes sv_hap_context.csv + sv_hap_context_genes.csv.
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from founder_genotype import build_genotype

CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
PANEL = "panel/arch3"
MAC_MIN = 3; nF = 231; NPERM = 2000
HAPGEA = "results/grenenet_gea/hapfreq_clq90/pipelineB_varlen/hap_gea.csv"
# TAIR10 centromere midpoints (bp)
CEN = {"Chr1": 15_086_000, "Chr2": 3_607_000, "Chr3": 13_588_000,
       "Chr4": 3_956_000, "Chr5": 11_726_000}
GFF = os.path.expanduser("~/ara_key_files/TAIR10_GFF3_genes_transposons.gff")
CLENS = {"Chr1": 30_427_671, "Chr2": 19_698_289, "Chr3": 23_459_830,
         "Chr4": 18_585_056, "Chr5": 26_975_502}


def r2_cols(a, B):
    a = a - a.mean(); Bc = B - B.mean(0)
    num = (a[:, None] * Bc).sum(0) ** 2; den = (a @ a) * (Bc ** 2).sum(0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def sv_r2_for_mac(Gp, regp, sv_mac, want_hits=False):
    M = Gp.shape[1]; sv_r2 = np.zeros(M)
    hits = []                                    # (chrom, pos, ref_len, alt_len, best_r2) for want_hits
    for ch in CHROMS:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        dl = np.abs(al - rl)
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        svj = np.where((dl > 50) & (na >= sv_mac) & (na <= nF - sv_mac) & (nc / nF >= 0.9))[0]
        rc = regp[regp.chrom == ch]
        for u, g in rc.groupby("unit"):
            s0, e0 = g.start.iloc[0], g.end.iloc[0]
            inb = svj[(pos[svj] >= s0) & (pos[svj] <= e0)]
            if len(inb) == 0:
                continue
            cols = g.index.to_numpy(); Bc = Gp[:, cols]
            for j in inb:
                a = vp[:, j].toarray().ravel().astype(float)
                rr = r2_cols(a, Bc)
                sv_r2[cols] = np.maximum(sv_r2[cols], rr)
                if want_hits:
                    hits.append((ch, int(pos[j]), int(rl[j]), int(al[j]), float(rr.max())))
    return (sv_r2, pd.DataFrame(hits, columns=["chrom", "pos", "ref_len", "alt_len", "best_r2"])) if want_hits else sv_r2


def coverage_frac(intervals, s, e, clen):
    """fraction of [s,e] covered by a set of (start,end) intervals, via bp cumsum."""
    cov = np.zeros(clen + 2, np.int8)
    for a, b in intervals:
        cov[max(a, 0):min(b, clen) + 1] = 1
    c = np.concatenate([[0], np.cumsum(cov)])
    return (c[np.minimum(e, clen) + 1] - c[np.maximum(s, 0)]) / np.maximum(e - s, 1)


def qbin(x, nb):
    ed = np.unique(np.quantile(x, np.linspace(0, 1, nb + 1)))
    return np.clip(np.digitize(x, ed[1:-1]), 0, len(ed) - 2)


def strata_null(tagvec, sel, strata, rng, nperm=NPERM):
    tot = np.zeros(nperm)
    idx = np.arange(len(tagvec))
    pools = {}
    for i in sel:
        s = strata[i]
        if s not in pools:
            pools[s] = idx[strata == s]
    for i in sel:
        pool = pools[strata[i]]
        tot += tagvec[pool[rng.integers(0, len(pool), nperm)]]
    return tot / len(sel)


def main():
    os.chdir("/global/scratch/users/tbellg/kmate")
    G, founders, reg = build_genotype()
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    reg["cluster"] = reg.groupby("block").cumcount()
    cnt = G.sum(0); blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk)); sbk = blk[order]
    is_ref = np.zeros(len(cnt), bool)
    is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
    poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
    Gp = G[:, poly].astype(float); regp = reg[poly].reset_index(drop=True)
    raw = lib.multisite_gwas_raw("clq90_pc1")
    assert list(regp.unit.values) == list(raw["unit"])
    p_joint = raw["p_joint"]; mac = raw["mac"].astype(int); M = len(p_joint)
    joint_order = np.argsort(p_joint)

    # ---- per-block covariates -> per-cluster ----
    print("annotating blocks (TE / gene / pericentromere)...")
    gff = pd.read_csv(GFF, sep="\t", header=None, comment="#",
                      names=["chrom", "src", "feat", "start", "end", "sc", "st", "fr", "attr"])
    gff = gff[gff.chrom.isin(CHROMS)]
    genes = gff[gff.feat == "gene"]; tes = gff[gff.feat == "transposable_element"]
    L = pd.read_csv("analysis/grenenet_selection/sv_adaptive/sv_landscape_clq0.9.csv")
    L["unit"] = L.block_id if "block_id" in L else (L.chrom + ":" + L.start.astype(str) + "-" + L.end.astype(str))
    nkept = L.set_index("unit").n_kept.to_dict() if "n_kept" in L else {}

    blocks = regp.drop_duplicates("unit")[["chrom", "start", "end", "unit"]].reset_index(drop=True)
    te_frac = np.zeros(len(blocks)); gene_frac = np.zeros(len(blocks)); gene_n = np.zeros(len(blocks))
    pcen = np.zeros(len(blocks))
    for ch in CHROMS:
        clen = CLENS[ch]
        te_iv = list(zip(tes[tes.chrom == ch].start.to_numpy(), tes[tes.chrom == ch].end.to_numpy()))
        ge = genes[genes.chrom == ch]
        gstart = np.sort(ge.start.to_numpy()); gend = np.sort(ge.end.to_numpy())
        ge_iv = list(zip(ge.start.to_numpy(), ge.end.to_numpy()))
        sel = np.where(blocks.chrom.to_numpy() == ch)[0]
        # precompute cumulative coverage once per chrom
        covte = np.zeros(clen + 2, np.int8)
        for a, b in te_iv:
            covte[max(a, 0):min(b, clen) + 1] = 1
        cte = np.concatenate([[0], np.cumsum(covte)])
        covg = np.zeros(clen + 2, np.int8)
        for a, b in ge_iv:
            covg[max(a, 0):min(b, clen) + 1] = 1
        cg = np.concatenate([[0], np.cumsum(covg)])
        for i in sel:
            s, e = int(blocks.start.iat[i]), int(blocks.end.iat[i])
            ln = max(e - s, 1)
            te_frac[i] = (cte[min(e, clen) + 1] - cte[max(s, 0)]) / ln
            gene_frac[i] = (cg[min(e, clen) + 1] - cg[max(s, 0)]) / ln
            gene_n[i] = (np.searchsorted(gstart, e, "right") - np.searchsorted(gend, s, "left"))
            pcen[i] = abs((s + e) / 2 - CEN[ch])
    blocks["te_frac"] = te_frac; blocks["gene_frac"] = gene_frac
    blocks["gene_n"] = gene_n; blocks["pcen"] = pcen
    blocks["nkept"] = blocks.unit.map(nkept).fillna(0).to_numpy()
    bm = blocks.set_index("unit")
    cu = regp.unit.to_numpy()
    C = {k: bm.loc[cu, k].to_numpy() for k in ["te_frac", "gene_frac", "gene_n", "pcen", "nkept"]}

    n1 = int(round(0.01 * M))
    macbin = qbin(mac, 8)

    rng = np.random.default_rng(0)
    rows = []
    for svmac in (24, 46):
        sv_r2 = sv_r2_for_mac(Gp, regp, svmac)
        sel = joint_order[:n1]
        obs = sv_r2[sel].mean()
        print(f"\n===== JOINT top-1%, SV_MAC={svmac} (obs mean-r²={obs:.4f}) =====")
        # covariate contrast selected vs all
        print("  covariate (selected vs all): "
              + " | ".join(f"{k} {C[k][sel].mean():.3g} vs {C[k].mean():.3g}"
                           for k in ["te_frac", "gene_frac", "gene_n", "pcen", "nkept"]))
        controls = {
            "mac_only":        macbin,
            "mac+size(nkept)": macbin * 10 + qbin(C["nkept"], 5),
            "mac+pericentro":  macbin * 10 + qbin(C["pcen"], 5),
            "mac+TE_frac":     macbin * 10 + qbin(C["te_frac"], 5),
            "mac+gene_n":      macbin * 10 + qbin(C["gene_n"], 5),
        }
        for name, strata in controls.items():
            nul = strata_null(sv_r2, sel, strata, rng)
            med = max(np.median(nul), 1e-12)
            fold = obs / med; p = (1 + (nul >= obs).sum()) / (NPERM + 1)
            print(f"  {name:>16}: x{fold:.2f}  p={p:.4f}")
            rows.append(dict(sv_mac=svmac, control=name, obs=round(obs, 4),
                             fold=round(fold, 3), p_perm=round(p, 4)))
        # arms-only (drop pericentromeric clusters within 3 Mb of centromere)
        arm = C["pcen"] >= 3e6
        selA = joint_order[arm[joint_order]][:n1]
        obsA = sv_r2[selA].mean()
        strata_arm = macbin.copy()
        nul = strata_null(sv_r2, selA, np.where(arm, strata_arm, -1), rng)
        med = max(np.median(nul), 1e-12); foldA = obsA / med
        pA = (1 + (nul >= obsA).sum()) / (NPERM + 1)
        print(f"  {'arms_only':>16}: x{foldA:.2f}  p={pA:.4f}  ({arm[sel].mean()*100:.0f}% of top were on arms)")
        rows.append(dict(sv_mac=svmac, control="arms_only", obs=round(obsA, 4),
                         fold=round(foldA, 3), p_perm=round(pA, 4)))
    pd.DataFrame(rows).to_csv(f"{lib.GEA}/sv_adaptive/sv_hap_context.csv", index=False)
    print(f"\n[wrote] {lib.GEA}/sv_adaptive/sv_hap_context.csv")

    # ---- gene-content readout: SVs tagging a selected haplotype (MAC24, r²>=0.5) ----
    print("\n=== gene content of SVs that tag a fitness-selected haplotype (MAC24, r²>=0.5) ===")
    sv_r2_24, hits = sv_r2_for_mac(Gp, regp, 24, want_hits=True)
    # SVs whose best-tagged cluster is a selected (top-1% JOINT) cluster: recompute per-SV best cluster
    sel_set = set(joint_order[:n1].tolist())
    # tag SVs that reach r²>=0.5 with ANY cluster; annotate genes on those in selected blocks
    hi = hits[hits.best_r2 >= 0.5].copy()
    hi = lib.annotate_svs(hi, flank=0)
    genic = (hi.n_genes > 0).mean()
    print(f"  {len(hi)} common SVs tag a haplotype at r²>=0.5 | {genic*100:.0f}% overlap a gene")
    top = hi[hi.n_genes > 0].sort_values("best_r2", ascending=False).head(20)
    print(top[["chrom", "pos", "ref_len", "alt_len", "best_r2", "gene", "gene_name"]].to_string(index=False))
    hi.to_csv(f"{lib.GEA}/sv_adaptive/sv_hap_context_genes.csv", index=False)
    print(f"[wrote] {lib.GEA}/sv_adaptive/sv_hap_context_genes.csv")


if __name__ == "__main__":
    main()
