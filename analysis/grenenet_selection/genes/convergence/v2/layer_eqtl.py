#!/usr/bin/env python
"""v2 layer: eQTL tag-SNP readout, keyed on `unit`.  Port of eqtl_tag_snps.py.

Tag SNP per unit (first that applies):
  shortread_r2  best_r2_shortread_pos (independent GrENE-Net short-read SNP LD), if testable
  panel_r2      best_r2_snp_pos (panel SNP LD in the founders), if testable
  pool_af       SNP within +-WIN (25 kb) of the unit with MAF(p0)>0.05 whose gen9 ALT-frequency vector over the 352 pools
                correlates best (|r|) with the unit's. Unit vector = min(sum of member-record frequencies, 1), exactly as
                event_gea/build_event_matrix.py. Not LD.
eQTL readout = old rule: cis_eqtl_all.csv rows of the TARGET gene (gemma -lm 4, gene body +-10 kb, 1001g SNPs); exact tag SNP if
tested, else the nearest tested SNP (status says which); nlp, beta, pct = % of the gene's tested SNPs with lower nlp.
Where the tag lies outside the scan window, no readout (status tag_outside_scan_window). Genes without a scan: status
no_eqtl_scan_for_gene (the eQTL scan only exists for the genes in functional_track_candidates.csv).
env: kmate, compute node.  Writes results/candidate_list_v2_events_cr50/v2_layer_eqtl.csv
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import units as UN
lib = UN.lib
CONV = UN.CONV
OUT = f"{CONV}/results/candidate_list_v2_events_cr50"
WIN = 25_000


def main():
    F = pd.read_csv(f"{OUT}/v2_functional_units.csv", usecols=["unit", "target_gene"])
    V = pd.read_csv(f"{OUT}/v2_layer_snp_visibility.csv").set_index("unit")
    I = UN.info(F.unit.values).set_index("unit")
    D = F.set_index("unit").join(I[["chrom", "pos", "ref_len", "alt_len", "dlen", "n_members"]]).join(V, how="left")
    rows = {}
    need = []
    for u, r in D.iterrows():
        sr, pr = r.best_r2_shortread_pos, r.best_r2_snp_pos
        if (not bool(r.tag_shortread_untestable)) and np.isfinite(sr) and sr > 0 and np.isfinite(r.best_r2_shortread):
            rows[u] = (int(sr), float(r.best_r2_shortread), "shortread_r2")
        elif (not bool(r.tag_untestable)) and np.isfinite(pr) and pr > 0 and np.isfinite(r.best_r2_snp):
            rows[u] = (int(pr), float(r.best_r2_snp), "panel_r2")
        else:
            need.append(u)
    if os.environ.get("LAYER_EQTL_TEST_POOLAF"):          # smoke test of the fallback only; output goes to the scratchpad
        k = int(os.environ["LAYER_EQTL_TEST_POOLAF"]); need = list(D.index[:k]); [rows.pop(u) for u in need]
    print(f"{len(D)} units; shortread {sum(v[2]=='shortread_r2' for v in rows.values())}, "
          f"panel {sum(v[2]=='panel_r2' for v in rows.values())}, pool_af needed {len(need)}", flush=True)

    if need:
        Mn = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
        Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy", mmap_mode="r")
        idx_s = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
        p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
        cs, ps_ = idx_s["chrom"].astype(str), idx_s["pos"]
        maf0 = np.minimum(p0s, 1 - p0s)
        for u in need:
            r = D.loc[u]
            mem = np.sort(UN.member_rows(u))
            v = np.minimum(np.asarray(Mn[:, mem], np.float64).sum(1), 1.0)      # sum of members, capped at 1
            sel = np.where((cs == r.chrom) & (ps_ >= r.pos - WIN) & (ps_ <= r.pos + WIN) & (maf0 > 0.05))[0]
            if not len(sel) or not np.isfinite(v).all() or np.ptp(v) == 0:
                rows[u] = (None, np.nan, "none"); continue
            S = np.asarray(Ms[:, sel], np.float64)
            vz = v - v.mean(); Sz = S - np.nanmean(S, 0)
            with np.errstate(all="ignore"):
                rr = np.nansum(Sz * vz[:, None], 0) / np.sqrt(np.nansum(Sz ** 2, 0) * np.sum(vz ** 2))
            if not np.isfinite(rr).any():
                rows[u] = (None, np.nan, "none"); continue
            j = int(np.nanargmax(np.abs(rr)))
            rows[u] = (int(ps_[sel][j]), float(rr[j]), "pool_af")

    Q = pd.read_csv(f"{CONV}/results/eqtl/cis_eqtl_all.csv", usecols=["gene", "ps", "beta", "se", "p_wald", "af", "nlp"])
    Qg = {g: d.reset_index(drop=True) for g, d in Q.groupby("gene")}
    out = []
    for u, r in D.iterrows():
        tp, st, meth = rows[u]
        o = dict(unit=u, chrom=r.chrom, pos=int(r.pos), ref_len=int(r.ref_len), alt_len=int(r.alt_len), n_members=int(r.n_members),
                 target_gene=r.target_gene, tag_pos=tp if tp is not None else np.nan, tag_method=meth,
                 tag_stat=round(st, 4) if np.isfinite(st) else np.nan,
                 tag_dist_bp=(tp - int(r.pos)) if tp is not None else np.nan)
        q = Qg.get(r.target_gene)
        if tp is None:
            o["eqtl_status"] = "no_tag_snp"
        elif q is None:
            o["eqtl_status"] = "no_eqtl_scan_for_gene"
        else:
            o["eqtl_gene_n_snps"] = len(q); o["eqtl_gene_peak_nlp"] = round(q.nlp.max(), 3)
            h = q[q.ps == tp]
            if len(h):
                h = h.iloc[0]; o["eqtl_status"] = "exact_tag_snp"; o["eqtl_snp_pos"] = int(h.ps)
            elif q.ps.min() <= tp <= q.ps.max():
                h = q.loc[(q.ps - tp).abs().idxmin()]; o["eqtl_status"] = "nearest_tested_snp"; o["eqtl_snp_pos"] = int(h.ps)
            else:
                h = None; o["eqtl_status"] = "tag_outside_scan_window"
            if h is not None:
                o.update(eqtl_nlp=round(h.nlp, 3), eqtl_beta=h.beta, eqtl_se=h.se, eqtl_p=h.p_wald, eqtl_snp_af=h.af,
                         eqtl_pct=round(100 * (q.nlp < h.nlp).mean(), 1), eqtl_snp_dist_to_tag=int(h.ps) - tp)
        out.append(o)
    T = pd.DataFrame(out)
    T.to_csv(os.environ.get("LAYER_EQTL_OUT", f"{OUT}/v2_layer_eqtl.csv"), index=False)
    print(T.tag_method.value_counts().to_string()); print(T.eqtl_status.value_counts().to_string())


if __name__ == "__main__":
    main()
