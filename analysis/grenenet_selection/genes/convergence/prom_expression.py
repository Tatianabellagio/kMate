#!/usr/bin/env python
"""Which gene does a shared-promoter variant act on? Test every one of them.

126 of 909 promoter variants sit in 2-4 promoters at once -- divergent (head-to-head) gene
pairs share their upstream DNA, and in Arabidopsis that stretch often acts as a
bidirectional promoter driving both genes. `classify()` keeps the first match by coordinate
as the primary `gene`, which is where the array starts, not a biological choice; a
nearest-TSS tie-break would move the primary for 58 of them. Either rule assumes one true
target, which is the assumption in question.

So nothing is picked. For every gene in `prom_genes` this runs the carrier-vs-non-carrier
expression test exactly as evidence_matrix.py stage C does (1001T rosette raw counts,
Mann-Whitney, the same 300-gene lineage control on the same carrier split) and lets the
data say which gene moves:

  primary_only      the listed gene responds, its neighbour does not
  alternative_only  the neighbour responds and the listed gene does not -- the call names
                    the wrong gene
  both              a bidirectional effect
  neither           unresolved, as before

"Responds" is the L_expression criterion (BH q < 0.10 and lineage-control p <= 0.05), with
BH taken over all gene tests in this run, on a gene that is not silent (median raw >= 1).
ATAC and motif turnover are deliberately not used: they sit in the shared DNA and so
describe both genes equally.

Writes results/prom_expression.csv (one row per variant x promoter gene) and
results/prom_expression_variants.csv (one row per variant). env: kmate. Compute node.
"""
import os
import sys
import numpy as np
import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import carrier_expression_test as CT                             # noqa: E402

OUT = f"{HERE}/results"
PROJ = CT.PROJ
KEY = ["chrom", "pos", "ref_len", "alt_len"]


def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); ok = np.isfinite(p)
    if ok.sum():
        v = p[ok]; o = np.argsort(v); m = len(v)
        r = np.minimum.accumulate((v[o] * m / np.arange(1, m + 1))[::-1])[::-1]
        out = np.empty(m); out[o] = np.minimum(r, 1); q[ok] = out
    return q


def main():
    A = pd.read_csv(f"{OUT}/evidence_matrix.csv")
    V = pd.read_csv(f"{OUT}/variants_classified.csv", low_memory=False)
    J = A[KEY + ["store_row", "target_gene", "mode", "n_lines"]].merge(
        V[KEY + ["prom_genes", "n_prom_genes", "prom_gene_nearest_tss"]].drop_duplicates(KEY),
        on=KEY, how="left")
    J = J[J.n_prom_genes.fillna(0) > 1].reset_index(drop=True)
    print(f"{len(J)} matrix variants sit in more than one promoter", flush=True)

    Gx = pd.read_csv(f"{CT.EQ}/expr_genes.csv")
    Cx = pd.read_csv(f"{CT.EQ}/expr_control.csv")
    graw = Gx.pivot_table(index="acc", columns="id", values="raw", aggfunc="mean")
    ctl = Cx.pivot_table(index="acc", columns="id", values="value", aggfunc="mean")
    med = graw.median()

    vfs, rows = {}, []
    for v in J.itertuples():
        f = f"{PROJ}/panel/arch3/{v.chrom.lower()}/merged_231_{v.chrom.lower()}_final.vcf.gz"
        vf = vfs.setdefault(f, pysam.VariantFile(f))
        car, ref = CT.carriers(vf, v.chrom, int(v.pos), int(v.ref_len), int(v.alt_len))
        genes = str(v.prom_genes).split(",")            # nearest-TSS first
        for rank, g in enumerate(genes, 1):
            o = dict(store_row=v.store_row, chrom=v.chrom, pos=int(v.pos),
                     ref_len=int(v.ref_len), alt_len=int(v.alt_len), gene=g,
                     tss_rank=rank, is_primary=(g == v.target_gene), n_prom_genes=len(genes))
            if g not in graw.columns:
                o["expr_status"] = "not in 1001T"; rows.append(o); continue
            r = graw[g]
            cc = [a for a in r.index if a in car and np.isfinite(r[a])]
            rr = [a for a in r.index if a in ref and np.isfinite(r[a])]
            o.update(gene_median_raw=float(med[g]), n_car_expr=len(cc), n_ref_expr=len(rr),
                     med_raw_car=r.reindex(cc).median(), med_raw_ref=r.reindex(rr).median())
            o["expr_fold"] = ((o["med_raw_car"] + 1) / (o["med_raw_ref"] + 1)
                              if cc and rr else np.nan)
            p = CT.mw(r.reindex(cc).values, r.reindex(rr).values)
            o["p_expr"] = p
            if np.isfinite(p):
                pc = np.array([CT.mw(ctl[c].reindex(cc).values, ctl[c].reindex(rr).values)
                               for c in ctl.columns])
                pc = pc[np.isfinite(pc)]
                o["p_expr_emp"] = float((np.sum(pc <= p) + 1) / (len(pc) + 1))
                o["expr_status"] = "tested"
            else:
                o["expr_status"] = f"too few expressed carriers ({len(cc)})"
            rows.append(o)

    T = pd.DataFrame(rows)
    T["q_expr"] = bh(T.get("p_expr", pd.Series(np.nan, index=T.index)))
    # silent genes cannot "respond": on median raw counts < 1 a carrier/non-carrier split is
    # a contrast between zeros and near-zeros, and it clears the lineage control anyway. The
    # first run called three such genes responders -- AT3G10750 (median 0.12), which would
    # have reassigned ASD1's variant, and AT4G27890 (0.00), which did the same for all 14
    # records of the AT4G27885 locus. Same threshold as score_evidence's gene_silent flag.
    T["gene_silent"] = T.gene_median_raw < 1
    T["responds"] = (T.q_expr < 0.10) & (T.p_expr_emp <= 0.05) & ~T.gene_silent
    T.to_csv(f"{OUT}/prom_expression.csv", index=False)

    def verdict(d):
        p = d[d.is_primary].responds.any(); a = d[~d.is_primary].responds.any()
        if p and a:
            return "both"
        if p:
            return "primary_only"
        if a:
            return "alternative_only"
        tested = (d.expr_status == "tested").sum()
        return "neither" if tested >= 2 else "neither (partly untestable)"
    W = (T.groupby("store_row").apply(lambda d: pd.Series(dict(
            genes=",".join(d.gene), primary=",".join(d[d.is_primary].gene),
            responders=",".join(d[d.responds].gene), verdict=verdict(d))),
            include_groups=False).reset_index())
    W = W.merge(J[["store_row", "target_gene", "mode", "n_lines"]], on="store_row", how="left")
    W.to_csv(f"{OUT}/prom_expression_variants.csv", index=False)

    print(f"\n{len(T)} variant x gene tests; tested {int((T.expr_status == 'tested').sum())}, "
          f"not in 1001T {int((T.expr_status == 'not in 1001T').sum())}, "
          f"too few carriers {int(T.expr_status.astype(str).str.startswith('too few').sum())}")
    print("verdicts:", W.verdict.value_counts().to_dict())
    k = W[W.verdict.isin(["alternative_only", "both", "primary_only"])]
    print(k.sort_values(["verdict", "n_lines"], ascending=[True, False]).to_string(index=False))


if __name__ == "__main__":
    main()
