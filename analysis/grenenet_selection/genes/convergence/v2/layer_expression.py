#!/usr/bin/env python
"""v2 layer -- does carrying the unit change its target gene's expression? (carrier vs non-carrier, 1001 Transcriptomes)

Same test as carrier_expression_test.py, on UNITS: carriers = founders carrying ANY member allele of the unit (exact union from var_pa, via
v2/units.py -- no VCF size-tolerance guess), non-carriers = founders called at the unit (any member) that do not carry it. Expression is
the accession's raw count / d_log2_batch of the unit's TARGET gene (the gene the functional handle is about) from data/eqtl_v2/ (exported by
export_expression_v2.R for the 946 ranked genes).
  gene level       Mann-Whitney, carriers vs called non-carriers, on raw counts (silent accessions stay in); medians; fraction silent
  isoform level    for multi-isoform genes, logit(isoform / sum) per accession, bottom 5% of total expression masked; best p, Bonferroni within gene
  lineage control  the SAME split tested on 300 random expressed genes; p_emp = fraction of control genes differing at least as strongly
                   (carriers share ancestry, and ancestry shifts expression genome-wide; only a difference that beats the control is attributable)
Tested units: every unit with a functional handle (v2_functional_units.csv). Output: results/candidate_list_v2_events_cr50/v2_layer_expression.csv
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from scipy.special import logit
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import units as UN                                              # noqa: E402
V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")
EQ = os.path.abspath(os.path.join(CONV, "..", "..", "..", "..", "data", "eqtl_v2"))
MIN_N = 3


def mw(a, b):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < MIN_N or len(b) < MIN_N:
        return np.nan
    return mannwhitneyu(a, b).pvalue


def main():
    F = pd.read_csv(f"{V2}/v2_functional_units.csv")[["unit", "target_gene", "vclass", "ftier"]].drop_duplicates()
    print(f"{len(F):,} (unit, target gene) pairs; loading expression ...", flush=True)
    Gx = pd.read_csv(f"{EQ}/expr_genes.csv"); Ix = pd.read_csv(f"{EQ}/expr_isoforms.csv"); Cx = pd.read_csv(f"{EQ}/expr_control.csv")
    graw = Gx.pivot_table(index="acc", columns="id", values="raw", aggfunc="mean")      # replicates of one accession averaged
    iraw = Ix.pivot_table(index="acc", columns="id", values="value", aggfunc="mean")
    ctl = Cx.pivot_table(index="acc", columns="id", values="value", aggfunc="mean")
    print(f"expression: {graw.shape[1]} genes x {graw.shape[0]} accessions; control {ctl.shape[1]} genes", flush=True)

    founders, car, cal = UN.carrier_matrix(F.unit.values)
    fid = np.array([int(x) if str(x).isdigit() else -1 for x in founders])
    rows = []
    for j, r in enumerate(F.itertuples()):
        out = dict(unit=r.unit, gene=r.target_gene, n_carriers_panel=int(car[:, j].sum()), n_called_noncarrier_panel=int((cal[:, j] & ~car[:, j]).sum()))
        if r.target_gene not in graw.columns:
            out["status"] = "gene not in 1001T"; rows.append(out); continue
        cset = set(fid[car[:, j]]); rset = set(fid[cal[:, j] & ~car[:, j]])
        rg = graw[r.target_gene]
        cc = [a for a in rg.index if a in cset and np.isfinite(rg[a])]
        rr = [a for a in rg.index if a in rset and np.isfinite(rg[a])]
        med_c, med_r = rg.reindex(cc).median(), rg.reindex(rr).median()
        out.update(n_car_expr=len(cc), n_ref_expr=len(rr), med_raw_car=med_c, med_raw_ref=med_r,
                   frac_car_silent=float((rg.reindex(cc) < 1).mean()) if cc else np.nan,
                   frac_ref_silent=float((rg.reindex(rr) < 1).mean()) if rr else np.nan,
                   fold=(med_c + 1) / (med_r + 1) if cc and rr else np.nan)
        out["p_gene"] = mw(rg.reindex(cc).values, rg.reindex(rr).values)
        if len(cc) >= MIN_N and len(rr) >= MIN_N:
            pc = np.array([mw(ctl[c].reindex(cc).values, ctl[c].reindex(rr).values) for c in ctl.columns]); pc = pc[np.isfinite(pc)]
            out["ctl_frac_p05"] = float((pc < 0.05).mean())
            if np.isfinite(out["p_gene"]):
                out["p_emp"] = float((np.sum(pc <= out["p_gene"]) + 1) / (len(pc) + 1))
        iso = [c for c in iraw.columns if c.rsplit(".", 1)[0] == r.target_gene]
        if len(iso) >= 2 and len(cc) >= MIN_N:
            X = iraw[iso]; tot = X.sum(1); X = X[tot > tot.quantile(0.05)]
            R = logit((X.div(X.sum(1), axis=0)).clip(0.01, 0.99))
            ps = {c: mw(R[c].reindex(cc).values, R[c].reindex(rr).values) for c in iso}; ps = {k: p for k, p in ps.items() if np.isfinite(p)}
            if ps:
                k = min(ps, key=ps.get); out.update(n_isoforms=len(iso), iso_best=k, p_iso=min(1.0, ps[k] * len(iso)))
        out["status"] = "tested" if np.isfinite(out.get("p_gene", np.nan)) else f"too few carriers with expression ({len(cc)})"
        rows.append(out)
        if (j + 1) % 300 == 0:
            print(f"  {j + 1:,}/{len(F):,}", flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(f"{V2}/v2_layer_expression.csv", index=False)
    print(D.status.value_counts().to_string())
    t = D[D.status == "tested"]
    print(f"\ntested {len(t):,} | p_emp <= 0.05: {int((t.p_emp <= 0.05).sum())} | fold < 0.5 (carriers lower): {int((t.fold < 0.5).sum())}")
    print(f"wrote {V2}/v2_layer_expression.csv")


if __name__ == "__main__":
    main()
