#!/usr/bin/env python
"""Does carrying the variant change the gene's expression or splicing? -- tested on the
variant itself, in the founders that have both a panel genotype and 1001T expression.

The SNP cis-eQTL reads expression at a tagging SNP, which fails exactly where our variants
are interesting (GPX6: no founder r2; PGIP1: the eQTL sits on another haplotype). For
large-effect variants the direct contrast works even with few carriers -- RLM3's deletion
carriers have ~0 counts against a median of 167.

Per variant (split panel records merged, as everywhere else):
  gene level     Mann-Whitney on d_log2_batch, carriers vs called non-carriers; median
                 raw counts in each group; fraction of carriers with ~no expression
  isoform level  for multi-isoform genes, logit(isoform / sum of isoforms) per accession,
                 bottom 5% of total expression masked (as the upstream sQTL); min p over
                 isoforms, Bonferroni-corrected within the gene
  lineage control  the SAME carrier split tested on 300 random expressed genes. Carriers
                 often share ancestry, and ancestry alone shifts expression genome-wide;
                 `p_emp` = fraction of control genes differing at least as strongly. Only
                 a difference that beats the control is attributable to the variant.

Writes results/carrier_expression_test.csv.
env: kmate. Compute node, a few minutes.
"""
import os, sys
import numpy as np, pandas as pd, pysam
from scipy.stats import mannwhitneyu
from scipy.special import logit
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = f"{HERE}/results"
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
EQ = f"{PROJ}/data/eqtl"
MIN_N = 3


def carriers(vf, ch, pos, rl, al):
    size = al - rl
    tol = max(10, int(abs(size) * 0.05)) if abs(size) > 50 else 0
    car, cal = set(), set()
    for r in vf.fetch(ch, pos - 1, pos):
        if r.pos != pos:
            continue
        for a in r.alts or ():
            s = len(a) - len(r.ref)
            if np.sign(s) != np.sign(size) or abs(s - size) > tol:
                continue
            for smp in vf.header.samples:
                g = r.samples[smp].get("GT", (None,))[0]
                if g is not None:
                    cal.add(int(smp))
                if g == 1:
                    car.add(int(smp))
    return car, cal - car


def mw(a, b):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < MIN_N or len(b) < MIN_N:
        return np.nan
    return mannwhitneyu(a, b).pvalue


def main():
    V = pd.read_csv(f"{OUT}/expression_test_variants.csv")
    Gx = pd.read_csv(f"{EQ}/expr_genes.csv")
    Ix = pd.read_csv(f"{EQ}/expr_isoforms.csv")
    Cx = pd.read_csv(f"{EQ}/expr_control.csv")
    # replicate samples of one accession averaged
    graw = Gx.pivot_table(index="acc", columns="id", values="raw", aggfunc="mean")
    gdl = Gx.pivot_table(index="acc", columns="id", values="dlog", aggfunc="mean")
    iraw = Ix.pivot_table(index="acc", columns="id", values="value", aggfunc="mean")
    ctl = Cx.pivot_table(index="acc", columns="id", values="value", aggfunc="mean")
    vfs, rows = {}, []
    for v in V.itertuples():
        f = f"{PROJ}/panel/arch3/{v.chrom.lower()}/merged_231_{v.chrom.lower()}_final.vcf.gz"
        vf = vfs.setdefault(f, pysam.VariantFile(f))
        car, ref = carriers(vf, v.chrom, int(v.pos), int(v.ref_len), int(v.alt_len))
        out = dict(sym=v.sym, gene=v.gene, chrom=v.chrom, pos=v.pos, ref_len=v.ref_len,
                   alt_len=v.alt_len, src=v.src, mode=v.mode,
                   carriers_panel=len(car), ref_panel=len(ref))
        if v.gene not in gdl.columns:
            out["status"] = "gene not in 1001T"; rows.append(out); continue
        e = gdl[v.gene]; r = graw[v.gene]
        # summaries from RAW counts: d_log2 turns silent accessions into NaN, which dropped
        # every RLM3 deletion carrier from the medians
        cc = [a for a in r.index if a in car and np.isfinite(r[a])]
        rr = [a for a in r.index if a in ref and np.isfinite(r[a])]
        out.update(n_car_expr=len(cc), n_ref_expr=len(rr),
                   med_raw_car=r.reindex(cc).median(), med_raw_ref=r.reindex(rr).median(),
                   frac_car_silent=float((r.reindex(cc) < 1).mean()) if cc else np.nan,
                   frac_ref_silent=float((r.reindex(rr) < 1).mean()) if rr else np.nan,
                   fold=(r.reindex(cc).median() + 1) / (r.reindex(rr).median() + 1))
        out["p_gene"] = mw(r.reindex(cc).values, r.reindex(rr).values)
        # lineage control on the same split
        if len(cc) >= MIN_N and len(rr) >= MIN_N:
            pc = np.array([mw(ctl[c].reindex(cc).values, ctl[c].reindex(rr).values)
                           for c in ctl.columns])
            pc = pc[np.isfinite(pc)]
            out["ctl_frac_p05"] = float((pc < 0.05).mean())
            out["p_emp"] = float((np.sum(pc <= out["p_gene"]) + 1) / (len(pc) + 1)) \
                if np.isfinite(out["p_gene"]) else np.nan
        # isoform ratios
        iso = [c for c in iraw.columns if c.rsplit(".", 1)[0] == v.gene]
        if len(iso) >= 2 and len(cc) >= MIN_N:
            X = iraw[iso]; tot = X.sum(1)
            X = X[tot > tot.quantile(0.05)]
            R = logit((X.div(X.sum(1), axis=0)).clip(0.01, 0.99))
            ps = {c: mw(R[c].reindex(cc).values, R[c].reindex(rr).values) for c in iso}
            ps = {k: p for k, p in ps.items() if np.isfinite(p)}
            if ps:
                k = min(ps, key=ps.get)
                out.update(n_isoforms=len(iso), iso_best=k,
                           p_iso=min(1.0, ps[k] * len(iso)))
        out["status"] = "tested" if np.isfinite(out.get("p_gene", np.nan)) else \
            f"too few carriers with expression ({len(cc)})"
        rows.append(out)
    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/carrier_expression_test.csv", index=False)
    print(D.status.value_counts().to_string())


if __name__ == "__main__":
    main()
