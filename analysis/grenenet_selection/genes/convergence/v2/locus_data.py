#!/usr/bin/env python
"""v2 locus data -- everything the locus figures need for one candidate gene, extracted once into a small folder.

For each gene (AGI id) of the NEW list: results/candidate_list_v2_events_cr50/locus/<GENE>/
  units.csv          the gene's units (v2_unit_evidence rows) + their GEA axes (sig_axes) from gea_pool
  zoom_<axis>.csv    tile_cr50 LFMM p in the window for every unit/SNP: class (snp | nonsnp), pos, ref_len, alt_len, MAF, pval, unit
  genes_window.csv   gene models in the window: gene, strand, feature (gene / CDS / exon / five_prime_UTR / three_prime_UTR), start, end
  te_window.csv      reference TEs in the window
  atac_window.csv    multi-tissue ATAC peaks in the window (start, end, n_tissues, flower/leaf/root/shoot)
  dp_pools.csv       per pool: site, plot, the climate axes of the unit's GEA, and the unit's delta p (gen 9 frequency - founding frequency)
  expr_<unit>.csv    per accession with RNA-seq: expression of the target gene (raw, d_log2) and carrier status of the unit (carrier / noncarrier / uncalled)
Window = gene +-WIN_KB (default 30 kb), extended to cover every unit. env: kmate. Compute node.
usage: python locus_data.py AT3G04610 AT2G22540
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
GEA = os.path.abspath(os.path.join(CONV, "..", ".."))
sys.path.insert(0, HERE); sys.path.insert(0, CONV); sys.path.insert(0, GEA); sys.path.insert(0, os.path.join(GEA, "genes", "dissection"))
import units as UN                                              # noqa: E402
import atac_overlap as AO                                       # noqa: E402
import screen_sig_blocks as SB                                  # noqa: E402
V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")
TILE = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/tile_cr50"
CM = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices_cr50"
EQ = os.path.abspath(os.path.join(GEA, "..", "..", "data", "eqtl_v2"))
WIN_KB = 30


def main(genes):
    U = pd.read_csv(f"{V2}/v2_unit_evidence.csv", low_memory=False)
    P = pd.read_csv(f"{V2}/gea_pool.csv")[["unit", "sig_axes", "best_axis", "MAF"]]
    pools = pd.read_csv(f"{CM}/gen9.pools.csv")
    units_csv = pd.read_csv(f"{GEA}/genes/convergence/event_gea/cr50/nonsnp_gen9_units.csv", usecols=["unit"])
    col_of = pd.Series(np.arange(len(units_csv)), index=units_csv.unit.values)
    dp = np.load(f"{GEA}/genes/convergence/event_gea/cr50/nonsnp_gen9_dp.npy", mmap_mode="r")
    G, sub, TE = SB.load_gff(); peaks = AO.load_peaks()
    Gx = pd.read_csv(f"{EQ}/expr_genes.csv"); graw = Gx.pivot_table(index="acc", columns="id", values="raw", aggfunc="mean")
    gdl = Gx.pivot_table(index="acc", columns="id", values="dlog", aggfunc="mean")
    for g in genes:
        out = f"{V2}/locus/{g}"; os.makedirs(out, exist_ok=True)
        X = U[U.target_gene == g].merge(P, on="unit", how="left")
        if X.empty:
            print(f"{g}: no units in v2_unit_evidence"); continue
        X.to_csv(f"{out}/units.csv", index=False)
        ge = G[G.gene == g].iloc[0]; ch = ge.chrom
        lo = int(min(ge.start, X.pos.min())) - WIN_KB * 1000; hi = int(max(ge.end, X.pos.max() + X.ref_len.max())) + WIN_KB * 1000
        axes = sorted({a for s in X.sig_axes.dropna() for a in str(s).split(",") if a} | set(X.best_axis.dropna()))
        for ax in axes:
            rows = []
            for cls, f in (("snp", f"{TILE}/lfmm_snp_gen9_{ax}.csv"), ("nonsnp", f"{TILE}/lfmm_nonsnp_gen9_{ax}.csv")):
                d = pd.read_csv(f, usecols=lambda c: c in ("chrom", "pos", "ref_len", "alt_len", "MAF", "pval", "unit"))
                d = d[(d.chrom == ch) & d.pos.between(lo, hi)].copy(); d["cls"] = cls; d["n_tested"] = int(pd.read_csv(f, usecols=["pos"]).shape[0]); rows.append(d)
            pd.concat(rows, ignore_index=True).to_csv(f"{out}/zoom_{ax}.csv", index=False)
        # gene models in the window
        gw = G[(G.chrom == ch) & (G.end >= lo) & (G.start <= hi)]
        recs = [dict(gene=r.gene, strand=r.strand, feature="gene", start=r.start, end=r.end) for r in gw.itertuples()]
        for k, nm in (("cds", "CDS"), ("exon", "exon"), ("u5", "five_prime_UTR"), ("u3", "three_prime_UTR")):
            s = sub[k]; s = s[(s.chrom == ch) & s.gene.isin(gw.gene) & (s.end >= lo) & (s.start <= hi)]
            recs += [dict(gene=r.gene, strand=gw.set_index("gene").strand[r.gene], feature=nm, start=r.start, end=r.end) for r in s.itertuples()]
        pd.DataFrame(recs).to_csv(f"{out}/genes_window.csv", index=False)
        t = TE[(TE.chrom == ch) & (TE.end >= lo) & (TE.start <= hi)]; t.to_csv(f"{out}/te_window.csv", index=False)
        ps, pe, num, tm = peaks[ch]; m = (pe >= lo) & (ps <= hi)
        pd.DataFrame({"start": ps[m], "end": pe[m], "n_tissues": num[m], **{n: tm[m, i] for i, n in enumerate(AO.TISSUES)}}).to_csv(f"{out}/atac_window.csv", index=False)
        # per-pool delta p of each unit + climate of the unit's axes
        D = pools[["pool", "site", "plot"] + sorted(set(axes))].copy()
        for u in X.unit:
            D[f"dp_unit{u}"] = np.asarray(dp[:, int(col_of[u])], float)
        D.to_csv(f"{out}/dp_pools.csv", index=False)
        # expression of the target gene by carrier status
        founders, car, cal = UN.carrier_matrix(X.unit.values)
        fid = np.array([int(x) if str(x).isdigit() else -1 for x in founders])
        for j, u in enumerate(X.unit):
            rows = []
            for a in graw.index:
                k = np.where(fid == a)[0]
                st = "uncalled" if not len(k) or not cal[k[0], j] else ("carrier" if car[k[0], j] else "noncarrier")
                rows.append(dict(acc=a, status=st, raw=graw[g].get(a, np.nan) if g in graw.columns else np.nan, dlog=gdl[g].get(a, np.nan) if g in gdl.columns else np.nan))
            pd.DataFrame(rows).to_csv(f"{out}/expr_{u}.csv", index=False)
        print(f"{g}: {len(X)} units, window {ch}:{lo:,}-{hi:,}, axes {axes} -> {out}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
