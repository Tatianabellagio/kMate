#!/usr/bin/env python
"""Are SVs enriched in the temporally-SELECTED clq0.9 haploblocks at one site? (site 4)

Selection is called at the r2>=0.9 haploblock level (blocks_mcf90/chrN_clq0.9_blocks),
by EMPIRICAL OUTLIER on the per-variant GLOBAL-mode temporal selection coefficients:

  block score  = |median_{SNP in block} s|      (SNP-based, so the block's SELECTION is
                 defined by the haplotype-tagging SNPs; SVs are NOT used to call selection
                 -> the enrichment test is not circular)
  SELECTED     = block score in the top --top-frac genome-wide tail (drift = the bulk)

Then: are SELECTED blocks more likely to CARRY an SV than non-selected blocks, controlling
for block size (bigger blocks trivially carry more variants)? Two size-aware tests:
  (1) size-matched null: for each selected block draw non-selected blocks of matched
      variant count; compare the SV-bearing fraction -> fold + empirical p.
  (2) logistic: selected ~ has_SV + log(n_var); the has_SV coefficient + Wald p.

PART A (per-variant, class-level, p0-matched): SVs start rarer (lower founding p0) so a
  rare variant makes a bigger logit swing; we bin by p0 and compare per-class median |s|
  within bins -> a fair "do SVs move faster per generation" readout, independent of blocks.

INPUT: site_temporal/site{S}_scoef_{snp,smallindel,sv}.npz  (build: site_variant_temporal_scoef.py)
       blocks_mcf90/chr{1..5}_clq0.9_blocks_clq0.9.tsv
OUTPUT (--out, default site_temporal):
       site{S}_clq90_blocks.parquet     per-block table (composition + score + selected)
       site{S}_sv_enrichment.json       summary numbers
       site{S}_sv_enrichment.png        figure

CAVEAT: GLOBAL-mode AF is projected through the intact founder panel, so variants in
perfect founder-LD share one trajectory. Scoring selection from a block's SNPs and asking
whether the block also carries an SV = "do selected founder haplotypes carry SVs", which is
the honest, LD-safe form of the question (not per-SV independent selection).

FOUNDER-PANEL MAC FLOOR (--min-mac): the raw per-sample/pool variant catalog behind
site_variant_temporal_scoef.py carries no founder-panel MAF floor of its own (only the
p0 pool-reachability check below, which is noisy and does not reliably track founder-level
rarity) -- genome-wide, 51% of its SV-class calls are literal founder singletons (MAC<=1),
vs the founder-panel MAF>=0.05 floor already used to define the LD blocks and sv_landscape.
Composition counts (n_snp/n_indel/n_sv/has_sv) and the SNPs used for block_score are both
floored to lib.founder_panel_keep() before scoring, so neither side of the enrichment test
runs on singleton noise. See KINSHIP_TEMPORAL_METHODS_RESEARCH.md context.

  PY=<plotting env>; $PY site_sv_enrichment.py --site 4 --top-frac 0.01 --min-mac 2
"""
from __future__ import annotations
import argparse, glob, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from site_variant_temporal_scoef import site_scoef, STORE

CLASSES = ("snp", "smallindel", "sv")
SV_MIN_BP = 50
BD = f"{lib.GEA}/blocks_mcf90"


def load_variants(site, sdir, min_p0=0.02, min_mac=2):
    """Per-variant temporal s + class + reachability. Uses cached per-class npz if
    present (site 4), else computes in-memory via site_scoef (no giant npz written).

    min_mac floors every variant (SNP too, not just SV/indel) to lib.founder_panel_keep
    -- founder-panel MAC>=min_mac & called-frac>=0.9 -- on top of the pool-level p0
    reachability check, so block_score and composition are scored from the same
    non-singleton set."""
    snp_npz = f"{sdir}/site{site}_scoef_snp.npz"
    if os.path.exists(snp_npz):
        frames = []
        for c in CLASSES:
            z = np.load(f"{sdir}/site{site}_scoef_{c}.npz", allow_pickle=True)
            k = z["keep"]
            frames.append(pd.DataFrame(dict(
                chrom=z["chrom"][k].astype(str), pos=z["pos"][k].astype(np.int64),
                s=z["s"][k].astype(float), p0=z["p0"][k].astype(float),
                size=z["size"][k].astype(np.int64), cls=c)))
        V = pd.concat(frames, ignore_index=True)
    else:
        idx_snp = np.load(f"{STORE}/index_snp.npz")
        idx_non = np.load(f"{STORE}/index_nonsnp.npz")
        Rsnp = site_scoef(site, "snp"); Rnon = site_scoef(site, "nonsnp")
        dlen = np.abs(idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64))
        frames = []
        specs = [("snp", idx_snp, Rsnp, np.ones(Rsnp["s"].shape[0], bool)),
                 ("smallindel", idx_non, Rnon, dlen <= SV_MIN_BP),
                 ("sv", idx_non, Rnon, dlen > SV_MIN_BP)]
        for name, idx, R, sel in specs:
            frames.append(pd.DataFrame(dict(
                chrom=idx["chrom"].astype("U5")[sel], pos=idx["pos"][sel].astype(np.int64),
                s=R["s"][sel].astype(float), p0=R["p0"][sel].astype(float),
                size=np.abs(idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64))[sel],
                cls=name)))
        V = pd.concat(frames, ignore_index=True)
    # reachability + finite filter + founder-panel MAC floor (no singletons)
    reach = (V.p0 >= min_p0) & (V.p0 <= 1 - min_p0) & np.isfinite(V.s)
    fp = lib.founder_panel_keep(V.chrom.to_numpy(), V.pos.to_numpy(), min_mac=min_mac)
    keep = (reach & fp).to_numpy()
    print(f"  reachable {int(reach.sum()):,} -> founder-panel MAC>={min_mac} keeps "
          f"{int(keep.sum()):,} ({', '.join(f'{c} {int((keep & (V.cls==c).to_numpy()).sum()):,}' for c in CLASSES)})")
    V = V[keep].reset_index(drop=True)
    V["abs_s"] = V.s.abs()
    return V


def load_blocks():
    b = pd.concat([pd.read_csv(f, sep="\t")
                   for f in sorted(glob.glob(f"{BD}/chr*_clq0.9_blocks_clq0.9.tsv"))])
    return b.reset_index(drop=True)


def assign_block(V, blocks):
    bid = np.full(len(V), -1, np.int64)
    for c, g in blocks.groupby("chrom"):
        g = g.sort_values("start_pos")
        st = g.start_pos.to_numpy(); en = g.end_pos.to_numpy(); idx = g.index.to_numpy()
        vi = np.where(V.chrom.values == c)[0]
        p = V.pos.values[vi]
        j = np.searchsorted(st, p, "right") - 1
        ok = (j >= 0)
        inb = ok & (p <= en[np.clip(j, 0, len(en) - 1)])
        bid[vi[inb]] = idx[j[inb]]
    return bid


def part_A(V):
    q = np.quantile(V.p0, np.linspace(0, 1, 11)); q[-1] += 1e-6
    V["p0_bin"] = np.clip(np.digitize(V.p0, q[1:-1]), 0, 9)
    med = V.groupby(["p0_bin", "cls"]).abs_s.median().unstack()
    binmid = V.groupby("p0_bin").p0.median()
    wsv = V[V.cls == "sv"].groupby("p0_bin").size().reindex(range(10), fill_value=0)
    ratio = (med["sv"] / med["snp"]).reindex(range(10))
    matched = float(np.nansum(ratio * wsv) / wsv[np.isfinite(ratio)].sum())
    raw = float(V[V.cls == "sv"].abs_s.median() / V[V.cls == "snp"].abs_s.median())
    return med, binmid, raw, matched


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=int, default=4)
    ap.add_argument("--out", default=f"{lib.GEA}/site_temporal")
    ap.add_argument("--top-frac", type=float, default=0.01, help="selected = top tail of block score")
    ap.add_argument("--min-snp", type=int, default=3, help="min SNPs to score a block")
    ap.add_argument("--min-mac", type=int, default=12,
                    help="founder-panel MAC floor for every variant (12 = MAF>=5% common SVs, "
                         "matching build_sv_landscape.py / the LD blocks; decision 2026-07-01 "
                         "to rely on common SVs since low-freq SV calls are more likely spurious)")
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    sdir = args.out; os.makedirs(args.out, exist_ok=True)
    V = load_variants(args.site, sdir, min_mac=args.min_mac)
    print(f"reachable variants: {len(V):,} | " +
          " ".join(f"{c}={int((V.cls==c).sum()):,}" for c in CLASSES))

    med, binmid, rawA, matchedA = part_A(V)
    print(f"[A] |s| SV/SNP  raw={rawA:.3f}  p0-matched={matchedA:.3f}")

    # ---- block assignment + per-block table -----------------------------------
    blocks = load_blocks()
    V["bid"] = assign_block(V, blocks)
    inb = V[V.bid >= 0]
    print(f"[B] variants in a clq0.9 block: {len(inb):,} ({len(inb)/len(V):.1%}); "
          f"rest fall in inter-block (low-LD) gaps")
    g = inb.groupby("bid")
    comp = g.cls.value_counts().unstack(fill_value=0)
    for c in CLASSES:
        if c not in comp:
            comp[c] = 0
    tab = pd.DataFrame(index=comp.index)
    tab["n_snp"] = comp["snp"]; tab["n_smallindel"] = comp["smallindel"]; tab["n_sv"] = comp["sv"]
    tab["n_var"] = tab[["n_snp", "n_smallindel", "n_sv"]].sum(1)
    # SV-INDEPENDENT size covariate for matching: non-SV variant count = "opportunity
    # to contain an SV". Matching on n_var would partly match on SV presence itself.
    tab["n_nonsv"] = tab["n_snp"] + tab["n_smallindel"]
    tab["has_sv"] = tab.n_sv > 0
    # SNP-based block selection score (signed median SNP s, and its magnitude)
    snp_s = inb[inb.cls == "snp"].groupby("bid").s.median()
    tab["block_s"] = snp_s.reindex(tab.index)
    tab["block_score"] = tab["block_s"].abs()
    tab = tab.join(blocks[["chrom", "start_pos", "end_pos"]], on="bid")

    scored = tab[tab.n_snp >= args.min_snp].copy()
    thr = scored.block_score.quantile(1 - args.top_frac)
    scored["selected"] = scored.block_score >= thr
    tab["selected"] = tab.index.isin(scored.index[scored.selected])
    nsel = int(scored.selected.sum())
    print(f"[B] scored blocks (>={args.min_snp} SNP): {len(scored):,} | "
          f"selected (top {args.top_frac:.1%}, |median SNP s|>={thr:.3f}): {nsel}")

    # ---- enrichment: has_SV in selected vs non-selected -----------------------
    sel = scored[scored.selected]; non = scored[~scored.selected]
    f_sel = float(sel.has_sv.mean()); f_non = float(non.has_sv.mean())
    # (1) size-matched null: match on n_nonsv (SV-independent size), nearest by count
    rng = np.random.default_rng(args.seed)
    non_by_size = {}
    for nv, sub in non.groupby("n_nonsv"):
        non_by_size[nv] = sub.index.to_numpy()
    non_sizes = np.array(sorted(non_by_size))
    null_frac = np.empty(args.n_perm)
    sel_sizes = sel.n_nonsv.to_numpy()
    for b in range(args.n_perm):
        picks = []
        for nv in sel_sizes:
            k = non_sizes[np.argmin(np.abs(non_sizes - nv))]
            picks.append(rng.choice(non_by_size[k]))
        null_frac[b] = non.loc[picks].has_sv.mean()
    p_emp = float((np.sum(null_frac >= f_sel) + 1) / (args.n_perm + 1))
    fold = f_sel / np.median(null_frac)
    # matcher balance: mean non-SV size of selected vs a matched draw (should agree)
    picks0 = [rng.choice(non_by_size[non_sizes[np.argmin(np.abs(non_sizes - nv))]])
              for nv in sel_sizes]
    bal_sel = float(sel.n_nonsv.mean()); bal_null = float(non.loc[picks0].n_nonsv.mean())
    print(f"[B] size balance (n_nonsv): selected mean={bal_sel:.1f}  matched-null mean={bal_null:.1f}")
    # (2) raw (size-UNcontrolled) 2x2 Fisher, for reference vs the matched null
    from scipy.stats import fisher_exact
    a = int(sel.has_sv.sum()); b = int((~sel.has_sv).sum())
    c = int(non.has_sv.sum()); d = int((~non.has_sv).sum())
    or_raw, p_raw = fisher_exact([[a, b], [c, d]])
    print(f"[B] has_SV  selected={f_sel:.3f}  matched-null(med)={np.median(null_frac):.3f}  "
          f"non-sel={f_non:.3f}  fold={fold:.2f}  p_emp={p_emp:.4f}")
    print(f"[B] raw 2x2 (size-uncontrolled) OR={or_raw:.2f}  p={p_raw:.3g}  "
          f"(selected {a}/{a+b} SV-bearing vs non-sel {c}/{c+d})")

    # ---- outputs --------------------------------------------------------------
    # CSV (not parquet: the kmate env used by the sbatch array has no pyarrow)
    tab.reset_index().to_csv(f"{args.out}/site{args.site}_clq90_blocks.csv.gz",
                             index=False, compression="gzip")
    summ = dict(site=args.site, block_set="blocks_mcf90 clq0.9", n_blocks=int(len(tab)),
                min_mac=args.min_mac, n_reachable_variants=int(len(V)),
                frac_variants_in_block=float(len(inb)/len(V)),
                A_raw_ratio_sv_snp=rawA, A_p0matched_ratio_sv_snp=matchedA,
                A_median_abs_s={c: float(V[V.cls == c].abs_s.median()) for c in CLASSES},
                B_top_frac=args.top_frac, B_min_snp=args.min_snp, B_score_thr=float(thr),
                B_n_scored=int(len(scored)), B_n_selected=nsel,
                B_hasSV_selected=f_sel, B_hasSV_null_median=float(np.median(null_frac)),
                B_hasSV_nonselected=f_non, B_fold=float(fold), B_p_emp=p_emp,
                B_raw_OR_hasSV=float(or_raw), B_raw_fisher_p=float(p_raw))
    json.dump(summ, open(f"{args.out}/site{args.site}_sv_enrichment.json", "w"), indent=2)

    # top selected SV-bearing blocks, gene-annotated
    topsv = sel[sel.has_sv].sort_values("block_score", ascending=False).head(30).copy()
    if len(topsv):
        ann = lib.annotate_svs(topsv.rename(columns={"start_pos": "pos"}).assign(
            ref_len=1), genes=None)
        topsv["genes"] = ann["genes_all"].values
        topsv.to_csv(f"{args.out}/site{args.site}_top_selected_sv_blocks.csv")

    # ---- plot data (rendered separately; matplotlib is broken in this env) -----
    scored2 = scored.copy()
    dec = pd.qcut(scored2.block_score, 10, labels=False, duplicates="drop")
    sv_by_dec = scored2.groupby(dec).has_sv.mean()
    np.savez(f"{args.out}/site{args.site}_sv_enrichment_plotdata.npz",
             A_binmid=binmid.values, A_med_snp=med["snp"].values,
             A_med_smallindel=med["smallindel"].values, A_med_sv=med["sv"].values,
             B_decile=np.arange(len(sv_by_dec)), B_sv_bearing=sv_by_dec.values,
             C_null_frac=null_frac, C_f_sel=f_sel, C_f_non=f_non,
             fold=fold, p_emp=p_emp)
    print(f"-> {args.out}/site{args.site}_sv_enrichment.json + _clq90_blocks.parquet "
          f"+ _plotdata.npz  (render figure separately)")


if __name__ == "__main__":
    main()
