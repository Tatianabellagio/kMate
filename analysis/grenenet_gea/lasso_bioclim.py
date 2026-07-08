#!/usr/bin/env python
"""Design A, step 2 -- multilocus LASSO/elastic-net analog of Zhou 2022 Fig. 3 (allelic
heterogeneity / incomplete LD), for founder ORIGIN CLIMATE (bio1-19).

The GRM/GBLUP pass (varexp_bioclim.py) answers the POLYGENIC question: does the SV / non-SNP
layer raise heritability & CV prediction beyond SNPs? Expected null in a selfer (K_snp~=K_all).
LASSO answers the SPARSE question the GRM cannot: can a multilocus model pick a SPECIFIC non-SNP
(indel/SV) marker that carries climate-adaptation signal SNPs cannot locally tag -- Zhou's
incomplete-LD / allelic-heterogeneity mechanism (their cis-LASSO that found causal SVs the MLM
missed). This is where a genuine pangenome win could still appear even when the aggregate GRM does
not move.

Design (per bioclim trait):
  * Feature arms, standardized founder x marker Z (0/1 haploid presence, uncalled->col freq p;
    MAC>=5, call>=0.9 -- IDENTICAL marker universe to class_grms.npz):
      - SNP-only     : LD-pruned SNP backbone (r^2<PRUNE within +-WINDOW_KB) -- the "linear tag set"
      - SNP+nonSNP   : same SNP backbone + ALL non-SNP (indel+SV) markers (the layer under test,
                       kept unpruned so an untagged SV cannot be pruned away)
      - SNP+SV       : SNP backbone + SV-only markers
  * ElasticNetCV (l1_ratio configurable; LASSO at l1_ratio=1), repeated k-fold, predictive R^2 per
    arm -> GAIN = R2(SNP+nonSNP) - R2(SNP).  Fairness: SNP backbone identical across arms, so any
    gain is the non-SNP layer, not marker count (cross-checked against the GRM snp_matched control).
  * For every non-SNP marker with a non-zero coefficient in the winning model: report its
    SNP-tagging status (best_r2 to any panel SNP within +-50kb, from nonsnp_tagging_chr*.npz) and
    gene annotation (lib.annotate). An untagged (r^2<0.2) selected marker = the incomplete-LD
    smoking gun.

Reuses build_class_grms._load_chrom / _classes for the exact marker filter. Env: kmate (sklearn).
Heavy (genome-wide markers x 19 traits x CV) -> run via sbatch (lasso_bioclim.sbatch), big mem.
Writes lasso_bioclim.csv (per trait x arm R^2 + gain) + lasso_bioclim_selected.csv (nonzero
non-SNP markers, coef, tagging, gene) + lasso_bioclim_summary.json.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
import build_class_grms as bcg           # _load_chrom, _classes, MIN_MAC, CALL_MIN, N

OUT = f"{lib.GEA}/varexp"
BIOCLIM_CSV = ("/global/scratch/users/tbellg/gea_grene-net/key_files/"
               "1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
BIOS = [f"bio{i}" for i in range(1, 20)]
CHROMS = [f"chr{i}" for i in range(1, 6)]


def build_Z(cl_want, prune_r2=None, window_kb=50, cap=None):
    """Streamed founder x marker standardized Z for the requested variant classes.

    cl_want: set of class codes {0 SNP,1 indel,2 SV}. prune_r2: if set, greedy LD-prune within
    +-window_kb keeping markers with max r^2 < prune_r2 (SNP backbone). Returns (Z[N,M], meta df).
    """
    N = bcg.N
    Zs, chrom_l, pos_l, cls_l = [], [], [], []
    for cl in CHROMS:
        rl, al, vp, vc, n_alt, n_cal = bcg._load_chrom(cl)
        klass = bcg._classes(rl, al)
        mac = np.minimum(n_alt, n_cal - n_alt)
        keep = (mac >= bcg.MIN_MAC) & (n_cal >= bcg.CALL_MIN * N) & np.isin(klass, list(cl_want))
        cols = np.where(keep)[0]
        meta = np.load(f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.meta.npz",
                       allow_pickle=True)
        pos = meta["pos"].astype(int)
        # standardized Z for the kept cols
        g = vp[:, cols].toarray().astype(np.float64)
        called = vc[:, cols].toarray().astype(bool)
        p = np.where(n_cal[cols] > 0, n_alt[cols] / np.maximum(n_cal[cols], 1), 0.0)
        miss = ~called
        if miss.any():
            g[miss] = np.take(p, np.where(miss)[1])
        sd = np.sqrt(np.maximum(p * (1 - p), 1e-6))
        Z = ((g - p) / sd).astype(np.float32)                       # N x len(cols)
        cpos = pos[cols]
        if prune_r2 is not None:
            kmask = _ld_prune(Z, cpos, prune_r2, window_kb)
            Z = Z[:, kmask]; cpos = cpos[kmask]; cols = cols[kmask]
        Zs.append(Z); pos_l.append(cpos)
        chrom_l.append(np.full(len(cpos), cl)); cls_l.append(klass[cols])
        print(f"  {cl}: kept {Z.shape[1]} markers (classes {sorted(cl_want)})", flush=True)
    Z = np.concatenate(Zs, axis=1)
    m = pd.DataFrame(dict(chrom=np.concatenate(chrom_l), pos=np.concatenate(pos_l),
                          cls=np.concatenate(cls_l)))
    if cap is not None and Z.shape[1] > cap:                        # keep highest-variance markers
        v = Z.var(0); idx = np.argsort(v)[::-1][:cap]
        Z = Z[:, idx]; m = m.iloc[idx].reset_index(drop=True)
    return Z, m


def _ld_prune(Z, pos, r2_thresh, window_kb):
    """Greedy position-ordered LD prune: keep a marker iff its max r^2 with already-kept markers
    within +-window_kb is < r2_thresh. Z columns assumed standardized (r = corr = ZiZj/(N-1))."""
    N = Z.shape[0]; order = np.argsort(pos, kind="stable")
    win = window_kb * 1000
    kept_idx, kept_pos, kept_Z = [], [], []
    Zn = Z / np.sqrt(np.maximum((Z ** 2).sum(0), 1e-12))            # unit-norm cols -> dot = r
    for j in order:
        pj = pos[j]
        if kept_pos:
            lo = np.searchsorted(kept_pos, pj - win)
            if lo < len(kept_pos):
                block = np.stack(kept_Z[lo:], axis=1)                # N x nrecent
                r = Zn[:, j] @ block
                if (r * r).max() >= r2_thresh:
                    continue
        kept_idx.append(j); kept_pos.append(pj); kept_Z.append(Zn[:, j])
    mask = np.zeros(Z.shape[1], bool); mask[kept_idx] = True
    return mask


def load_tagging():
    """Per non-SNP marker best_r2 to any SNP within +-50kb, keyed (chrom,pos)."""
    d = {}
    for i, cl in enumerate(CHROMS, 1):
        t = np.load(f"{OUT}/nonsnp_tagging_chr{i}.npz", allow_pickle=True)
        for pos, r2 in zip(t["pos"].astype(int), t["best_r2"].astype(float)):
            d[(cl, int(pos))] = float(r2)
    return d


def cv_r2(Z, y, l1_ratio, nfold, nrep, seed):
    """Repeated k-fold predictive R^2 of ElasticNetCV (alpha chosen by inner CV per fold)."""
    N = len(y); preds = np.zeros(N); cnt = np.zeros(N)
    for rep in range(nrep):
        kf = KFold(nfold, shuffle=True, random_state=seed + rep)
        for tr, te in kf.split(Z):
            en = ElasticNetCV(l1_ratio=l1_ratio, n_alphas=50, cv=4,
                              max_iter=5000, tol=1e-3, n_jobs=-1, selection="random",
                              random_state=seed)
            en.fit(Z[tr], y[tr])
            preds[te] += en.predict(Z[te]); cnt[te] += 1
    yhat = preds / np.maximum(cnt, 1)
    return 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prune-r2", type=float, default=0.9)
    ap.add_argument("--window-kb", type=int, default=50)
    ap.add_argument("--snp-cap", type=int, default=150000)
    ap.add_argument("--l1-ratio", type=float, default=1.0)      # 1.0 = pure LASSO (Zhou)
    ap.add_argument("--nfold", type=int, default=5)
    ap.add_argument("--nrep", type=int, default=5)
    ap.add_argument("--traits", nargs="*", default=BIOS)
    args = ap.parse_args()

    print("building SNP backbone (LD-pruned) ...", flush=True)
    Z_snp, m_snp = build_Z({0}, prune_r2=args.prune_r2, window_kb=args.window_kb, cap=args.snp_cap)
    print("building non-SNP (indel+SV, unpruned) ...", flush=True)
    Z_ns, m_ns = build_Z({1, 2})
    is_sv = (m_ns["cls"].values == 2)
    print(f"SNP backbone {Z_snp.shape[1]} | non-SNP {Z_ns.shape[1]} "
          f"(SV {int(is_sv.sum())}, indel {int((~is_sv).sum())})", flush=True)

    Z_snp_sv = np.concatenate([Z_snp, Z_ns[:, is_sv]], axis=1)
    Z_snp_ns = np.concatenate([Z_snp, Z_ns], axis=1)
    ns_off = Z_snp.shape[1]                                     # non-SNP feature offset in SNP+nonSNP

    bio = pd.read_csv(BIOCLIM_CSV).set_index("ecotypeid")
    bio.index = bio.index.astype(str)
    founders = np.load(f"{lib.PROJ}/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz",
                       allow_pickle=True)["founders"].astype(str)
    assert list(founders) == list(bio.loc[founders].index)      # order sanity
    tags = load_tagging()

    rows, selrows = [], []
    for b in args.traits:
        y = bio.loc[founders, b].astype(float).values
        y = (y - y.mean()) / y.std()
        r2_snp = cv_r2(Z_snp, y, args.l1_ratio, args.nfold, args.nrep, 0)
        r2_ssv = cv_r2(Z_snp_sv, y, args.l1_ratio, args.nfold, args.nrep, 0)
        r2_snn = cv_r2(Z_snp_ns, y, args.l1_ratio, args.nfold, args.nrep, 0)
        rows += [dict(trait=b, arm="snp", r2=r2_snp),
                 dict(trait=b, arm="snp+sv", r2=r2_ssv),
                 dict(trait=b, arm="snp+nonsnp", r2=r2_snn),
                 dict(trait=b, arm="gain_sv", r2=r2_ssv - r2_snp),
                 dict(trait=b, arm="gain_nonsnp", r2=r2_snn - r2_snp)]
        # which non-SNP markers does the full model select? (fit once on all data)
        en = ElasticNetCV(l1_ratio=args.l1_ratio, n_alphas=50, cv=5, max_iter=8000,
                          tol=1e-3, n_jobs=-1, selection="random", random_state=0).fit(Z_snp_ns, y)
        coef_ns = en.coef_[ns_off:]
        nz = np.where(np.abs(coef_ns) > 1e-8)[0]
        for j in nz:
            chrom, pos, cls = m_ns.iloc[j][["chrom", "pos", "cls"]]
            r2tag = tags.get((chrom, int(pos)), np.nan)
            selrows.append(dict(trait=b, chrom=chrom, pos=int(pos),
                                cls={1: "indel", 2: "sv"}[int(cls)], coef=float(coef_ns[j]),
                                snp_tag_r2=r2tag, untagged=bool(r2tag < 0.2) if r2tag == r2tag else None))
        n_unt = sum(1 for r in selrows if r["trait"] == b and r["untagged"])
        print(f"  {b:6s} r2_snp={r2_snp:+.3f} gain_sv={r2_ssv-r2_snp:+.3f} "
              f"gain_nonsnp={r2_snn-r2_snp:+.3f} | nonSNP selected={len(nz)} untagged={n_unt}",
              flush=True)

    pd.DataFrame(rows).to_csv(f"{OUT}/lasso_bioclim.csv", index=False)
    sel = pd.DataFrame(selrows)
    sel.to_csv(f"{OUT}/lasso_bioclim_selected.csv", index=False)
    g = pd.DataFrame(rows)
    gsv = g[g.arm == "gain_sv"].r2.values; gns = g[g.arm == "gain_nonsnp"].r2.values
    summary = dict(
        n_snp_backbone=int(Z_snp.shape[1]), n_nonsnp=int(Z_ns.shape[1]), n_sv=int(is_sv.sum()),
        prune_r2=args.prune_r2, l1_ratio=args.l1_ratio,
        gain_sv={"median": float(np.median(gsv)), "max": float(np.max(gsv)),
                 "n_pos": int((gsv > 0).sum()), "n": len(gsv)},
        gain_nonsnp={"median": float(np.median(gns)), "max": float(np.max(gns)),
                     "n_pos": int((gns > 0).sum()), "n": len(gns)},
        n_nonsnp_selected=int(len(sel)),
        n_untagged_selected=int(sel["untagged"].sum()) if len(sel) else 0)
    with open(f"{OUT}/lasso_bioclim_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\n=== LASSO SUMMARY ===")
    print(json.dumps(summary, indent=2))
    print(f"wrote {OUT}/lasso_bioclim.csv + _selected.csv + _summary.json")


if __name__ == "__main__":
    main()
