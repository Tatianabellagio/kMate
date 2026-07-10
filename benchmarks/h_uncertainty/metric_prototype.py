#!/usr/bin/env python3
"""Two clean per-region certainty metrics, validated against known targets.

Motivated by the user's decomposition:
  (A) RESOLVABILITY — "some regions are just hard from PANEL complexity"
      (founder collinearity / how many distinguishing k-mers exist). A property of
      the panel K alone, a priori, INDEPENDENT of coverage. It is the floor that
      persists at infinite coverage (the global 0.0485 bias; per-window err_noiseless).
  (B) COVERAGE — "given the k-mers we actually realized, how much EXTRA error".
      A variance that shrinks with depth. Target = coverage_excess = err_cov - err_noiseless.

The old per-window eff_rank (window_certainty.py) was a WEAK resolvability predictor
(r=-0.09) because it was computed on the COUNT-WEIGHTED, h-conditioned Fisher J — i.e.
contaminated by coverage. Here we build PANEL-ONLY statistics (no counts, no h) and ask
which best predicts err_noiseless, and a clean COVERAGE statistic for the variance part.
Both are reported in BOTH modes (per-window track + whole-chromosome global).

Targets are reused from results/window_certainty.tsv (err_noiseless, err_cov*, nzk*,
old eff_rank/fisherSE) and results/convergence_diag.tsv (global floor 0.0485). We only
recompute deterministic panel quantities here — no EM reruns.

Run via sbatch (loads the ~2.5 GB dense K). p80 Chr1, 10 kb windows (matched to the tsv).
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, str(ROOT / "src"))
from kmate.block_em import define_windows, assign_kmers_to_blocks   # noqa
from kmate.h_uncertainty import fisher_information_h, _tangent_pinv_on_support  # noqa
PRE = ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2"
HU  = ROOT / "benchmarks/h_uncertainty"
COVS = [3.0, 10.0, 30.0]
WINDOW_BP = 10000


def load_all():
    kp = sparse.load_npz(PRE / "kmer_pa_Chr1.kmer_pa.npz")
    K = np.asarray(kp.todense(), dtype=np.float32) if sparse.issparse(kp) else np.asarray(kp, np.float32)
    m = np.load(PRE / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
    fnd = np.asarray(m["founders"]).astype(str); bid = np.asarray(m["bubble_id"])
    omega = (1.0 / np.bincount(bid)[bid]).astype(np.float32)
    meta = dict(bid=bid, chrom=np.asarray(m["bubble_chrom"]).astype(str),
                start=np.asarray(m["bubble_start"]), end=np.asarray(m["bubble_end"]))
    return K, fnd, omega, meta


def truth(fnd, pool="cov10_n50_g0_s42_hotspots_p80_chr1"):
    w = pd.read_csv(ROOT / f"benchmarks/p80/sims/{pool}/pool_weights.tsv", sep="\t")
    wm = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wm.get(str(f), 0.0) for f in fnd], dtype=np.float64); return h / h.sum()


def eff_rank_gram(G):
    """Effective rank = exp(spectral entropy) of a (small) PSD Gram matrix."""
    ev = np.clip(np.linalg.eigvalsh(G), 0, None)
    pos = ev[ev > ev.max() * 1e-12] if ev.size and ev.max() > 0 else ev[:0]
    if pos.size == 0:
        return 0.0
    p = pos / pos.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def n_distinct_signatures(Ksub):
    """# of distinct binary founder rows (founders with identical presence collapse)."""
    if Ksub.shape[0] == 0:
        return 0
    return int(np.unique(Ksub > 0, axis=0).shape[0])


def main():
    t0 = time.time()
    K, fnd, omega, meta = load_all()
    F, Kn = K.shape
    h_true = truth(fnd); supp = np.flatnonzero(h_true > 0)
    mu = (h_true @ K).astype(np.float64)                       # noiseless rate per k-mer
    print(f"p80 F={F} K={Kn:,} support={supp.size}  load {time.time()-t0:.0f}s", flush=True)

    # ---- window assignment, IDENTICAL to window_certainty.py ----
    blocks = define_windows(meta["chrom"], meta["start"], meta["end"], window_bp=WINDOW_BP)
    kb = assign_kmers_to_blocks(meta["bid"], meta["chrom"], meta["start"], meta["end"], blocks)
    order = np.argsort(kb, kind="stable"); kbs = kb[order]
    st = np.searchsorted(kbs, np.arange(len(blocks)), "left")
    en = np.searchsorted(kbs, np.arange(len(blocks)), "right")
    print(f"{len(blocks)} windows  ({time.time()-t0:.0f}s)", flush=True)

    Ksupp = K[supp]                                            # present founders only
    rows = []
    for wi in range(len(blocks)):
        idx = order[st[wi]:en[wi]]; nk = idx.size
        if nk == 0:
            continue
        Kw = K[:, idx]; Kws = Ksupp[:, idx]; omw = omega[idx]; muw = mu[idx]
        csa = Kw.sum(0); css = Kws.sum(0)
        rec = dict(window=wi, n_kmers=nk)
        # ---- RESOLVABILITY (panel-only, no counts/h) ----
        # supply: # distinguishing k-mers (vary across founders)
        rec["ndist_apriori"] = int(((csa > 0) & (csa < F)).sum())
        rec["ndist_cond"]    = int(((css > 0) & (css < supp.size)).sum())
        # identity-resolvable founders: distinct presence signatures
        rec["nsig_apriori"]  = n_distinct_signatures(Kw)
        rec["nsig_cond"]     = n_distinct_signatures(Kws)
        # spectral: eff_rank of the omega-weighted PANEL Gram  G = (K*w) Kᵀ  (no counts, no h)
        rec["gram_eff_apriori"] = eff_rank_gram((Kw * omw) @ Kw.T)
        rec["gram_eff_cond"]    = eff_rank_gram((Kws * omw) @ Kws.T)
        # founders with at least one present k-mer here (panel analog of n_local_truefounders)
        rec["nfnd_present"] = int((Kws.sum(1) > 0).sum())
        # CRLB: predicted sqrt(total variance) = sqrt(trace of tangent-pinv of the
        # NOISELESS Fisher J = K diag(ω/μ) Kᵀ on the present support. This is the
        # PRINCIPLED identifiability->error predictor (eff_rank/cond are crude proxies).
        # Counts=muw => k-mers absent under h_true (μ=0) drop out, matching the EM.
        Jnl = fisher_information_h(h_true, Kw, muw, omega=omw)
        spres = supp[(Kws.sum(1) > 0)]                        # present founders seen here
        for rc in (1e-2, 1e-3, 1e-6):
            Sig = _tangent_pinv_on_support(Jnl, spres, rcond=rc)
            rec[f"crlb_sqrt_rc{rc:g}"] = float(np.sqrt(np.clip(np.trace(Sig), 0, None)))
        # ---- COVERAGE (analytic, deterministic): expected realized informative support ----
        dmask = (css > 0) & (css < supp.size)                 # distinguishing among present
        for lam in COVS:
            # E[# distinguishing k-mers with >=1 count] = Σ (1 - e^{-λ μ_k})
            rec[f"exp_realized_dist_cov{lam:g}"] = float((1.0 - np.exp(-lam * muw[dmask])).sum())
        rows.append(rec)
        if wi % 500 == 0:
            print(f"  window {wi}/{len(blocks)}  {time.time()-t0:.0f}s", flush=True)

    panel = pd.DataFrame(rows)
    tgt = pd.read_csv(HU / "results/window_certainty.tsv", sep="\t")
    df = panel.merge(tgt, on="window", suffixes=("", "_tgt"))
    df = df[df.n_kmers >= 200].copy()
    for lam in COVS:
        df[f"excess_cov{lam:g}"] = df[f"err_cov{lam:g}"] - df["err_noiseless"]
    df.to_csv(HU / "results/metric_prototype_windows.tsv", sep="\t", index=False)

    def corr(a, b):
        d = df[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
        return np.corrcoef(d[a], d[b])[0, 1] if len(d) > 5 else np.nan
    def corrlog(a, b):
        d = df[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
        d = d[d[a] > 0]
        return np.corrcoef(np.log10(d[a]), d[b])[0, 1] if len(d) > 5 else np.nan

    print(f"\n{'='*70}\nWINDOW MODE  (n={len(df)} windows >=200 k-mers)\n{'='*70}", flush=True)
    print("(A) RESOLVABILITY — predict the noiseless floor err_noiseless:", flush=True)
    print(f"   [baselines from window_certainty.py]", flush=True)
    print(f"     corr(eff_rank[count-weighted],   floor) = {corr('eff_rank','err_noiseless'):+.3f}", flush=True)
    print(f"     corr(n_local_truefounders,       floor) = {corr('n_local_truefounders','err_noiseless'):+.3f}", flush=True)
    print(f"   [new PANEL-ONLY candidates]", flush=True)
    for c in ["ndist_apriori","ndist_cond","nsig_apriori","nsig_cond",
              "gram_eff_apriori","gram_eff_cond","nfnd_present"]:
        print(f"     corr({c:18s}, floor) = {corr(c,'err_noiseless'):+.3f}   "
              f"corr(log10, floor) = {corrlog(c,'err_noiseless'):+.3f}", flush=True)
    print(f"   [CRLB — principled identifiability predictor, vs floor]", flush=True)
    for rc in (1e-2, 1e-3, 1e-6):
        c = f"crlb_sqrt_rc{rc:g}"
        print(f"     corr({c:18s}, floor) = {corr(c,'err_noiseless'):+.3f}   "
              f"corr(log10, floor) = {corrlog(c,'err_noiseless'):+.3f}", flush=True)
    print("\n(B) COVERAGE — predict coverage_excess = err_cov - err_noiseless:", flush=True)
    for lam in COVS:
        ex = f"excess_cov{lam:g}"
        print(f"   cov {lam:g}:  corr(fisherSEtot, excess) = {corr(f'fisherSEtot_cov{lam:g}', ex):+.3f}   "
              f"corr(nzk, excess) = {corr(f'nzk_cov{lam:g}', ex):+.3f}   "
              f"corr(exp_realized_dist, excess) = {corr(f'exp_realized_dist_cov{lam:g}', ex):+.3f}", flush=True)

    # ---------------- GLOBAL MODE (whole chromosome = one region) ----------------
    print(f"\n{'='*70}\nGLOBAL MODE  (whole Chr1, one region)\n{'='*70}", flush=True)
    # omega-weighted panel Gram over ALL k-mers, chunked to stay in memory
    def gram_chunked(rowsel):
        Ksel = K[rowsel]; G = np.zeros((Ksel.shape[0], Ksel.shape[0]), dtype=np.float64)
        step = 1_000_000
        for s0 in range(0, Kn, step):
            sl = slice(s0, min(s0 + step, Kn))
            Kc = Ksel[:, sl].astype(np.float64) * omega[sl]
            G += Kc @ Ksel[:, sl].astype(np.float64).T
        return G
    g_ap = eff_rank_gram(gram_chunked(np.arange(F)))
    g_cd = eff_rank_gram(gram_chunked(supp))
    print(f"   gram_eff_apriori (all {F})      = {g_ap:.2f}", flush=True)
    print(f"   gram_eff_cond    (present {supp.size}) = {g_cd:.2f}", flush=True)
    print(f"   [ref] count-weighted Fisher eff_rank ~ 20 of 62 (window_certainty/summary)", flush=True)
    print(f"   [ref] global noiseless floor err = 0.0485 (coverage-flat); n_hat_supp 52 vs 50", flush=True)
    print(f"   GLOBAL coverage metric: at cohort depth >>30x, realized support ~ full -> excess ~ 0", flush=True)
    print(f"\nwrote results/metric_prototype_windows.tsv   ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
