#!/usr/bin/env python
"""Build a registry of every founder-burden model fitted under results/sv_adaptive/final_models*/.

Retrospective collector: the fits do not record their own specification, but the specification is fully
encoded in the output directory name plus the filename suffix that `_sv_burden_final_models.R` and
`_sv_burden_pooled_models.R` build from their environment variables. This walks the coefficient files,
decodes that specification, and emits two tidy tables:

  model_specs.csv     one row per fitted model -- the human-readable index of "what have we run"
  model_registry.csv  one row per (model x term) -- estimates, SEs, p-values, and the LRT/Bonferroni
                      for the interaction terms; the table a specification curve is drawn from

Both are written to results/sv_adaptive/. Nothing is refitted and nothing existing is modified.

Suffix grammar, in the order the R scripts build it:
    coef_<model><_x_EXTRA><_SUBSET><_reRE_DROP>.csv
EXTRA and RE_DROP are alphanumeric-squashed by the R scripts (gsub("[^A-Za-z0-9]+","")), so decoding them
is lossy; known tokens are mapped back through EXTRA_MAP / REDROP_MAP and anything unrecognised is carried
through verbatim and flagged in `spec_decoded`, never silently guessed. `basic` env (pandas, scipy).
"""
import os, re, sys, glob
import pandas as pd, numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HERE, "results", "sv_adaptive")

# ---- what each output directory means -------------------------------------------------------------
# unit  : whether insertions/deletions are separate terms (split) or summed per context (pooled)
# origin: which variants build the burdens -- `shared` (carried by >=1 cold AND >=1 warm founder) or `all`
# te    : sv6 = reference TE is its own burden term; sv4 = no TE term, TE-context variants enter no burden
DIRS = {
    "final_models":                 dict(unit="split",  origin="shared", te="sv6"),
    "final_models_sv4_noTE":        dict(unit="split",  origin="shared", te="sv4"),
    "final_models_sv4_all":         dict(unit="split",  origin="all",    te="sv4"),
    "final_models_sv6_all":         dict(unit="split",  origin="all",    te="sv6"),
    "final_models_pooled":          dict(unit="pooled", origin="shared", te="sv6"),
    "final_models_pooled_all":      dict(unit="pooled", origin="all",    te="sv6"),
    "final_models_pooled_all_noTE": dict(unit="pooled", origin="all",    te="sv4"),
}
SPLIT_TERMS  = ["ins_CDS", "ins_TE", "ins_promoter", "del_CDS", "del_TE", "del_promoter"]
POOLED_TERMS = ["sv_CDS", "sv_TE", "sv_promoter"]
# the interaction set actually tested in each (unit, te) combination -- this is the Bonferroni denominator
TESTED = {
    ("split",  "sv6"): SPLIT_TERMS,
    ("split",  "sv4"): ["ins_CDS", "ins_promoter", "del_CDS", "del_promoter"],
    ("pooled", "sv6"): POOLED_TERMS,
    ("pooled", "sv4"): ["sv_CDS", "sv_promoter"],
}
# EXTRA_FIXED / RE_DROP are alphanumeric-squashed in the filename; map the tokens we have actually run
EXTRA_MAP = {
    "": "",
    "obio1cobio1cbio1c": "obio1c + obio1c:bio1c",
    "obio1cbio1c":       "obio1c:bio1c",
    "cactusbio1c":       "cactus:bio1c",
}
REDROP_MAP = {
    "": [],
    "fid":            ["fid"],
    "fidsvCDS":       ["fid", "sv_CDS"],
    "svCDSsvTE":      ["sv_CDS", "sv_TE"],
    "fiddelCDSinsTE": ["fid", "del_CDS", "ins_TE"],
}
SUBSET_N = {"": 231, "cactus": 80, "shortread": 151}
FREQ_BIN_NAMES = ["private", "rare", "common",          # frequency bins from _sv_burden_strata.py STRATUM=freq
                  "small", "medium", "large"]          # size bins from STRATUM=size
N_GARDENS = 31


def parse_name(stem, unit):
    """Split `<model><suffix>` into its parts. Returns (model, extra_tok, subset, redrop_tok) or None."""
    terms = POOLED_TERMS if unit == "pooled" else SPLIT_TERMS
    # BURDEN_STRATIFY=freq renames a burden sv_<ctx> -> sv_<ctx>_<bin>, so the drop-refits for those
    # carry the stratified name; include them or the LRT files for stratified fits go unparsed
    terms = terms + [f"{t}_{b}" for t in terms for b in FREQ_BIN_NAMES]
    p = "p" if unit == "pooled" else ""
    bases = [f"M1{p}_reml", f"M2{p}_reml", f"M2{p}_ml"] + [f"M2{p}_ml_drop_{t}" for t in terms]
    # longest first so "M2p_ml_drop_sv_CDS" wins over "M2p_ml"
    for b in sorted(bases, key=len, reverse=True):
        if stem == b or stem.startswith(b + "_"):
            suf = stem[len(b):]
            m = re.fullmatch(r"(?:_x_([A-Za-z0-9]+))?(?:_(cactus|shortread))?(?:_re([A-Za-z0-9]+))?"
                             r"(?:_freq([A-Za-z0-9]+))?", suf)
            if m is None:
                return None
            return b, m.group(1) or "", m.group(2) or "", m.group(3) or "", m.group(4) or ""
    return None


def main():
    spec_rows, coef_rows, unknown = [], [], []
    for d, meta in DIRS.items():
        path = os.path.join(G, d)
        if not os.path.isdir(path):
            continue
        for f in sorted(glob.glob(os.path.join(path, "coef_*.csv"))):
            stem = os.path.basename(f)[len("coef_"):-len(".csv")]
            parsed = parse_name(stem, meta["unit"])
            if parsed is None:
                unknown.append(os.path.join(d, os.path.basename(f)))
                continue
            model, extra_tok, subset, redrop_tok, freq_tok = parsed
            decoded = (extra_tok in EXTRA_MAP) and (redrop_tok in REDROP_MAP)
            extra  = EXTRA_MAP.get(extra_tok, f"<UNDECODED:{extra_tok}>")
            redrop = REDROP_MAP.get(redrop_tok, [f"<UNDECODED:{redrop_tok}>"])
            suffix = stem[len(model):]
            tested = TESTED[(meta["unit"], meta["te"])]
            if freq_tok:   # burden names differ; take the tested set from the coefficient file itself
                tested = [t[:-len(":bio1c")] for t in pd.read_csv(f)["term"] if t.endswith(":bio1c")
                          and not t.startswith(("obio1c", "bio1c", "cactus"))]
            n_f = SUBSET_N[subset]

            ll_f = os.path.join(path, f"loglik_{stem}.csv")
            ll = pd.read_csv(ll_f).iloc[0].to_dict() if os.path.exists(ll_f) else {}
            co = pd.read_csv(f).set_index("term")

            # A fit written after the spec-writer was added carries its own specification row. Prefer it:
            # it holds the literal formula and the true founder/row counts, so it cannot be wrong the way
            # a decoded filename can. Older fits fall back to the decoding above.
            sp_f = os.path.join(path, f"spec_{stem}.csv")
            sp = pd.read_csv(sp_f).iloc[0].to_dict() if os.path.exists(sp_f) else {}
            if sp:
                # R writes an empty string for an absent EXTRA_FIXED / RE_DROP, which pandas reads back as
                # NaN -- coerce those to "" rather than letting the string "nan" leak into the registry
                txt = lambda k: "" if not isinstance(sp.get(k), str) else sp[k]
                extra   = txt("extra_fixed")
                redrop  = [x for x in txt("re_drop").split(",") if x]
                subset  = "" if txt("subset") in ("", "all_founders") else txt("subset")
                n_f     = int(sp.get("n_founders", n_f))
                decoded = True

            spec_rows.append(dict(
                model_id=f"{d}/{stem}", directory=d, model=model, suffix=suffix,
                unit=meta["unit"], origin=meta["origin"], te_term=meta["te"],
                local_adaptation=("obio1c" in extra), extra_fixed=extra, stratify=freq_tok,
                subset=subset or "all_founders", re_drop=",".join(redrop),
                n_founders=n_f, n_gardens=N_GARDENS, n_rows=int(sp.get("n_rows", n_f * N_GARDENS)),
                formula=sp.get("formula", ""), fixed=sp.get("fixed", ""), random=sp.get("random", ""),
                onestage_dir=sp.get("onestage_dir", ""), fitted_at=sp.get("fitted_at", ""),
                git_commit=sp.get("git_commit", ""), self_documented=bool(sp),
                estimator="REML" if model.endswith("_reml") else "ML",
                has_interactions=model.startswith(("M2", "M2p")),
                is_lrt_refit="_drop_" in model,
                logLik=ll.get("logLik", np.nan), df=ll.get("df", np.nan),
                singular=ll.get("singular", np.nan),
                n_interactions_tested=len(tested), spec_decoded=decoded))

            # LRT + Bonferroni attach to the REML fit whose ML twin has the drop refits
            lrt = {}
            if model in (f"M2{'p' if meta['unit']=='pooled' else ''}_reml",):
                mlb = model.replace("_reml", "_ml")
                base = os.path.join(path, f"loglik_{mlb}{suffix}.csv")
                if os.path.exists(base):
                    b = pd.read_csv(base).logLik[0]
                    for t in tested:
                        dr = os.path.join(path, f"loglik_{mlb}_drop_{t}{suffix}.csv")
                        if os.path.exists(dr):
                            lrt[t] = stats.chi2.sf(2 * (b - pd.read_csv(dr).logLik[0]), 1)
            for t, r in co.iterrows():
                base_t = t.split(":")[0]
                p_lrt = lrt.get(base_t, np.nan) if ":bio1c" in t else np.nan
                coef_rows.append(dict(
                    model_id=f"{d}/{stem}", directory=d, model=model, suffix=suffix,
                    unit=meta["unit"], origin=meta["origin"], te_term=meta["te"],
                    local_adaptation=("obio1c" in extra), extra_fixed=extra, stratify=freq_tok,
                    subset=subset or "all_founders", re_drop=",".join(redrop),
                    n_founders=n_f, n_rows=n_f * N_GARDENS,
                    estimator="REML" if model.endswith("_reml") else "ML",
                    singular=ll.get("singular", np.nan),
                    term=t,
                    term_kind=("interaction" if ":" in t else
                               "intercept" if t == "(Intercept)" else
                               "local_adaptation" if "obio1c" in t else
                               "burden" if t.startswith(("sv_", "ins_", "del_", "snp_")) else "covariate"),
                    estimate=r["estimate"], se=r["se"], z=r.get("z", np.nan), p_wald=r["p"],
                    lrt_p=p_lrt,
                    bonferroni=min(p_lrt * len(tested), 1.0) if p_lrt == p_lrt else np.nan))

    S = pd.DataFrame(spec_rows).sort_values(["directory", "model", "suffix"]).reset_index(drop=True)
    C = pd.DataFrame(coef_rows).sort_values(["directory", "model", "suffix", "term"]).reset_index(drop=True)
    S.to_csv(os.path.join(G, "model_specs.csv"), index=False)
    C.to_csv(os.path.join(G, "model_registry.csv"), index=False)

    print(f"[wrote] {G}/model_specs.csv     {len(S):>4} fitted models")
    print(f"[wrote] {G}/model_registry.csv  {len(C):>4} model x term rows\n")
    if unknown:
        print(f"WARNING: {len(unknown)} coefficient file(s) could not be parsed:")
        for u in unknown: print(f"  {u}")
    bad = S[~S.spec_decoded]
    if len(bad):
        print(f"\nWARNING: {len(bad)} model(s) have an undecoded EXTRA_FIXED / RE_DROP token "
              f"(add it to EXTRA_MAP / REDROP_MAP):")
        for r in bad.itertuples(): print(f"  {r.model_id}  extra={r.extra_fixed} re_drop={r.re_drop}")
    print("\nDISTINCT SPECIFICATIONS (REML fits; M1 = no interactions, M2 = with):")
    rep = S[(S.estimator == "REML") & (~S.is_lrt_refit)].copy()
    rep["models"] = rep.model
    g = (rep.groupby(["unit", "origin", "te_term", "extra_fixed", "subset", "re_drop", "stratify"])
            .agg(fits=("model", lambda v: "+".join(sorted(set(v)))),
                 n=("model", "size"), singular=("singular", lambda v: bool(np.any(v == True))))
            .reset_index())
    print(g.to_string(index=False))
    print(f"\n{len(g)} distinct specifications, {len(rep)} reported fits, {len(S)} fits in total "
          f"(the rest are ML refits used only for the LRTs).")

    # self-check: the numbers quoted for the current main model must come back out of the registry
    q = C[(C.model_id == "final_models_pooled_all/M2p_reml_x_obio1cobio1cbio1c")]
    if len(q):
        print("\nSELF-CHECK — the current main model, read back from the registry:")
        for t in ("sv_CDS", "sv_CDS:bio1c", "bio1c:obio1c"):
            r = q[q.term == t]
            if len(r):
                r = r.iloc[0]
                extra = f"   LRT p {r.lrt_p:.4g}, Bonferroni {r.bonferroni:.4g}" if r.lrt_p == r.lrt_p else ""
                print(f"  {t:<16} {r.estimate:+.5f} +- {r.se:.5f}   Wald p {r.p_wald:.4g}{extra}")
    return S, C


if __name__ == "__main__":
    main()
