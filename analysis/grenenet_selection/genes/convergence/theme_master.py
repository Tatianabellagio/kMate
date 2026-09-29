#!/usr/bin/env python
"""Tag master-table candidates with the biology of interest: stress / temperature /
circadian-light / flowering.

The theme regexes are **imported** from `dissection/theme_filter.py`, not copied -- one
definition of what counts as a flowering gene, in one place.

What is different here: `theme_filter.py` matches against `symbol` + `protein_name` +
`categories` only. The master table also carries UniProt free-text **FUNCTION**,
UniProt **keywords**, and **GO biological process** for most genes, and those are where
the functional evidence actually lives -- a gene called "Uncharacterized protein" can
still have a curated FUNCTION mentioning cold acclimation. Matching the richer text is
the whole reason CLAUDE.md mandates the TAIR+UniProt annotator over mygene alone, so
not using those fields would waste the annotation.

`theme_fields` records which field produced the hit, so a match on free text can be
told apart from a match on the gene's name.

Outputs -> results/
  master_candidate_genes.csv    themes/theme_fields columns added in place
  themed_candidates.csv         themed genes only, ranked by evidence then mechanism

env: kmate.
"""
from __future__ import annotations
import os
import re
import sys
import importlib.util as ilu
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = f"{HERE}/results"

# the text fields to search, richest first
FIELDS = ["symbol", "protein_name", "uniprot_function", "uniprot_keywords",
          "go_bp", "categories"]


def _themes_mod():
    p = os.path.join(os.path.dirname(HERE), "dissection", "theme_filter.py")
    spec = ilu.spec_from_file_location("theme_filter", p)
    mod = ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tag(M: pd.DataFrame, THEMES: dict, CAT_MAP: dict) -> pd.DataFrame:
    themes, fields = [], []
    for r in M.itertuples():
        hits, where = set(), set()
        for f in FIELDS:
            txt = str(getattr(r, f, "") or "")
            if not txt or txt == "nan":
                continue
            for th, pat in THEMES.items():
                if re.search(pat, txt, re.I):
                    hits.add(th)
                    where.add(f)
        for c in str(getattr(r, "categories", "") or "").split(","):
            c = c.strip()
            if c in CAT_MAP:
                hits.add(CAT_MAP[c])
                where.add("categories")
        themes.append(",".join(sorted(hits)))
        fields.append(",".join(sorted(where)))
    M = M.copy()
    M["themes"] = themes
    M["theme_fields"] = fields
    M["n_themes"] = M.themes.apply(lambda s: len(s.split(",")) if s else 0)
    return M


def main():
    tf = _themes_mod()
    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv")
    M = M.drop(columns=[c for c in ("themes", "theme_fields", "n_themes")
                        if c in M.columns])
    M = tag(M, tf.THEMES, tf.CAT_MAP)
    M.to_csv(f"{OUT}/master_candidate_genes.csv", index=False)

    T = M[M.themes.ne("")].copy()
    T = T.sort_values(["n_independent", "gea_n_clusters", "gwas_n_gardens"],
                      ascending=False)
    T.to_csv(f"{OUT}/themed_candidates.csv", index=False)

    print(f"{len(T):,} of {len(M):,} candidate genes carry a theme\n")
    for th in tf.THEMES:
        sub = M[M.themes.str.contains(th, na=False)]
        short = sub[sub.n_independent >= 1]
        print(f"  {th:16s} {len(sub):>4} genes   "
              f"({len(short)} with >=1 gated line of evidence)")

    # how much the extra fields bought
    namey = M[M.themes.ne("") & M.theme_fields.str.contains("symbol|protein_name")]
    textonly = M[M.themes.ne("") & ~M.theme_fields.str.contains("symbol|protein_name")]
    print(f"\nmatched on name/symbol: {len(namey):,}   "
          f"ONLY on free text (FUNCTION/keywords/GO): {len(textonly):,}")
    print(f"\nwrote {OUT}/themed_candidates.csv, themes added to "
          f"master_candidate_genes.csv")


if __name__ == "__main__":
    main()
