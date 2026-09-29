#!/usr/bin/env python
"""From the functional candidates (SV/indel in CDS/UTR/promoter), pull the genes
with a role in CLIMATE / FLOWERING / CIRCADIAN / STRESS -- the biology of interest.

Themes matched on annotated protein_name + symbol + the TAIR/UniProt `categories`
tag. Restricted to vclass in {sv, smallindel} (the variant itself is an SV or indel,
placed in CDS/UTR/promoter, so the link is physical, not LD-tagged).
Writes themed_candidates.csv, prints the hits per theme.
env: kmate.
"""
import os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
THEMES = {
    "flowering": r"flowering|florig|\bMADS\b|\bAGL\d|\bFLC\b|CONSTANS|\bSOC1\b|\bSVP\b|"
                 r"vernaliz|floral|photoperiod|\bFT\b|GIGANTEA|FRIGIDA|APETALA|LEAFY|"
                 r"\bMAF\d|meristem identity|AGAMOUS",
    "circadian_light": r"circadian|clock|\bCCA1\b|\bLHY\b|\bTOC1\b|pseudo-response regulat|"
                       r"\bPRR\d|EARLY FLOWERING|\bELF\d|\bLUX\b|zeitlupe|\bZTL\b|"
                       r"phytochrome|cryptochrome|photoreceptor|light.?harvest|"
                       r"light signal|\bCOP1\b|\bPIF\d|phototropin|blue.?light",
    "climate_temp": r"\bcold\b|\bheat\b|temperature|freezing|chilling|thermo|"
                    r"heat shock|\bHSP\d|\bHSF|\bCBF\d|\bDREB|C-repeat|"
                    r"dehydration-responsive|frost",
    "stress": r"\bstress\b|drought|abscisic|\bABA\b|osmotic|\bsalt\b|salin|dehydr|"
              r"desiccation|\bLEA\b|late embryogenesis|\bRD\d\d|oxidative|"
              r"reactive oxygen|glutathione|peroxidas|catalase|superoxide|"
              r"aquaporin|water channel|proline|osmo|redox|thioredoxin",
}
CAT_MAP = {"flowering": "flowering", "temperature": "climate_temp",
           "water": "stress", "oxidative": "stress", "light": "circadian_light"}


def themes_for(row):
    text = f"{row.get('symbol','')} {row.get('protein_name','')}".lower()
    hits = set()
    for th, pat in THEMES.items():
        if re.search(pat, text, re.I):
            hits.add(th)
    for c in str(row.get("categories", "")).split(","):
        if c in CAT_MAP:
            hits.add(CAT_MAP[c])
    return ",".join(sorted(hits))


def main():
    d = pd.read_csv(f"{HERE}/functional_candidates.csv")
    d = d[d.vclass.isin(["sv", "smallindel"])].copy()      # variant IS an SV or indel
    d["themes"] = d.apply(themes_for, axis=1)
    T = d[d.themes != ""].copy().sort_values(
        ["themes", "tier", "best_nlp"], ascending=[True, True, False])
    cols = ["chrom", "pos", "size", "vclass", "tier", "region", "gene", "symbol",
            "protein_name", "themes", "best_axis", "best_nlp", "n_axes", "MAF",
            "snp_cosig_2kb", "block_gene_mismatch"]
    T[cols].to_csv(f"{HERE}/themed_candidates.csv", index=False)
    print(f"themed SV/indel candidates (CDS/UTR/promoter): {len(T)} of {len(d)} "
          f"functional SV/indel hits\n")
    for th in THEMES:
        sub = T[T.themes.str.contains(th)].sort_values("best_nlp", ascending=False)
        print(f"===== {th.upper()}  (n={len(sub)}) =====")
        show = sub[["chrom", "pos", "size", "vclass", "tier", "gene", "symbol",
                    "protein_name", "best_axis", "best_nlp", "n_axes",
                    "snp_cosig_2kb"]].head(14)
        with pd.option_context("display.width", 250, "display.max_colwidth", 46):
            print(show.to_string(index=False))
        print()


if __name__ == "__main__":
    main()
