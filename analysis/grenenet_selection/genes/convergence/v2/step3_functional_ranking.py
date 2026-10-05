#!/usr/bin/env python
"""v2 step 3 -- functional handle, themes and ranking for the NEW candidate list (GEA only; events + call rate 0.50).

RANKING PRINCIPLE (user, 2026-10-02): what matters most is whether an SV (or indel) plausibly has a FUNCTION -- it lands in a gene,
in a promoter, or beside one. So genes are ranked first by their best functional handle, then by the GEA line, then by theme.

Functional handle per unit (SV events and small indels; MNPs dropped -- they are substitutions, not length changes). Strongest first:
  F1_CDS_disruptive  in a CDS and (a small indel with size % 3 != 0, i.e. a frameshift) OR an SV event
  F2_CDS_inframe     in a CDS, small indel, size % 3 == 0
  F3_UTR             5'/3' UTR or non-coding exon
  F4_promoter        <= 1 kb upstream of the TSS (strand-aware; PROMOTER TAKES PRECEDENCE OVER TE for this list)
  F5_downstream      <= 1 kb past a gene's 3' end (only for units the classifier left non-genic)
  F6_intron
  F7_proximal        none of the above, but the nearest gene is within 5 kb (target = that gene)
Within a tier an SV event ranks ahead of a small indel (rank value = tier + 0 for SV, + 0.5 for indel).

Themes (user interests): flowering, circadian, oxidative_stress, and the climate axes temperature / water / light, matched on the gene's
protein name + GO biological-process terms + UniProt keywords (TAIR GO + UniProt annotator). Annotation coverage is partial -- an
unannotated gene cannot be flagged -- so `n_annot_terms` is carried and the flags are a net, not a classification.

Gene-level GEA evidence is recomputed over the units that have a handle on the gene (target_gene): E1 recurrence over >= 2 of 7 climate
clusters, E2 a clean axis (lambda < 2.0), L_gea = E1 & E2 (one gated line). Raw uncalibrated LFMM p.
Inputs : results/candidate_list_v2_events_cr50/{v2_variants_classified.csv, gea_pool.csv}
Outputs: v2_functional_units.csv (one row per unit with a handle), v2_candidate_ranking.csv (one row per gene, ranked)
env: kmate. Compute node (annotator needs outbound HTTPS).
"""
from __future__ import annotations
import os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, CONV); sys.path.insert(0, os.path.abspath(os.path.join(CONV, "..", "..")))
import lib                                                      # noqa: E402
import build_convergence as bc                                  # noqa: E402
import build_functional_variants as fv                          # noqa: E402  (downstream_flank, FLANK_BP)

V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")
PROXIMAL_BP = 5000
TIERS = ["F1_CDS_disruptive", "F2_CDS_inframe", "F3_UTR", "F4_promoter", "F5_downstream", "F6_intron", "F7_proximal"]

THEMES = {
    "flowering": ["flower", "floral", "flowering time", "vernaliz", "photoperiod", "inflorescence", "meristem identity",
                  "reproductive development", "transition to flowering", "flc", "short vegetative phase"],
    "circadian": ["circadian", "rhythm", "clock", "entrainment"],
    "oxidative_stress": ["oxidative stress", "reactive oxygen", "hydrogen peroxide", "superoxide", "glutathione", "ascorbate",
                         "peroxidase", "peroxiredoxin", "thioredoxin", "redox", "antioxidant"],
    "temperature": ["cold", "freezing", "chilling", "frost", "heat", "thermo", "temperature", "heat shock", "cbf", "dreb",
                    "cold acclimation", "vernalization"],
    "water": ["drought", "water deprivation", "dehydrat", "desiccation", "osmotic", "salt stress", "salinity", "abscisic",
              "stomat", "hyperosmotic"],
    "light": ["response to light", "photomorphogenesis", "red light", "blue light", "far-red", "ultraviolet", "response to uv",
              "uv-b", "shade avoidance", "photoreceptor", "phytochrome", "cryptochrome"],
}
THEME_ORDER = list(THEMES)


def functional_units(V: pd.DataFrame) -> pd.DataFrame:
    n0 = len(V)
    V = V[V["size"] > 0].copy()                                  # SV events + small indels only
    print(f"{n0:,} units -> {len(V):,} length-changing (SV events + small indels)")
    nong = ~V.tier.isin(bc.GENIC)
    V["downstream_gene"] = ""; V["downstream_dist"] = np.nan
    flk = fv.downstream_flank(V[nong])
    V.loc[flk.index, "downstream_gene"] = flk.downstream_gene; V.loc[flk.index, "downstream_dist"] = flk.downstream_dist
    is_cds = V.tier.eq("1_CDS")
    V["frameshift"] = is_cds & V.vclass.eq("smallindel") & V["size"].mod(3).ne(0)
    V["cds_sv"] = is_cds & V.vclass.eq("sv")
    f = pd.Series("", index=V.index, dtype=object)
    f[V.dist_to_gene.le(PROXIMAL_BP) & nong] = "F7_proximal"
    f[V.tier.eq("5_intron")] = "F6_intron"
    f[V.downstream_gene.ne("") & f.isin(["", "F7_proximal"])] = "F5_downstream"
    f[V.tier.eq("3_promoter")] = "F4_promoter"
    f[V.tier.isin(["2_UTR", "4_exon_noncoding"])] = "F3_UTR"
    f[is_cds] = "F2_CDS_inframe"
    f[V.frameshift | V.cds_sv] = "F1_CDS_disruptive"
    V["ftier"] = f
    V["target_gene"] = np.select([V.ftier.eq("F5_downstream"), V.ftier.eq("F7_proximal")],
                                 [V.downstream_gene, V.nearest_gene], default=V.gene)
    V = V[V.ftier.ne("") & V.target_gene.ne("")].copy()
    V["tier_rank"] = V.ftier.map({t: i + 1 for i, t in enumerate(TIERS)})
    V["rank_value"] = V.tier_rank + np.where(V.vclass.eq("sv"), 0.0, 0.5)
    return V


def theme_flags(ann: pd.DataFrame) -> pd.DataFrame:
    text = (ann.symbol.fillna("") + " " + ann.protein_name.fillna("") + " " + ann.go_bp.fillna("") + " " + ann.uniprot_keywords.fillna("")).str.lower()
    out = pd.DataFrame(index=ann.index)
    for th, words in THEMES.items():
        pat = re.compile("|".join(re.escape(w) for w in words))
        out[f"theme_{th}"] = text.map(lambda t: bool(pat.search(t)))
    out["themes"] = out.apply(lambda r: ",".join(th for th in THEME_ORDER if r[f"theme_{th}"]), axis=1)
    out["theme_any"] = out.themes.ne("")
    out["n_annot_terms"] = ann.go_bp.fillna("").map(lambda s: len([x for x in s.split(";") if x]))
    return out


def main():
    V = pd.read_csv(f"{V2}/v2_variants_classified.csv")
    F = functional_units(V)
    print(F.ftier.value_counts().reindex(TIERS).fillna(0).astype(int).to_string())
    print("\nby tier x class:"); print(pd.crosstab(F.ftier, F.vclass).reindex(TIERS).fillna(0).astype(int).to_string())

    # ---- gene-level table over the units with a handle ----
    def agg(d: pd.DataFrame) -> pd.Series:
        best = d.sort_values(["rank_value", "gea_n_clusters", "gea_nlp"], ascending=[True, False, False]).iloc[0]
        return pd.Series({
            "chrom": d.chrom.iat[0], "best_ftier": best.ftier, "best_rank_value": best.rank_value, "best_unit_class": best.vclass,
            "best_unit": int(best.unit), "best_unit_size": int(best["size"]), "best_unit_te_overlap": bool(best.te_overlap),
            "n_units": len(d), "n_sv_events": int(d.vclass.eq("sv").sum()), "n_smallindel": int(d.vclass.eq("smallindel").sum()),
            "n_sv_with_handle_F1_F4": int((d.vclass.eq("sv") & d.tier_rank.le(4)).sum()),
            "tiers": ",".join(sorted(set(d.ftier))),
            "gea_nlp": d.gea_nlp.max(), "gea_n_clusters": int(d.gea_n_clusters.max()), "gea_n_axes": int(d.gea_n_axes.max()),
            "gea_min_lam": d.gea_min_lam.min(),
            "gea_sig_clusters": ",".join(sorted({c for s in d.gea_sig_clusters for c in str(s).split(",")})),
            "units": ",".join(map(str, sorted(d.unit))),
        })
    A = F.groupby("target_gene").apply(agg, include_groups=False).reset_index().rename(columns={"target_gene": "gene"})
    A["E1_gea_multi_cluster"] = A.gea_n_clusters >= 2
    A["E2_gea_clean_axis"] = A.gea_min_lam < 2.0
    A["L_gea"] = A.E1_gea_multi_cluster & A.E2_gea_clean_axis
    A = bc.add_loci(A)

    ann = bc.annotate(A.gene.tolist()).set_index("gene")
    for col in ("symbol", "protein_name", "uniprot_function", "uniprot_keywords", "go_bp", "categories"):
        A[col] = A.gene.map(ann[col]) if col in ann.columns else ""
    A = A.reset_index(drop=True)                                 # add_loci sorts by position and keeps the old labels
    T = theme_flags(A)                                           # same index as A -> join by label, never by position
    A = A.join(T)
    assert A.theme_any.eq(A.themes.ne("")).all()

    A = A.sort_values(["best_rank_value", "L_gea", "theme_any", "gea_n_clusters", "gea_nlp"],
                      ascending=[True, False, False, False, False]).reset_index(drop=True)
    A.insert(0, "rank", np.arange(1, len(A) + 1))
    F.to_csv(f"{V2}/v2_functional_units.csv", index=False)
    A.to_csv(f"{V2}/v2_candidate_ranking.csv", index=False)

    print(f"\ngenes with a functional handle: {len(A):,} | L_gea {int(A.L_gea.sum())} | any theme {int(A.theme_any.sum())}")
    print("best tier per gene:", A.best_ftier.value_counts().reindex(TIERS).fillna(0).astype(int).to_dict())
    print("themes:", {th: int(A[f'theme_{th}'].sum()) for th in THEME_ORDER})
    print(f"wrote {V2}/v2_functional_units.csv, v2_candidate_ranking.csv")


if __name__ == "__main__":
    main()
