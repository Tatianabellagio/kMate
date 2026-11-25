#!/usr/bin/env python
"""v2 step 4 -- merge every evidence layer into per-unit evidence and the final ranked gene list (GEA only; events + call rate 0.50; NO GWAS).

Unit level (one row per unit with a functional handle, keyed on `unit`), layers merged on `unit`:
  v2_functional_units.csv   handle tier (F1..F7), target gene, GEA recurrence (the unit's own gea_* columns)
  v2_layer_atac.csv         open chromatin overlap (multi-tissue ATAC)
  v2_layer_snp_visibility.csv   SNP shadowing (Bonferroni SNP within 2 kb) and independent short-read SNP tagging r2
  v2_layer_coding.csv       effect on the protein (consequence, fraction lost, Pfam domain, protein_grade)
  v2_layer_expression.csv   carrier vs non-carrier expression of the target gene (1001 Transcriptomes), lineage-controlled
  v2_layer_tfbs.csv         REF-vs-ALT TF-motif turnover (if present)
  v2_layer_eqtl.csv         cis-eQTL readout at the tag SNP (INFORMATIONAL only: the scan covers 22 genes)

Lines of evidence beyond the functional handle (unit level, independent of each other as far as the data allow):
  U_gea       gated GEA recurrence of THIS unit: >= 2 of 7 climate clusters AND a clean axis (lambda < 2.0)
  U_protein   protein_grade >= 2 (gene loss / NMD / >= 10% of the protein lost / coding SV / in-frame >= 10 aa / domain hit)
  U_expr      carriers differ from non-carriers in the target gene's expression AND beat the lineage control (p_emp <= 0.05)
  U_regulatory  overlaps open chromatin (ATAC) OR a clean (non-repeat) TF-motif gain/loss
  U_snp_blind SNP-unique (no Bonferroni SNP within 2 kb) AND not tagged by short-read SNPs (best r2 < 0.2): invisible to a SNP study
n_extra_lines = U_gea + U_expr + U_regulatory + U_snp_blind (U_protein is part of "function", which leads the ranking).

RANKING (user, 2026-10-02): what matters most is that the SV / indel plausibly has a FUNCTION. Genes are sorted by
  1 best functional handle (tier, SV ahead of indel within a tier)   2 protein grade of that unit   3 n_extra_lines of the best unit
  4 gated GEA line of the gene   5 any user theme (flowering, circadian, oxidative stress, temperature, water, light)   6 GEA clusters, p.
Raw uncalibrated LFMM p; a candidate net, not a calibrated hit list.
Outputs: v2_unit_evidence.csv, v2_candidate_list.csv.   env: kmate. Compute node.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, CONV)
V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")


def rd(name, required=True):
    p = f"{V2}/{name}"
    if not os.path.exists(p):
        if required:
            raise SystemExit(f"missing {name}")
        print(f"  (optional layer absent: {name})"); return None
    return pd.read_csv(p, low_memory=False)


def main():
    F = rd("v2_functional_units.csv")
    U = F.copy()
    layers = {
        "atac": (rd("v2_layer_atac.csv"), ["in_atac", "atac_n_tissues"]),
        "snp": (rd("v2_layer_snp_visibility.csv"), ["snp_cosig_2kb", "snp_cosig_dist", "best_r2_snp", "best_r2_shortread", "tag_shortread_untestable",
                                                     "snp_blind_shortread", "tag_join_ambiguous"]),
        "coding": (rd("v2_layer_coding.csv"), ["site_feature", "protein_effect", "consequence", "nmd", "protein_frac_lost", "domain_hit", "domains_hit",
                                               "protein_grade", "repeat_copy", "effect_ambiguous"]),
        "expr": (rd("v2_layer_expression.csv"), ["status", "n_car_expr", "n_ref_expr", "med_raw_car", "med_raw_ref", "fold", "p_gene", "p_emp"]),
        "tfbs": (rd("v2_layer_tfbs.csv", required=False), None),
        "eqtl": (rd("v2_layer_eqtl.csv"), ["eqtl_status", "eqtl_nlp", "eqtl_beta", "eqtl_pct", "tag_method", "tag_stat"]),
    }
    for name, (L, cols) in layers.items():
        if L is None:
            continue
        keep = ["unit"] + (cols if cols else [c for c in L.columns if c not in U.columns])
        keep = [c for c in keep if c in L.columns]
        if name in ("expr", "tfbs", "eqtl"):                      # these are per (unit, target gene); the unit's handle gene is the one tested
            L = L.copy()
        ren = {c: f"{name}_{c}" for c in keep if c != "unit" and c in U.columns}
        U = U.merge(L[keep].rename(columns=ren).drop_duplicates("unit"), on="unit", how="left")
    assert U.unit.is_unique

    f = lambda c: U[c].fillna(False).astype(bool) if c in U.columns else pd.Series(False, index=U.index)
    U["U_gea"] = (U.gea_n_clusters >= 2) & (U.gea_min_lam < 2.0)
    U["U_protein"] = U.get("protein_grade", pd.Series(0, index=U.index)).fillna(0) >= 2
    U["U_expr"] = (U.get("expr_status", U.get("status", pd.Series("", index=U.index))).astype(str).eq("tested")) & (U.get("p_emp", pd.Series(np.nan, index=U.index)) <= 0.05)
    # clean TF-motif turnover: a motif lost or gained (FIMO p<1e-4, q<0.005, low-complexity dropped) that is not flagged as repeat-derived
    # (tfbs_repeat covers redundant and microsatellite matches -- the AT4G13200 lesson)
    tf_clean = pd.Series(False, index=U.index)
    if "tfbs_lost" in U.columns:
        turn = (pd.to_numeric(U.tfbs_lost, errors="coerce").fillna(0) > 0) | (pd.to_numeric(U.tfbs_gained, errors="coerce").fillna(0) > 0)
        rep = U.tfbs_repeat.fillna(False).astype(str).str.lower().isin(["true", "1", "1.0"])
        tf_clean = turn & ~rep
        U["tfbs_clean_turnover"] = tf_clean
    U["U_regulatory"] = f("in_atac") | tf_clean
    U["U_snp_blind"] = (~f("snp_cosig_2kb")) & f("snp_blind_shortread")
    U["n_extra_lines"] = U[["U_gea", "U_expr", "U_regulatory", "U_snp_blind"]].sum(axis=1)
    U["protein_grade"] = U.get("protein_grade", pd.Series(0, index=U.index)).fillna(0).astype(int)
    # tier refinement: the classifier places a unit in the coding tiers by POSITION; keep it there only if the transcript-level effect confirms a
    # protein-changing CDS consequence (not synonymous, not a splice-region / UTR call). Others drop to F3 (still a gene-body variant).
    sf = U.get("site_feature", pd.Series("", index=U.index)).fillna("").astype(str)
    cons = U.get("consequence", pd.Series("", index=U.index)).fillna("").astype(str)
    coding_tier = U.ftier.isin(["F1_CDS_disruptive", "F2_CDS_inframe"])
    demote = coding_tier & ~(sf.str.startswith("CDS") & cons.ne("synonymous"))
    U["ftier_classifier"] = U.ftier; U["rank_value_classifier"] = U.rank_value
    U.loc[demote, "ftier"] = "F3_UTR"
    U["tier_rank"] = U.ftier.map({"F1_CDS_disruptive": 1, "F2_CDS_inframe": 2, "F3_UTR": 3, "F4_promoter": 4, "F5_downstream": 5, "F6_intron": 6, "F7_proximal": 7})
    U["rank_value"] = U.tier_rank + np.where(U.vclass.eq("sv"), 0.0, 0.5)
    print(f"tier refinement: {int(demote.sum())} of {int(coding_tier.sum())} coding-tier units demoted to F3 (no protein-changing CDS effect)")
    U.to_csv(f"{V2}/v2_unit_evidence.csv", index=False)
    print(f"unit evidence: {len(U):,} units | " + " | ".join(f"{c} {int(U[c].sum()):,}" for c in ("U_gea", "U_protein", "U_expr", "U_regulatory", "U_snp_blind")))

    # ---- gene level ----
    R = rd("v2_candidate_ranking.csv")                           # step 3: annotation, themes, gene-level GEA evidence, loci
    best = (U.sort_values(["target_gene", "rank_value", "protein_grade", "n_extra_lines", "gea_nlp"], ascending=[True, True, False, False, False])
             .groupby("target_gene").head(1).set_index("target_gene"))
    G = R.set_index("gene").copy()
    for c in ("unit", "ftier", "vclass", "protein_grade", "n_extra_lines", "rank_value", "consequence", "site_feature", "protein_effect", "domains_hit",
              "U_gea", "U_protein", "U_expr", "U_regulatory", "U_snp_blind", "in_atac", "snp_cosig_2kb", "best_r2_shortread", "expr_fold", "fold", "p_emp", "size"):
        if c in best.columns:
            G[f"best_{c}" if not c.startswith("best_") else c] = best[c]
    G["n_units_regulatory"] = U.groupby("target_gene").U_regulatory.sum()
    G["n_units_expr"] = U.groupby("target_gene").U_expr.sum()
    G["n_units_snp_blind"] = U.groupby("target_gene").U_snp_blind.sum()
    G["n_units_protein"] = U.groupby("target_gene").U_protein.sum()
    G = G.reset_index()
    G["best_n_extra_lines"] = G.best_n_extra_lines.fillna(0).astype(int)
    G = G.sort_values(["best_rank_value", "best_protein_grade", "best_n_extra_lines", "L_gea", "theme_any", "gea_n_clusters", "gea_nlp"],
                      ascending=[True, False, False, False, False, False, False]).reset_index(drop=True)
    G["rank"] = np.arange(1, len(G) + 1)
    G.insert(1, "list_version", "v2_events_cr50_GEA_only")
    G.to_csv(f"{V2}/v2_candidate_list.csv", index=False)
    print(f"\nranked genes: {len(G):,} | best handle: {G.best_ftier.value_counts().to_dict()}")
    print(f"wrote {V2}/v2_unit_evidence.csv, v2_candidate_list.csv")


if __name__ == "__main__":
    main()
