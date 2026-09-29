#!/usr/bin/env python
"""Gene-by-gene review of the themed shortlist: one locus figure + one dossier each.

For every gene it renders `dissection/plot_locus_combined.py`'s 3-panel figure, parses
the **LD-confirm** the plotter prints, and writes a dossier section pulling together
everything the master table knows: how it was found, on which axes/gardens, its blocks,
its variants and their mechanism, whether SNPs shadow or tag it, and the curated
UniProt function.

The verdict follows the triage rule this project arrived at independently in
`dissection/results/loci/CANDIDATE_VERDICTS.md`: **recurrence + founder frequency +
placement, never p-value.** Three criteria, counted:

  C1 recurrent   >=3 climate clusters or >=3 gardens, OR found by both scans
  C2 common      >=20 of 231 founders carry the lead allele
  C3 direct      the strongest mechanism is CDS / UTR / promoter (not intron-only
                 or downstream-flank)

  ROBUST 3 met | MODERATE 2 | WEAK 1 | FRAGILE 0

**Why LD is NOT the verdict here.** An earlier version of this script ranked on the
LD-confirm (founder r2 between the lead and the gene's own variants), which is the right
test for a *block-attributed* candidate -- it is what caught the CARK case, where the
lead sat 29 kb from the gene it was attributed to. It does not transfer to this
pipeline: every gene here is assigned at its variant's own position, so all 23 leads are
physically inside their target gene and detachment is impossible by construction. Worse,
the LD rule inverts the evidence for promoter variants -- it graded GPX6 WEAK on
r2=0.24, where `CANDIDATE_VERDICTS.md` graded the same variant at the same r2 ROBUST,
correctly, because a promoter variant sitting on its own haplotype background is a
BETTER-localized causal candidate, not a worse one.

So the LD numbers are kept as an annotation, not a verdict:

  r2_gene    lead vs the gene's own variants -- high means the lead tags the gene's
             common haplotype, so the signal could belong to any variant on it
  r2_local   the best r2 the lead reaches anywhere in the window
  ld_ambiguous  r2_local >= 0.5 and more than double r2_gene: a tightly-linked
             alternative in the window could carry the signal instead. This is about
             WHICH VARIANT, not which gene.

These sort the list. They do not replace reading the figure.

Outputs -> results/
  GENE_DOSSIERS.md          one section per gene, verdict + evidence + function
  gene_review.csv           the same as a table (LD numbers + verdict), for sorting
  plots/loci/<sym>_combined.{png,pdf}

env: kmate.  Compute node -- needs bcftools and the founder panel VCFs. ~20 s per gene.
"""
from __future__ import annotations
import os
import re
import sys
import json
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = f"{HERE}/results"
PLOTTER = os.path.join(os.path.dirname(HERE), "dissection", "plot_locus_combined.py")
SRC_LOCI = os.path.join(os.path.dirname(HERE), "dissection", "loci")
DST_LOCI = f"{OUT}/plots/loci"
PY = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/python"

LD_RE = re.compile(
    r"lead carriers=(\d+).*?max=([\d.]+)\s+med=([\d.]+).*?max r² local=([\d.]+)")


DIRECT_FTIERS = {"F1_CDS_frameshift", "F2_CDS_inframe", "F3_UTR", "F4_promoter"}


def verdict(r, carriers, r2_gene, r2_local):
    """Recurrence + frequency + placement. LD is annotation only -- see module docstring."""
    recurrent = (max(int(r.gea_n_clusters or 0), int(r.gwas_n_gardens or 0)) >= 3
                 or r.found_by == "GEA+GWAS")
    common = carriers >= 20
    direct = str(r.best_ftier) in DIRECT_FTIERS
    met = int(recurrent) + int(common) + int(direct)
    verd = {3: "ROBUST", 2: "MODERATE", 1: "WEAK", 0: "FRAGILE"}[met]
    why = ", ".join([("recurrent" if recurrent else "not recurrent"),
                     (f"common ({carriers} carriers)" if common
                      else f"rare ({carriers} carriers)"),
                     (f"direct placement ({r.best_ftier})" if direct
                      else f"indirect placement ({r.best_ftier})")])
    amb = (np.isfinite(r2_local) and np.isfinite(r2_gene)
           and r2_local >= 0.5 and r2_gene > 0 and r2_local > 2 * r2_gene)
    return verd, why, bool(amb)


def config_for(r, V):
    ch, pos = r.lead_variant.split(":")
    pos = int(pos)
    v = V[(V.chrom == ch) & (V.pos == pos)]
    if not len(v):
        return None
    v = v.iloc[0]
    sym = r.symbol if isinstance(r.symbol, str) and r.symbol else r.gene
    axis = r.gea_best_axis if isinstance(r.gea_best_axis, str) and r.gea_best_axis \
        else "bio1"
    return {"gene": r.gene, "sym": re.sub(r"[^A-Za-z0-9_.-]", "_", sym),
            "chrom": ch, "gstart": int(r.gene_start), "gend": int(r.gene_end),
            "vpos": pos, "ref_len": int(v.ref_len), "alt_len": int(v.alt_len),
            "axis": axis, "pad": 8000,
            "region": str(v.region) if isinstance(v.region, str) else "intergenic"}


def fmt(x, nd=2, dash="—"):
    return dash if x is None or (isinstance(x, float) and not np.isfinite(x)) \
        else f"{x:.{nd}f}"


def dossier(r, cfg, ld, verd, why, amb, vars_here):
    sym = cfg["sym"]
    L = []
    L.append(f"## {sym} ({r.gene}) — **{verd}**\n")
    L.append(f"*{why}.*\n")
    prot = r.protein_name if isinstance(r.protein_name, str) else "—"
    L.append(f"- **Protein** — {prot}")
    L.append(f"- **Themes** — {r.themes}")
    L.append(f"- **Location** — {r.chrom}:{int(r.gene_start):,}–{int(r.gene_end):,} "
             f"({r.strand}) · locus `{r.locus}` holding **{int(r.locus_n_genes)}** "
             f"candidate genes")
    L.append(f"- **Found by** — {r.found_by}"
             + (f" · GEA {int(r.gea_n_clusters)} climate clusters / "
                f"{int(r.gea_n_axes)} axes ({r.gea_sig_clusters}), "
                f"best {r.gea_best_axis} nlp {fmt(r.gea_nlp, 1)}, λ_min "
                f"{fmt(r.gea_min_lam)}" if r.gea_n_clusters else "")
             + (f" · GWAS {int(r.gwas_n_gardens)} gardens ({r.gwas_gardens}), "
                f"nlp {fmt(r.gwas_nlp, 1)}, MAC {fmt(r.gwas_min_mac, 0)}"
                if r.gwas_n_gardens else ""))
    L.append(f"- **Evidence** — {int(r.n_independent)} independent line(s): "
             + ", ".join([n for n, f in (("GEA recurrence", r.L_gea),
                                         ("GWAS gardens", r.L_gwas),
                                         ("both scans", r.L_cross)) if f] or ["—"]))
    L.append(f"- **Blocks** — {r.blocks}"
             + ("  ⚠ block attribution would name a different gene"
                if r.block_gene_mismatch else ""))
    L.append(f"- **Variants** — {int(r.n_variants)} "
             f"({int(r.n_sv)} SV / {int(r.n_smallindel)} indel / {int(r.n_mnp)} MNP), "
             f"max {int(r.max_size)} bp, in {r.regions}; "
             f"strongest mechanism **{r.best_ftier}**"
             + (" · **frameshift**" if r.has_frameshift else "")
             + (" · SV in CDS" if r.has_cds_sv else ""))
    L.append(f"- **Lead** — {r.lead_variant} ({r.lead_vclass}, {r.lead_size} bp, "
             f"{cfg['region']}), {ld['carriers']} founder carriers")
    L.append(f"- **LD (annotation, not the verdict)** — r²(lead↔gene) max "
             f"**{fmt(ld['r2_gene'])}** / med {fmt(ld['r2_med'])} · best local r² "
             f"{fmt(ld['r2_local'])}"
             + ("  ⚠ **LD-ambiguous**: a tightly-linked alternative in the window "
                "could carry the signal instead (about which VARIANT, not which gene)"
                if amb else ""))
    L.append(f"- **SNPs** — "
             + ("a Bonferroni SNP within 2 kb" if r.snp_cosig_any
                else "no Bonferroni SNP within 2 kb (SNP-unique)")
             + " · tagging r² "
             + ("not testable" if r.tag_untestable_all
                else f"{fmt(r.best_r2_snp)}"
                     + (" (**SNP-blind**)" if r.snp_blind_any else "")))
    if r.size_inferred_any:
        L.append("- ⚠ **Allele inferred** — a variant here sits at a multiallelic "
                 "position; its length was guessed as the largest record, so class "
                 "and any frameshift call are unreliable.")
    fn = r.uniprot_function
    if isinstance(fn, str) and fn.strip() and fn.lower() != "nan":
        L.append(f"- **Function (UniProt)** — {fn.strip()[:600]}")
    kw = r.uniprot_keywords
    if isinstance(kw, str) and kw.strip() and kw.lower() != "nan":
        L.append(f"- **Keywords** — {kw.strip()[:300]}")
    L.append(f"\n![{sym}](plots/loci/{sym}_combined.png)\n")
    if len(vars_here):
        L.append("| variant | class | size | region | mechanism | GEA nlp | "
                 "clusters | gardens |")
        L.append("|---|---|---|---|---|---|---|---|")
        for v in vars_here.rename(columns={"size": "vsize"}).itertuples():
            L.append(f"| {v.chrom}:{int(v.pos):,} | {v.vclass} | {int(v.vsize)} | "
                     f"{v.region} | {v.ftier if isinstance(v.ftier,str) else '—'} | "
                     f"{fmt(v.gea_nlp,1)} | {fmt(v.gea_n_clusters,0)} | "
                     f"{fmt(v.gwas_n_gardens,0)} |")
        L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-shortlist", action="store_true",
                    help="all 77 gated genes instead of the themed 23")
    ap.add_argument("--reuse-ld", action="store_true",
                    help="reuse LD numbers and figures from a previous run "
                         "(skips re-rendering; ~20 s/gene saved)")
    a = ap.parse_args()

    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv")
    V = pd.read_csv(f"{OUT}/variants_classified.csv")
    F = pd.read_csv(f"{OUT}/functional_variants.csv")

    sel = M[M.n_independent >= 1]
    if not a.all_shortlist:
        sel = sel[sel.themes.notna() & sel.themes.ne("")]
    sel = sel.sort_values(["n_independent", "gea_n_clusters", "gwas_n_gardens"],
                          ascending=False)
    os.makedirs(DST_LOCI, exist_ok=True)
    CACHE = pd.DataFrame()
    if a.reuse_ld and os.path.exists(f"{OUT}/gene_review.csv"):
        CACHE = pd.read_csv(f"{OUT}/gene_review.csv").set_index("gene")
        print(f"reusing LD for {len(CACHE)} genes from gene_review.csv")
    print(f"reviewing {len(sel)} genes\n")

    rows, secs = [], []
    for i, r in enumerate(sel.itertuples(), 1):
        cfg = config_for(r, V)
        if cfg is None:
            print(f"[{i}/{len(sel)}] {r.gene}: no lead variant, skipped")
            continue
        print(f"[{i}/{len(sel)}] {cfg['sym']} ({r.gene}) ...", flush=True)
        if a.reuse_ld and r.gene in CACHE.index:
            c = CACHE.loc[r.gene]
            ld = {"carriers": int(c.lead_carriers), "r2_gene": float(c.r2_gene_max),
                  "r2_med": float(c.r2_gene_med), "r2_local": float(c.r2_local_max)}
        else:
            p = subprocess.run([PY, PLOTTER, json.dumps(cfg)],
                               capture_output=True, text=True)
            m = LD_RE.search(p.stdout or "")
            ld = {"carriers": int(m.group(1)) if m else 0,
                  "r2_gene": float(m.group(2)) if m else np.nan,
                  "r2_med": float(m.group(3)) if m else np.nan,
                  "r2_local": float(m.group(4)) if m else np.nan}
            if not m:
                print(f"    no LD-confirm parsed; stderr: {(p.stderr or '')[-300:]}")
            for ext in (".png", ".pdf"):
                src = f"{SRC_LOCI}/{cfg['sym']}_combined{ext}"
                if os.path.exists(src):
                    shutil.copy2(src, f"{DST_LOCI}/{cfg['sym']}_combined{ext}")
        verd, why, amb = verdict(r, ld["carriers"], ld["r2_gene"], ld["r2_local"])
        print(f"    carriers={ld['carriers']} r2_gene={fmt(ld['r2_gene'])} "
              f"local={fmt(ld['r2_local'])} -> {verd}"
              + ("  [LD-ambiguous]" if amb else ""))

        vh = F[F.target_gene.eq(r.gene)].sort_values("gea_nlp", ascending=False)
        secs.append(dossier(r, cfg, ld, verd, why, amb, vh.head(12)))
        rows.append({"gene": r.gene, "symbol": cfg["sym"], "verdict": verd,
                     "reason": why, "themes": r.themes, "found_by": r.found_by,
                     "r2_gene_max": ld["r2_gene"], "r2_gene_med": ld["r2_med"],
                     "r2_local_max": ld["r2_local"], "lead_carriers": ld["carriers"],
                     "ld_ambiguous": amb,
                     "n_independent": r.n_independent,
                     "gea_n_clusters": r.gea_n_clusters,
                     "gwas_n_gardens": r.gwas_n_gardens,
                     "best_ftier": r.best_ftier, "lead_variant": r.lead_variant,
                     "locus": r.locus, "locus_n_genes": r.locus_n_genes,
                     "block_gene_mismatch": r.block_gene_mismatch})

    R = pd.DataFrame(rows)
    ORD = {"ROBUST": 0, "MODERATE": 1, "WEAK": 2, "FRAGILE": 3}
    R["_o"] = R.verdict.map(ORD)
    R = R.sort_values(["_o", "r2_gene_max"], ascending=[True, False]).drop(columns="_o")
    R.to_csv(f"{OUT}/gene_review.csv", index=False)

    order = {g: i for i, g in enumerate(R.gene)}
    secs = [s for _, s in sorted(zip([order.get(x.gene, 99) for x in sel.itertuples()
                                      if config_for(x, V)], secs),
                                 key=lambda t: t[0])]
    with open(f"{OUT}/GENE_DOSSIERS.md", "w") as fh:
        fh.write("# Gene-by-gene review — themed shortlist\n\n")
        fh.write(
            f"{len(R)} genes, one dossier each.\n\n"
            "⚠ **The `verdict` column is a 3-criterion count, not a judgement.** It "
            "counts whether the gene is recurrent (>=3 climate clusters or gardens, or "
            "found by both scans), common (>=20 of 231 founders carry the lead) and "
            "directly placed (CDS/UTR/promoter). Every threshold in that is arbitrary, "
            "and it visibly mis-ranks: CRK18 -- cross-scan, frameshift, r2=0.87 to its "
            "own gene, in the one locus that survived honest WZA recalibration -- is "
            "demoted to MODERATE purely by having 15 carriers instead of 20, while "
            "genes with r2~0 to their gene are promoted for being common. Use the "
            "component columns and the figures; treat the label as a sort key.\n\n"
            "**LD is deliberately not in the verdict.** Every gene here is assigned at "
            "its variant's own position, so all leads sit physically inside their "
            "target gene and CARK-type detachment is impossible by construction. An "
            "earlier version of this script ranked on LD and graded GPX6 WEAK at "
            "r2=0.24, where `dissection/results/loci/CANDIDATE_VERDICTS.md` graded the "
            "same variant at the same r2 ROBUST -- correctly, because a promoter "
            "variant on its own haplotype background is a better-localized causal "
            "candidate, not a worse one. `ld_ambiguous` flags where a tightly-linked "
            "alternative could carry the signal instead: that is about which "
            "VARIANT, not which gene.\n\n")
        fh.write("| gene | criteria met | why | r²(lead↔gene) | best local r² | "
                 "carriers | themes | found by | locus (n genes) |\n"
                 "|---|---|---|---|---|---|---|---|---|\n")
        for r in R.itertuples():
            fh.write(f"| {r.symbol} | **{r.verdict}** | {r.reason} | "
                     f"{fmt(r.r2_gene_max)}"
                     + (" ⚠" if r.ld_ambiguous else "")
                     + f" | {fmt(r.r2_local_max)} | {r.lead_carriers} | {r.themes} "
                     f"| {r.found_by} | {r.locus} ({int(r.locus_n_genes)}) |\n")
        fh.write("\n---\n\n")
        fh.write("\n---\n\n".join(secs))

    print(f"\n{R.verdict.value_counts().to_string()}")
    print(f"\nwrote {OUT}/GENE_DOSSIERS.md and {OUT}/gene_review.csv")


if __name__ == "__main__":
    main()
