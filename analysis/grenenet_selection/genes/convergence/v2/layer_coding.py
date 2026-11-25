#!/usr/bin/env python
"""v2 layer -- what does each GEA unit do to the gene it lands in? (site effect, repeat-copy change, protein consequence, Pfam domain, grade)

Per UNIT, from its REPRESENTATIVE record (v2/units.py: the member with most carriers; for an SV event that is the allele the event is named
after). Reuses the old machinery unchanged -- build_repeat_variants.repeat_call (tandem repeat-copy change), build_sv_hit_dossier.site_effect /
tx_effect (effect on EVERY transcript of every gene the REF span touches), coding_evidence.consequence / nmd_call / pfam_for (premature stop and
nonsense-mediated decay call, fraction of the protein lost, Pfam domain overlap) -- but keyed on `unit`, with the allele taken from the unit,
not looked up by (chrom,pos,ref_len,alt_len).

protein_grade (as coding_evidence.py): 3 gene loss / PTC with NMD / >= 50% of protein lost; 2 >= 10% lost, a C-terminal change hitting a domain,
a coding SV, or an in-frame change >= 10 aa; 1 an in-frame / C-terminal change inside a domain; 0 else.  Flags, not filters: repeat_copy,
effect_ambiguous (isoforms disagree), uniprot_len_mismatch.

Writes results/candidate_list_v2_events_cr50/v2_layer_coding.csv (one row per pooled unit; coding fields empty where the unit touches no CDS) and
its own Pfam cache (v2/ results copy of the old cache, extended) -- the old coding_pfam_cache.json is read, never written.
env: kmate. Compute node; outbound HTTPS for UniProt / InterPro on genes not yet cached.
"""
from __future__ import annotations
import os, sys, json, shutil
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, CONV); sys.path.insert(0, os.path.abspath(os.path.join(CONV, "..", "..")))
import units as UN                                              # noqa: E402
import build_repeat_variants as BRV                             # noqa: E402  (repeat_call; .SC.load_ref; .DS = build_sv_hit_dossier)
import coding_evidence as CE                                    # noqa: E402
V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")
OLD_CACHE = os.path.join(CONV, "results", "coding_pfam_cache.json")
CE.CACHE = os.path.join(V2, "coding_pfam_cache.json")           # never touch the old cache file
if not os.path.exists(CE.CACHE) and os.path.exists(OLD_CACHE):
    shutil.copy(OLD_CACHE, CE.CACHE)
DS = BRV.DS


def main():
    P = pd.read_csv(f"{V2}/gea_pool.csv")[["unit", "vclass"]]
    A = UN.rep_alleles(P.unit.values).merge(P, on="unit")
    A["REF"] = A.ref.str.upper(); A["ALT"] = A.alt.str.upper()
    print(f"{len(A):,} units (rep alleles loaded); loading reference + TAIR10 features ...", flush=True)
    ref = BRV.SC.load_ref()
    iv, tx, by_gene = DS.load_tair_features()

    # repeat-copy change (indels/SVs only; an MNP changes no length)
    rc = []
    for r in A.itertuples():
        rc.append(BRV.repeat_call(r, ref) if r.vclass != "mnp" else {})
    RC = pd.DataFrame(rc, index=A.index)
    EF = pd.DataFrame([DS.site_effect(r, r.ALT, ref, iv, tx, by_gene) for r in A.itertuples()], index=A.index)
    D = pd.concat([A.drop(columns=["ref", "alt"]), RC, EF], axis=1)
    print(f"site effect on a transcript: {int((D.site_feature != '').sum()):,} of {len(D):,} units; "
          f"CDS-affecting: {int(D.site_feature.astype(str).str.startswith('CDS').sum()):,}", flush=True)

    cod = D[D.site_feature.astype(str).str.startswith("CDS")].copy()
    C = pd.DataFrame([CE.consequence(r, tx.get(r.site_tx)) for r in cod.itertuples()], index=cod.index)
    cod = pd.concat([cod, C], axis=1)
    genes = sorted({g for g in cod.site_effect_gene.dropna() if str(g).startswith("AT")})
    Pf = CE.pfam_for(genes)
    hit, names, mism = [], [], []
    for r in cod.itertuples():
        p = Pf.get(r.site_effect_gene)
        if not p:
            hit.append(np.nan); names.append(""); mism.append(np.nan); continue
        wt = pd.to_numeric(r.wt_aa, errors="coerce")
        mism.append(bool(np.isfinite(wt) and int(wt) != p["length"]))
        if not np.isfinite(r.aff_lo):
            hit.append(np.nan); names.append(""); continue
        h = [d for d in p["domains"] if d[3] >= r.aff_lo and d[2] <= r.aff_hi]
        hit.append(bool(h)); names.append(";".join(sorted({d[1] for d in h})))
    cod["domain_hit"], cod["domains_hit"], cod["uniprot_len_mismatch"] = hit, names, mism
    fl, dh, cq = cod.protein_frac_lost.fillna(0), cod.domain_hit.fillna(False).astype(bool), cod.consequence
    inframe_aa = (cod.alt_len - cod.ref_len).abs() / 3
    cod["protein_grade"] = np.select([cq.isin(["gene_loss", "ptc_nmd"]) | (fl >= 0.5),
                                      (fl >= 0.10) | ((cq == "cterm") & dh) | (cq == "coding_sv") | ((cq == "inframe") & (inframe_aa >= 10)),
                                      dh & cq.isin(["inframe", "cterm"])], [3, 2, 1], 0)
    for c in ("consequence", "nmd", "nmd_note", "protein_frac_lost", "aff_lo", "aff_hi", "domain_hit", "domains_hit", "uniprot_len_mismatch", "protein_grade"):
        D[c] = cod[c].reindex(D.index)
    D["repeat_copy"] = D.get("repeat_copy_change", pd.Series(False, index=D.index)).astype(str).str.lower().eq("true")
    D["effect_ambiguous"] = D.isoforms_disagree.astype(str).str.lower().eq("true")
    keep = ["unit", "chrom", "pos", "ref_len", "alt_len", "vclass", "n_members", "n_carriers", "site_effect_gene", "site_gene_type", "site_tx", "site_feature",
            "protein_effect", "protein_pos", "wt_aa", "mut_aa", "n_isoforms", "isoforms_disagree", "n_genes_hit", "isoform_effects", "repeat_copy",
            "effect_ambiguous", "consequence", "nmd", "nmd_note", "protein_frac_lost", "domain_hit", "domains_hit", "uniprot_len_mismatch", "protein_grade"]
    keep += [c for c in ("unit_len", "n_units_changed", "repeat_copy_change") if c in D.columns]
    D[[c for c in keep if c in D.columns]].to_csv(f"{V2}/v2_layer_coding.csv", index=False)
    print("\nconsequence:", cod.consequence.replace("", "(none)").value_counts().to_dict())
    print("protein_grade:", cod.protein_grade.value_counts().sort_index().to_dict())
    print("by class (CDS-affecting):", cod.vclass.value_counts().to_dict(), "| grade >= 2 by class:", cod[cod.protein_grade >= 2].vclass.value_counts().to_dict())
    print(f"wrote {V2}/v2_layer_coding.csv")


if __name__ == "__main__":
    main()
