#!/usr/bin/env python
"""v2 layer -- is each GEA unit visible to SNPs? Two separate questions, both keyed on `unit`.

1. snp_cosig_2kb  -- a genome-wide Bonferroni SNP hit (any of the 22 axes, tile_cr50 snp scan, MAF > 0.05) within 2 kb of the unit's
   representative position: statistical shadowing (the locus would also have been found by a SNP study). Same method as build_snp_tagging.py.
2. best_r2_snp    -- the best founder-panel r2 between the unit and any panel SNP within +-50 kb (masked, MAC >= 2): genetic detectability,
   independent of significance. snp_blind = best r2 < 0.2; NaN r2 = untestable, which is NOT the same as untagged.
   Two SNP sources: the arch3 PANEL SNPs (best_r2_snp) and the GrENE-Net SHORT-READ SNPs (best_r2_shortread). The panel SNPs near an SV are genotyped
   from the same graph haplotypes as the SV, so a perfect panel r2 can be partly built in; the short-read r2 is the independent measure of whether a
   SNP study would see the unit (snp_blind_shortread).
   * SV EVENTS: r3_persite_gwas/results/sv_snp_ld_v2/tagging_svevent_panel_{Chr}.npz (build_tagging_events.py: anchor = union over member
     alleles). Joined EXACTLY on (chrom, rep_col) -- the panel column of the event's representative record.
   * short indels / MNPs: the record-level tagging_indel_panel_{Chr}.npz (one record = one unit), joined on (chrom, pos, |dlen|) because that
     file carries no record id; tag_join_ambiguous flags positions where several records collapse onto that key.
Output: results/candidate_list_v2_events_cr50/v2_layer_snp_visibility.csv
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, CONV); sys.path.insert(0, os.path.abspath(os.path.join(CONV, "..", "..")))
import units as UN                                              # noqa: E402
import build_snp_tagging as bs                                  # noqa: E402  (sig_snps, add_cosig: pure functions)
import lib                                                      # noqa: E402
V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")
TAGDIR = f"{lib.GEA}/r3_persite_gwas/results/sv_snp_ld_v2"

P = pd.read_csv(f"{V2}/gea_pool.csv")[["unit", "chrom", "pos", "vclass"]]
d = UN.info(P.unit.values)                                      # rep_col, dlen, ...
d = d.merge(P[["unit", "vclass"]], on="unit")
print("collecting Bonferroni SNP hits, 22 axes ...", flush=True)
S = bs.sig_snps()
d = bs.add_cosig(d, S)

def tag_tables(mode):
    ev, rec = [], []
    for ch in UN.CHROMS:
        z = np.load(f"{TAGDIR}/tagging_svevent_{mode}_{ch}.npz", allow_pickle=True)
        ev.append(pd.DataFrame({"chrom": ch, "rep_col": z["rep_col"].astype(np.int64), "best_r2_snp": z["best_r2"], "best_r2_snp_pos": z["best_snp_pos"],
                                "tag_n_snp_window": z["n_snp_window"], "tag_an": z["an"], "tag_maf": z["maf"]}))
        z = np.load(f"{TAGDIR}/tagging_indel_{mode}_{ch}.npz", allow_pickle=True)
        rec.append(pd.DataFrame({"chrom": ch, "pos": z["pos"].astype(np.int64), "size": z["size"].astype(np.int64), "best_r2_snp": z["best_r2"],
                                 "best_r2_snp_pos": z["best_snp_pos"], "tag_n_snp_window": z["n_snp_window"], "tag_an": z["an"], "tag_maf": z["maf"]}))
    EV = pd.concat(ev, ignore_index=True); assert not EV.duplicated(["chrom", "rep_col"]).any()
    RC = pd.concat(rec, ignore_index=True).sort_values("best_r2_snp", ascending=False)
    n = RC.groupby(["chrom", "pos", "size"]).size().rename("tag_n_records")
    RC = RC.groupby(["chrom", "pos", "size"]).first().join(n).reset_index(); RC["tag_join_ambiguous"] = RC.tag_n_records > 1
    is_sv = d.vclass.eq("sv")
    A = d[is_sv].merge(EV, on=["chrom", "rep_col"], how="left"); A["tag_source"] = "svevent rep_col"; A["tag_join_ambiguous"] = False
    B = d[~is_sv].assign(size=lambda x: x.dlen.abs()).merge(RC, on=["chrom", "pos", "size"], how="left"); B["tag_source"] = "indel record pos+size"
    keep = ["unit", "snp_cosig_2kb", "snp_cosig_dist", "snp_cosig_axes", "best_r2_snp", "best_r2_snp_pos", "tag_n_snp_window", "tag_an", "tag_maf",
            "tag_join_ambiguous", "tag_source"]
    return pd.concat([A[keep], B[keep]], ignore_index=True).sort_values("unit").reset_index(drop=True)


L = tag_tables("panel")
Sh = tag_tables("shortread")
assert (L.unit.values == Sh.unit.values).all()
L["best_r2_shortread"] = Sh.best_r2_snp.values
L["best_r2_shortread_pos"] = Sh.best_r2_snp_pos.values
L["tag_shortread_untestable"] = Sh.best_r2_snp.isna().values
L["snp_blind_shortread"] = (Sh.best_r2_snp < 0.2).values
L["tag_untestable"] = L.best_r2_snp.isna()
L["snp_blind"] = L.best_r2_snp < 0.2
assert L.unit.is_unique and len(L) == len(P)
L.to_csv(f"{V2}/v2_layer_snp_visibility.csv", index=False)

Lj = L.merge(P[["unit", "vclass"]], on="unit")
print(f"\nwrote {V2}/v2_layer_snp_visibility.csv ({len(L):,} units)")
for c in ("sv", "smallindel", "mnp"):
    x = Lj[Lj.vclass == c]
    print(f"  {c:11s} n={len(x):4d} | SNP-shadowed (cosig<=2kb) {100*x.snp_cosig_2kb.mean():5.1f}% | untestable tagging {100*x.tag_untestable.mean():5.1f}% "
          f"| median r2 panel {x.best_r2_snp.median():.2f} / short-read {x.best_r2_shortread.median():.2f} | SNP-blind panel {int(x.snp_blind.sum())} / short-read {int(x.snp_blind_shortread.sum())} | ambiguous join {int(x.tag_join_ambiguous.sum())}")
t = Lj[Lj.vclass == "sv"]
print("\nSV events: shadowed x tagged/blind")
print(pd.crosstab(t.snp_cosig_2kb.map({True: "shadowed", False: "SNP-unique"}), t.snp_blind_shortread.map({True: "SNP-blind (short-read r2<0.2)", False: "taggable"}).where(~t.tag_shortread_untestable, "untestable")).to_string())
