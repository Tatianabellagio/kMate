#!/usr/bin/env python
"""Cargo as a line of evidence: what an INSERTION brings with it.

Every other line here is reference-based and therefore describes the landing site. For an
insertion that is half the question -- the inserted sequence is absent from TAIR10, so no
ATAC peak, no 1001T eQTL and no gene model covers it. The `cargo/` tree characterises it
from the sequence itself (`cargo/tfbs_insertion_payload.py`: FIMO sites in each inserted
sequence, each locus against its OWN dinucleotide-shuffled null, 100 shuffles). This script
reads those results, applies the repeat diagnostics the promoter-turnover work needed, and
turns the survivors into one column the scorer can use.

Read-only across the boundary: reads `cargo/results/tfbs_payload/*`, writes only
`results/cargo_evidence.csv`.

What the line requires, and why each clause is there:

  1. `emp_p <= 0.05` against the locus' own shuffled null. The group result is an excess
     with NO locus surviving multiple testing (24 of 193 at p <= 0.05 vs 9.7 expected, BH
     q <= 0.20 retains none), so this is a candidate filter, not a test -- same standing as
     every other line here.
  2. Not repeat-driven. `hits_per_motif` (FIMO hits / distinct motifs) and the tandem
     fraction of the inserted sequence, the two measures that separated AT4G11800's
     307-sites-from-27-motifs artefact from GPX6's real turnover (repeat_context.py). The
     enriched cargo loci pass comfortably -- 1.0-3.7 hits per motif against AT4G11800's
     11.4 -- so this clause currently removes nothing. It is kept because it is the check
     that would catch the next one.
  3. Not carried by an artefact family. BBR-BPC enrichment is the GA-repeat artefact
     (documented in METHODS_FUNCTIONAL_TRACKS.md). The telomeric 7-mer artefact is already
     removed upstream by the payload's `long_repeat_tract` filter.
  4. Lands somewhere a cis-element could act: promoter / 5'UTR / UTR, or an ATAC peak.
     Cargo enrichment in a gene desert is a sequence fact with no route to a phenotype.

`cargo_note` records which families drive it, so ERF/GATA-led loci stay visible: the
dinucleotide shuffle controls composition but not clustered GC, and the GCC-box is GC-rich
in an AT-rich genome, so an ERF-led enrichment is weaker evidence than a bZIP-led one.

env: kmate. Compute node.
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tfbs_turnover as TF                                       # noqa: E402
import sv_content as SC                                          # noqa: E402
import repeat_context as RC                                      # noqa: E402

OUT = f"{HERE}/results"
CARGO = f"{os.path.dirname(os.path.dirname(HERE))}/cargo/results/tfbs_payload"
KEY = ["chrom", "pos", "ref_len", "alt_len"]          # never chrom:pos alone (multiallelic)

ARTEFACT_FAMILIES = {"BBR-BPC"}                       # GA repeats
CAUTION_FAMILIES = {"ERF", "GATA"}                    # clustered GC, not controlled by the shuffle
HITS_PER_MOTIF_MAX = 4.0                              # repeat_context.SPM_MAX
REGULATORY = ("promoter", "UTR")


def main():
    L = pd.read_csv(f"{CARGO}/payload_loci.csv")
    S = pd.read_csv(f"{CARGO}/payload_sites.csv")
    F = pd.read_csv(f"{CARGO}/payload_families.csv")
    A = pd.read_csv(f"{OUT}/evidence_matrix.csv")

    rows = []
    for r in L.itertuples():
        d = S[S.locus_id == r.locus_id]
        hpm = len(d) / d.motif_id.nunique() if len(d) and d.motif_id.nunique() else np.nan
        ref, alt = TF.fetch_allele(r.chrom, int(r.pos), int(r.ref_len), int(r.alt_len))
        # strip the prefix AND suffix shared with REF: a naive alt[len(ref):] is wrong for a
        # record that is not left-anchored (Chr3:10,204,785 is REF=TG -> ALT starting AA, so
        # the whole 187 bp ALT is new, not the last 185)
        ins = SC.inserted_part(ref, alt) if ref else ""
        tandem = RC.tandem_fraction(ins) if len(ins) > 20 else 0.0
        fam = F[F.locus_id == r.locus_id].sort_values("enrich", ascending=False)
        top = fam.tf_family.iloc[0] if len(fam) else ""
        lead = ", ".join(f"{x.tf_family} {x.enrich:.1f}x" for x in fam.head(3).itertuples())
        rows.append(dict(
            chrom=r.chrom, pos=int(r.pos), ref_len=int(r.ref_len), alt_len=int(r.alt_len),
            locus_id=r.locus_id, cargo_bp=int(r.ins_bp), cargo_class=r.content_class,
            cargo_repeat_family=r.repeat_family, cargo_sites=int(r.n_sites),
            cargo_motifs=int(r.n_motifs), cargo_enrich=round(float(r.enrich), 2),
            cargo_p=float(r.emp_p), cargo_hits_per_motif=round(hpm, 2),
            cargo_tandem=round(tandem, 2), cargo_top_family=top, cargo_note=lead,
            cargo_site_cat=r.site_cat, cargo_atac=bool(r.atac_peak)))
    C = pd.DataFrame(rows)

    C["cargo_repeat_driven"] = ((C.cargo_hits_per_motif >= HITS_PER_MOTIF_MAX)
                               | (C.cargo_tandem >= RC.TANDEM_MIN))
    C["cargo_artefact_family"] = C.cargo_top_family.isin(ARTEFACT_FAMILIES)
    C["cargo_gc_caution"] = C.cargo_top_family.isin(CAUTION_FAMILIES)
    landed = C.cargo_site_cat.str.startswith(REGULATORY) | C.cargo_atac
    # WITHDRAWN 2026-09-20. Not "suspended pending recomputation" -- the recomputation ran,
    # and then a length-matched control killed the test itself. Three steps:
    #
    #  1. The payload statistic counts FIMO HITS, which are not independent: at
    #     Chr5:19,636,028 one 44 bp ABRE draws 59 bZIP hits from 23 motifs, and a
    #     dinucleotide shuffle destroys clustered elements so it rarely piles hits that way.
    #  2. `cargo_sites_null.py` redid all 193 loci on merged SITES. 25/193 reached
    #     p <= 0.05 family-merged -- but the rate is a function of LENGTH, not biology:
    #     4.9% below 200 bp (nominal) against 63% above 2 kb, at enrichments of only
    #     1.14-1.19x with p pinned to the 1/101 floor.
    #  3. `cargo_null_control.py` ran the identical statistic on 1,158 SIGNAL-FREE
    #     length-matched windows (3 random genomic + 3 TE-overlapping per locus). Any
    #     p <= 0.05 there is a false positive by construction, and the false-positive rate
    #     is 9.5% -> 70.2% across the same length bins on fam_sites_p (11.1% -> 64.9% on
    #     the original hits_p). The dinucleotide shuffle preserves composition and nothing
    #     else, so its bias is per-base and roughly constant while its sampling noise
    #     shrinks with length -- a fixed ~15% excess is unbeatable at 5 kb and invisible at
    #     150 bp. The 25 survivors are that, not cargo.
    #
    # Head-to-head, real insertions are at or BELOW the false-positive rate of signal-free
    # sequence in every length bin of every statistic (pooled fam_sites_p: real 25/193 =
    # 13.0% vs control 16.5%; no bin has a one-sided Fisher p < 0.43). Only all-motif-merged
    # `any_sites_p` is calibrated at all lengths (control FPR 2.8%, if anything
    # conservative) -- and on that one the real loci give 5/193 = 2.6% against a 2.8%
    # control rate. So there is no cargo TF-binding-site enrichment to score, at any length,
    # on any of the three counting conventions.
    #
    # This line therefore scores 0 everywhere and is NOT conditional on the null files
    # existing -- an earlier version re-enabled itself as soon as cargo_sites_null.csv
    # appeared, which would have scored candidates on the statistic the control invalidated.
    # The cargo_* description columns stay: the sequence facts are still facts, and the
    # merged-site columns are carried for auditing. What would revive this line is a
    # different control, not a different threshold -- the honest null for non-reference
    # insertion sequence is the OTHER panel SV insertions that are not GEA/GWAS hits, i.e.
    # signal-free sequence of the same class, rather than reference windows.
    C["cargo_hits_p"] = C.cargo_p
    C["L_cargo"] = 0
    _ = landed                     # kept: the regulatory-landing clause, if the line revives
    for _f, _cols in ((f"{OUT}/cargo_sites_null.csv",
                       ["fam_sites_enrich", "fam_sites_p", "any_sites_enrich", "any_sites_p"]),):
        if os.path.exists(_f):     # descriptive only -- never feeds L_cargo
            C = C.merge(pd.read_csv(_f)[["locus_id"] + _cols], on="locus_id", how="left")

    # store_row on the full four-part key: one position can carry several records
    C = C.merge(A[KEY + ["store_row", "target_gene", "mode", "n_lines"]], on=KEY, how="left")
    C.to_csv(f"{OUT}/cargo_evidence.csv", index=False)

    print(f"{len(C)} insertions with cargo; {int(C.store_row.notna().sum())} join the matrix")
    print(f"  enriched (p<=0.05): {int((C.cargo_p <= 0.05).sum())}"
          f"  repeat-driven {int(C.cargo_repeat_driven.sum())}"
          f"  artefact-family {int(C.cargo_artefact_family.sum())}"
          f"  GC-caution {int(C.cargo_gc_caution.sum())}")
    print(f"  L_cargo = 1: {int(C.L_cargo.sum())}")
    k = C[C.L_cargo == 1].sort_values("cargo_enrich", ascending=False)
    cols = ["locus_id", "target_gene", "mode", "cargo_bp", "cargo_class", "cargo_enrich",
            "cargo_p", "cargo_atac", "n_lines", "cargo_note"]
    print(k[cols].to_string(index=False))


if __name__ == "__main__":
    main()
