#!/usr/bin/env python
"""Repeat / TE classification of every SV record (>= 50 bp) in the arch3 panel, and of every SV event.

Sources -- one masking protocol for all 227,851 records:
  graffite/raw/{chrom}/       GraffiTE's own repeatmask_VCF output for the 155,986 records it masked
  graffite/skipped/{chrom}/   the same steps (_submit_graffite_skipped.sh, verified row-identical to
                              GraffiTE on 300 of its own sequences) for the 71,865 records GraffiTE's
                              `bcftools view --types indels` skipped -- HPRC annotate_vcf COMPLEX
                              records and INS/DEL-labelled records that also swap a reference base
Both are RepeatMasker -s with combined_tair10_repbase.rm.lib, low-complexity masking ON (GraffiTE's
setting; it warns against -nolow), on the whole allele: ALT for insertions, REF for deletions.

Spans are computed from the RAW RepeatMasker .out, as the union of hit intervals over the allele:
  span_any      any repeat
  span_te       TE classes only (everything except Simple_repeat / Low_complexity / Satellite)
  span_simple   Simple_repeat + Low_complexity + Satellite
  span_graffite GraffiTE's own total_match_span (union of OneCode copies), kept for reference only.
                OneCodeToFindThemAll emits no RC/Helitron copies at all (41,555 raw Helitron hits on the
                GraffiTE-masked records, 0 in its output), so this span -- and GraffiTE's
                `pangenome.vcf`, which filters on it at 0.80 -- is blind to Helitrons.
`cargo` uses the thresholds of _rm_cargo_classify.py on span_te (te_full >= 0.80, te_most >= 0.50,
te_partial >= 0.20, te_trace > 0, non_te), so it is comparable to rm_cargo_allins.csv -- which differs in
having been run with -nolow -norna at default speed, on insertions only.

Not used: GraffiTE's TSD calls. Since upstream commit 3f062c1 (2025-03-26) prepTSD.sh names the
reference flanks by region (samtools faidx -r) while TSD_Match_v2.sh looks them up as `{id}__L/__R`, so
every TSD search ran on the SV's own unmasked edges with no flanking sequence (395 PASS of 42,837).

Replacement records (hprc_type COMPLEX, and INS/DEL with unshared_bp > 0): the classified sequence is the
whole ALT (insertion) or REF (deletion), which includes sequence the other allele replaces. For the
36,733 with > 50 unshared bp on both alleles, the class describes the longer allele only.

Events: arch3_event_map.npz (truvari collapse defaults). The representative record is the first by
(position, record) -- truvari --keep first. The 1,383 records of exactly 50 bp are outside the event map
(it uses > 50 bp) and have no event.

Env: basic. Writes results/sv_adaptive/graffite/sv_repeat_class_records.csv, sv_repeat_class_events.csv.
"""
import gzip, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = "/global/scratch/projects/fc_moilab/tbellg/kmate/panel/arch3"
OUT = f"{HERE}/results/sv_adaptive"
G = f"{OUT}/graffite"
CHROMS = [f"Chr{i}" for i in range(1, 6)]
NON_TE = {"Simple_repeat", "Low_complexity", "Satellite"}


def read_out(path, onecode=False):
    """RepeatMasker .out (or GraffiTE's OneCode .onecode.out) -> DataFrame of hits."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue                                    # header / blank lines
            rows.append((f[4], int(f[5]), int(f[6]), f[9], f[10], float(f[1])))
    return pd.DataFrame(rows, columns=["id", "qs", "qe", "family", "rmclass", "div"])


def union_bp(g):
    """total length of the union of 1-based closed intervals, per id"""
    s = g.sort_values(["id", "qs"])
    out = {}
    for k, d in s.groupby("id", sort=False):
        tot, cs, ce = 0, None, None
        for a, b in zip(d.qs.to_numpy(), d.qe.to_numpy()):
            if cs is None:
                cs, ce = a, b
            elif a <= ce + 1:
                ce = max(ce, b)
            else:
                tot += ce - cs + 1
                cs, ce = a, b
        out[k] = tot + (ce - cs + 1)
    return pd.Series(out, dtype=np.int64)


def unshared(r, a):
    """bases of the shorter allele not covered by the prefix/suffix it shares with the longer one"""
    s, L = (r, a) if len(r) < len(a) else (a, r)
    p = 0
    while p < len(s) and s[p] == L[p]:
        p += 1
    q = 0
    while q < len(s) - p and s[-1 - q] == L[-1 - q]:
        q += 1
    return len(s) - p - q


def main():
    # --- records: the exact set and alleles GraffiTE was given --------------------------------------
    # the HPRC annotate_vcf type (INS / DEL / COMPLEX) is the third field of the panel record ID; re-deriving
    # it from allele lengths disagrees on a few dozen records, so it is read from the panel itself
    vt = {}
    for i, ch in enumerate(CHROMS, 1):
        n = 0
        with gzip.open(f"{PANEL}/chr{i}/merged_231_chr{i}_final.vcf.gz", "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 5)
                if "," not in f[4] and "<" not in f[4] and abs(len(f[4]) - len(f[3])) >= 50:
                    vt[f"{ch}_{n}"] = f[2].split("-")[2]
                n += 1
    recs = []
    for line in open(f"{G}/sv_sites.vcf"):
        if line[0] == "#":
            continue
        c, p, i, r, a = line.split("\t", 5)[:5]
        recs.append((i, c, int(p), "ins" if len(a) > len(r) else "del", abs(len(a) - len(r)),
                     max(len(r), len(a)), vt[i], unshared(r.upper(), a.upper())))
    assert len(recs) == len(vt), (len(recs), len(vt))
    R = pd.DataFrame(recs, columns=["id", "chrom", "pos", "kind", "size_bp", "seq_len", "hprc_type",
                                    "unshared_bp"])
    R["key"] = R.id.str.replace("_", "|", n=1)
    R["rec"] = R.id.str.split("_").str[1].astype(int)
    print(f"{len(R):,} SV records; hprc_type {R.hprc_type.value_counts().to_dict()}", flush=True)

    # --- masking output from both sources ------------------------------------------------------------
    H, OC, src = [], [], {}
    for ch in CHROMS:
        for tag, d in (("graffite", f"{G}/raw/{ch}"), ("graffite_steps_rerun", f"{G}/skipped/{ch}")):
            ids = [l[1:].strip() for l in open(f"{d}/indels.fa") if l[0] == ">"]
            assert os.path.exists(f"{d}/repeatmasker_dir/indels.fa.onecode.out"), f"{d} not finished"
            for k in ids:
                assert k not in src, f"{k} masked twice"
                src[k] = tag
            H.append(read_out(f"{d}/repeatmasker_dir/indels.fa.out"))
            OC.append(read_out(f"{d}/repeatmasker_dir/indels.fa.onecode.out"))
    H = pd.concat(H, ignore_index=True); OC = pd.concat(OC, ignore_index=True)
    R["annotated_by"] = R.id.map(src)
    assert R.annotated_by.notna().all(), f"{R.annotated_by.isna().sum()} records never masked"
    assert len(src) == len(R), (len(src), len(R))
    print(f"masked: {R.annotated_by.value_counts().to_dict()}; {len(H):,} RM hits, "
          f"{len(OC):,} OneCode copies", flush=True)
    print("RM classes:", H.rmclass.value_counts().to_dict(), flush=True)

    H["is_te"] = ~H.rmclass.isin(NON_TE)
    H["superfamily"] = H.rmclass
    H["len"] = H.qe - H.qs + 1
    te = H[H.is_te]
    R = R.set_index("id")
    R["span_any"] = union_bp(H).reindex(R.index).fillna(0) / R.seq_len
    R["te_bp"] = union_bp(te).reindex(R.index).fillna(0).astype(int)
    R["span_te"] = R.te_bp / R.seq_len
    R["span_simple"] = union_bp(H[~H.is_te]).reindex(R.index).fillna(0) / R.seq_len
    R["helitron_bp"] = union_bp(te[te.rmclass == "RC/Helitron"]).reindex(R.index).fillna(0).astype(int)
    R["span_graffite"] = union_bp(OC).reindex(R.index).fillna(0) / R.seq_len
    R["onecode_copies"] = OC.groupby("id").size().reindex(R.index).fillna(0).astype(int)

    # the TE family covering the most of the allele (union bp per family), and its superfamily
    fb = union_bp(te.assign(id=te.id + "\t" + te.family)).rename("fam_bp").reset_index()
    fb[["id", "family"]] = fb["index"].str.split("\t", expand=True)
    top = fb.sort_values(["id", "fam_bp", "family"], ascending=[True, False, True]) \
        .drop_duplicates("id").set_index("id")
    R["top_family"] = top.family.reindex(R.index).fillna("")
    R["top_family_frac"] = top.fam_bp.reindex(R.index).fillna(0) / R.seq_len
    sf = te.groupby(["id", "family"]).agg(superfamily=("superfamily", "first"), div=("div", "min"))
    ix = list(zip(top.index, top.family))
    R["top_superfamily"] = pd.Series(sf.superfamily.reindex(ix).to_numpy(), index=top.index) \
        .reindex(R.index).fillna("")
    R["top_div"] = pd.Series(sf["div"].reindex(ix).to_numpy(), index=top.index).reindex(R.index)
    R["n_te_families"] = te.groupby("id").family.nunique().reindex(R.index).fillna(0).astype(int)
    for c in ("span_any", "span_te", "span_simple", "span_graffite", "top_family_frac"):
        R[c] = R[c].clip(upper=1.0)

    def lab(f):
        if f >= 0.80: return "te_full"
        if f >= 0.50: return "te_most"
        if f >= 0.20: return "te_partial"
        if f > 0:     return "te_trace"
        return "non_te"
    R["cargo"] = R.span_te.map(lab)
    # what covers most of the allele: TE, simple/low-complexity/satellite, or nothing in the library
    R["dominant"] = np.select([(R.span_te >= R.span_simple) & (R.span_te >= 1 - R.span_any),
                               R.span_simple >= 1 - R.span_any], ["te", "simple"], "unmasked")

    # --- events --------------------------------------------------------------------------------------
    E = np.load(f"{OUT}/arch3_event_map.npz")
    em = pd.concat([pd.DataFrame(dict(chrom=ch, rec=E[f"{ch}_rec"], event=E[f"{ch}_bubble"]))
                    for ch in CHROMS], ignore_index=True)
    R = R.reset_index().merge(em, on=["chrom", "rec"], how="left")
    R["event_key"] = np.where(R.event.notna(), R.chrom + ":" + R.event.fillna(-1).astype(int).astype(str), "")
    print(f"records without an event: {(R.event_key == '').sum():,} "
          f"(sizes {sorted(R.loc[R.event_key == '', 'size_bp'].unique())[:5]})", flush=True)
    R = R.sort_values(["chrom", "pos", "rec"])
    R["is_rep"] = False
    first = R[R.event_key != ""].drop_duplicates("event_key").index
    R.loc[first, "is_rep"] = True

    cols = ["key", "id", "chrom", "pos", "kind", "size_bp", "seq_len", "hprc_type", "unshared_bp",
            "annotated_by", "cargo", "dominant", "span_te", "span_simple", "span_any", "span_graffite",
            "te_bp", "helitron_bp", "top_family", "top_superfamily", "top_family_frac", "top_div",
            "n_te_families", "onecode_copies", "event_key", "is_rep"]
    R = R[cols]
    R.to_csv(f"{G}/sv_repeat_class_records.csv", index=False)

    ev = R[R.event_key != ""]
    g = ev.groupby("event_key")
    rep = ev[ev.is_rep].set_index("event_key")
    EV = rep[["chrom", "pos", "kind", "key", "size_bp", "hprc_type", "cargo", "dominant", "span_te",
              "span_simple", "top_family", "top_superfamily", "top_family_frac", "top_div"]] \
        .rename(columns={"key": "rep_key"})
    EV["n_records"] = g.size()
    EV["n_records_te_full"] = g.cargo.apply(lambda s: int((s == "te_full").sum()))
    EV["frac_records_same_cargo"] = ev.assign(rc=ev.event_key.map(rep.cargo)) \
        .assign(same=lambda d: d.cargo == d.rc).groupby("event_key").same.mean()
    EV["frac_records_same_superfamily"] = ev.assign(rs=ev.event_key.map(rep.top_superfamily)) \
        .assign(same=lambda d: d.top_superfamily == d.rs).groupby("event_key").same.mean()
    EV["any_replacement"] = ev.assign(x=(ev.hprc_type == "COMPLEX") | (ev.unshared_bp > 0)) \
        .groupby("event_key").x.any()
    EV = EV.reset_index()
    EV.to_csv(f"{G}/sv_repeat_class_events.csv", index=False)

    # --- report --------------------------------------------------------------------------------------
    print(f"\nwrote {G}/sv_repeat_class_records.csv ({len(R):,} records)")
    print(f"wrote {G}/sv_repeat_class_events.csv ({len(EV):,} events)\n")
    print("RECORDS: cargo x kind x annotated_by")
    print(pd.crosstab([R.kind, R.annotated_by], R.cargo, margins=True).to_string(), "\n")
    print("EVENTS: cargo x kind")
    print(pd.crosstab(EV.kind, EV.cargo, margins=True).to_string(), "\n")
    print("EVENTS: dominant x kind")
    print(pd.crosstab(EV.kind, EV.dominant, margins=True).to_string(), "\n")
    print("EVENTS te_full: top superfamilies")
    print(EV[EV.cargo == "te_full"].groupby("kind").top_superfamily.value_counts().groupby(level=0)
          .head(8).to_string(), "\n")
    print(f"events whose records disagree on cargo class: "
          f"{(EV.frac_records_same_cargo < 1).sum():,} of {(EV.n_records > 1).sum():,} multi-record events")
    hel = R[R.cargo == "te_full"]
    print(f"te_full records whose top superfamily is RC/Helitron: {(hel.top_superfamily == 'RC/Helitron').sum():,} "
          f"of {len(hel):,}; of those, span_graffite > 0.80 (in pangenome.vcf): "
          f"{((hel.top_superfamily == 'RC/Helitron') & (hel.span_graffite > 0.8)).sum():,}")


if __name__ == "__main__":
    main()
