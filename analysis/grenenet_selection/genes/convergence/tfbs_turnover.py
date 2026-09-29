#!/usr/bin/env python
"""Allele-specific TF-binding turnover at a candidate indel/SV (FIMO, PlantTFDB).

Adapted from the MOI-LAB `NARROW_MYB6_tfbs_fimo.py` (see METHODS_FUNCTIONAL_TRACKS.md).
Their version substitutes SNV alleles, so REF and ALT sequences have equal length and a
FIMO hit maps back to the genome with `win_start + start - 1`. Our candidates are indels
and SVs, where that mapping is wrong for everything downstream of the variant: after a
`d = len(REF) - len(ALT)` bp deletion the ALT sequence is shifted by `d`, so the naive
mapping pairs motif hits at DIFFERENT places and invents gains and losses.

So here the ALT coordinates are mapped back explicitly:

    p_alt <= off                    -> same reference position   (upstream of the variant)
    off < p_alt <= off + len(ALT)   -> inside the ALT allele     (no reference position;
                                       any hit overlapping it is allele-specific by
                                       construction and is reported as `in_alt_allele`)
    p_alt >  off + len(ALT)         -> p_alt + len(REF) - len(ALT)   (downstream)

A motif hit is keyed on (motif, strand, reference start, does it span the variant) and
compared between the two sequences. The end coordinate is deliberately NOT in the key: a
hit that SPANS an indel necessarily ends `len(REF)-len(ALT)` bp apart on the two alleles,
so keying on it splits one site into a LOST + a GAINED (measured: 28 + 29 spurious calls
on a 2 bp insertion at AT4G13200). Spanning hits are matched on their start and judged by
score, with the length difference reported.

    ref only -> LOST        the candidate allele destroys the site
    alt only -> GAINED      the candidate allele creates it
    both     -> WEAKENED / STRENGTHENED / unchanged, by FIMO score

Their thresholds are kept: FIMO --thresh 1e-4, keep q < 0.005, and drop low-complexity
matches (homopolymer run >= 5, or <= 2 distinct bases -- "the CAM5 lesson"). The
asymmetry, not the absolute q, is the signal: in a window this short the per-window q
inflates for both sequences equally.

Sequence sources (same coordinates as everything else in this tree):
  reference  TAIR10.chr.iupacN.fa -- the pangenome backbone, so positions match our panel
  alleles    panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz (REF/ALT as called)

Outputs -> results/tfbs/
  <sym>_hits.csv       every FIMO hit kept, both sequences, in reference coordinates
  <sym>_turnover.csv   one row per motif site that differs between the alleles
  turnover_summary.csv one row per candidate (counts + the TF families gained/lost)

Usage:
  python tfbs_turnover.py --genes FUS3,GPX6,AT4G13200 [--pad 400] [--thresh 1e-4]
env: kmate for this script; FIMO comes from the `meme` env (--fimo to override).
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile

import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402

OUT = f"{HERE}/results/tfbs"
REF_FA = ("/global/scratch/users/tbellg/pang/pang_1001gplus/"
          "20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa")
PANEL = f"{lib.PROJ}/panel/arch3"
MOTIFS = f"{lib.PROJ}/data/motifs/Ath_TF_binding_motifs.meme"
FAMILY = f"{lib.PROJ}/data/motifs/Ath_TF_binding_motifs_information.txt"
FIMO = "/global/home/users/tbellg/miniforge3/envs/meme/bin/fimo"
Q_MAX = 0.005


def low_complexity(s: str) -> bool:
    s = (s or "").upper()
    if not s:
        return True
    if re.search(r"(A{5,}|C{5,}|G{5,}|T{5,})", s):
        return True
    return len(set(s)) <= 2


def repeat_tract(s: str, min_units: int = 5) -> bool:
    """Is the match mostly a di/tri-nucleotide repeat tract?

    The homopolymer rule above does not catch (CT)n / (GA)n microsatellites: they have 3-4
    distinct bases once flanking bases are included, so they pass. Measured at AT4G13200,
    where a 2 bp insertion sits inside a CT tract and 53 of 57 spanning hits are BBR-BPC
    (GA-repeat binders) -- motif turnover there is repeat-length variation, not a
    transcription-factor story. Flagged rather than dropped: repeat-tract binding can be
    real (BPC proteins do bind GA repeats), but it cannot be told from noise by FIMO alone.
    """
    s = (s or "").upper()
    for k in (2, 3):
        for i in range(k):
            unit = s[i:i + k]
            if len(unit) < k:
                continue
            run = 0
            for j in range(i, len(s) - k + 1, k):
                run = run + 1 if s[j:j + k] == unit else 0
                if run >= min_units:
                    return True
    return False


def candidates(genes: list[str]) -> pd.DataFrame:
    """Resolve symbols/AGIs against the current review tables -> one row per candidate."""
    R = f"{HERE}/results"
    frames = []
    T = pd.read_csv(f"{R}/screen_top_loci_ownaxis.csv")
    frames.append(T.assign(set="GEA_ownaxis")[
        ["target_gene", "symbol", "chrom", "pos", "ref_len", "alt_len", "ftier", "set"]])
    G = pd.read_csv(f"{R}/screen_gwas_rescreen_v2.csv")
    frames.append(G.assign(set="GWAS")[
        ["target_gene", "symbol", "chrom", "pos", "ref_len", "alt_len", "ftier", "set"]])
    S = pd.read_csv(f"{R}/screen_3criteria.csv")
    V = pd.read_csv(f"{R}/screen_visual_review.csv")
    r1 = V.merge(S[["store_row", "ftier"]], on="store_row", how="left", suffixes=("", "_s"))
    frames.append(r1.assign(set="round1", symbol=r1.sym)[
        ["target_gene", "symbol", "chrom", "pos", "ref_len", "alt_len", "ftier", "set"]])
    C = pd.concat(frames, ignore_index=True).drop_duplicates(
        ["target_gene", "chrom", "pos", "ref_len", "alt_len"])
    want = {g.strip().upper() for g in genes}
    sel = C[C.symbol.astype(str).str.upper().isin(want) | C.target_gene.str.upper().isin(want)]
    missing = want - set(sel.symbol.astype(str).str.upper()) - set(sel.target_gene.str.upper())
    if missing:
        print(f"  not found in the candidate tables: {sorted(missing)}")
    return sel


def fetch_allele(chrom: str, pos: int, ref_len: int, alt_len: int):
    """REF/ALT sequences of the panel record with these lengths at this position."""
    vcf = f"{PANEL}/{chrom.lower()}/merged_231_{chrom.lower()}_final.vcf.gz"
    with pysam.VariantFile(vcf) as vf:
        for rec in vf.fetch(chrom, pos - 1, pos):
            if rec.pos != pos:
                continue
            for alt in rec.alts or ():
                if len(rec.ref) == ref_len and len(alt) == alt_len:
                    return rec.ref.upper(), alt.upper()
    return None, None


def fimo_hits(fa: str, thresh: float, fimo_bin: str) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as td:
        cmd = [fimo_bin, "--oc", td, "--thresh", str(thresh), "--verbosity", "1",
               MOTIFS, fa]
        p = subprocess.run(cmd, capture_output=True, text=True)
        tsv = f"{td}/fimo.tsv"
        if p.returncode != 0 or not os.path.exists(tsv):
            raise SystemExit(f"FIMO failed: {(p.stderr or '')[-400:]}")
        d = pd.read_csv(tsv, sep="\t", comment="#")
    d = d.dropna(subset=["motif_id", "sequence_name"])
    return d.rename(columns={"p-value": "pval", "q-value": "qval"})


def map_to_ref(p: int, off: int, len_ref: int, len_alt: int, seq: str):
    """ALT-sequence coordinate (1-based) -> reference coordinate.

    `off` is the 0-based index of the variant start, so the variant occupies 1-based
    positions off+1 .. off+len_alt on the ALT sequence. Returns None inside that stretch:
    those bases have no reference position.
    """
    if seq == "ref":
        return p
    if p <= off:
        return p
    if p <= off + len_alt:
        return None
    return p + len_ref - len_alt


def run_one(r, pad: int, thresh: float, fimo_bin: str, fa_ref) -> dict | None:
    chrom, pos = r.chrom, int(r.pos)
    ref_a, alt_a = fetch_allele(chrom, pos, int(r.ref_len), int(r.alt_len))
    if ref_a is None:
        print(f"  {r.symbol}: no panel record with these allele lengths -- skipped")
        return None
    lo, hi = max(1, pos - pad), pos + pad + len(ref_a)
    window = fa_ref.fetch(chrom, lo - 1, hi).upper()
    off = pos - lo                      # 0-based index of the variant start in `window`
    if window[off:off + len(ref_a)] != ref_a:
        print(f"  {r.symbol}: REF does not match the genome at {chrom}:{pos} "
              f"({window[off:off+len(ref_a)][:20]} vs {ref_a[:20]}) -- skipped")
        return None
    alt_window = window[:off] + alt_a + window[off + len(ref_a):]

    sym = str(r.symbol).replace("/", "_").replace(" ", "_")
    os.makedirs(OUT, exist_ok=True)
    fa = f"{OUT}/{sym}_locus.fa"
    with open(fa, "w") as fh:
        fh.write(f">ref\n{window}\n>alt\n{alt_window}\n")
    d = fimo_hits(fa, thresh, fimo_bin)
    os.remove(fa)

    fam = {}
    if os.path.exists(FAMILY):
        f = pd.read_csv(FAMILY, sep="\t")
        gid = "Gene_id" if "Gene_id" in f.columns else f.columns[0]
        fam = dict(zip(f[gid].astype(str), f["Family"].astype(str)))

    rows = []
    for h in d.itertuples():
        keep_q = float(h.qval) < Q_MAX
        matched = getattr(h, "matched_sequence", "")
        lowc = low_complexity(matched)
        rept = repeat_tract(matched)
        n_var = len(ref_a) if h.sequence_name == "ref" else len(alt_a)
        spans = int(h.start) <= off + n_var and int(h.stop) > off
        s_ref = map_to_ref(int(h.start), off, len(ref_a), len(alt_a), h.sequence_name)
        e_ref = map_to_ref(int(h.stop), off, len(ref_a), len(alt_a), h.sequence_name)
        # a hit lying wholly inside the ALT allele has no reference position at all
        in_allele = (h.sequence_name == "alt") and s_ref is None and e_ref is None
        rows.append(dict(
            symbol=r.symbol, gene=r.target_gene, seq=h.sequence_name, motif=h.motif_id,
            family=fam.get(str(h.motif_id), "?"), strand=h.strand,
            gstart=lo + (s_ref - 1) if s_ref else pos,
            gend=lo + (e_ref - 1) if e_ref else pos + len(ref_a) - 1,
            span_variant=bool(spans), in_alt_allele=in_allele, length=int(h.stop) - int(h.start) + 1,
            score=float(h.score), pval=float(h.pval),
            qval=float(h.qval), matched=getattr(h, "matched_sequence", ""),
            low_complexity=lowc, repeat_tract=rept,
            kept=bool(keep_q and not lowc)))
    H = pd.DataFrame(rows)
    H.sort_values(["gstart", "pval"]).to_csv(f"{OUT}/{sym}_hits.csv", index=False)

    K = H[H.kept].copy()
    # end coordinate excluded on purpose (see the module docstring)
    K["key"] = list(zip(K.motif, K.strand, K.gstart, K.span_variant, K.in_alt_allele))
    ref_m = {k: v for k, v in zip(K[K.seq == "ref"].key, K[K.seq == "ref"].itertuples())}
    alt_m = {k: v for k, v in zip(K[K.seq == "alt"].key, K[K.seq == "alt"].itertuples())}
    turn = []
    for k in sorted(set(ref_m) | set(alt_m), key=lambda x: (x[2], x[0])):
        rr, aa = ref_m.get(k), alt_m.get(k)
        if rr and not aa:
            eff = "LOST"
        elif aa and not rr:
            eff = "GAINED"
        elif aa.score < rr.score:
            eff = "WEAKENED"
        elif aa.score > rr.score:
            eff = "STRENGTHENED"
        else:
            eff = "unchanged"
        x = rr or aa
        turn.append(dict(symbol=r.symbol, gene=r.target_gene, motif=k[0], family=x.family,
                         strand=k[1], gstart=k[2], span_variant=k[3], in_alt_allele=k[4],
                         ref_len_bp=rr.length if rr else None,
                         alt_len_bp=aa.length if aa else None,
                         ref_score=round(rr.score, 1) if rr else None,
                         alt_score=round(aa.score, 1) if aa else None, effect=eff,
                         best_q=round(min(v.qval for v in (rr, aa) if v), 4)))
    T = pd.DataFrame(turn)
    if len(T):
        T[T.effect != "unchanged"].to_csv(f"{OUT}/{sym}_turnover.csv", index=False)

    n = {e: int((T.effect == e).sum()) if len(T) else 0
         for e in ("LOST", "GAINED", "WEAKENED", "STRENGTHENED")}
    n_rep = int(K.repeat_tract.sum())
    fams = lambda e: ",".join(sorted({str(f) for f in T[T.effect == e].family})) if len(T) else ""
    kind = "deletion" if len(ref_a) > len(alt_a) else "insertion"
    print(f"  {str(r.symbol):12s} {chrom}:{pos:>9,} {abs(len(alt_a)-len(ref_a)):>5} bp {kind}"
          f"  kept {int(K.seq.eq('ref').sum()):>3} ref / {int(K.seq.eq('alt').sum()):>3} alt hits"
          f"  LOST {n['LOST']:>2}  GAINED {n['GAINED']:>2}  "
          f"WEAK {n['WEAKENED']:>2}  STRONG {n['STRENGTHENED']:>2}"
          + (f"   [{n_rep}/{len(K)} hits are repeat-tract matches]" if n_rep else ""))
    return dict(symbol=r.symbol, gene=r.target_gene, set=r.set, ftier=r.ftier,
                chrom=chrom, pos=pos, size=abs(len(alt_a) - len(ref_a)), kind=kind,
                window=f"{lo}-{hi}", n_ref_hits=int(K.seq.eq("ref").sum()),
                n_alt_hits=int(K.seq.eq("alt").sum()), n_repeat_tract_hits=n_rep,
                **{f"n_{k.lower()}": v for k, v in n.items()},
                families_lost=fams("LOST"), families_gained=fams("GAINED"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="comma-separated symbols or AGIs")
    ap.add_argument("--pad", type=int, default=400, help="bp of reference either side (400)")
    ap.add_argument("--thresh", type=float, default=1e-4, help="FIMO p threshold (1e-4)")
    ap.add_argument("--fimo", default=FIMO)
    a = ap.parse_args()

    for f in (REF_FA, MOTIFS):
        if not os.path.exists(f):
            sys.exit(f"missing input: {f}")
    if not os.path.exists(a.fimo):
        sys.exit(f"fimo not found at {a.fimo} (mamba create -n meme -c bioconda meme)")

    C = candidates(a.genes.split(","))
    if not len(C):
        sys.exit("no candidates resolved")
    os.makedirs(OUT, exist_ok=True)
    fa_ref = pysam.FastaFile(REF_FA)
    print(f"TFBS turnover, FIMO p<{a.thresh}, q<{Q_MAX}, +-{a.pad} bp:")
    rows = [x for x in (run_one(r, a.pad, a.thresh, a.fimo, fa_ref) for r in C.itertuples())
            if x]
    if rows:
        S = pd.DataFrame(rows)
        out = f"{OUT}/turnover_summary.csv"
        if os.path.exists(out):
            old = pd.read_csv(out)
            S = pd.concat([old[~old.gene.isin(S.gene)], S], ignore_index=True)
        S.sort_values(["set", "symbol"]).to_csv(out, index=False)
        print(f"\nwrote {out} and per-gene hits/turnover CSVs in {OUT}")


if __name__ == "__main__":
    main()
