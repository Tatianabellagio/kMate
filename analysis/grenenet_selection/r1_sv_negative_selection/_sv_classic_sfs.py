#!/usr/bin/env python
"""Classic site-frequency-spectrum data: SNPs at 4-fold degenerate sites (synonymous -- the putatively
neutral reference), SNPs at 0-fold degenerate sites (nonsynonymous -- the known-selected comparator),
SNPs by genomic context, short indels and SV events.

SYNONYMOUS / NONSYNONYMOUS (the default neutral reference / selected comparator): every SNP in CDS is
called from the reference codon and the codon carrying the ALT base (standard code, representative
`.1` model) -- the SnpEff-style annotation. 4-fold / 0-fold classes are kept as the stricter version.

DEGENERACY. TAIR10 reference (Ensembl Plants release, same assembly / coordinates as the panel) and
TAIR10 GFF CDS of the representative gene model (`.1`). For every codon position the number of the
three alternative bases that leave the amino acid unchanged (standard code): 3 -> 4-fold, 0 -> 0-fold.
Sites covered by CDS of two genes with conflicting degeneracy are dropped; models whose CDS length is
not a multiple of 3 are skipped.

FREQUENCIES. Panel presence/absence in 231 inbred founders (haploid sample). Per variant: carriers c
among k called founders; kept if k >= 90% of 231 and 1 <= c <= k-1 (singletons INCLUDED -- this is a
spectrum). SNPs and indels per record; SVs per Truvari event (sv_event_catalog.csv).
Missing data: every variant is projected to a common sample of M = 200 by hypergeometric sampling
(the standard down-projection), giving the expected SFS xi_j, j = 1..M-1, per class.

Env: kmate. Writes results/sv_adaptive/classic_sfs.npz:
  {class}_c, {class}_k  per-variant carriers and called;  {class}_sfs  projected unfolded SFS (M = 200)
"""
import os, sys, gzip
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.stats import hypergeom
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_H = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("sfsctx", os.path.join(_H, "_sv_sfs_context.py"))
sc = importlib.util.module_from_spec(_s); _s.loader.exec_module(sc)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
FA = f"{lib.GEA}/r1_sv_negative_selection/results/ref/TAIR10.fa.gz"
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9))   # founder call-rate floor; != 0.9 -> outputs tagged _cr{pct}
CRTAG = "" if CALL_MIN == 0.9 else f"_cr{int(round(CALL_MIN * 100))}"
SV_BP = 50
M = min(200, int(np.floor(CALL_MIN * 231 / 10) * 10))   # projection size: must not exceed the smallest k kept (0.9 -> 200, 0.7 -> 160)
CODE = dict(zip([a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
COMP = str.maketrans("ACGT", "TGCA")


def read_fasta():
    seq, name, buf = {}, None, []
    with gzip.open(FA, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if name: seq[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip().upper())
    seq[name] = "".join(buf)
    return {f"Chr{i}": seq[str(i)] for i in range(1, 6)}


def degeneracy(seq):
    g = pd.read_csv(lib.TAIR10_GENES_TE, sep="\t", header=None, comment="#",
                    names=["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"])
    g = g[(g.feat == "CDS") & g.chrom.isin(lib.CHROMS)].copy()
    g["tx"] = g.attr.str.extract(r"Parent=([^;,]+)")
    g = g[g.tx.str.endswith(".1")]
    deg = {c: np.full(len(seq[c]) + 1, -1, np.int8) for c in lib.CHROMS}
    # for synonymous / nonsynonymous calls: codon (0..63, TCAG order, coding strand), position in codon, strand
    cod = {c: np.full(len(seq[c]) + 1, -1, np.int8) for c in lib.CHROMS}
    cps = {c: np.zeros(len(seq[c]) + 1, np.int8) for c in lib.CHROMS}
    neg = {c: np.zeros(len(seq[c]) + 1, bool) for c in lib.CHROMS}
    CIDX = {k: i for i, k in enumerate(CODE)}
    nt = nskip = 0
    for (c, tx), d in g.groupby(["chrom", "tx"]):
        st = d.strand.iloc[0]; d = d.sort_values("start")
        posn = np.concatenate([np.arange(s, e + 1) for s, e in zip(d.start, d.end)])
        if st == "-":
            posn = posn[::-1]
        if posn.size % 3:
            nskip += 1; continue
        s = seq[c]; bases = "".join(s[p - 1] for p in posn)
        if st == "-":
            bases = bases.translate(COMP)
        if set(bases) - set("ACGT"):
            nskip += 1; continue
        nt += 1
        for i in range(0, posn.size, 3):
            cod_ = bases[i:i + 3]; aa = CODE[cod_]
            for j in range(3):
                same = sum(CODE[cod_[:j] + b + cod_[j + 1:]] == aa for b in "ACGT" if b != cod_[j])
                v = 4 if same == 3 else 0 if same == 0 else 2
                p = posn[i + j]
                if deg[c][p] == -1:
                    cod[c][p] = CIDX[cod_]; cps[c][p] = j; neg[c][p] = st == "-"
                deg[c][p] = v if deg[c][p] in (-1, v) else -2          # -2 = conflicting models
    print(f"degeneracy: {nt:,} gene models used, {nskip:,} skipped; "
          f"4-fold sites {sum((d == 4).sum() for d in deg.values()):,}, 0-fold {sum((d == 0).sum() for d in deg.values()):,}")
    return deg, cod, cps, neg


def project(c, k):
    """expected unfolded SFS in a sample of M, from per-variant carriers c of k called."""
    sfs = np.zeros(M + 1)
    u, cnt = np.unique(np.c_[c, k], axis=0, return_counts=True)
    j = np.arange(M + 1)
    for (ci, ki), n in zip(u, cnt):
        sfs += n * hypergeom.pmf(j, ki, ci, M)
    return sfs[1:M]                                                   # segregating in the subsample


def main():
    seq = read_fasta(); deg, cod, cps, neg = degeneracy(seq); del seq
    # SYN[codon, position, alt base (ACGT, coding strand)] = 1 synonymous, 0 nonsynonymous (incl. stop gained / lost)
    CL = list(CODE); SYN = np.zeros((64, 3, 4), np.int8)
    for ci, cdn in enumerate(CL):
        for j in range(3):
            for bi, b in enumerate("ACGT"):
                SYN[ci, j, bi] = CODE[cdn[:j] + b + cdn[j + 1:]] == CODE[cdn]
    BIDX = {"A": 0, "C": 1, "G": 2, "T": 3}; CMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
    IV = sc.layers(); N = 231
    acc = {}
    def add(name, c, k):
        a = acc.setdefault(name, [[], []]); a[0].append(c); a[1].append(k)
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        c = np.asarray((vp > 0).sum(0)).ravel(); k = np.asarray((vc > 0).sum(0)).ravel()
        ok = (k >= CALL_MIN * N) & (c >= 1) & (c <= k - 1)
        ctx = sc.context(IV, np.full(pos.size, C, dtype=object), pos, pos + np.maximum(rl - 1, 0))
        snp = (rl == 1) & (al == 1) & ok
        d = deg[C][np.minimum(pos, deg[C].size - 1)]
        add("snp_4fold", c[snp & (d == 4)], k[snp & (d == 4)]); add("snp_0fold", c[snp & (d == 0)], k[snp & (d == 0)])
        # synonymous / nonsynonymous: reference codon vs codon carrying the ALT base (unambiguous CDS sites only)
        pp = np.minimum(pos, deg[C].size - 1); cds = snp & (d >= 0)
        ii = np.where(cds)[0]; alt = meta["alt"][ii]
        ab = np.array([BIDX.get(CMP[a] if ng else a, -1) if a in BIDX else -1 for a, ng in zip(alt, neg[C][pp[ii]])])
        okb = ab >= 0; ii = ii[okb]; syn = SYN[cod[C][pp[ii]], cps[C][pp[ii]], ab[okb]].astype(bool)
        add("snp_syn", c[ii[syn]], k[ii[syn]]); add("snp_nonsyn", c[ii[~syn]], k[ii[~syn]])
        for i, x in enumerate(sc.CTX):
            add(f"snp_{x}", c[snp & (ctx == i)], k[snp & (ctx == i)])
        ind = (np.abs(al - rl) >= 1) & (np.abs(al - rl) <= SV_BP) & ~((rl == 1) & (al == 1)) & ok
        for i, x in enumerate(sc.CTX):
            add(f"indel_{x}", c[ind & (ctx == i)], k[ind & (ctx == i)])
        print(f"[{C}] SNPs {int(snp.sum()):,} (4-fold {int((snp & (d == 4)).sum()):,}, 0-fold {int((snp & (d == 0)).sum()):,}; "
              f"synonymous {int(syn.sum()):,}, nonsynonymous {int((~syn).sum()):,})", flush=True)
    E = pd.read_csv(f"{OUT}/sv_event_catalog.csv")
    E = E[(E.n_called >= CALL_MIN * N) & (E.n_carriers >= 1) & (E.n_carriers <= E.n_called - 1)]
    for kind in ("ins", "del"):
        for x in sc.CTX:
            d = E[(E.kind == kind) & (E.context == x)]; add(f"{kind}_{x}", d.n_carriers.to_numpy(), d.n_called.to_numpy())
    for cls in ("te_derived", "gene_dup", "novel", "dispersed_dup", "local_dup"):
        d = E[(E.kind == "ins") & (E.context == "CDS") & (E.seq_class == cls)]
        add(f"insCDS_{cls}", d.n_carriers.to_numpy(), d.n_called.to_numpy())
    S = {}
    for name, (cc, kk) in acc.items():
        c = np.concatenate(cc).astype(np.int64); k = np.concatenate(kk).astype(np.int64)
        S[f"{name}_c"] = c.astype(np.int16); S[f"{name}_k"] = k.astype(np.int16); S[f"{name}_sfs"] = project(c, k)
        print(f"  {name:<22} n = {c.size:>9,}", flush=True)
    np.savez_compressed(f"{OUT}/classic_sfs{CRTAG}.npz", M=M, **S)
    print(f"[wrote] {OUT}/classic_sfs{CRTAG}.npz (M = {M}, call rate >= {CALL_MIN})")


if __name__ == "__main__":
    main()
