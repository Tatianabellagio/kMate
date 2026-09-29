#!/usr/bin/env python
"""Map every SV record of the arch3 panel to an SV EVENT: alleles of one pangenome bubble that are the
same insertion (or deletion), merged.

WHY. The allele (panel record) is too strict: arch3 splits one insertion into one record per
assembly path, so the same event appears as several alleles with small sequence differences
(73% of passing insertion alleles share their bubble with another passing allele; 70% have one
within 10% of their length). The bubble is too loose: it can hold genuinely different insertions
at the same site. The event sits between the two.

MODE=truvari (default) = `truvari collapse` DEFAULTS (v5; truvari/collapse.py argparse):
  within a bubble, same type; a pair matches if
    size similarity  min(len)/max(len)                       >= 0.95  (--pctsize 0.95)
    seq similarity   (lenA + lenB - editdist) / (lenA + lenB) >= 0.95  (--pctseq 0.95; edlib, on the
                     inserted / deleted sequence, B rotated by the position offset = unroll)
    |posA - posB| <= 500 bp                                            (--refdist 500)
  greedy, no chaining (--chain off), representative = --keep first (earliest position; ties -> panel
  order): each allele joins the first kept representative it matches, else becomes a new one.
  NOTE: 0.7/0.7 are the `truvari bench` defaults (matching calls to a truth set), NOT the merge
  defaults; that more permissive version is kept as arch3_event_map_bench07.npz (PCTSEQ=PCTSIZE=0.7,
  KEEP=common) for comparison only.
MODE=strict (sensitivity), RULE (within a bubble, same type only, single linkage):
  insertions : lengths within 10% (min/max >= 0.9) AND inserted-sequence 15-mer Jaccard >= J_MIN
               (default 0.5, ~ <= 2-3% divergence). Long sequences are compared on a 1/4 k-mer
               subsample (hash mod 4 == 0) -- an unbiased Jaccard estimate.
  deletions  : lengths within 10% AND reciprocal overlap of the deleted reference spans >= 90%.
Records keep their var_pa column index; output has the same layout as arch3_bubble_map.npz
({C}_rec, {C}_bubble = EVENT code), so downstream bubble code paths run on events unchanged.
Also reports how many events a J_MIN of 0.3 / 0.8 would give (sensitivity).

Env: kmate. Writes results/sv_adaptive/arch3_event_map.npz (collapse defaults), _bench07 (0.7/0.7), _strict.
"""
import os, sys, gzip
import numpy as np
import scipy.sparse as sp
import edlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PANEL = "/global/scratch/users/tbellg/kmate/panel/arch3"
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP, K, LEN_TOL, RO_MIN = 50, 15, 0.9, 0.9
J_MIN = float(os.environ.get("EVENT_J", 0.5))
J_SENS = (0.3, 0.8)
MODE = os.environ.get("EVENT_MODE", "truvari")
PCTSIZE = float(os.environ.get("PCTSIZE", 0.95)); PCTSEQ = float(os.environ.get("PCTSEQ", 0.95)); REFDIST = 500
KEEP = os.environ.get("KEEP", "first")             # first (truvari default) | common


def seqsim(a, b, p=0):
    """Truvari seqsim with unroll: (lenA+lenB-ED)/(lenA+lenB); ED bounded for speed."""
    if p and len(b):
        f = p % len(b)
        if f:
            b = b[-f:] + b[:-f]
    tot = len(a) + len(b)
    kmax = int((1 - PCTSEQ) * tot)
    ed = edlib.align(a.upper(), b.upper(), task="distance", k=kmax)["editDistance"]
    return 0.0 if ed < 0 else (tot - ed) / tot


def truvari_events(g, R, cols, ncar):
    """greedy collapse of the alleles g (indices into cols) of one bubble & type -> list of labels."""
    seq = {x: R[cols[x]][3] for x in g}
    pos = {x: R[cols[x]][0] for x in g}
    order = (sorted(g, key=lambda x: (-ncar[x], x)) if KEEP == "common"
             else sorted(g, key=lambda x: (pos[x], x)))
    reps, lab = [], {}
    for x in order:
        hit = None
        for r in reps:
            la, lb = len(seq[x]), len(seq[r])
            if min(la, lb) / max(la, lb, 1) < PCTSIZE or abs(pos[x] - pos[r]) > REFDIST:
                continue
            if seqsim(seq[r], seq[x], pos[x] - pos[r]) >= PCTSEQ:
                hit = r; break
        if hit is None:
            reps.append(x); lab[x] = x
        else:
            lab[x] = lab[hit]
    return [lab[x] for x in g]


def kmers(s):
    if len(s) < K:
        return {s}
    ks = {s[i:i + K] for i in range(len(s) - K + 1)}
    if len(s) > 1000:
        ks = {x for x in ks if hash(x) % 4 == 0}
    return ks


def jaccard(a, b):
    return len(a & b) / max(1, len(a | b))


class UF:
    def __init__(self, n): self.p = list(range(n))
    def find(self, i):
        while self.p[i] != i:
            self.p[i] = self.p[self.p[i]]; i = self.p[i]
        return i
    def union(self, i, j): self.p[self.find(i)] = self.find(j)


def sv_records(ch, want):
    """Stream the merged panel VCF; for SV columns in `want` return pos, ref, inserted seq."""
    out = {}
    i = 0
    with gzip.open(f"{PANEL}/{ch.lower()}/merged_231_{ch.lower()}_final.vcf.gz", "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            if i in want:
                f = line.split("\t", 6)
                ref, alt = f[3], f[4]
                if len(alt) > len(ref):
                    s = alt[len(ref):] if alt.startswith(ref) else alt        # inserted sequence
                else:
                    s = ref[len(alt):] if ref.startswith(alt) else ref        # deleted sequence
                out[i] = (int(f[1]), len(ref), len(alt), s)
            i += 1
    return out


def main():
    os.environ.setdefault("PYTHONHASHSEED", "0")
    BM = np.load(f"{OUT}/arch3_bubble_map.npz")
    acc = {}
    for ch in lib.CHROMS:
        cols, bub = BM[f"{ch}_rec"], BM[f"{ch}_bubble"]
        R = sv_records(ch, set(cols.tolist()))
        vp = sp.load_npz(f"{lib.PROJ}/panel/arch3/{ch.lower()}/var_pa_231_arch3_{ch.lower()}.var_pa.npz").tocsc()
        ncar = np.asarray((vp[:, cols] > 0).sum(0)).ravel()
        ev = np.full(cols.size, -1, np.int64); nxt = 0
        sens = {j: 0 for j in J_SENS}
        order = np.argsort(bub, kind="stable")
        bounds = np.flatnonzero(np.diff(bub[order])) + 1
        for grp in np.split(order, bounds):
            for is_ins in (True, False):
                g = [x for x in grp if (R[cols[x]][2] - R[cols[x]][1] > 0) == is_ins]
                if not g:
                    continue
                n = len(g)
                if MODE == "truvari":
                    labs = truvari_events(g, R, cols, ncar)
                    rc = {r: nxt + k for k, r in enumerate(sorted(set(labs)))}
                    for a in range(n):
                        ev[g[a]] = rc[labs[a]]
                    nxt += len(rc)
                    continue
                L = [abs(R[cols[x]][2] - R[cols[x]][1]) for x in g]
                if is_ins:
                    KS = [kmers(R[cols[x]][3]) for x in g]
                else:
                    iv = [(R[cols[x]][0], R[cols[x]][0] + R[cols[x]][1] - 1) for x in g]
                pairs = []
                for a in range(n):
                    for b in range(a + 1, n):
                        if min(L[a], L[b]) / max(L[a], L[b]) < LEN_TOL:
                            continue
                        if is_ins:
                            s = jaccard(KS[a], KS[b])
                        else:
                            ov = min(iv[a][1], iv[b][1]) - max(iv[a][0], iv[b][0]) + 1
                            s = 1.0 if ov / max(iv[a][1] - iv[a][0] + 1, iv[b][1] - iv[b][0] + 1) >= RO_MIN else 0.0
                        pairs.append((a, b, s))
                for thr in (J_MIN,) + J_SENS:
                    uf = UF(n)
                    for a, b, s in pairs:
                        if s >= thr:
                            uf.union(a, b)
                    roots = {uf.find(a) for a in range(n)}
                    if thr == J_MIN:
                        rc = {r: nxt + k for k, r in enumerate(sorted(roots))}
                        for a in range(n):
                            ev[g[a]] = rc[uf.find(a)]
                        nxt += len(roots)
                    else:
                        sens[thr] += len(roots)
        assert (ev >= 0).all()
        acc[f"{ch}_rec"] = cols; acc[f"{ch}_bubble"] = ev
        isin = np.array([R[c][2] > R[c][1] for c in cols])
        evi = np.unique(ev[isin]).size; evd = np.unique(ev[~isin]).size
        print(f"[{ch}] {cols.size:,} SV records ({isin.sum():,} ins) -> {len(np.unique(bub)):,} bubbles -> "
              f"{len(np.unique(ev)):,} events ({MODE}; ins {evi:,}, del {evd:,})"
              + (f"; J>=0.3: {sens[0.3]:,}, J>=0.8: {sens[0.8]:,}" if MODE == "strict" else ""), flush=True)
    suf = ("_strict" if MODE == "strict" else
           "" if (PCTSEQ, PCTSIZE, KEEP) == (0.95, 0.95, "first") else
           "_bench07" if (PCTSEQ, PCTSIZE) == (0.7, 0.7) else f"_seq{PCTSEQ}_size{PCTSIZE}_{KEEP}")
    fn = f"{OUT}/arch3_event_map{suf}.npz"
    np.savez_compressed(fn, chroms=np.array(lib.CHROMS, dtype=str), mode=MODE, j_min=J_MIN, **acc)
    print(f"[wrote] {fn}")


if __name__ == "__main__":
    main()
