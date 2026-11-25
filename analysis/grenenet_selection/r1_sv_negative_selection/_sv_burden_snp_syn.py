#!/usr/bin/env python
"""Per-founder SYNONYMOUS and NONSYNONYMOUS coding-SNP burden, for the full burden model.

Same variant filter as the SNP burden in `_sv_burden_build.py`: SNP records called in >= 90% of the 231
founders, >= 2 carriers, SHARED origin (carried by >= 1 cold- and >= 1 warm-origin founder; labels from
sv_origin_partition_founders.csv). Synonymous / nonsynonymous exactly as `_sv_classic_sfs.py` (imported,
not copied): TAIR10 reference, representative (.1) gene models, sites with conflicting models dropped;
a SNP is synonymous if the codon carrying the ALT base codes the same amino acid (stop gained / lost counts
as nonsynonymous). Per founder: variants carried and variants called, so the model uses the same rate
(carried / called) as every other burden term.

Written to a SEPARATE file so `founder_burden{TAG}.csv` and the kinship matrix are untouched.
Env: kmate. SV_CALL_MIN (default 0.9). Writes results/sv_adaptive/founder_burden_snpsyn{CRTAG}.csv.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_H = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("csfs", os.path.join(_H, "_sv_classic_sfs.py"))
cs = importlib.util.module_from_spec(_s); _s.loader.exec_module(cs)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9)); MIN_CAR = 2
CRTAG = "" if CALL_MIN == 0.9 else f"_cr{int(round(CALL_MIN * 100))}"


def main():
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy(); N = len(FO)
    seq = cs.read_fasta(); deg, cod, cps, neg = cs.degeneracy(seq); del seq
    SYN = np.zeros((64, 3, 4), np.int8)
    for ci, cdn in enumerate(cs.CODE):
        for j in range(3):
            for bi, b in enumerate("ACGT"):
                SYN[ci, j, bi] = cs.CODE[cdn[:j] + b + cdn[j + 1:]] == cs.CODE[cdn]
    BIDX = {"A": 0, "C": 1, "G": 2, "T": 3}; CMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
    acc = {f"{k}_shared_CDS_all": np.zeros(N) for k in ("snpsyn", "snpnon")}
    acc.update({f"called_{k}": np.zeros(N) for k in list(acc)})
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        d = deg[C][np.minimum(pos, deg[C].size - 1)]
        cand = np.where((rl == 1) & (al == 1) & (d >= 0))[0]                   # SNPs at unambiguous coding sites
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()[:, cand]; vc = sp.load_npz(f"{base}.var_called.npz").tocsc()[:, cand]
        P = (vp > 0).tocsc(); Q = (vc > 0).tocsc(); assert P.shape[0] == N
        c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
        nc = np.asarray(P[GC].sum(0)).ravel(); nw = np.asarray(P[GW].sum(0)).ravel()
        keep = (k >= CALL_MIN * N) & (c >= MIN_CAR) & (nc > 0) & (nw > 0)
        pp = np.minimum(pos[cand], deg[C].size - 1); alt = meta["alt"][cand]
        ab = np.array([BIDX.get(CMP[a] if ng else a, -1) if a in BIDX else -1 for a, ng in zip(alt, neg[C][pp])])
        keep &= ab >= 0
        syn = np.zeros(cand.size, bool)
        syn[keep] = SYN[cod[C][pp[keep]], cps[C][pp[keep]], ab[keep]].astype(bool)
        for lab, m in (("snpsyn", keep & syn), ("snpnon", keep & ~syn)):
            s = np.where(m)[0]
            acc[f"{lab}_shared_CDS_all"] += np.asarray(P[:, s].sum(1)).ravel()
            acc[f"called_{lab}_shared_CDS_all"] += np.asarray(Q[:, s].sum(1)).ravel()
        print(f"[{C}] shared coding SNPs kept: synonymous {int((keep & syn).sum()):,}, nonsynonymous {int((keep & ~syn).sum()):,}", flush=True)
    D = pd.DataFrame(acc); D.insert(0, "founder", FO.founder.to_numpy())
    f = f"{OUT}/founder_burden_snpsyn{CRTAG}.csv"; D.to_csv(f, index=False)
    r = lambda k: D[f"{k}_shared_CDS_all"] / D[f"called_{k}_shared_CDS_all"]
    print(f"median carried / called: synonymous {D.snpsyn_shared_CDS_all.median():.0f} / {D.called_snpsyn_shared_CDS_all.median():.0f}, "
          f"nonsynonymous {D.snpnon_shared_CDS_all.median():.0f} / {D.called_snpnon_shared_CDS_all.median():.0f}")
    print(f"correlation of the two rates across founders: r = {np.corrcoef(r('snpsyn'), r('snpnon'))[0, 1]:.3f}")
    print(f"[wrote] {f}")


if __name__ == "__main__":
    main()
