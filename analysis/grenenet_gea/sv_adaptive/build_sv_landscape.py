#!/usr/bin/env python3
"""Per-haploblock structural-variant landscape on the clq0.9 (r2>=0.9) blocks.

For each clq0.9 haploblock, count the SNP / small-indel / SV records it contains,
on the SAME founder variant set the blocks were built from (founder MAF>=0.05,
called-frac>=0.9). SV = |alt_len - ref_len| > 50 bp; small-indel = other non-SNP.

Emits one row per block with class counts, length, SV density (per kb) and SV
fraction, plus block midpoint + distance to the (peri)centromere — the covariates
the enrichment test will need to control for (block size, genomic context).

Run in the `kmate` env (numpy 2 / scipy). Output CSV -> results/grenenet_gea/sv_adaptive/.
"""
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp

PANEL = Path("panel/arch3")
BLKDIR = Path("results/grenenet_gea/blocks_mcf90")
OUT = Path("results/grenenet_gea/sv_adaptive"); OUT.mkdir(parents=True, exist_ok=True)
CHROMS = [f"Chr{i}" for i in range(1, 6)]
# COMMON-SV founder floor: MAF>=0.05 (MAC>=12 among the 231 founders). Decision
# (2026-07-01, reverting the brief MAC>=2 no-singleton trial): rely on COMMON SVs only,
# because low-frequency SV calls are more likely spurious (kMate/panel calling error near
# the detection threshold) -- so the conservative, reliable set is MAF>=0.05, even though
# it keeps only ~3.5% of SV records. This matches the MAF floor the LD blocks themselves
# were built at (blocks_mcf90). For ONE consistent SV definition across instruments, the
# per-site temporal (site_sv_enrichment.py) should use the same floor (--min-mac 12).
MAC_MIN, CALLED_MIN, SV_BP = 12, 0.9, 50
# approx TAIR10 centromere midpoints (bp) for a genomic-context covariate
CEN = {"Chr1": 15_086_000, "Chr2": 3_607_000, "Chr3": 13_799_000,
       "Chr4": 3_956_000, "Chr5": 11_725_000}


def per_chrom(ch):
    meta = np.load(PANEL / ch.lower() / f"var_pa_231_arch3_{ch.lower()}.meta.npz", allow_pickle=True)
    pos = meta["pos"].astype(np.int64)
    dl = (meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
    is_snp = (meta["ref_len"] == 1) & (meta["alt_len"] == 1)
    vp = sp.load_npz(PANEL / ch.lower() / f"var_pa_231_arch3_{ch.lower()}.var_pa.npz")
    vc = sp.load_npz(PANEL / ch.lower() / f"var_pa_231_arch3_{ch.lower()}.var_called.npz")
    nF = vp.shape[0]
    n_alt = np.asarray(vp.sum(0)).ravel().astype(np.float64)
    n_cal = np.asarray(vc.sum(0)).ravel().astype(np.float64)
    # MAC on the full founder panel (identical to lib.founder_panel_keep): no singletons
    # on either allele + call-rate floor. NOT carrier-freq among-called, so it matches the
    # per-site temporal floor byte-for-byte.
    keep = (n_alt >= MAC_MIN) & (n_alt <= nF - MAC_MIN) & (n_cal / nF >= CALLED_MIN)
    # class of each KEPT variant
    snp = is_snp[keep]
    sv = np.abs(dl[keep]) > SV_BP
    indel = (~snp) & (~sv)
    kpos = pos[keep]
    cls = np.where(snp, 0, np.where(sv, 2, 1))  # 0=snp 1=indel 2=sv
    print(f"{ch}: {keep.sum():,}/{len(pos):,} variants kept "
          f"(snp {int(snp.sum()):,} / indel {int(indel.sum()):,} / sv {int(sv.sum()):,})", flush=True)

    blk = pd.read_csv(BLKDIR / f"{ch.lower()}_clq0.9_blocks_clq0.9.tsv", sep="\t")
    starts = blk.start_pos.to_numpy(); ends = blk.end_pos.to_numpy()
    order = np.argsort(starts); starts, ends = starts[order], ends[order]
    blk = blk.iloc[order].reset_index(drop=True)
    # assign each kept variant to a block (blocks non-overlapping, sorted)
    j = np.searchsorted(starts, kpos, side="right") - 1
    ok = (j >= 0) & (kpos <= ends[np.clip(j, 0, len(ends) - 1)])
    j = j[ok]; cj = cls[ok]
    nb = len(blk)
    cnt = np.zeros((nb, 3), np.int64)
    np.add.at(cnt, (j, cj), 1)
    out = blk.copy()
    out["chrom"] = ch
    out["n_snp"], out["n_indel"], out["n_sv"] = cnt[:, 0], cnt[:, 1], cnt[:, 2]
    out["n_kept"] = cnt.sum(1)
    out["len_bp"] = (out.end_pos - out.start_pos + 1).clip(lower=1)
    out["mid"] = (out.start_pos + out.end_pos) / 2
    out["dist_cen"] = (out["mid"] - CEN[ch]).abs()
    out["sv_density_kb"] = out.n_sv / (out.len_bp / 1000.0)
    out["sv_frac"] = np.divide(out.n_sv, out.n_kept, out=np.zeros(nb), where=out.n_kept.to_numpy() > 0)
    out["has_sv"] = (out.n_sv > 0).astype(int)
    return out


parts = [per_chrom(ch) for ch in CHROMS]
df = pd.concat(parts, ignore_index=True)
df["block_id"] = df.chrom + ":" + df.start_pos.astype(str) + "-" + df.end_pos.astype(str)
cols = ["block_id", "chrom", "start_pos", "end_pos", "len_bp", "mid", "dist_cen",
        "n_variants", "n_kept", "n_snp", "n_indel", "n_sv", "sv_density_kb", "sv_frac", "has_sv"]
df = df[cols]
df.to_csv(OUT / "sv_landscape_clq0.9.csv", index=False)

print(f"\n== SV landscape over {len(df):,} clq0.9 blocks ==")
print(f"blocks with >=1 SV: {int(df.has_sv.sum()):,} ({100*df.has_sv.mean():.1f}%)")
print(f"total: snp {df.n_snp.sum():,} | indel {df.n_indel.sum():,} | sv {df.n_sv.sum():,} "
      f"(sv = {100*df.n_sv.sum()/df.n_kept.sum():.2f}% of kept records)")
print(f"per-block n_kept: median {df.n_kept.median():.0f}  (min {df.n_kept.min()}, max {df.n_kept.max()})")
print(f"per-block n_sv:   median {df.n_sv.median():.0f}  mean {df.n_sv.mean():.2f}  max {df.n_sv.max()}")
print(f"sv_frac (SV-containing blocks): median {df.loc[df.has_sv==1,'sv_frac'].median():.3f}")
print("\nby chrom (blocks | %with SV | total SV):")
for ch, g in df.groupby("chrom"):
    print(f"  {ch}: {len(g):>6,} | {100*g.has_sv.mean():5.1f}% | {int(g.n_sv.sum()):>6,}")
print(f"\n-> {OUT/'sv_landscape_clq0.9.csv'}")
