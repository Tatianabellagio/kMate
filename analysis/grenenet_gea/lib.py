"""Shared helpers for the GrENE-Net SV-GEA analysis (kMate outputs).

The goal (mirroring the phase-1 SNP paper, but for SVs): find structural
variants whose frequency rises in hot sites and falls in cold ones — climate-
dependent directional selection — and ask whether SVs add adaptive signal SNPs
miss. See README.md.

Data:
  - kMate per-sample AF TSVs: results/grenenet_kmate_arch3/<MLFH...>.tsv
      cols: chrom pos ref_len alt_len alt_freq info n_called se
  - founding p0 (gen 0): the SEEDMIX kMate outputs (mean over 8 reps).
  - sample metadata (site/plot/generation/coverage): Table_S5.
"""
from __future__ import annotations
import os, glob
import numpy as np
import pandas as pd

PROJ = "/global/scratch/users/tbellg/kmate"
OUT = f"{PROJ}/results/grenenet_kmate_arch3"
SEEDMIX = f"{PROJ}/results/seedmix_kmate_arch3"
GEA = f"{PROJ}/results/grenenet_gea"
T5 = ("/global/scratch/users/tbellg/pang/grenenet_reads/"
      "Table_S5_sample_collection_sequencing_library.csv")

# Pilot sites (extend as the cohort grows). site -> (label, role)
SITE_CLIMATE = {4: "hot", 54: "cold"}   # 4=Cadiz/Madrid region Spain, 54=Cologne DE


def list_samples(base: str = OUT) -> list[str]:
    """Genome-wide per-sample TSVs present in `base` (excludes per-chrom files)."""
    return sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{base}/*.tsv")
                  if "_Chr" not in os.path.basename(p))


def sample_map(samples: list[str] | None = None) -> pd.DataFrame:
    """Map sample_id -> site, plot, date, year, generation, coverage, climate.

    Indexed by sample_id. `coverage` is Table_S5 weighted_mean_coverage.
    """
    t5 = pd.read_csv(T5).set_index("sampleid")
    if samples is None:
        samples = list_samples()
    keep = [s for s in samples if s in t5.index]
    # `coverage` = per-SAMPLE; `weighted_mean_coverage` = per-PLOT aggregate
    # (constant across a plot's generations) — keep both, default to per-sample.
    m = t5.loc[keep, ["site", "plot", "date", "year", "generation",
                      "coverage", "weighted_mean_coverage"]].copy()
    m = m.rename(columns={"weighted_mean_coverage": "wmean_cov_plot"})
    m["climate"] = m["site"].map(SITE_CLIMATE)
    return m


def load_af(sample: str, base: str = OUT) -> pd.DataFrame:
    """Load one sample's per-record AF TSV."""
    return pd.read_csv(f"{base}/{sample}.tsv", sep="\t")


def rec_key(df: pd.DataFrame) -> pd.Series:
    """Stable per-record key (chrom:pos:ref_len:alt_len) for joining samples.

    kMate TSVs carry allele *lengths*, not sequences, so this keys on
    (locus, ref_len, alt_len) — unique per panel record within a chrom/pos.
    """
    return (df.chrom.astype(str) + ":" + df.pos.astype(str) + ":"
            + df.ref_len.astype(str) + ":" + df.alt_len.astype(str))


SV_MIN_BP = 50   # SV threshold: length change must exceed this many bp


def is_snp(df: pd.DataFrame) -> pd.Series:
    return (df.ref_len == 1) & (df.alt_len == 1)


def is_indel(df: pd.DataFrame) -> pd.Series:
    """Any non-SNP record (ref or alt not a single base) — includes small indels."""
    return (df.ref_len != 1) | (df.alt_len != 1)


def is_sv(df: pd.DataFrame, min_bp: int = SV_MIN_BP) -> pd.Series:
    """SV = large indel: length change |alt_len - ref_len| > min_bp (default 50 bp)."""
    return (df.alt_len - df.ref_len).abs() > min_bp


def sv_size(df: pd.DataFrame) -> pd.Series:
    """SV length proxy = |alt_len - ref_len| (0 for SNPs/MNPs of equal length)."""
    return (df.alt_len - df.ref_len).abs()


def build_p0(sv_only: bool = False, cache: bool = True) -> pd.Series:
    """Founding (gen-0) per-record AF = mean alt_freq over the 8 SEEDMIX reps.

    Returns a Series indexed by rec_key. Cached to results/grenenet_gea/.
    """
    tag = "sv" if sv_only else "all"
    cache_path = f"{GEA}/p0_seedmix_{tag}.pkl"   # pickle: no pyarrow dependency
    if cache and os.path.exists(cache_path):
        return pd.read_pickle(cache_path)
    fs = [f for f in glob.glob(f"{SEEDMIX}/SEEDMIX_S*.tsv")
          if "_Chr" not in os.path.basename(f)]
    if not fs:
        raise FileNotFoundError(f"no SEEDMIX TSVs in {SEEDMIX}")
    tot = None
    n = np.zeros(0)
    for f in fs:
        d = pd.read_csv(f, sep="\t")
        if sv_only:
            d = d[is_sv(d)]
        s = pd.Series(d.alt_freq.to_numpy(), index=rec_key(d))
        if tot is None:
            tot = s.fillna(0.0)
            n = s.notna().astype(float)
        else:
            tot = tot.add(s.fillna(0.0), fill_value=0.0)
            n = n.add(s.notna().astype(float), fill_value=0.0)
    p0 = (tot / n.replace(0, np.nan)).rename("p0")
    if cache:
        os.makedirs(GEA, exist_ok=True)
        p0.to_pickle(cache_path)
    return p0


def build_qc_table(samples: list[str] | None = None, base: str = OUT,
                   out_csv: str | None = None) -> pd.DataFrame:
    """One-time per-sample QC summary (reads each TSV's lean columns once).

    Per sample: record/SV/SNP counts, finite-AF fraction, AF means (all + SV),
    called-mask stats, + joined coverage/generation/climate + eff_n_founders.
    Cache to results/grenenet_gea/pilot_qc.csv so notebooks just load it.
    """
    if samples is None:
        samples = list_samples(base)
    m = sample_map(samples)
    rows = []
    for s in samples:
        d = pd.read_csv(f"{base}/{s}.tsv", sep="\t",
                        usecols=["ref_len", "alt_len", "alt_freq", "n_called"])
        sv = (d.ref_len != 1) | (d.alt_len != 1)
        fin = d.alt_freq.notna()
        enf = eff_n_founders(s, base)
        rows.append(dict(
            sample=s,
            n_records=len(d), n_sv=int(sv.sum()), n_snp=int((~sv).sum()),
            finite_af_frac=float(fin.mean()),
            af_mean=float(d.alt_freq.mean()),
            sv_af_mean=float(d.loc[sv, "alt_freq"].mean()),
            median_n_called=float(d.n_called.median()),
            frac_ncalled_lt50=float((d.n_called < 50).mean()),
            eff_n_founders_mean=float(np.nanmean(list(enf.values()))) if enf else np.nan,
        ))
    qc = pd.DataFrame(rows).set_index("sample")
    qc = qc.join(m[["site", "plot", "generation", "coverage", "climate"]])
    if out_csv is None:
        out_csv = f"{GEA}/pilot_qc.csv"
    os.makedirs(GEA, exist_ok=True)
    qc.to_csv(out_csv)
    return qc


def _af_array(path: str) -> np.ndarray:
    """Just the alt_freq column as a float array (records are panel-ordered)."""
    return pd.read_csv(path, sep="\t", usecols=["alt_freq"]).alt_freq.to_numpy(dtype=float)


def build_group_means(samples: list[str] | None = None, base: str = OUT,
                      cache: bool = True) -> pd.DataFrame:
    """Per-record mean alt_freq for each (climate, generation) group + founding p0.

    All TSVs (evolved + SEEDMIX) share the SAME record order, so we aggregate by
    POSITION (fast numpy), not string keys. Returns record meta
    (chrom,pos,ref_len,alt_len) + p0 (mean over SEEDMIX reps) + one column per
    group `<climate>_g<gen>` (NaN-aware mean) + `<...>_n`. This is the Δp
    foundation: Δp_group = group_col - p0. Cached to results/grenenet_gea/.
    """
    # numpy .npz cache: readable across pandas versions (the pickle cache broke
    # when written by pandas 3.x and read by 2.x — StringDtype pickle mismatch).
    cache_path = f"{GEA}/group_means.npz"
    if cache and os.path.exists(cache_path):
        z = np.load(cache_path, allow_pickle=False)
        names = [str(c) for c in z["_columns"]]
        df = pd.DataFrame({c: z[c] for c in names})
        df["chrom"] = df["chrom"].astype(str)
        return df
    if samples is None:
        samples = list_samples(base)
    m = sample_map(samples)
    meta = pd.read_csv(f"{base}/{samples[0]}.tsv", sep="\t",
                       usecols=["chrom", "pos", "ref_len", "alt_len"])
    N = len(meta)

    # founding p0 = NaN-aware mean over the 8 SEEDMIX reps
    sm = [f for f in glob.glob(f"{SEEDMIX}/SEEDMIX_S*.tsv")
          if "_Chr" not in os.path.basename(f)]
    s_sum = np.zeros(N); s_cnt = np.zeros(N)
    for f in sm:
        a = _af_array(f); ok = np.isfinite(a)
        s_sum[ok] += a[ok]; s_cnt += ok
    out = meta.copy()
    out["p0"] = s_sum / np.where(s_cnt > 0, s_cnt, np.nan)

    # one accumulator per (climate, generation); each sample read once
    keys = list(m.groupby(["climate", "generation"]).groups.keys())
    acc = {k: [np.zeros(N), np.zeros(N), 0] for k in keys}
    for s in samples:
        k = (m.loc[s, "climate"], m.loc[s, "generation"])
        a = _af_array(f"{base}/{s}.tsv"); ok = np.isfinite(a)
        acc[k][0][ok] += a[ok]; acc[k][1] += ok; acc[k][2] += 1
    for (clim, gen), (gsum, gcnt, n) in acc.items():
        out[f"{clim}_g{gen}"] = gsum / np.where(gcnt > 0, gcnt, np.nan)
        out[f"{clim}_g{gen}_n"] = n

    if cache:
        os.makedirs(GEA, exist_ok=True)
        arrs = {}
        for c in out.columns:
            a = out[c].to_numpy()
            arrs[c] = a.astype("U10") if a.dtype == object else a
        np.savez(cache_path, _columns=np.array(list(out.columns)), **arrs)
    return out


def eff_n_founders(sample: str, base: str = OUT) -> dict:
    """Per-chrom effective #founders = 1/sum(h^2) from the saved h vectors.

    The per-chrom driver writes one file per chrom: <sample>_Chr{N}.h_per_chrom.npz
    (global mode), each holding that chrom's h under key 'Chr{N}'.
    """
    out = {}
    for p in sorted(glob.glob(f"{base}/{sample}_Chr*.h_per_chrom.npz")):
        z = np.load(p, allow_pickle=True)
        for k in z.files:
            if k.startswith("Chr"):
                h = z[k]
                out[k] = float(1.0 / np.sum(h ** 2)) if h.size else np.nan
    return out
