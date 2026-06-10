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
# Authoritative GrENE-Net sample table (has usesample, flowerscollected,
# fix_57_generation) — the phase-1 pooling keys. site_gen_plot pools built from
# its usesample==True rows reproduce the phase-1 merged_hapFIRE columns exactly.
SAMPLES_DATA = ("/global/scratch/projects/fc_moilab/projects/grenenet-phase1/"
                "frequency/hapFIRE_frequencies/samples_data_fix57.csv")
# per-site ERA5 bioclim (bio1-19), all 31 cohort sites (+ more); col `site`.
BIOCLIM = ("/global/scratch/projects/fc_moilab/projects/grenenet-phase1/"
           "drive_zenodo/data-intermediate/bioclimvars_experimental_sites_era5.csv")
AF_STORE = f"{GEA}/af_store"             # compact per-sample NPY store (build_af_store.py)
AF_SCALE, AF_NAN = 10000, 65535          # uint16 AF encoding (4-decimal + NaN sentinel)
# TAIR10 gene annotation (Chr1..Chr5, matches our SV `chrom` exactly): one row per
# gene. Sibling *_genes_transposons.gff adds TEs (the adaptive-SV class). Canonical
# copy lives in ~/ara_key_files (stable home location).
ARA_KEYS = "/global/home/users/tbellg/ara_key_files"
TAIR10_GENES = f"{ARA_KEYS}/TAIR10_GFF3_genes.gff"
TAIR10_GENES_TE = f"{ARA_KEYS}/TAIR10_GFF3_genes_transposons.gff"
# phase-1 hapFIRE LD haploblocks: 1.05M SNPs -> 16,674 blocks (id `chrom_idx`,
# e.g. '1_0'), TAIR10 Chr coords. The 12 kendall_* files share one SNP->block map.
LD_BLOCKS = ("/global/scratch/users/tbellg/gea_grene-net/ARCHIVE/"
             "linages_wza_picmin/kendall_0_w_id_n_blocks.csv")

# Pilot sites (extend as the cohort grows). site -> (label, role)
SITE_CLIMATE = {4: "hot", 54: "cold"}   # 4=Cadiz/Madrid region Spain, 54=Cologne DE
BIO_COLS = [f"bio{i}" for i in range(1, 20)]


def load_climate() -> pd.DataFrame:
    """Per-site bio1-19 (ERA5), indexed by int site. The full-GEA climate axis."""
    c = pd.read_csv(BIOCLIM)
    c["site"] = c["site"].astype(int)
    return c.set_index("site")[BIO_COLS]


def decode_af(u: np.ndarray) -> np.ndarray:
    """uint16 AF store -> float32 alt_freq in [0,1]; AF_NAN sentinel -> NaN."""
    f = u.astype(np.float32) / AF_SCALE
    f[u == AF_NAN] = np.nan
    return f


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


def cohort_meta(samples: list[str]) -> pd.DataFrame:
    """Full-GEA per-sample metadata: site/plot/generation/coverage + bio1-19.

    Indexed by sample_id, in the given `samples` order. Joins Table_S5 to the
    per-site ERA5 bioclim — the continuous climate axis for the 31-site GEA
    (supersedes the binary SITE_CLIMATE used by the 2-site pilot).
    """
    t5 = pd.read_csv(T5).set_index("sampleid")
    keep = [s for s in samples if s in t5.index]
    m = t5.loc[keep, ["site", "plot", "date", "year", "generation",
                      "coverage", "weighted_mean_coverage"]].copy()
    m = m.rename(columns={"weighted_mean_coverage": "wmean_cov_plot"})
    m["site"] = m["site"].astype(int)
    return m.join(load_climate(), on="site")


def pool_table(store: str = AF_STORE) -> pd.DataFrame:
    """Timepoint -> site_gen_plot pool map with flower weights, for samples in the
    store. One row per timepoint sample; columns: sample_id, pool, site, plot,
    generation (=fix_57_generation), flowerscollected, coverage. `pool` =
    "{site}_{generation}_{plot}" — the phase-1 merge unit (745 pools)."""
    sd = pd.read_csv(SAMPLES_DATA)
    sd = sd[sd["usesample"]].copy()
    have = {str(s) for s in np.load(f"{store}/samples.npy", allow_pickle=True)}
    sd = sd[sd.sampleid.isin(have)]
    sd["generation"] = sd["fix_57_generation"].astype(int)
    sd["pool"] = (sd.site.astype(str) + "_" + sd.generation.astype(str)
                  + "_" + sd["plot"].astype(str))
    return sd[["sampleid", "pool", "site", "plot", "generation",
               "flowerscollected", "coverage"]].reset_index(drop=True)


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


def load_genes() -> pd.DataFrame:
    """TAIR10 protein-coding genes: chrom, start, end, gene (AT-ID), name.

    chrom is 'Chr1'..'Chr5' (matches the SV records). ChrC/ChrM dropped.
    """
    g = pd.read_csv(TAIR10_GENES, sep="\t", header=None,
                    names=["chrom", "src", "feat", "start", "end",
                           "score", "strand", "frame", "attr"])
    g = g[g.chrom.isin([f"Chr{i}" for i in range(1, 6)])].copy()
    g["gene"] = g.attr.str.extract(r"ID=([^;]+)")
    g["name"] = g.attr.str.extract(r"Name=([^;]+)")
    return g[["chrom", "start", "end", "strand", "gene", "name"]].reset_index(drop=True)


def annotate_svs(df: pd.DataFrame, flank: int = 0, genes: pd.DataFrame | None = None
                 ) -> pd.DataFrame:
    """Tag each SV row (needs `chrom`,`pos`) with the gene(s) it overlaps.

    SV span = [pos, pos + ref_len) when ref_len present, else the point `pos`.
    `flank` bp widens the SV span on each side (e.g. 2000 for nearby/promoter
    hits). Adds columns: gene (first hit AT-ID), gene_name, n_genes (overlaps),
    genes_all (';'-joined). Vectorized per chrom via sorted-interval search.
    """
    if genes is None:
        genes = load_genes()
    out_gene = np.full(len(df), "", dtype=object)
    out_name = np.full(len(df), "", dtype=object)
    out_n = np.zeros(len(df), dtype=int)
    out_all = np.full(len(df), "", dtype=object)
    pos = df["pos"].to_numpy()
    rl = df["ref_len"].to_numpy() if "ref_len" in df else np.ones(len(df), int)
    lo = pos - flank
    hi = pos + np.where(rl > 1, rl, 1) + flank
    idx = df.index.to_numpy()
    for c, gc in genes.groupby("chrom"):
        gs = gc.start.to_numpy(); ge = gc.end.to_numpy()
        gid = gc.gene.to_numpy(); gnm = gc.name.fillna("").to_numpy()
        sel = np.where(df["chrom"].to_numpy() == c)[0]
        for j in sel:
            ov = np.where((gs <= hi[j]) & (ge >= lo[j]))[0]
            if ov.size:
                out_gene[j] = gid[ov[0]]; out_name[j] = gnm[ov[0]]
                out_n[j] = ov.size; out_all[j] = ";".join(gid[ov])
    res = df.copy()
    res["gene"] = out_gene; res["gene_name"] = out_name
    res["n_genes"] = out_n; res["genes_all"] = out_all
    return res


def assign_ld_blocks(chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Assign each record (chrom='Chr1'..'Chr5', pos) to its phase-1 hapFIRE LD
    block, by inheriting the block of the NEAREST genotyped SNP (per chrom).

    Returns a string array of block ids ('chrom_idx'); '' where the chrom has no
    SNP map. SNP density ~1/114 bp, so the nearest-SNP block is reliable.
    """
    b = pd.read_csv(LD_BLOCKS, usecols=["pos", "chrom", "block"])
    out = np.full(len(pos), "", dtype=object)
    pos = np.asarray(pos, dtype=np.int64)
    for ci, g in b.groupby("chrom"):
        m = chrom == f"Chr{ci}"
        if not m.any():
            continue
        g = g.sort_values("pos")
        sp = g.pos.to_numpy(); bl = g.block.to_numpy()
        p = pos[m]
        j = np.clip(np.searchsorted(sp, p), 0, len(sp) - 1)
        jm1 = np.clip(j - 1, 0, len(sp) - 1)
        # pick whichever flanking SNP is closer in bp
        use_prev = (j > 0) & (np.abs(p - sp[jm1]) <= np.abs(sp[j] - p))
        nearest = np.where(use_prev, jm1, j)
        out[np.where(m)[0]] = bl[nearest]
    return out


def collapse_to_blocks(df: pd.DataFrame, stat: str = "z_emp",
                       block_col: str = "block") -> pd.DataFrame:
    """Collapse a per-SV GEA frame to one LEAD SV per LD block (max |stat|).

    Adds n_sv_block (SVs in the block). Use for candidate dedup / Manhattan so a
    single low-recomb block (up to ~900 correlated SVs) can't dominate the list.
    NOTE: lead-by-|stat| is selection-biased — for unbiased GIF/QQ use a
    stat-independent representative instead (e.g. first by pos, or .sample()).
    """
    d = df.copy()
    d["_abs"] = d[stat].abs()
    d["n_sv_block"] = d.groupby(block_col)[stat].transform("size")
    lead = d.sort_values("_abs", ascending=False).groupby(block_col, as_index=False).first()
    return lead.drop(columns="_abs")


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
