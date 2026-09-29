"""Shared helpers for the GrENE-Net SV-GEA analysis (kMate outputs).

The goal (mirroring the phase-1 SNP paper, but for SVs): find structural
variants whose frequency rises in hot sites and falls in cold ones — climate-
dependent directional selection — and ask whether SVs add adaptive signal SNPs
miss. See README.md.

Data:
  - kMate per-sample AF TSVs: analysis/grenenet_selection/common/rerun_kfw_hb/evolved/<MLFH...>.tsv
      (--unit chrom + full-panel Kf_w; old-panel grenenet_kmate_arch3 deleted 2026-07-08)
      cols: chrom pos ref_len alt_len alt_freq info n_called se
  - founding p0 (gen 0): the SEEDMIX kMate outputs (mean over 8 reps).
  - sample metadata (site/plot/generation/coverage): Table_S5.
"""
from __future__ import annotations
import os, glob, json
import numpy as np
import pandas as pd
from scipy import stats

# PROJ derived from this file's own location (analysis/grenenet_selection/lib.py), not
# hardcoded, so the repo tolerates future scratch moves without a repeat of the
# 2026-07-22 break (user-scratch kmate tree migrated to project/fc_moilab scratch;
# the old hardcoded /global/scratch/users/tbellg/kmate silently stopped existing).
PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Repointed 2026-07-07 to the full-panel-Kf_w + haploblock (--unit chrom) rerun.
# Prev: rerun_perfounder (pre-Kf_w-fullpanel, 2026-07-06 — STALE). See
# analysis/grenenet_selection/common/rerun_kfw_hb/README.md.
OUT = f"{PROJ}/analysis/grenenet_selection/common/rerun_kfw_hb/evolved"
SEEDMIX = f"{PROJ}/analysis/grenenet_selection/common/rerun_kfw_hb/seedmix"
GEA = f"{PROJ}/analysis/grenenet_selection"
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
AF_STORE = f"{GEA}/common/results/af_store"             # compact per-sample NPY store (build_af_store.py)
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
# Finer LD blocks from BigLD on the 231-founder panel (kMate mcf90 recompute):
# per-chrom TSVs (chrom, start_pos, end_pos, n_variants), one r2 threshold each.
# clq0.9 = 58,376 blocks (vs 16,674 hapFIRE) -> the finer window definition for the
# phase-1 GEA re-run. Interval-based (unlike hapFIRE's nearest-SNP map), so they do
# NOT tile the genome: variants in inter-block gaps are unassigned (return '').
CLQ_BLOCKS_DIR = f"{GEA}/blocks/results/blocks_mcf90"

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


# QC exclusion: samples with too little USABLE panel data (Chr1 nonzero-k-mer
# fraction < 0.10) — dead/contaminated libraries whose h/AF are garbage. NOT a
# sequencing-depth cut (depth doesn't isolate them; corr(depth,nzfrac)~0.57). See
# analysis/grenenet_selection/qc/results/qc_coverage_audit.csv + analysis/grenenet_selection/notebooks/qc_coverage_audit.ipynb. Applied globally so no pool/
# site/analysis sees them. Set 2026-07-07.
QC_EXCLUDE_FILE = f"{os.path.dirname(os.path.abspath(__file__))}/../../data/qc_lowcov_exclude.txt"
def qc_excluded() -> set[str]:
    """Sample IDs excluded by the low-usable-data QC (see QC_EXCLUDE_FILE)."""
    try:
        with open(QC_EXCLUDE_FILE) as fh:
            return {ln.split("\t")[0].strip() for ln in fh
                    if ln.strip() and not ln.startswith("#")}
    except FileNotFoundError:
        return set()


CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def genome_h(samp: str, base: str) -> np.ndarray | None:
    """genome-wide founder h for one sample = mean over chroms of the per-chrom h.

    Reads the --unit chrom cohort format `{base}/{samp}_{ch}.h_per_chrom.npz`, key
    `{ch}` (a single 231-vector per chrom).
    """
    gs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_per_chrom.npz"
        if not os.path.exists(f):
            return None
        gs.append(np.load(f, allow_pickle=True)[ch].astype(np.float64))
    return np.mean(gs, 0)


def list_samples(base: str = OUT) -> list[str]:
    """Genome-wide per-sample TSVs present in `base` (excludes per-chrom files and
    QC-excluded low-usable-data samples)."""
    excl = qc_excluded()
    return sorted(s for s in (os.path.basename(p)[:-4] for p in glob.glob(f"{base}/*.tsv")
                              if "_Chr" not in os.path.basename(p))
                  if s not in excl)


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
    have -= qc_excluded()                       # drop low-usable-data / contaminated libraries
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

    Returns a Series indexed by rec_key. Cached to analysis/grenenet_selection/.
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
    Cache to analysis/grenenet_selection/pilot_qc.csv so notebooks just load it.
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
    foundation: Δp_group = group_col - p0. Cached to analysis/grenenet_selection/.
    """
    # numpy .npz cache: readable across pandas versions (the pickle cache broke
    # when written by pandas 3.x and read by 2.x — StringDtype pickle mismatch).
    cache_path = f"{GEA}/common/results/group_means.npz"
    if cache and os.path.exists(cache_path):
        z = np.load(cache_path, allow_pickle=False)
        names = [str(c) for c in z["_columns"]]
        df = pd.DataFrame({c: z[c] for c in names})
        df["chrom"] = df["chrom"].astype(str)
        return df
    if samples is None:
        samples = list_samples(base)
    m = sample_map(samples)
    # group_means is the hot/cold 2-site pilot aggregation (the Δp foundation): only
    # climate-classified samples (SITE_CLIMATE: hot/cold sites) contribute. Every other
    # site has climate=NaN and would KeyError the (climate,gen) accumulator below — and
    # groupby() silently drops NaN groups, so they were never meant to be aggregated.
    # (Old cache built clean only because list_samples then returned the pilot alone;
    # it now returns all 2168.) Drop NaN-climate samples up front.
    samples = [s for s in samples if s in m.index and pd.notna(m.loc[s, "climate"])]
    if not samples:
        raise ValueError("group_means: no climate-classified (hot/cold) samples found")
    m = m.loc[samples]
    meta = pd.read_csv(f"{base}/{samples[0]}.tsv", sep="\t",
                       usecols=["chrom", "pos", "ref_len", "alt_len"])
    N = len(meta)

    # founding p0 = NaN-aware mean over the 8 SEEDMIX reps
    sm = [f for f in glob.glob(f"{SEEDMIX}/SEEDMIX_S*.tsv")
          if "_Chr" not in os.path.basename(f)]
    s_sum = np.zeros(N); s_cnt = np.zeros(N)
    for f in sm:
        a = _af_array(f)
        # GUARD: by-position aggregation REQUIRES every TSV on the SAME segregating
        # panel. The stale 06-02 cache was built on the old 10.33M panel (a superset
        # carrying AC=0 monomorphic records) — a record-count mismatch here means an
        # old-panel TSV slipped in and would silently misalign the sum. Fail loudly.
        if len(a) != N:
            raise ValueError(
                f"panel mismatch: {f} has {len(a):,} records, expected {N:,} "
                f"(all TSVs must share the segregating panel order)")
        ok = np.isfinite(a)
        s_sum[ok] += a[ok]; s_cnt += ok
    out = meta.copy()
    out["p0"] = s_sum / np.where(s_cnt > 0, s_cnt, np.nan)

    # one accumulator per (climate, generation); each sample read once
    keys = list(m.groupby(["climate", "generation"]).groups.keys())
    acc = {k: [np.zeros(N), np.zeros(N), 0] for k in keys}
    for s in samples:
        k = (m.loc[s, "climate"], m.loc[s, "generation"])
        a = _af_array(f"{base}/{s}.tsv")
        if len(a) != N:                                    # same panel-alignment guard
            raise ValueError(
                f"panel mismatch: {s} has {len(a):,} records, expected {N:,} "
                f"(segregating panel); refusing to misalign by-position aggregation")
        ok = np.isfinite(a)
        acc[k][0][ok] += a[ok]; acc[k][1] += ok; acc[k][2] += 1
    for (clim, gen), (gsum, gcnt, n) in acc.items():
        out[f"{clim}_g{gen}"] = gsum / np.where(gcnt > 0, gcnt, np.nan)
        out[f"{clim}_g{gen}_n"] = n

    # NOTE on monomorphic sites: the AC=0 panel scaffolding that contaminated the old
    # cache came from the PRE-segregating 10.33M panel. It is excluded simply by reading
    # the segregating-panel TSVs (the input record set has no AC=0 records) — guaranteed
    # upstream, not here. We deliberately do NOT drop "all-absent" rows: on the segregating
    # panel those are founder-segregating variants merely unobserved in the 175-sample
    # hot/cold pilot (p0=0, dp=0) — legitimate panel records, kept so group_means stays
    # one-row-per-segregating-record and any key-join resolves. The alignment assertions
    # above are the real guard against a stale/mixed-panel input.

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


def assign_clq_blocks(chrom: np.ndarray, pos: np.ndarray, r2: float = 0.9) -> np.ndarray:
    """Assign each record (chrom='Chr1'..'Chr5', pos) to its BigLD clq{r2} block.

    Reads the per-chrom interval TSVs in CLQ_BLOCKS_DIR (chrom, start_pos,
    end_pos, n_variants). Unlike the hapFIRE nearest-SNP map, these are LD islands
    that do NOT tile the genome, so a variant landing in an inter-block gap gets ''
    (drop it upstream of WZA). Block ids are 'Chr{n}_{idx}', idx = the block's rank
    by start_pos within its chrom (stable, unique). Returns a string object array.
    """
    out = np.full(len(pos), "", dtype=object)
    pos = np.asarray(pos, dtype=np.int64)
    chrom = np.asarray(chrom, dtype=str)
    tag = f"clq{r2}"
    for ci in range(1, 6):
        c = f"Chr{ci}"
        f = f"{CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv"
        if not os.path.exists(f):
            continue
        g = pd.read_csv(f, sep="\t").sort_values("start_pos").reset_index(drop=True)
        st = g["start_pos"].to_numpy(np.int64)
        en = g["end_pos"].to_numpy(np.int64)
        m = np.where(chrom == c)[0]
        if len(m) == 0:
            continue
        p = pos[m]
        j = np.searchsorted(st, p, "right") - 1          # candidate block (last start <= p)
        inb = (j >= 0) & (p <= en[np.clip(j, 0, len(en) - 1)])
        idx = np.where(inb)[0]
        out[m[idx]] = [f"{c}_{int(j[k])}" for k in idx]
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


def bh(pv: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    m = len(pv); o = pv.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((pv[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def lamgc(pv: np.ndarray) -> float:
    """Genomic-control lambda from a p-value vector (median chi2_1 / expected median)."""
    return float(np.median(stats.chi2.isf(np.clip(pv, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1))


def multisite_gwas_raw(tag: str = "clq90_pc1") -> dict:
    """Load the multisite founder-GWAS Z/C (+recomputed JOINT/GLOBAL/CLIMATE) for further
    C-whitened contrasts (e.g. climate-cluster, locality-index) without a new GWAS run.

    Recomputes chi2_joint/p_joint/z_global/z_clim from the raw Z (M markers x S sites,
    genomic-controlled) and cross-site null covariance C in
    analysis/grenenet_selection/archive/window_hapfreq_retired/hapfreq/multisite_founder_gwas_{tag}.npz, using the identical
    formulas as founder_gwas_multisite.py. Self-checks the recomputed JOINT lambda_GC
    against the saved _meta.json so a silent Z/C misread can't drift unnoticed (the npz
    doesn't store the derived per-marker stats, and the sibling .csv is p_joint-sorted so
    it can't be joined back to Z positionally).
    """
    h = f"{GEA}/hapfreq"
    z = np.load(f"{h}/multisite_founder_gwas_{tag}.npz", allow_pickle=True)
    meta = json.load(open(f"{h}/multisite_founder_gwas_{tag}_meta.json"))
    Z, C, sites, bio1 = z["Z"], z["C"], z["sites"], z["bio1"]
    chrom, start, end, mac = z["chrom"], z["start"], z["end"], z["mac"]
    M, S = Z.shape
    Cinv = np.linalg.pinv(C)
    one = np.ones(S)
    dg = float(one @ Cinv @ one)
    z_global = (Z @ Cinv @ one) / np.sqrt(dg)
    c0 = (bio1 - bio1.mean()) / bio1.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one
    z_clim = (Z @ Cinv @ c) / np.sqrt(float(c @ Cinv @ c))
    chi2_joint = np.einsum("mi,ij,mj->m", Z, Cinv, Z)
    p_joint = stats.chi2.sf(chi2_joint, S)
    lam = lamgc(p_joint)
    assert abs(lam - meta["tests"]["JOINT"]["lam"]) < 1e-3, (
        f"recomputed JOINT lambda {lam:.4f} != saved meta {meta['tests']['JOINT']['lam']:.4f} "
        "-- Z/C recompute drifted from the original GWAS run")
    return dict(Z=Z, C=C, Cinv=Cinv, sites=sites, bio1=bio1, one=one, dg=dg,
                unit=np.array([f"{c_}:{s}-{e}" for c_, s, e in zip(chrom, start, end)]),
                chrom=chrom, mac=mac, chi2_joint=chi2_joint, p_joint=p_joint, q_joint=bh(p_joint),
                z_global=z_global, z_clim=z_clim, M=M, S=S, meta=meta)


# founder MAC count matching build_sv_landscape.py's MAF>=0.05 floor (~5.2% of 231
# founders) -- the bar that already defines the LD blocks and the SV landscape. Default for
# founder_panel_keep below; independent of maf_filter's own min_maf (used by the
# climate-cluster / locality-index marker ranking, currently 1% per project decision).
COMMON_MAC = 12


def maf_filter(raw: dict, min_maf: float = 0.01) -> dict:
    """Restrict a lib.multisite_gwas_raw() dict to markers with founder MAF >= min_maf.

    Rare-founder-allele markers have noisier kinship-corrected z-scores, so any NEW
    per-marker contrast (not already lambda-validated the way JOINT/GLOBAL/CLIMATE were)
    should be ranked only among markers clearing a standard MAF floor (default 1%, the
    common GWAS QC convention). min_maf is converted to a founder MAC via the cohort's
    n_founders (in raw["meta"]) so it's independent of the founder count. Per-site fields
    (S-length: sites/bio1/one/C/Cinv) are untouched; per-marker fields (M-length) are
    masked and M is updated.
    """
    min_mac = int(np.ceil(min_maf * raw["meta"]["n_founders"]))
    keep = raw["mac"] >= min_mac
    out = dict(raw)
    for k in ("Z", "unit", "chrom", "mac", "chi2_joint", "p_joint", "q_joint", "z_global", "z_clim"):
        out[k] = raw[k][keep]
    out["M"] = int(keep.sum())
    return out


def founder_panel_keep(chrom: np.ndarray, pos: np.ndarray, min_mac: int = COMMON_MAC,
                       called_min: float = 0.9) -> np.ndarray:
    """Founder-panel MAC/call-rate mask for arbitrary (chrom,pos) records.

    Same floor as COMMON_MAC / build_sv_landscape.py (MAF>=0.05, called-frac>=0.9 among
    the 231 founders) -- the bar that already defines the LD blocks and the SV landscape.
    Use to floor any OTHER per-variant catalog (e.g. the raw per-sample/pool kMate variant
    calls behind site_variant_temporal_scoef.py) before counting SV/indel/SNP composition
    or scoring selection on it: unlike the founder panel, that catalog carries no MAF floor
    of its own (only a noisy pool-level p0 reachability check), so >50% of its SV-class
    calls are literal founder singletons -- see KINSHIP_TEMPORAL_METHODS_RESEARCH context.
    Positions absent from the founder panel entirely do not pass (conservative: no MAC to
    vouch for them). Returns a bool array aligned to (chrom,pos), True where at least one
    panel record at that exact position passes.
    """
    import scipy.sparse as sp
    out = np.zeros(len(pos), bool)
    chrom = np.asarray(chrom); pos = np.asarray(pos, dtype=np.int64)
    for ch in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        m = chrom == ch
        if not m.any():
            continue
        cl = ch.lower()
        base = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        lpos = meta["pos"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz"); vc = sp.load_npz(f"{base}.var_called.npz")
        n_alt = np.asarray(vp.sum(0)).ravel(); n_cal = np.asarray(vc.sum(0)).ravel()
        passmaf = (n_alt >= min_mac) & (n_alt <= vp.shape[0] - min_mac) & (n_cal >= called_min * vp.shape[0])
        order = np.argsort(lpos); lpos_s = lpos[order]; pass_s = passmaf[order]
        p = pos[m]
        j = np.searchsorted(lpos_s, p)
        found = (j < len(lpos_s)) & (lpos_s[np.clip(j, 0, len(lpos_s) - 1)] == p)
        res = np.zeros(len(p), bool)
        res[found] = pass_s[np.clip(j[found], 0, len(lpos_s) - 1)]
        out[m] = res
    return out


def matched_perm_test(values: np.ndarray, bins: np.ndarray, idx: np.ndarray,
                      nperm: int = 10000, seed: int = 0) -> tuple[float, float, float, float]:
    """Size-matched permutation test (the design used throughout sv_adaptive/): is
    mean(values[idx]) higher than expected from len(idx) draws matched on `bins` (e.g.
    block-size bin)?

    Each idx item's null draws come from its OWN bin, `nperm` times, so bigger/smaller
    items are compared to same-size controls. Returns (observed, null_mean, obs/null,
    one-sided p = P(null_mean_over_draws >= observed)).
    """
    rng = np.random.default_rng(seed)
    bin_members = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    obs = values[idx].mean()
    null = np.empty((len(idx), nperm))
    for k, i in enumerate(idx):
        null[k] = rng.choice(values[bin_members[bins[i]]], nperm)
    nm = null.mean(0)
    p = (1 + (nm >= obs).sum()) / (nperm + 1)
    ratio = obs / nm.mean() if nm.mean() != 0 else np.nan
    return float(obs), float(nm.mean()), float(ratio), float(p)


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
