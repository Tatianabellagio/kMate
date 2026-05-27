"""Generate recombinant mosaic founder FASTAs for the recomb pool simulation.

Each "individual" is a haploid recombinant of multiple founder genomes after
n_generations of meiotic recombination at a Poisson-distributed crossover rate
(default 4 cM/Mb, the A. thaliana literature consensus). Output:

  out_dir/
    pool_weights.tsv             (founder, count, weight) — for record-keeping
    haps/
      s_ind001/h1.fa             (1.fai)  — mosaic FASTA + index for VISOR
      s_ind002/h1.fa
      ...
    ancestry.tsv                 (ind_id, chrom, start, end, founder)

Design notes:
  - Initial gen 0: n_indiv haploids, founder identity sampled by pool_weights
  - Each generation: pair every individual with a random other one; offspring
    is a recombinant haploid (Poisson crossover positions).
  - After n_generations, n_indiv haploid mosaics remain.
  - Mosaic FASTA is built by stitching founder segments at recombination
    boundaries. We use the founder FASTAs in `cactus_dir/<founder>.chr.fa`.

Truth (computed in companion script): per-record AF =
  sum_i cn_var[founder_at_(chrom,pos)_in_ind_i, record] × weight_i / sum_i weight_i

This collapses to the same expression cactus_em projects through
(founder_freq × cn_var), but with founder_freq varying per genomic window.
"""
from __future__ import annotations
import argparse, os, random, subprocess
from pathlib import Path
import numpy as np
import pandas as pd


# ----- hand-rolled .fai-indexed FASTA reader (no pyfaidx dependency) -------
def load_fai(fai_path):
    """Return dict {seqname: (length, offset, line_blen, line_len)} from .fai."""
    out = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            name, length, offset, lblen, llen = parts[:5]
            out[name] = (int(length), int(offset), int(lblen), int(llen))
    return out


def fa_fetch(fa_path, fai, name, start, end):
    """Fetch [start, end] (1-based inclusive) from FASTA at fa_path using .fai."""
    length, offset, lblen, llen = fai[name]
    s, e = start - 1, end  # 0-based [s, e)
    s_off = offset + (s // lblen) * llen + (s % lblen)
    e_off = offset + (e // lblen) * llen + (e % lblen)
    with open(fa_path, 'rb') as f:
        f.seek(s_off)
        raw = f.read(e_off - s_off)
    return raw.decode('ascii').replace('\n', '').replace('\r', '').upper()


class IndexedFasta:
    """Minimal pyfaidx replacement: __getitem__(name) returns a slicer, .close()."""
    def __init__(self, fa_path):
        self.path = fa_path
        fai_path = fa_path + ".fai"
        if not os.path.exists(fai_path):
            raise FileNotFoundError(f"missing .fai for {fa_path}")
        self.fai = load_fai(fai_path)

    def fetch(self, name, start, end):
        return fa_fetch(self.path, self.fai, name, start, end)
# ---------------------------------------------------------------------------

# A. thaliana literature: ~3.6-4 cM/Mb ≈ 4e-8 crossover/bp/meiosis
DEFAULT_RECOMB_RATE = 4e-8

CHROM_LENGTHS = {
    "Chr1": 30427671,
    "Chr2": 19698289,
    "Chr3": 23459830,
    "Chr4": 18585056,
    "Chr5": 26975502,
}


def sample_crossovers(chrom_len, rate, rng, allowed_positions=None):
    """Sample crossover positions on one chromosome via Poisson process.

    If `allowed_positions` is None: positions sampled uniformly along the chrom
    (the original behaviour, useful as a stress-test).

    If `allowed_positions` is an array of bp positions (e.g., LD-block
    boundaries from BigLD on the panel): the Poisson-distributed *number* of
    crossovers per chromosome is preserved, but each crossover is placed at a
    randomly chosen allowed position. This models real meiosis where crossover
    hotspots define LD block boundaries — within-block ancestry is preserved.
    """
    n = rng.poisson(chrom_len * rate)
    if n == 0:
        return np.array([], dtype=np.int64)
    if allowed_positions is None:
        pos = np.sort(rng.integers(1, chrom_len, size=n))
    else:
        if len(allowed_positions) == 0:
            return np.array([], dtype=np.int64)
        # Sample with replacement from allowed positions, then sort + dedup.
        # Replacement is fine: at our rates each chrom has ~1-4 crossovers and
        # ~1000s of allowed positions, so collisions are vanishingly rare.
        chosen = rng.choice(allowed_positions, size=n, replace=True)
        pos = np.unique(np.sort(chosen))
    return pos


def load_ld_block_boundaries(path):
    """Load LD-block boundaries from a hapfire_block_index.npz (build_hapfire_block_index.py output).

    Returns dict {chrom: sorted np.array of bp positions where crossovers can fire}.
    Boundaries are taken as block_pos_end of each block (= the last SNP in the
    block, just before the next block starts). Chrom names are converted from
    panel VCF style ('1','2',...) to FASTA style ('Chr1','Chr2',...).
    """
    npz = np.load(path, allow_pickle=True)
    block_chrom = np.asarray(npz['block_chrom']).astype(str)
    block_pos_end = np.asarray(npz['block_pos_end']).astype(np.int64)
    out = {}
    for raw_chrom in np.unique(block_chrom):
        m = block_chrom == raw_chrom
        # Convert '1' → 'Chr1' to match CHROM_LENGTHS keys
        target = f"Chr{raw_chrom}" if not str(raw_chrom).startswith("Chr") else str(raw_chrom)
        # Filter to positions strictly inside the chrom (not at ends — those don't
        # produce ancestry switches anyway)
        L = CHROM_LENGTHS.get(target, None)
        ends = block_pos_end[m]
        if L is not None:
            ends = ends[(ends > 0) & (ends < L)]
        out[target] = np.sort(np.unique(ends))
    return out


def make_one_mosaic(parent_a_id, parent_b_id, rate, rng, allowed_per_chrom=None):
    """Build segment list for a 1-gen recombinant of (parent_a, parent_b).

    Returns dict {chrom: [(start, end, founder_id), ...]}.
    """
    out = {}
    for chrom, L in CHROM_LENGTHS.items():
        chrom_allowed = allowed_per_chrom.get(chrom) if allowed_per_chrom else None
        cuts = sample_crossovers(L, rate, rng, allowed_positions=chrom_allowed)
        # Walk left-to-right alternating founders (random starting parent)
        which = rng.integers(0, 2)
        boundaries = [0] + list(cuts) + [L]
        segs = []
        for i in range(len(boundaries) - 1):
            s, e = boundaries[i], boundaries[i + 1]
            f = parent_a_id if which == 0 else parent_b_id
            segs.append((s + 1, e, f))     # 1-based [start, end] inclusive
            which = 1 - which
        out[chrom] = segs
    return out


def make_individuals(n_indiv, n_generations, weights_dict, rate, rng, allowed_per_chrom=None):
    """Build n_indiv recombinant haploids after n_generations.

    Returns list of length n_indiv, each entry is dict {chrom: [(start,end,founder),...]}.
    """
    founders = list(weights_dict.keys())
    fweights = np.array([weights_dict[f] for f in founders], dtype=np.float64)
    fweights /= fweights.sum()

    # Generation 0: each individual is one founder (sampled by weight)
    pop = []
    for _ in range(n_indiv):
        f = rng.choice(founders, p=fweights)
        # whole-chromosome single-founder
        pop.append({c: [(1, L, f)] for c, L in CHROM_LENGTHS.items()})

    # n_generations of pairing + recombination
    for gen in range(n_generations):
        new_pop = []
        for _ in range(n_indiv):
            # pick two parents at random from pop
            i, j = rng.integers(0, len(pop)), rng.integers(0, len(pop))
            child = recombine_segments(pop[i], pop[j], rate, rng, allowed_per_chrom)
            new_pop.append(child)
        pop = new_pop

    return pop


def recombine_segments(parent_a, parent_b, rate, rng, allowed_per_chrom=None):
    """Recombine two segment-track parents into one haploid offspring."""
    out = {}
    for chrom, L in CHROM_LENGTHS.items():
        chrom_allowed = allowed_per_chrom.get(chrom) if allowed_per_chrom else None
        cuts = sample_crossovers(L, rate, rng, allowed_positions=chrom_allowed)
        which = rng.integers(0, 2)
        boundaries = [0] + list(cuts) + [L]
        segs = []
        for i in range(len(boundaries) - 1):
            s, e = boundaries[i] + 1, boundaries[i + 1]
            parent = parent_a if which == 0 else parent_b
            # pull all parent segments overlapping [s,e]
            for ps, pe, f in parent[chrom]:
                ov_s = max(s, ps); ov_e = min(e, pe)
                if ov_s <= ov_e:
                    segs.append((ov_s, ov_e, f))
            which = 1 - which
        # Merge adjacent same-founder segments
        merged = []
        for s, e, f in segs:
            if merged and merged[-1][2] == f and merged[-1][1] + 1 == s:
                merged[-1] = (merged[-1][0], e, f)
            else:
                merged.append((s, e, f))
        out[chrom] = merged
    return out


def write_ancestry_tsv(pop, out_path):
    rows = []
    for i, ind in enumerate(pop):
        ind_id = f"ind{i+1:03d}"
        for chrom in sorted(ind):
            for s, e, f in ind[chrom]:
                rows.append((ind_id, chrom, s, e, f))
    df = pd.DataFrame(rows, columns=["ind_id", "chrom", "start", "end", "founder"])
    df.to_csv(out_path, sep="\t", index=False)
    return df


def write_mosaic_fasta(ind_segs, founder_fastas, out_fa, samtools=None):
    """Stitch founder FASTA segments into one mosaic FASTA, indexed via samtools faidx."""
    with open(out_fa, "w") as fout:
        for chrom in sorted(ind_segs):
            fout.write(f">{chrom}\n")
            seqs = []
            for s, e, f in ind_segs[chrom]:
                seq = founder_fastas[f].fetch(chrom, s, e)
                seqs.append(seq)
            full = "".join(seqs)
            for i in range(0, len(full), 80):
                fout.write(full[i:i+80] + "\n")
    if samtools:
        rc = subprocess.run([samtools, "faidx", out_fa], capture_output=True)
        if rc.returncode != 0:
            raise RuntimeError(f"samtools faidx failed:\n{rc.stderr.decode()}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-indiv", type=int, default=50)
    ap.add_argument("--n-generations", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--recomb-rate", type=float, default=DEFAULT_RECOMB_RATE)
    ap.add_argument("--cactus-dir", required=True,
                    help="dir holding <founder>.chr.fa[.fai]")
    ap.add_argument("--founders-meta", required=True,
                    help="poolfreq cn_full meta.npz with 'founders' field")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samtools", default="/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/samtools",
                    help="samtools binary, used to index the mosaic FASTAs (default: shared install)")
    ap.add_argument("--crossovers-from-ld-blocks", default=None,
                    help="Path to a hapfire_block_index.npz (output of build_hapfire_block_index.py). "
                         "If given, sample crossover positions only at LD-block boundaries (= recombination "
                         "hotspots). If omitted, fall back to uniform-position Poisson sampling.")
    ap.add_argument("--chroms", default=None,
                    help="Space-separated chroms to include (e.g. 'Chr1'). Default: all 5. "
                         "Use this for Chr1-only sims when founder FASTAs only have Chr1.")
    ap.add_argument("--source-weights", default=None,
                    help="Optional TSV (founder, prob) defining the source-population founder "
                         "probabilities. The gen-0 pool is then a Stage-1 multinomial draw of "
                         "N_indiv individuals from these probabilities. Default = uniform 1/F "
                         "(neutral source). Use to model differential founder fitness ('selection') "
                         "by passing skewed probabilities. Probabilities are renormalized to sum to 1.")
    args = ap.parse_args()

    # Filter CHROM_LENGTHS in place so make_individuals / recombine_segments /
    # write_mosaic_fasta all see the restricted set.
    if args.chroms is not None:
        keep = set(args.chroms.split())
        bad = keep - set(CHROM_LENGTHS)
        if bad:
            raise ValueError(f"unknown chroms: {bad}; valid keys: {list(CHROM_LENGTHS)}")
        for c in list(CHROM_LENGTHS):
            if c not in keep:
                del CHROM_LENGTHS[c]
        print(f"  --chroms {args.chroms} -> using only: {list(CHROM_LENGTHS)}", flush=True)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    haps_dir = out / "haps"
    haps_dir.mkdir(exist_ok=True)

    meta = np.load(args.founders_meta, allow_pickle=True)
    founders = list(meta["founders"])
    rng = np.random.default_rng(args.seed)

    # Two-stage pool-seq model (2026-05-21):
    #   Source population: defined by founder probabilities p_f (uniform by default;
    #                      can be skewed via --source-weights to model selection).
    #   Pool (Stage 1):    Multinomial(N_indiv, p_f) draw → realized founder counts.
    # The source weights are saved so compute_recomb_truth.py can produce the
    # "infinite-pool" expectation alongside the realized-pool truth.
    if args.source_weights is not None:
        sw = pd.read_csv(args.source_weights, sep="\t")
        if not {"founder", "prob"}.issubset(sw.columns):
            raise ValueError(
                f"--source-weights TSV must have columns 'founder' and 'prob'; "
                f"got {list(sw.columns)}"
            )
        prob_lookup = dict(zip(sw["founder"].astype(str), sw["prob"].astype(float)))
        source_probs = np.array([prob_lookup.get(str(f), 0.0) for f in founders],
                                dtype=np.float64)
        missing = [str(f) for f in founders if str(f) not in prob_lookup]
        if missing:
            print(f"  WARN: {len(missing)} founders not in --source-weights "
                  f"(setting p=0): {missing[:5]}{'...' if len(missing)>5 else ''}",
                  flush=True)
        if source_probs.sum() <= 0:
            raise ValueError("--source-weights probabilities sum to 0")
        source_probs = source_probs / source_probs.sum()
        print(f"  source population: --source-weights {args.source_weights}  "
              f"(max p={source_probs.max():.4f}, min p={source_probs[source_probs>0].min():.6f}, "
              f"n founders with p>0: {int((source_probs>0).sum())}/{len(founders)})",
              flush=True)
    else:
        source_probs = np.full(len(founders), 1.0 / len(founders), dtype=np.float64)
        print(f"  source population: uniform 1/{len(founders)} (neutral, default)",
              flush=True)

    # Save the source-population definition (founder, prob). Constant across records,
    # this is the population parameter we'd estimate with an infinite pool.
    source_df = pd.DataFrame({
        "founder": [str(f) for f in founders],
        "prob":    source_probs,
    })
    source_df.to_csv(out / "source_weights.tsv", sep="\t", index=False)

    # Stage 1: multinomial draw of N_indiv individuals from source population.
    counts = rng.multinomial(args.n_indiv, source_probs)
    total = int(counts.sum())
    weights_dict = {str(f): int(c) for f, c in zip(founders, counts) if c > 0}
    rows = [(str(f), int(c), float(c) / total) for f, c in zip(founders, counts)]
    weights_df = pd.DataFrame.from_records(rows, columns=["founder", "count", "weight"])
    weights_df.to_csv(out / "pool_weights.tsv", sep="\t", index=False)
    n_used = (counts > 0).sum()
    print(f"  pool (Stage 1): {args.n_indiv} individuals drawn from source, seed={args.seed}", flush=True)
    print(f"  n founders with count>0 in pool: {n_used} ({n_used/len(founders)*100:.1f}%)", flush=True)
    print(f"  recomb rate: {args.recomb_rate:.2e} per bp ({args.n_generations} generation{'s' if args.n_generations != 1 else ''})", flush=True)

    # Optional: load LD-block boundaries to constrain crossover positions
    allowed_per_chrom = None
    if args.crossovers_from_ld_blocks:
        allowed_per_chrom = load_ld_block_boundaries(args.crossovers_from_ld_blocks)
        n_allowed = sum(len(v) for v in allowed_per_chrom.values())
        print(f"  crossover model: hotspot-aligned (LD-block boundaries from "
              f"{args.crossovers_from_ld_blocks})", flush=True)
        print(f"    n_allowed_positions: {n_allowed:,} across "
              f"{len(allowed_per_chrom)} chroms; per-chrom counts: "
              f"{ {c: len(v) for c, v in allowed_per_chrom.items()} }", flush=True)
    else:
        print(f"  crossover model: uniform-position Poisson", flush=True)

    # Build mosaic individuals
    pop = make_individuals(args.n_indiv, args.n_generations, weights_dict,
                           args.recomb_rate, rng,
                           allowed_per_chrom=allowed_per_chrom)
    write_ancestry_tsv(pop, out / "ancestry.tsv")
    print(f"  wrote ancestry.tsv: {sum(len(p[c]) for p in pop for c in p):,} segments", flush=True)

    # Pre-load founder FASTAs once (only those actually used)
    used_founders = set(f for ind in pop for chrom in ind for _, _, f in ind[chrom])
    founder_fastas = {}
    cactus = Path(args.cactus_dir)
    for f in used_founders:
        fa_path = cactus / f"{f}.chr.fa"
        if not fa_path.exists():
            raise FileNotFoundError(fa_path)
        founder_fastas[f] = IndexedFasta(str(fa_path))
    print(f"  loaded {len(founder_fastas)} founder FASTAs", flush=True)

    # Stitch mosaics
    for i, ind in enumerate(pop):
        ind_id = f"ind{i+1:03d}"
        hap_dir = haps_dir / f"s_{ind_id}"
        hap_dir.mkdir(exist_ok=True)
        out_fa = hap_dir / "h1.fa"
        if out_fa.exists() and (hap_dir / "h1.fa.fai").exists():
            continue
        write_mosaic_fasta(ind, founder_fastas, str(out_fa), samtools=args.samtools)
        if (i + 1) % 10 == 0:
            print(f"  built {i+1}/{args.n_indiv} mosaic FASTAs", flush=True)
    print(f"  done — {args.n_indiv} mosaics in {haps_dir}/", flush=True)


if __name__ == "__main__":
    main()
