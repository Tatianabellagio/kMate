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
    These are consensus sequences whose coordinates drift from TAIR10, so each
    segment (TAIR10 coordinates, as in ancestry.tsv and the truth) is cut at the
    founder's own consensus coordinates via `chain_dir/<founder>.chain`; mosaic
    lengths therefore differ from TAIR10 (a pure founder = its whole consensus).

Truth (computed in companion script): per-record AF =
  sum_i var_pa[founder_at_(chrom,pos)_in_ind_i, record] × weight_i / sum_i weight_i

This collapses to the same expression cactus_em projects through
(founder_freq × var_pa), but with founder_freq varying per genomic window.
"""
from __future__ import annotations
import argparse, os, subprocess, sys, shutil
from pathlib import Path
import numpy as np
import pandas as pd


def _default_samtools():
    """Resolve samtools beside the running interpreter (so the conda env that
    runs this script supplies its own samtools), else fall back to $PATH. Avoids
    hard-coding an env path that breaks on cluster/env migrations."""
    cand = os.path.join(os.path.dirname(sys.executable), "samtools")
    if os.path.exists(cand):
        return cand
    return shutil.which("samtools") or "samtools"


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
    """Fetch [start, end] (1-based inclusive, within the sequence) from FASTA at fa_path using .fai."""
    length, offset, lblen, llen = fai[name]
    s = start - 1                      # 0-based inclusive start
    assert 0 <= s < end <= length, (fa_path, name, start, end, length)
    s_off = offset + (s // lblen) * llen + (s % lblen)
    e_off = offset + (end // lblen) * llen + (end % lblen)
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


class ConsensusChain:
    """TAIR10 -> founder-consensus coordinates for one chromosome, from the chain that `bcftools consensus -c` writes
    alongside the founder FASTA (t = TAIR10, q = consensus; ungapped blocks, 0-based).

    Founder FASTAs are consensus sequences (indels and SVs applied), so their coordinates drift from TAIR10 by up to
    ~0.4 Mb. Ancestry tracks and the truth are in TAIR10 coordinates, so each segment must be cut from the founder at
    the consensus bases that correspond to its TAIR10 span (sims/AUDIT_2026-06-18.md #3; benchmarks/p80/audit_sims/).
    """
    def __init__(self, chain_path, chrom, t_len, q_len):
        ts, qs, sz, found, on = [], [], [], False, False
        with open(chain_path) as fh:
            for line in fh:
                p = line.split()
                if not p:
                    continue
                if p[0] == "chain":
                    on = p[2] == chrom
                    if on:
                        if found:
                            raise ValueError(f"{chain_path}: more than one chain for {chrom}")
                        found = True
                        if int(p[3]) != t_len or int(p[8]) != q_len or p[4] != "+" or p[9] != "+":
                            raise ValueError(f"{chain_path} {chrom}: chain sizes {p[3]}/{p[8]} do not match TAIR10 "
                                             f"{t_len} / FASTA {q_len} (or not + strand)")
                        t, q = int(p[5]), int(p[10])
                    continue
                if not on:
                    continue
                ts.append(t); qs.append(q); sz.append(int(p[0]))
                if len(p) == 3:
                    t += int(p[0]) + int(p[1]); q += int(p[0]) + int(p[2])
        if not found:
            raise ValueError(f"{chain_path}: no chain for {chrom}")
        self.ts, self.qs, self.sz = np.array(ts), np.array(qs), np.array(sz)
        self.t_len, self.q_len = t_len, q_len

    def cut(self, t):
        """0-based consensus index where a piece starting at TAIR10 base t (1-based; t_len+1 = past the end) begins.
        Consecutive pieces of one founder tile its consensus exactly: bases inserted between TAIR10 t-1 and t, and a
        TAIR10 base t deleted in the founder, both resolve to the start of the next aligned block, so inserted bases
        stay with the piece on the left."""
        if t > self.t_len:
            return self.q_len
        t0 = t - 1
        i = int(np.searchsorted(self.ts, t0, side="right")) - 1
        if i < 0:
            return 0
        if t0 < self.ts[i] + self.sz[i]:
            return int(self.qs[i] + (t0 - self.ts[i]))
        return int(self.qs[i + 1]) if i + 1 < len(self.qs) else self.q_len

    def piece(self, fasta, chrom, s, e):
        """Consensus bases of TAIR10 span [s, e] (1-based inclusive)."""
        a, b = self.cut(s), self.cut(e + 1)
        if b <= a:
            return ""
        assert b <= self.q_len
        return fasta.fetch(chrom, a + 1, b)

# A. thaliana literature: ~3.6-4 cM/Mb ≈ 4e-8 crossover/bp/meiosis
DEFAULT_RECOMB_RATE = 4e-8

CHROM_LENGTHS = {
    "Chr1": 30427671,
    "Chr2": 19698289,
    "Chr3": 23459830,
    "Chr4": 18585056,
    "Chr5": 26975502,
}


def sample_crossovers(chrom_len, rate, rng):
    """Sample crossover positions on one chromosome via a Poisson process.

    The number of crossovers ~ Poisson(chrom_len * rate); each crossover is
    placed at a uniformly random position along the chromosome.
    """
    n = rng.poisson(chrom_len * rate)
    if n == 0:
        return np.array([], dtype=np.int64)
    return np.sort(rng.integers(1, chrom_len, size=n))


def make_individuals(n_indiv, n_generations, weights_dict, rate, rng,
                     gen0_no_replace=False, selfing_rate=0.0):
    """Build n_indiv recombinant haploids after n_generations.

    gen0_no_replace: when True, draw gen-0 founders to maximize uniformity of
    pool composition (SEEDMIX mimicry):
      * n_indiv <= F: sample WITHOUT replacement (each founder ≤ 1 time);
        truth_h = uniform 1/F for the n_indiv chosen, 0 for the rest.
      * n_indiv >  F: BALANCED allocation — each founder gets count in
        {floor(n/F), ceil(n/F)} as evenly as possible (avoids Poisson(n/F)
        quantization spread). Truth_h then has only 1-2 distinct nonzero values.
    Default False = legacy multinomial with replacement.

    Returns list of length n_indiv, each entry is dict {chrom: [(start,end,founder),...]}.
    """
    founders = list(weights_dict.keys())
    fweights = np.array([weights_dict[f] for f in founders], dtype=np.float64)
    fweights /= fweights.sum()
    F = len(founders)

    # Generation 0: each individual is one founder.
    if gen0_no_replace:
        if n_indiv <= F:
            chosen = rng.choice(founders, size=n_indiv, replace=False, p=fweights)
        else:
            base = n_indiv // F            # floor
            extra = n_indiv - base * F     # this many founders get +1
            counts = np.full(F, base, dtype=int)
            extra_idx = rng.choice(F, size=extra, replace=False)
            counts[extra_idx] += 1
            chosen = np.repeat(founders, counts)
            rng.shuffle(chosen)            # randomize ind_id order
        pop = [{c: [(1, L, f)] for c, L in CHROM_LENGTHS.items()} for f in chosen]
    else:
        pop = []
        for _ in range(n_indiv):
            f = rng.choice(founders, p=fweights)
            # whole-chromosome single-founder
            pop.append({c: [(1, L, f)] for c, L in CHROM_LENGTHS.items()})

    # n_generations of reproduction. Each offspring either SELFS (prob
    # selfing_rate) or OUTCROSSES (prob 1 - selfing_rate).
    #   - Selfing: A. thaliana founders are inbred/homozygous and this framework
    #     carries one haploid mosaic per individual (no within-individual het;
    #     §2.3/caveat #3). A self-cross of a homozygote reproduces the parent
    #     exactly — recombination between two identical haplotypes is a no-op —
    #     so selfing = a clonal copy of one random parent's ancestry track.
    #   - Outcrossing: pair two random parents and recombine (legacy behaviour).
    # selfing_rate=0.0 (default) never draws the extra random, so existing sims
    # reproduce bit-identically.
    for gen in range(n_generations):
        new_pop = []
        for _ in range(n_indiv):
            if selfing_rate > 0.0 and rng.random() < selfing_rate:
                # SELF: clone one random parent (homozygous → identical offspring)
                i = rng.integers(0, len(pop))
                child = {c: list(segs) for c, segs in pop[i].items()}
            else:
                # OUTCROSS: pick two parents at random from pop
                i, j = rng.integers(0, len(pop)), rng.integers(0, len(pop))
                child = recombine_segments(pop[i], pop[j], rate, rng)
            new_pop.append(child)
        pop = new_pop

    return pop


def recombine_segments(parent_a, parent_b, rate, rng):
    """Recombine two segment-track parents into one haploid offspring."""
    out = {}
    for chrom, L in CHROM_LENGTHS.items():
        cuts = sample_crossovers(L, rate, rng)
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


def write_mosaic_fasta(ind_segs, founder_fastas, out_fa, samtools=None, chains=None):
    """Stitch founder FASTA segments into one mosaic FASTA, indexed via samtools faidx.

    chains: {(founder, chrom): ConsensusChain}. Each TAIR10 segment is cut from the founder at its own consensus
    coordinates, so a whole-chromosome founder contributes its whole consensus and mosaic lengths vary. None = founder
    genomes in TAIR10 coordinates (SNP-only consensus; main() checks every length equals TAIR10)."""
    with open(out_fa, "w") as fout:
        for chrom in sorted(ind_segs):
            fout.write(f">{chrom}\n")
            seqs = []
            for s, e, f in ind_segs[chrom]:
                if chains is None:
                    seq = founder_fastas[f].fetch(chrom, s, e)
                else:
                    seq = chains[(f, chrom)].piece(founder_fastas[f], chrom, s, e)
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
                    help="kmer_pa meta.npz with 'founders' field")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samtools", default=_default_samtools(),
                    help="samtools binary used to index the mosaic FASTAs "
                         "(default: samtools beside the running python, else $PATH)")
    ap.add_argument("--chroms", default=None,
                    help="Space-separated chroms to include (e.g. 'Chr1'). Default: all 5. "
                         "Use this for Chr1-only sims when founder FASTAs only have Chr1.")
    ap.add_argument("--gen0-no-replace", action="store_true",
                    help="SEEDMIX mimicry: gen-0 sampling without replacement (n<=F) "
                         "or balanced allocation (n>F) — each founder gets floor(n/F) "
                         "or ceil(n/F) individuals. Truth_h becomes uniform 1/F or "
                         "the closest possible. Default = legacy multinomial w/ replace.")
    ap.add_argument("--selfing-rate", type=float, default=0.0,
                    help="Per-offspring probability of SELFING (clonal copy of one parent) "
                         "vs OUTCROSSING (recombinant of two parents) at each generation. "
                         "A. thaliana is >97%% selfing (Exposito-Alonso et al. 2018), so "
                         "--selfing-rate 0.97 throttles effective recombination ~30x and "
                         "keeps most individuals pure single-founder even after g3. "
                         "Default 0.0 = forced outcrossing (legacy; existing sims reproduce "
                         "bit-identically). Only takes effect when --n-generations >= 1.")
    ap.add_argument("--chain-dir", default=None,
                    help="dir holding <founder>.chain (bcftools consensus -c, written with the founder FASTAs): segments "
                         "are cut at each founder's own consensus coordinates. Omit only for founder genomes in TAIR10 "
                         "coordinates (SNP-only consensus); their lengths are checked.")
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

    # Stage 1: founder count assignment.
    # Default = multinomial (legacy); when --gen0-no-replace is set, use
    # no-replacement (n<=F) or balanced allocation (n>F) so truth_h is uniform.
    F = len(founders)
    if args.gen0_no_replace:
        if args.n_indiv <= F:
            chosen = rng.choice(F, size=args.n_indiv, replace=False, p=source_probs)
            counts = np.zeros(F, dtype=int); counts[chosen] = 1
        else:
            base = args.n_indiv // F
            extra = args.n_indiv - base * F
            counts = np.full(F, base, dtype=int)
            extra_idx = rng.choice(F, size=extra, replace=False)
            counts[extra_idx] += 1
        print(f"  gen-0 sampling: NO-REPLACE/balanced ({'sampled n<=F' if args.n_indiv<=F else 'balanced n>F'})", flush=True)
    else:
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

    print(f"  crossover model: uniform-position Poisson (random)", flush=True)

    # Build mosaic individuals
    if not (0.0 <= args.selfing_rate <= 1.0):
        raise ValueError(f"--selfing-rate must be in [0, 1]; got {args.selfing_rate}")
    if args.selfing_rate > 0.0:
        print(f"  reproduction: SELFING rate {args.selfing_rate:.3f} "
              f"(outcross {1-args.selfing_rate:.3f}); selfing = clonal copy", flush=True)
    else:
        print(f"  reproduction: forced outcrossing (selfing rate 0.0, legacy)", flush=True)

    pop = make_individuals(args.n_indiv, args.n_generations, weights_dict,
                           args.recomb_rate, rng,
                           gen0_no_replace=args.gen0_no_replace,
                           selfing_rate=args.selfing_rate)

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

    chains = None
    if not args.chain_dir:
        for f in used_founders:
            for c, L in CHROM_LENGTHS.items():
                if founder_fastas[f].fai[c][0] != L:
                    raise ValueError(f"{f} {c}: {founder_fastas[f].fai[c][0]:,} bp, not the TAIR10 {L:,}: "
                                     f"a consensus with indels needs --chain-dir")
        print("  no --chain-dir: founder genomes in TAIR10 coordinates (lengths checked)", flush=True)
    else:
        chains = {}
        for f in used_founders:
            for c, L in CHROM_LENGTHS.items():
                chains[(f, c)] = ConsensusChain(Path(args.chain_dir) / f"{f}.chain", c, L,
                                                founder_fastas[f].fai[c][0])
        print(f"  loaded consensus chains for {len(used_founders)} founders "
              f"(segments cut at each founder's own coordinates)", flush=True)

    # Stitch mosaics
    for i, ind in enumerate(pop):
        ind_id = f"ind{i+1:03d}"
        hap_dir = haps_dir / f"s_{ind_id}"
        hap_dir.mkdir(exist_ok=True)
        out_fa = hap_dir / "h1.fa"
        if out_fa.exists() and (hap_dir / "h1.fa.fai").exists():
            continue
        write_mosaic_fasta(ind, founder_fastas, str(out_fa), samtools=args.samtools, chains=chains)
        if (i + 1) % 10 == 0:
            print(f"  built {i+1}/{args.n_indiv} mosaic FASTAs", flush=True)
    print(f"  done — {args.n_indiv} mosaics in {haps_dir}/", flush=True)


if __name__ == "__main__":
    main()
