"""Mosaic genomes for one simulated pool: founder draw, generations of reproduction, mosaic FASTAs.

n individuals are drawn without replacement from the founder list (`--founders-meta`, its order fixes the draw with the
seed; `--source-weights` restricts / weights the draw), each starting as one whole founder. Each of `--n-generations`
generations replaces every individual by an offspring that is, with probability `--selfing-rate`, a clonal copy of one
random parent (inbred founders: selfing reproduces the parent) and otherwise a recombinant of two random parents
(Poisson crossovers at `--recomb-rate`, uniform positions). Segments are in TAIR10 coordinates; each is cut from its
founder at the founder's own coordinates through the chain written with the founder genome (`--chain-dir`), so a whole
founder keeps its whole sequence and mosaic lengths differ. Without `--chain-dir` the genomes must be in TAIR10
coordinates (SNP-only consensus; lengths checked).

Writes to --out-dir: pool_weights.tsv (founder, count, weight), ancestry.tsv (ind_id, chrom, start, end, founder),
haps/s_ind###/h1.fa (+ .fai). Called by benchmarks/simulations/simulate.py. Env: kmate.
"""
import argparse, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd

RECOMB_RATE = 4e-8      # crossovers per bp per meiosis (~4 cM/Mb, A. thaliana)
CHROM_LENGTHS = {"Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830, "Chr4": 18585056, "Chr5": 26975502}
SAMTOOLS = os.path.join(os.path.dirname(sys.executable), "samtools")
SAMTOOLS = SAMTOOLS if os.path.exists(SAMTOOLS) else (shutil.which("samtools") or "samtools")


class IndexedFasta:
    """Fetch 1-based inclusive ranges from a .fai-indexed FASTA."""
    def __init__(self, path):
        self.path = path
        self.fai = {}
        for line in open(f"{path}.fai"):
            name, length, offset, lblen, llen = line.split("\t")[:5]
            self.fai[name] = (int(length), int(offset), int(lblen), int(llen))

    def fetch(self, name, start, end):
        length, offset, lblen, llen = self.fai[name]
        s = start - 1
        assert 0 <= s < end <= length, (self.path, name, start, end, length)
        a, b = offset + (s // lblen) * llen + s % lblen, offset + (end // lblen) * llen + end % lblen
        with open(self.path, "rb") as f:
            f.seek(a)
            return f.read(b - a).decode("ascii").replace("\n", "").replace("\r", "").upper()


class ConsensusChain:
    """TAIR10 -> founder-consensus coordinates for one chromosome, from the chain `bcftools consensus -c` writes beside
    the founder FASTA (t = TAIR10, q = consensus; ungapped blocks, 0-based)."""
    def __init__(self, chain_path, chrom, t_len, q_len):
        ts, qs, sz, found, on = [], [], [], False, False
        for line in open(chain_path):
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
            if on:
                ts.append(t); qs.append(q); sz.append(int(p[0]))
                if len(p) == 3:
                    t += int(p[0]) + int(p[1]); q += int(p[0]) + int(p[2])
        if not found:
            raise ValueError(f"{chain_path}: no chain for {chrom}")
        self.ts, self.qs, self.sz, self.t_len, self.q_len = np.array(ts), np.array(qs), np.array(sz), t_len, q_len

    def cut(self, t):
        """0-based consensus index where a piece starting at TAIR10 base t (1-based; t_len + 1 = past the end) begins.
        Consecutive pieces of one founder tile its consensus exactly: bases inserted between TAIR10 t-1 and t, and a
        TAIR10 base deleted in the founder, both resolve to the start of the next aligned block, so inserted bases stay
        with the piece on the left."""
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
        return fasta.fetch(chrom, a + 1, b) if b > a else ""


def sample_crossovers(chrom_len, rate, rng):
    n = rng.poisson(chrom_len * rate)
    return np.sort(rng.integers(1, chrom_len, size=n)) if n else np.array([], dtype=np.int64)


def recombine(parent_a, parent_b, rate, rng):
    """One haploid offspring of two segment-track parents: parents alternate between crossovers."""
    out = {}
    for chrom, L in CHROM_LENGTHS.items():
        bounds = [0] + list(sample_crossovers(L, rate, rng)) + [L]
        which = rng.integers(0, 2)
        segs = []
        for s0, e in zip(bounds[:-1], bounds[1:]):
            s = s0 + 1
            for ps, pe, f in (parent_a if which == 0 else parent_b)[chrom]:
                if max(s, ps) <= min(e, pe):
                    segs.append((max(s, ps), min(e, pe), f))
            which = 1 - which
        merged = []
        for s, e, f in segs:                                  # merge adjacent pieces of the same founder
            if merged and merged[-1][2] == f and merged[-1][1] + 1 == s:
                merged[-1] = (merged[-1][0], e, f)
            else:
                merged.append((s, e, f))
        out[chrom] = merged
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n-indiv", type=int, required=True)
    ap.add_argument("--n-generations", type=int, default=0)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--selfing-rate", type=float, default=0.0)
    ap.add_argument("--recomb-rate", type=float, default=RECOMB_RATE)
    ap.add_argument("--cactus-dir", required=True, help="founder genomes <founder>.chr.fa (+ .fai)")
    ap.add_argument("--chain-dir", default=None, help="<founder>.chain; omit only for genomes in TAIR10 coordinates")
    ap.add_argument("--founders-meta", required=True, help=".npz with 'founders': the founder list, in draw order")
    ap.add_argument("--source-weights", default=None, help="TSV (founder, prob): draw probabilities (default uniform)")
    ap.add_argument("--chroms", default=" ".join(CHROM_LENGTHS), help="space-separated chromosomes")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    for c in list(CHROM_LENGTHS):
        if c not in a.chroms.split():
            del CHROM_LENGTHS[c]
    assert 0.0 <= a.selfing_rate <= 1.0
    out = Path(a.out_dir)
    (out / "haps").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    # founder draw: n of the F founders without replacement (the same random calls as every earlier simulation)
    founders = [str(f) for f in np.load(a.founders_meta, allow_pickle=True)["founders"]]
    F = len(founders)
    p = np.full(F, 1.0 / F)
    if a.source_weights:
        sw = pd.read_csv(a.source_weights, sep="\t", dtype={"founder": str})
        p = sw.set_index("founder").prob.reindex(founders).fillna(0.0).to_numpy(np.float64)
        p = p / p.sum()
    if a.n_indiv > F:
        raise ValueError(f"{a.n_indiv} individuals from {F} founders: a pool holds each founder at most once")
    counts = np.zeros(F, dtype=int)
    counts[rng.choice(F, size=a.n_indiv, replace=False, p=p)] = 1
    pd.DataFrame({"founder": founders, "count": counts, "weight": counts / counts.sum()}).to_csv(
        out / "pool_weights.tsv", sep="\t", index=False)
    drawn = [f for f, c in zip(founders, counts) if c]
    order = rng.choice(drawn, size=a.n_indiv, replace=False, p=np.full(len(drawn), 1.0 / len(drawn)))   # individual order
    pop = [{c: [(1, L, str(f))] for c, L in CHROM_LENGTHS.items()} for f in order]

    for _ in range(a.n_generations):
        new = []
        for _ in range(a.n_indiv):
            if a.selfing_rate > 0.0 and rng.random() < a.selfing_rate:
                new.append({c: list(s) for c, s in pop[rng.integers(0, len(pop))].items()})        # self = clone
            else:
                i, j = rng.integers(0, len(pop)), rng.integers(0, len(pop))
                new.append(recombine(pop[i], pop[j], a.recomb_rate, rng))                          # outcross
        pop = new

    pd.DataFrame([(f"ind{i + 1:03d}", c, s, e, f) for i, ind in enumerate(pop) for c in sorted(ind) for s, e, f in ind[c]],
                 columns=["ind_id", "chrom", "start", "end", "founder"]).to_csv(out / "ancestry.tsv", sep="\t", index=False)

    used = sorted({f for ind in pop for c in ind for _, _, f in ind[c]})
    fasta = {f: IndexedFasta(f"{a.cactus_dir}/{f}.chr.fa") for f in used}
    if a.chain_dir:
        chains = {(f, c): ConsensusChain(f"{a.chain_dir}/{f}.chain", c, L, fasta[f].fai[c][0])
                  for f in used for c, L in CHROM_LENGTHS.items()}
    else:
        bad = [(f, c) for f in used for c, L in CHROM_LENGTHS.items() if fasta[f].fai[c][0] != L]
        assert not bad, f"genomes not in TAIR10 coordinates, pass --chain-dir: {bad[:3]}"
        chains = None
    for i, ind in enumerate(pop):
        fa = out / "haps" / f"s_ind{i + 1:03d}" / "h1.fa"
        fa.parent.mkdir(exist_ok=True)
        with open(fa, "w") as fh:
            for c in sorted(ind):
                seq = "".join(chains[(f, c)].piece(fasta[f], c, s, e) if chains else fasta[f].fetch(c, s, e)
                              for s, e, f in ind[c])
                fh.write(f">{c}\n" + "".join(seq[k:k + 80] + "\n" for k in range(0, len(seq), 80)))
        subprocess.run([SAMTOOLS, "faidx", str(fa)], check=True)
    print(f"{a.n_indiv} individuals ({len(used)} founders), {a.n_generations} generations, selfing {a.selfing_rate}: "
          f"{out}", flush=True)


if __name__ == "__main__":
    main()
