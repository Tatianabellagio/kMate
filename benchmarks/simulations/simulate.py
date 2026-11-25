#!/usr/bin/env python
"""Simulate one pool: mosaic genomes, paired-end reads, per-record truth.

Usage: python simulate.py POOL_ID      (POOL_ID as listed by runs.py; writes config.POOLS/POOL_ID/)

1. Mosaics (sims/scripts/make_recomb_mosaics.py). n individuals are drawn from the founder set without replacement
   (at most one per founder), then reproduce for `gen` generations, each offspring a clone of
   one parent with probability `selfing` and otherwise a recombinant of two (4 cM/Mb). Each TAIR10 segment is cut from
   its founder at the founder's own coordinates (chains), so a whole founder keeps its whole sequence (SNP-only genomes,
   for hapFIRE, are in TAIR10 coordinates and need no chain).
   Founder-type pools fix the 50 founders instead: 40 long-read + 10 short-read ("long") or 5 + 45 ("short").
2. Reads (wgsim, seeded). Individual i has pool share w_i (1/n each, or `dominant` for individual 1 and the rest split
   equally) and gets round(depth * w_i * L_i / (2 * 150)) read pairs, L_i its length without N, so every base of every
   individual is read at depth * w_i.
3. Truth (sims/scripts/compute_recomb_truth.py): each record's allele frequency in the pool, sum_i w_i * allele of the
   founder individual i carries there, over the founders called at the record; one table per variant panel.

Writes ancestry.tsv, pool_weights.tsv, shares.tsv (ind_id, share, length, non_N, read_pairs, seed), reads/r{1,2}.fq,
truth_<panel>.tsv.gz. The mosaic FASTAs are deleted once read (the pool is reproducible from its id). Env: kmate.
"""
import os, sys, shutil, subprocess, zlib
import numpy as np, pandas as pd, pysam
import config as C
from runs import all_pools


def run(*cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def founder_type_pool(pool, out):
    """Source weights putting the 50 founders of a founder-type pool at 1/50 (and every other founder at 0)."""
    order = [str(f) for f in np.load(C.FOUNDER_SETS[pool.founders]["order"], allow_pickle=True)["founders"]]
    long_read = set(open(C.LONG_READ).read().split())
    n_long = {"long": 40, "short": 5}[pool.founder_type]
    rng = np.random.default_rng(pool.seed)
    chosen = set(rng.choice(sorted(f for f in order if f in long_read), n_long, replace=False)) | \
        set(rng.choice(sorted(f for f in order if f not in long_read), pool.n - n_long, replace=False))
    path = f"{out}/founder_type_source.tsv"
    pd.DataFrame({"founder": order, "prob": [1.0 / pool.n if f in chosen else 0.0 for f in order]}).to_csv(path, sep="\t", index=False)
    return path


def non_n_length(fa):
    seq = pysam.FastaFile(fa).fetch(C.CHROM).upper()
    return len(seq), len(seq) - seq.count("N")


def main(pool_id):
    pool = all_pools()[pool_id]
    fs = C.FOUNDER_SETS[pool.founders]
    out = f"{C.POOLS}/{pool_id}"
    shutil.rmtree(out, ignore_errors=True)                         # always from scratch (a requeued job may have left a part)
    os.makedirs(f"{out}/reads")

    # 1. mosaics
    extra = ["--selfing-rate", pool.selfing] if pool.selfing else []
    if pool.founder_type:
        extra += ["--source-weights", founder_type_pool(pool, out)]
    run(C.PYTHON, C.MOSAICS, "--n-indiv", pool.n, "--n-generations", pool.gen, "--seed", pool.seed,
        "--cactus-dir", fs["genomes"], *(["--chain-dir", fs["genomes"]] if fs["chains"] else []), "--founders-meta", fs["order"],
        "--out-dir", out, "--chroms", C.CHROM, *extra)

    # 2. reads
    inds = sorted(os.listdir(f"{out}/haps"))                       # s_ind001, s_ind002, ...
    share = np.full(len(inds), 1.0 / len(inds))
    if pool.dominant:
        share[0], share[1:] = pool.dominant, (1 - pool.dominant) / (len(inds) - 1)
    rows = []
    with open(f"{out}/reads/r1.fq", "w") as r1, open(f"{out}/reads/r2.fq", "w") as r2:
        for ind, w in zip(inds, share):
            fa = f"{out}/haps/{ind}/h1.fa"
            length, non_n = non_n_length(fa)
            pairs = round(pool.depth * w * non_n / (2 * C.READ_LEN))
            seed = zlib.crc32(f"{pool_id}/{ind}".encode()) & 0x7FFFFFFF
            t1, t2 = f"{out}/reads/_{ind}.1.fq", f"{out}/reads/_{ind}.2.fq"
            subprocess.run([C.WGSIM, "-N", str(pairs), "-1", str(C.READ_LEN), "-2", str(C.READ_LEN),
                            "-d", str(C.INSERT), "-s", str(C.INSERT_SD), "-e", str(C.ERROR), "-r", str(C.MUTATION),
                            "-R", str(C.INDEL_FRAC), "-X", str(C.INDEL_EXT), "-A", str(C.MAX_N), "-S", str(seed),
                            fa, t1, t2], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for t, o in ((t1, r1), (t2, r2)):
                with open(t) as fh:
                    shutil.copyfileobj(fh, o)
                os.remove(t)
            rows.append(dict(ind_id=ind.removeprefix("s_"), share=w, length=length, non_N=non_n, read_pairs=pairs, seed=seed))
    S = pd.DataFrame(rows)
    S.to_csv(f"{out}/shares.tsv", sep="\t", index=False)
    print(f"{len(S)} individuals, {S.read_pairs.sum():,} read pairs; lengths {S.length.min():,}-{S.length.max():,}", flush=True)
    shutil.rmtree(f"{out}/haps")

    # 3. truth
    for panel in fs["truth"]:
        pre = C.VARIANT_PANELS[panel]["prefix"]
        run(C.PYTHON, C.TRUTH, "--ancestry", f"{out}/ancestry.tsv", "--shares", f"{out}/shares.tsv",
            "--var-pa", f"{pre}.var_pa.npz", "--var-called", f"{pre}.var_called.npz", "--var-meta", f"{pre}.meta.npz",
            "--out", f"{out}/truth_{panel}.tsv.gz")
    open(f"{out}/DONE", "w").close()


if __name__ == "__main__":
    main(sys.argv[1])
