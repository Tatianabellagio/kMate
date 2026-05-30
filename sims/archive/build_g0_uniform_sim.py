"""Build a g0 (no-recomb) uniform pool sim from a sampled founder set.

Sample n founders from the panel of 231 with a chosen cactus/PG composition,
wgsim per-founder reads at 1/n contribution, concatenate to pool r1.fq, r2.fq.

Saves:
    h_truth.tsv     — per-founder truth h (1/n if sampled, 0 otherwise)
    founders.txt    — sampled founder IDs
    reads/r1.fq, r2.fq
    config.json
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys
from pathlib import Path
import numpy as np

WGSIM = '/global/home/users/tbellg/miniforge3/envs/pang/bin/wgsim'
SAMTOOLS = '/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/samtools'
ROOT = Path(__file__).resolve().parents[2]
FOUNDER_DIR = ROOT / 'sims/visor_freqk/founder_fastas_231_v3'


def sample_founders(panel: list[str], side_map: dict[str,str],
                    n: int, n_cactus_target: int | None, seed: int) -> list[str]:
    """Sample n founders, optionally biased toward cactus_count = n_cactus_target."""
    rng = np.random.default_rng(seed)
    cactus = [f for f in panel if side_map[f] == 'cactus']
    pg     = [f for f in panel if side_map[f] == 'PG']
    if n_cactus_target is None:
        return sorted(rng.choice(panel, size=n, replace=False).tolist())
    n_pg = n - n_cactus_target
    if n_cactus_target > len(cactus) or n_pg > len(pg):
        raise ValueError(f'cannot sample {n_cactus_target} cactus + {n_pg} PG')
    chosen_c = rng.choice(cactus, size=n_cactus_target, replace=False).tolist()
    chosen_p = rng.choice(pg, size=n_pg, replace=False).tolist()
    return sorted(chosen_c + chosen_p)


def build_sim(out_dir: Path, panel: list[str], chosen: list[str],
              n_pairs_per_founder: int, seed: int, keep_intermediate: bool=False):
    out_dir.mkdir(parents=True, exist_ok=True)
    reads_dir = out_dir / 'reads'
    per_dir   = out_dir / 'per_founder'
    reads_dir.mkdir(exist_ok=True)
    per_dir.mkdir(exist_ok=True)

    r1_chunks, r2_chunks = [], []
    for i, fid in enumerate(chosen):
        chr1_fa = per_dir / f'{fid}.chr1.fa'
        if not chr1_fa.exists() or chr1_fa.stat().st_size == 0:
            subprocess.run([SAMTOOLS, 'faidx', str(FOUNDER_DIR / f'{fid}.chr.fa'), 'Chr1'],
                            stdout=open(chr1_fa, 'w'), check=True)
            subprocess.run([SAMTOOLS, 'faidx', str(chr1_fa)], check=True)
        r1 = per_dir / f'{fid}_1.fq'
        r2 = per_dir / f'{fid}_2.fq'
        # Use fid+seed for per-founder wgsim seed (stable per-founder, varies by replicate)
        per_seed = (int(fid) * 1000 + seed) % 2147483647
        subprocess.run([WGSIM, '-N', str(n_pairs_per_founder),
                        '-1', '150', '-2', '150', '-d', '500', '-s', '50',
                        '-e', '0.005', '-r', '0', '-R', '0', '-X', '0',
                        '-S', str(per_seed),
                        str(chr1_fa), str(r1), str(r2)],
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        r1_chunks.append(r1); r2_chunks.append(r2)
        if (i+1) % 20 == 0:
            print(f'  [{i+1}/{len(chosen)}] wgsim done', flush=True)

    print(f'  concatenating into pool reads...', flush=True)
    with open(reads_dir/'r1.fq','wb') as f:
        for r in r1_chunks: f.write(r.read_bytes())
    with open(reads_dir/'r2.fq','wb') as f:
        for r in r2_chunks: f.write(r.read_bytes())

    if not keep_intermediate:
        import shutil
        shutil.rmtree(per_dir, ignore_errors=True)

    # h_truth: 1/n for sampled founders, 0 for others
    h = 1.0/len(chosen)
    with open(out_dir/'h_truth.tsv','w') as f:
        f.write('founder\th_truth\n')
        for fid in panel:
            f.write(f'{fid}\t{h if fid in set(chosen) else 0.0:.6f}\n')
    with open(out_dir/'founders.txt','w') as f:
        for fid in chosen: f.write(fid+'\n')
    return reads_dir/'r1.fq', reads_dir/'r2.fq'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True, help='pool size (sampled founders)')
    ap.add_argument('--n-cactus', type=int, default=None, help='# cactus to include (default: random)')
    ap.add_argument('--rep', type=int, required=True, help='replicate index for naming')
    ap.add_argument('--seed', type=int, default=None, help='RNG seed; default = rep')
    ap.add_argument('--cov', type=int, default=10, help='total pool coverage (default 10×)')
    ap.add_argument('--out-base', type=str, required=True,
                    help='base dir; sim writes to {base}/g0_n{N}_rep{REP}_{tag}')
    args = ap.parse_args()
    seed = args.rep if args.seed is None else args.seed

    # Panel from existing founder FASTAs (231 founders)
    panel = sorted([p.stem.replace('.chr','') for p in FOUNDER_DIR.glob('*.chr.fa')])
    assert len(panel) == 231, f'expected 231 founders, got {len(panel)}'

    # Side map from data/founder_split_cactus_pg.json
    with open(ROOT/'data/founder_split_cactus_pg.json') as f:
        split = json.load(f)
    cactus_set = set(map(str, split['cactus']))
    side_map = {f: ('cactus' if f in cactus_set else 'PG') for f in panel}

    chosen = sample_founders(panel, side_map, args.n, args.n_cactus, seed)
    n_c = sum(1 for f in chosen if side_map[f]=='cactus')
    n_p = sum(1 for f in chosen if side_map[f]=='PG')
    print(f'sampled n={len(chosen)}: cactus={n_c}, PG={n_p}')

    # Compose tag
    if args.n_cactus is None:
        tag = 'rand'
    elif args.n_cactus >= 0.7*args.n: tag = 'cact'
    elif args.n_cactus <= 0.2*args.n: tag = 'pg'
    else: tag = 'bal'
    out_dir = Path(args.out_base) / f'g0_n{args.n}_rep{args.rep}_{tag}'
    print(f'out_dir = {out_dir}')

    # Per-founder pairs: cov × 30Mb / (2×150) / n
    GENOME_CHR1 = 30_003_408
    total_pairs = args.cov * GENOME_CHR1 // 300
    n_pairs = max(1, total_pairs // len(chosen))
    print(f'  per-founder N pairs = {n_pairs}  (cov {args.cov}×, total ≈ {total_pairs:,} pairs)')

    build_sim(out_dir, panel, chosen, n_pairs, seed)

    # Save config
    with open(out_dir/'config.json','w') as f:
        json.dump({'n':args.n,'rep':args.rep,'seed':seed,'cov':args.cov,
                   'n_cactus':n_c,'n_pg':n_p,
                   'n_pairs_per_founder':n_pairs,
                   'tag':tag}, f, indent=2)
    print(f'DONE  {out_dir}')


if __name__ == '__main__':
    main()
