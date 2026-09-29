# tests/ — correctness checks for the kMate pipeline

Two different things live here, and they are not run the same way:

1. **Unit-ish correctness tests** that run on a workstation or a compute node in seconds
   to minutes.
2. **SLURM validation harnesses** (`run_*.sh`) that build real matrices and run genome-wide
   EM. These are jobs, not tests — submit them with `sbatch`, don't run them inline.

Written 2026-08-24 during the root-level cleanup.

> ⚠ Per `CLAUDE.md`: never run any of this on a login node (`ln00X`). Get onto a compute
> node with `salloc` or submit with `sbatch` first.

## Tests

| File | What it checks |
|---|---|
| `test_kmer_count.py` | Validates the k-mer counter against a known reference and a known BAM. |
| `test_genomewide_validation.py` | Concatenates per-chromosome `kmer_pa` matrices, counts k-mers in each pool against the full set, runs EM, and compares the result to truth. |

## SLURM harnesses

| File | What it runs |
|---|---|
| `run_build_kmer_pa.sh` | Builds the `K_pa` (k-mer presence/absence) matrices. |
| `run_build_var_pa.sh` | Builds the `V_pa` (variant presence/absence) matrices. |
| `run_genomewide_validation.sh` | Drives `test_genomewide_validation.py` genome-wide. |

## Fixture generation

`make_selftest_fixture.py` builds the tiny deterministic end-to-end fixture under
`src/kmate/data/selftest/` that backs `kmate selftest`. **Run once; its outputs are
committed** so they ship in the wheel/conda package — `.gitignore` carries explicit
negations (`!src/kmate/data/selftest/`) to override the global bulk-data rules for that one
directory. Don't regenerate it casually.

## Subdirectories

- `archive/` — superseded one-off test and diagnostic scripts (`compare_h.py`,
  `af_scatter*.py`, the g0/n50 subset builders, …), kept per `archive/README.md`.
- `logs/` — SLURM output; git-ignored (`tests/logs/`), prunable.
