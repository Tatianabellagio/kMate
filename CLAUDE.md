<!-- ╔══════════════════════════════════════════════════════════════════════╗ -->
<!-- ║                                                                        ║ -->
<!-- ║   ⛔  NEVER RUN COMPUTE ON A SAVIO LOGIN NODE  (ln001 / ln002 / ...)    ║ -->
<!-- ║                                                                        ║ -->
<!-- ╚══════════════════════════════════════════════════════════════════════╝ -->

# ⛔ FIRST: check what node you are on

**Before running ANY command, run `hostname`.**

- Login nodes are named `ln00X` (e.g. `ln002`). **Do NOT run compute there** —
  python, mamba/conda, samtools/bcftools, jellyfish, snakemake, R, pipelines,
  notebook execution, or any real work. It chokes the shared login node, it is
  against cluster policy, and it makes sessions "chug" silently.
- Compute nodes look like `n0221.savio3`, `n0000.savio2`, etc. Those are fine.

**If you find yourself on a login node, STOP.** Get onto a compute node first:

```bash
salloc --partition=savio3 --nodes=1 --time=2:00:00 --pty bash    # interactive
# or submit batch work with sbatch / srun
```

Only light, instant operations (`git`, `ls`, `cat`, editing files, `squeue`,
`sbatch`, `salloc`, `srun`) are acceptable on a login node — and only to get
onto a compute node.

> A user-level PreToolUse hook (`~/.claude/hooks/block_login_node.sh`) hard-blocks
> compute commands when the hostname is a login node. If you hit that block, it is
> working as intended — hop to a compute node, don't try to route around it.

---

## Project notes

- **Install tools via mamba/conda, never pip** — envs are conda-managed.
- **No co-author trailers in commits** (no `Co-Authored-By` / "Generated with
  Claude") — project policy.
- Production pipeline is **arch3** (see the kMate docs / memory for details).
