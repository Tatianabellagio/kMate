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
- **Plotting convention: no chart titles, no subplot titles** (no
  `ax.set_title()` / `fig.suptitle()`) in any matplotlib/seaborn figure —
  notebook markdown headers + axis labels carry context instead. For small
  multiples where each panel needs identity, use an in-panel corner annotation
  (`ax.annotate`/`ax.text` in axes-fraction coords), not a title.
- **Gene annotation: use `phase1_replication/annotate_genes_tair_uniprot.py`**
  (TAIR GO + UniProt), *not* the older mygene.info-only
  `annotate_gene_function.py`. NCBI has no free-text summaries for Arabidopsis
  loci, so mygene alone left the `summary` column empty (0/73 genes on the
  non-SNP WZA hit list) and the functional categories were name-driven. The
  TAIR/UniProt version gives protein name + keywords 73/73, curated FUNCTION
  text 39/73, and lifts tagged categories 25→31/73. Gene *coordinates* stay
  local TAIR10 GFF (`lib.load_genes()`) — do not switch to Ensembl for those:
  Ensembl Plants serves the same TAIR10/Araport11 models and adds nothing.
  Both APIs are public/no-auth but need outbound HTTPS (works from savio4
  compute nodes) and an explicit `User-Agent` — `current.geneontology.org`
  returns 403 to urllib's default. TAIR's own `arabidopsis.org/download_files`
  URLs are 403/login-walled; the GO Consortium GAF mirror is the usable route.
