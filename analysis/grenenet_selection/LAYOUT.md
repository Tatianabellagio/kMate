# Where things go

Rules for adding work to this tree, so it stops sprawling. This file is about
**file placement only** — it makes no claims about results. For what the
analyses found, see the per-section docs and `README.md`.

Set 2026-08-24, when the tree was reorganised: 198 loose files and 40
undocumented directories at this top level became the eight sections below.

---

## The shape

```
analysis/
  grenenet_selection/          <- this tree: GrENE-Net selection analyses
      lib.py                   <- shared loaders. STAYS at this level.
      README.md  LAYOUT.md  GLOBAL_MODE_DECISION.md  EXPORT_MANIFEST.md
      logs/                    <- SLURM job logs (sbatch writes here)
      notebooks/               <- rendered .ipynb for the whole tree
      archive/                 <- retired work, never deleted
      <section>/
          *.py *.sh *.sbatch   <- the scripts, loose in the section
          results/             <- everything those scripts write
              plots/           <- figures. never loose beside the data
          <subproject>/        <- same scripts+results/ shape, one level down
  panel_qc/                    <- kMate panel/method QC (NOT GrENE-Net biology)
```

| section | what belongs in it | scripts | subprojects |
|---|---|---|---|
| `r1_sv_negative_selection/` | Result 1 — SVs vs SNPs matched on initial frequency; climate correlation | 32 | — |
| `r2_gea_nonsnp/` | Result 2 — GEA hits visible only in non-SNP data | 7 | `phase1_replication/`, `cam5_replication/` |
| `r3_persite_gwas/` | Result 3 — per-site GWAS, ecotype-selection coefficient as trait | 34 | — |
| `extras/` | side investigations that are **not** one of the three results | 23 | `driver_passenger/` |
| `common/` | inputs shared by more than one section (AF store, per-gen and pool matrices, p0, founder h) | 10 | `rerun_kfw_hb/` |
| `blocks/` | LD-block / analysis-unit definition (what a test unit *is*) | 39 | `hap_blocks/`, `bigld_env/` |
| `wza/` | the WZA block-aggregation method: shared `wza_script.py` + the investigation that settled the regime | 1 | `investigation/` |
| `qc/` | QC of *this* analysis (coverage, panel overlap, seed-mix identifiability) | 10 | `seedmix_validation/` |

`blocks/` and `wza/` are **method** sections, not results: `blocks/` defines the
unit, `wza/` aggregates per-variant p-values over it. `wza_script.py` lives there
rather than in a results section because it has ~19 referrers spanning
`phase1_replication/`, `wza/investigation/` and `extras/driver_passenger/`.

`r2_gea_nonsnp/` is small at the top level because its production pipeline is the
`phase1_replication/` subproject; the gen-3 SV pilot that used to sit beside it
was retired 2026-08-24 (gen 9 is the current generation).

`panel_qc/` is a sibling of `grenenet_selection/`, not part of it: it holds
kMate method/panel validation (k-mer index comparison, panel stats, seed-mix
duplicate rate). It used to be a second top-level `results/` directory at the
repo root, which is what made "analysis vs results" confusing.

---

## The rules

**1. Do not create a new directory at this top level.** New work goes inside
the section that owns the question. If it genuinely fits none of them, it is
almost certainly an `extras/` item.

**2. Output goes under `<section>/results/`, never beside the scripts and
never at the top level.** This is the rule that had been half-applied since
2026-07-10 and is the main reason the tree sprawled: `gea/`, `lfmm/`,
`varexp/`, `site_temporal/` and eleven others had accumulated as top-level
siblings of the code that wrote them, so `blocks/` (code) sat next to
`blocks_mcf90/` (its output) and read like a duplicate.

**3. A new output directory needs a `.gitignore` rule in the same commit.**
Rules here are literal paths, so they silently stop matching when a parent is
renamed. Always confirm with:

```bash
git check-ignore -v analysis/grenenet_selection/<section>/results/<newdir>
```

**4. Figures go in a `plots/` subdirectory. Never loose in a results dir.**
A results directory should hold data (`.npz`, `.csv`, `.tsv`) and one `plots/`
subdirectory — not 30 PNGs interleaved with the CSVs. Write figures as:

```python
plt.savefig(f"{OUT}/plots/my_figure.png", dpi=150, bbox_inches="tight")
```

If `OUT` is a new directory, create `OUT/plots` alongside it. 416 existing
figures were swept into 19 `plots/` dirs on 2026-08-24 and 83 `savefig` paths
repointed; **25 calls in 16 files still write loose** because they pass a
variable built at runtime rather than an inline path, so they need a look by
hand:

`blocks/plot_sv_landscape.py`, `blocks/plot_unit_distributions.py`,
`qc/plot_coverage_distribution.py`,
`r1_sv_negative_selection/{_build_sv_polarity_manhattan_nb,_render_site4_enrichment_fig}.py`,
`r2_gea_nonsnp/cam5_replication/{build_nb,build_nb_byclass}.py`,
`r2_gea_nonsnp/phase1_replication/{_build_manhattan_3x3_nb,_build_manhattan_clq90_nb,plot_clq90_manhattan}.py`,
`r3_persite_gwas/plot_class_gwas_pngs.py`,
`wza/investigation/_build_manhattan_nb.py`,
`notebooks/_recreate_density_snp_only.py`, and the three `genes_expl/plot_*.py`
(left untouched as uncommitted work).

**5. Scripts reach `lib.py` by counting directory levels.** The idiom is

```python
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
```

with **one `os.path.dirname()` per level between the script and `lib.py`** — two
for a script in `<section>/`, three in `<section>/<subproject>/`. Moving a
script to a different depth without fixing this breaks `import lib` at runtime,
not at move time, so it fails silently later.

Two related traps:
- Many scripts carry a **second** `sys.path.insert` for their *own* directory,
  so they can import sibling modules. Only the chain that resolves to `lib.py`
  should ever be re-levelled; deepening the other one breaks sibling imports.
- Scripts that hardcode `/global/scratch/users/tbellg/kmate/...` still work —
  that path is a symlink to the real repo root — but prefer the `__file__`-
  relative form.

**6. Retire, don't delete.** Superseded work moves to `archive/` with a note
saying what replaced it. Before moving anything, grep for inbound references:
several directories here look stale but are load-bearing (see below).

**7. Do not rewrite paths inside `archive/`.** Archived scripts were written
against the layout of their time; repointing them at today's paths would
falsify the record.

---

## Naming

- Sections are `r1_`/`r2_`/`r3_` + a short description of the *question*, so the
  three results sections sort together and are self-describing.
- Avoid names that describe a *filter or version* as a new top-level thing
  (`_v2`, `_maconly`, `_nofilter`). Those are variants of one analysis and
  belong together — the four `sv_snp_ld*` directories are why this rule exists.
- Do not name a directory after the method if the tree also holds other
  methods. This tree was called `grenenet_gea`, but only `r2` is GEA; `r1` is
  selection and `r3` is GWAS. It was renamed `grenenet_selection` 2026-08-24.

---

## Load-bearing despite appearances

Checked during the cleanup — do not retire these on a name-based hunch:

| looks stale | actually |
|---|---|
| `blocks/results/blocks_recompute/` | superseded as the production block map, but `blocks/_build_blockcoherence_nb.py` still globs its `panel_support_tag` / `missing_sensitivity` CSVs for two notebook panels |
| `r3_persite_gwas/results/sv_snp_ld/` | not an LD result at all — holds the GrENE-Net VCF and extracted short-read genotypes, i.e. the *input* to the `sv_snp_ld_v2*` tables |
| `r1_sv_negative_selection/_temporal_s_*.py`, `_compute_s_*.py` | audited 2026-07-15 and found to test genuinely different statistics. The docs around them contain *retractions* of a superseded label, not declarations of one |
| `r2_gea_nonsnp/phase1_replication/annotate_gene_function.py` | superseded as an annotator, but still imported for its `mygene_batch` helper |
| `qc/seedmix_validation/fix_seedmix/`, `fix_norm/` | labelled abandoned arms, but later arms still read their cached counts |

A compute script having "no in-repo reference" does **not** mean it is dead:
notebooks in this tree load results by *output filename*, not by importing the
producer. Check for the data file before concluding anything is orphaned.

---

## Known loose ends

- `notebooks/` is a single flat directory for the whole tree, while every
  builder script now lives in its section — so a notebook and its `_build_*.py`
  are no longer adjacent. Kept flat deliberately; revisit if it gets painful.
- A few uncommitted items remain at this top level (`genes_expl/` and several
  `_build_*_nb.py`), held back as a separate workstream.
- 25 `savefig` calls in 16 files still write figures loose (rule 4 above).
- Sections should not read each other's `results/`. The gene-layer scripts that
  did (they read *and wrote* r3's `varexp`) moved r2 -> r3 on 2026-08-24; if a
  new cross-section dependency appears, promote the shared input to `common/`
  rather than reaching across.
