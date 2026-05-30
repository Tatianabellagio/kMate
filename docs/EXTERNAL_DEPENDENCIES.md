# External (vendored) dependencies

The `external/` directory holds third-party tools the pipeline calls. It is
**gitignored** (not committed) — the tools are large upstream repositories with
their own history and licensing, so we do not vendor their source into this repo.
This file is the tracked record of *what* `external/` must contain and *where to
get it*, so a fresh clone can reconstruct it and reproduce the pipeline.

## What production actually uses

Only **two** external scripts are on the production path (the arch3 panel build,
`panel/arch3/`):

| Tool (dir) | Upstream | Pinned commit | Used by |
|---|---|---|---|
| `external/genotyping-pipelines` | https://github.com/eblerjana/genotyping-pipelines.git | `ea1e52b` (2026-05-07) | `panel/arch3/chr1/jobA1_annotate_chr1.sh` → `prepare-vcf-MC/workflow/scripts/annotate_vcf.py` |
| `external/pangenie-tools` | https://github.com/eblerjana/pangenie.git | `ebd227a` (2026-05-07) | `panel/arch3/chr1/jobA2_pg_chr1.sh` + `jobA3_cactus_chr1.sh` → `pipelines/run-from-callset/scripts/convert-to-biallelic.py` |
| `external/HapFIRE` | git@github.com:xingwu2/HapFIRE.git | `a573877` (2026-03-18) | **comparator only** — methods-comparison runs; *not* on the production path |

Within `genotyping-pipelines` only the `prepare-vcf-MC` sub-pipeline is used; the
other sub-pipelines (benchmarking, cohort-genotyping, …) come with the upstream
clone and are unused here. `pangenie-tools` is the `eblerjana/pangenie` repo
(vendored under the `pangenie-tools/` name).

## Reproducing `external/`

```bash
mkdir -p external && cd external
git clone https://github.com/eblerjana/genotyping-pipelines.git
git -C genotyping-pipelines checkout ea1e52b
git clone https://github.com/eblerjana/pangenie.git pangenie-tools
git -C pangenie-tools checkout ebd227a
# comparator only (skip unless reproducing the hapFIRE comparison column):
git clone git@github.com:xingwu2/HapFIRE.git
git -C HapFIRE checkout a573877
```

## Related runtime dependencies (not in `external/`)

- **PanGenie** binary + **Jellyfish 2** + **bcftools** — panel genotyping / k-mer
  counting (conda env `hapfm`; see `HANDOFF.md` Reproducibility).
- **VISOR SHORtS** — pool-seq read simulation (conda env `pang`; see
  `docs/SIMULATIONS_METHODS.md` §10).
- **HARP** binary — required on `PATH` for HapFIRE comparator runs only
  (`BACKGROUND.md` References).
