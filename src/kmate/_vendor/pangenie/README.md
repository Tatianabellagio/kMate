# Vendored: two scripts from PanGenie (Jana Ebler, MIT)

Unmodified copies, so `kmate decompose` works from a conda/PyPI install with no
second checkout. `kmate decompose` runs them as scripts; nothing imports them.

| file | upstream path | role |
|---|---|---|
| `annotate_vcf.py` | `pipelines/prepare-vcf-from-MC/workflow/scripts/annotate_vcf.py` | decompose each Minigraph-Cactus bubble into its nested variants; write the ID-annotated multi-allelic VCF and the biallelic catalog |
| `convert-to-biallelic.py` | `pipelines/run-from-callset/scripts/convert-to-biallelic.py` | move each sample's genotype from the multi-allelic records onto the catalog by variant ID |

- Source: https://github.com/eblerjana/pangenie, commit `ebd227a5d3ab24955b6578a60501b85930476e5a`
- License: MIT, `LICENSE.md` in this directory (Copyright (c) 2020 Jana Ebler)
- These are the same files the GrENE-Net arch3 panel was built with
  (`panel/arch3/README.md`).

Cite Ebler et al. (2022) *Nature Genetics* 54:518-525 and Liao et al. (2023)
*Nature* 617:312-324 if the decomposition matters to your results
(`kmate decompose --citation`).

To update: copy the two files from a newer upstream commit, record the commit
here, and rerun the decompose check in `tests/README.md`.
