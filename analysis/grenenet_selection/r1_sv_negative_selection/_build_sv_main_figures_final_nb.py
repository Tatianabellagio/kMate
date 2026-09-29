#!/usr/bin/env python
"""Build+execute `sv_main_figures_final.ipynb`: the two headline main-text figures --
Figure 1 (purifying selection on SVs in the founders) and Figure 2's lead panel (the same signal
acting in real time in GrENE-Net). Everything else -- Figure 2's full panel a-d, the walk-through,
the per-context/origin supplements, Figure 3 (which SVs) -- stays in `sv_main_fig{1,2,3}.ipynb` as
supplementary material. Reuses the exact cells from `_build_main_fig{1,2}_nb.py` (import only, no
duplication); those two scripts guard their own `build()` call behind `if __name__ == "__main__"` so
importing them here does not also rebuild `sv_main_fig{1,2}.ipynb`. Shared setup/data code in
`_main_fig_common.py`. `basic` env."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_fig_common import code_setup, code_grenenet, build, new_markdown_cell, new_code_cell
from _build_main_fig1_nb import md_fig1, code_fig1
from _build_main_fig2_nb import md_lead, code_lead, md_sfs, code_sfs, md_sfs_garden, code_sfs_garden

md_title = """# Selection on structural variants — main figures

The two headline results:

1. **Figure 1** — SVs, insertions particularly, carry the classic site-frequency signature of
   purifying selection in the 231-founder panel: shifted towards rare relative to synonymous SNPs,
   more negative Tajima's *D*, and depleted from coding sequence relative to genomic expectation.
2. **Figure 2** — the same signal acting in real time: in GrENE-Net, insertions (and, in CDS,
   deletions moving the other way) shift in frequency relative to matched same-context SNPs across
   the 31 gardens.

**Unit and filters** (as in `sv_main_fig1.ipynb`/`sv_main_fig2.ipynb`): the SV **event** (arch3
alleles of a pangenome bubble merged with `truvari collapse` defaults); called in ≥ 90% of the 231
founders throughout; ≥ 2 carrier founders for the GrENE-Net panel, singletons included in the
spectra. Context, first match wins: CDS > UTR > intron > reference TE > promoter (≤ 2 kb upstream)
> intergenic.

The last section re-expresses result 2 as a site-frequency spectrum (seed mix → final generation), which
is where the drift baseline and the limits of the CDS-specific claim are easiest to see.

**Otherwise this notebook is just the two main-text figures.** Full per-context panels, robustness checks, the
origin/climate walk-through, and every supplement are in `sv_main_fig1.ipynb` and
`sv_main_fig2.ipynb`; Figure 3 (which SVs are being purified) is in `sv_main_fig3.ipynb`."""

if __name__ == "__main__":
    build("sv_main_figures_final.ipynb", [
        new_markdown_cell(md_title), new_code_cell(code_setup),
        new_markdown_cell(md_fig1), new_code_cell(code_fig1),
        new_code_cell(code_grenenet),
        new_markdown_cell(md_lead), new_code_cell(code_lead),
        new_markdown_cell(md_sfs), new_code_cell(code_sfs),
        new_markdown_cell(md_sfs_garden), new_code_cell(code_sfs_garden),
    ])
