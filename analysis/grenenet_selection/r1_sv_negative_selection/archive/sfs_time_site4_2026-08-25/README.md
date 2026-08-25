# sfs_time_site4 (archived 2026-08-25)

AF spectrum across generations 0→3 at site 4 (the hot pilot site), by variant class.
Archived under the paper-streamlining pass: descriptive, single-site, and by its own
conclusion not an independent signal.

> All three classes start from a similar seedmix shape (median alt-AF ~0.02; SVs already
> slightly rarer, median 0.011 vs 0.021). Across gen 1→3 the bulk SNP/indel spectra barely
> move, while the **fraction lost creeps up fastest for SVs** — by gen 3, ~0.23% of SV
> records vs ~0.07–0.09% for SNP/indel.

The notebook's own reading is the reason it does not earn a slot: *"the raw SFS shift shown
here is consistent with [the de-trended null] rather than a new independent signal"* —
it is confounded with SVs' lower starting frequency, which the frequency-matched analyses
control for and this one does not.

**The same observation, done properly, is a live result.** `sfs_shift_by_site.ipynb` tests
loss with a frequency-matched null across all 31 sites and finds the SV extinction excess
robust to every aggregation tried. Cite that, not this.

## Contents

    sfs_time_site4.ipynb
    scripts/_build_sfs_time_site4_nb.py
    scripts/_compute_sfs_time_site4.py
