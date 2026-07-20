# Superseded — buggy k-mer→block assignment (2026-07-07)

These scripts/notebook assigned k-mers to blocks using **per-bubble** positions
(`bubble_start`, length = n_bubbles = 190,569) to slice the **per-k-mer** presence
matrix (10.9M columns). That mismapped each block to ~190k wrong columns
(~4,854 "k-mers"/block) instead of the true ~1.9M k-mers/block — a ~40x undercount
that spuriously deflated K_b to a fake median of **142** and produced a false
"haploblock beats production 18/18" result.

**Do not trust any K_b / TVD numbers from these.** The bug: `bpos = (bubble_start +
bubble_end)//2` (per-bubble) used to index the per-k-mer matrix. Correct mapping:
`kpos = (bubble_start[bubble_id] + bubble_end[bubble_id])//2` (per-k-mer, via bubble_id).

Corrected artifacts (in the parent dir):
- haploblock_r2010_eps0_validation.ipynb  (CORRECTED; K_b≈231)
- haploblock_r2020_eps0_validation.ipynb
- verify_kb_genomewide.py  + kb_genomewide_summary.json
- _build_validation_nb.py, _build_validation_nb_r2020.py

Kept (not deleted) for audit trail.
