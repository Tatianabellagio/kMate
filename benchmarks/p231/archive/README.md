# p231/archive/

Superseded p231 simulation sets. Not consumed by anything — the live sims are in
`../sims/`.

## `sims_buggy_preFix_2026-06-18/`

Six sim regimes (`cov10_n{50,231}_g{0,1}` ± `self97`, `hotspots_p231_chr1`) set aside
2026-06-18 when the simulation was fixed and every regime regenerated.

⚠️ **The specific defect is not recorded anywhere** — not in the git history, not in
any doc. The directory name is the entire provenance. Stated plainly rather than
guessed at: if you need to know what was wrong, it has to be re-derived by diffing
these truth tables against the regenerated ones in `../sims/`.

What *is* safe to say: any benchmark number computed before 2026-06-18 used these
sims. Since truth tables are the yardstick every arm is scored against, a pre-06-18
number is not comparable to a post-06-18 one even for the same estimator.
