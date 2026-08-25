# benchmarks/founder_h_uncertainty/

**Question:** how certain is the founder mixture `ĥ`, and where does that certainty
break down?

**Status: closed-out investigation.** Retained for its finding, not as a maintained
benchmark. Full log: [`SESSION_2026-06-29.md`](SESSION_2026-06-29.md).

## The finding (the reason this dir exists)

The error in `ĥ` is dominated by **identifiability / collinearity bias, not sampling
variance.** Consequences, which are what make this load-bearing rather than trivia:

- In `--unit chrom` (global) the bias is **coverage-independent** — the same at 3× as
  at infinite coverage. You cannot sequence your way out of it.
- In window mode it is ~3.5× larger and *does* respond to coverage, via local k-mer
  support.
- **Variance-based SEs cannot see it.** Fisher information and the parametric
  bootstrap both estimate sampling variance, so their nominal 95% intervals do not
  cover the true `h` at the advertised rate. Treat `--emit-af-se` / bootstrap output
  as a variance estimate only, never as a total-error bar on `ĥ`.

This is why the per-founder `h` vector is reported with the caveat that aggregates
(sums over multiple carriers) are reliable while per-ecotype values are noisy — see
`BACKGROUND.md` limitation 4.

## Contents

| script | what it probed |
|---|---|
| `bench_h_uncertainty.py` | analytic Fisher SE on `ĥ` vs parametric bootstrap (the gold standard) |
| `diag_convergence.py` | the EM convergence / bias floor seen in that benchmark |
| `diag_window.py` | the window-mode arm of the same question |
| `af_certainty.py` | whether the EM bias that kills h-certainty **cancels** in the AF projection |
| `af_calibrate_floor.py` | calibrates the AF identifiability floor `c` for the total AF SE |

`af_certainty.py` is the practically important one: AF is a *sum* over carriers, so the
per-founder bias partially cancels there — which is why AF accuracy stays high (R²>0.99)
even though per-founder `h` is noisy.

Related: `../../docs/EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md`,
`../../src/kmate/h_uncertainty.py`, `ALGORITHM.md`.
