# Open items — things known to be unfinished or untested

Deliberately short. Each entry says what is unresolved, why it matters, and what
"done" looks like. Delete an entry when it is closed; do not let it rot into a
list of things that are actually finished.

Last reviewed: 2026-09-29.

---

## 1. `--het split` has never run on a real phased VCF

**Status:** implemented, tested only on synthetic input.

`kmate decompose --haploidize --het split` turns each diploid sample into two
haplotype columns, which is what makes kMate usable on an outbred, phased panel
(HPRC-style) instead of only on inbred lines. It was verified on a hand-written
3-sample VCF: phase preserved, founder axis doubled, missing calls excluded from
`var_called`, unphased heterozygotes refused.

**Why it matters:** an outbred panel is the case where the default
(`--het missing`) silently discards half the data, so `split` is the path someone
outside this project is most likely to need — and it is the least exercised.

**Done when:** a real phased multi-sample VCF (an HPRC release is the obvious
candidate) has been run end to end — `decompose --het split` → `build-var-pa` →
`build-kmer-pa` → `kmate run` on one pool — and the founder axis, carrier counts
and AF projection have been checked against the source VCF. Watch for phasing
conventions this code does not expect (`|` vs `/` per record, partially phased
blocks, multi-sample phase-set `PS` tags, and ploidy ≠ 2 on sex chromosomes,
which `_gt_to_pair` does not handle).

## 2. Calibrated per-record AF uncertainty (future work)

**Status:** removed from the CLI; the idea is unfinished, not wrong.

`se` is the binomial SE from panel support, `sqrt(p(1-p)/n_called)`. It says how many
founder haplotypes had a genotype call at that record, and nothing about how well the
founder mixture itself was resolved.

An attempt at a calibrated alternative existed as `--emit-af-se`: a Fisher delta-method
SE on the projected AF, floored by a panel identifiability constant,
`sqrt(Fisher_SE^2 + c^2)`, on the reasoning that per-record AF error is bias-dominated so
the Fisher term alone under-covers. It also recorded `eff_rank` and `cond` of the
`h`-covariance as resolvability diagnostics. It was never validated against known truth,
so it was removed rather than left as an option that looks endorsed.

Worth knowing before picking this up: on evonet sample MEAJM013-38 the two are **not a
rescaling of each other**. Over 400k records they are uncorrelated (r = -0.000), with the
calibrated value a median 2.5x the binomial and 12.8x at the 90th percentile. So they
measure different things, and choosing between them needs truth data rather than a
plausibility argument.

**Done when:** a calibrated SE is validated on simulated pools where the true AF is known
(coverage of nominal intervals, across coverage and panel-support strata), and either
replaces `se` or ships as a clearly-labelled second column.

## 3. bioconda still serves the broken 0.1.0

**Status:** 0.1.2 published to PyPI 2026-09-29; bioconda PR pending.

PyPI serves 0.1.2 (verified: the hash matches the built artifact byte-for-byte, and it
passes `selftest` in a pristine env). bioconda still serves **0.1.0**, which cannot count
k-mers and predates the per-founder normalization fix, so the install docs route users to
PyPI and warn against `conda install kmate`.

**Done when:** the bioconda recipe update is merged, and the install docs collapse back to
a single `mamba create ... kmate` line with the warning removed.

## 4. Inherited gap in the decomposition catalog

**Status:** known, accepted for production.

`annotate_vcf.py` builds its atomic catalog with `vcfwave` internally, so ~0.8%
of atomic variants are absent from the catalog and cannot be recovered by symbolic-ID
matching. This is upstream, not ours.

**Done when:** either the Arch-1 graph-native path (`vg deconstruct -a -e` +
`resolve-nested-genotypes`) replaces it — the documented upgrade route in
`panel/arch3/README.md` — or the gap is measured on a current panel and shown to be
tolerable for the analyses that depend on it.

## 5. Docs link into untracked `archive/` trees

**Status:** cosmetic, affects outside readers only.

`archive/` and `archive_gea/` are no longer tracked (they carried a 104 MB file that
blocked every push). Many docs still reference `archive/...` paths, which resolve in a
local checkout but 404 on GitHub. The top-level README says so.

**Done when:** either the references are pruned from user-facing docs, or the
archive is published somewhere citable and the links point there.
