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

## 2. Release 0.1.2 is not published

**Status:** version bumped in the repo; nothing uploaded.

PyPI serves **0.1.1**; bioconda serves **0.1.0**, which is broken (wrong
`jellyfish` dependency, and a pre-July-2026 estimator with no `--normalize` /
`--unit` / `--emit-af-se`). The install docs route around this by telling users
to take 0.1.1 from PyPI and warning them off `conda install kmate`.

**Done when:** a 0.1.2 sdist is built and uploaded to PyPI, and the bioconda
recipe's `sha256` — currently flagged stale in `conda/bioconda/meta.yaml`, since
it is the 0.1.1 hash — is regenerated from what PyPI serves.

## 3. The bioconda PR is held on purpose

**Status:** recipe ready and verified; not submitted.

Deliberate: kMate is still changing, and submitting now means asking bioconda
reviewers to re-review each iteration. See the status banner in
`conda/bioconda/SUBMISSION.md`. Note this is an **update** to an existing recipe,
not a new submission.

**Done when:** the API has settled, 0.1.2 (or later) is on PyPI, and the PR is
opened per `SUBMISSION.md`.

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
