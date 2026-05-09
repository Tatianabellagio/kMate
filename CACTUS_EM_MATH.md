# cactus_em — algorithm with math notation

Companion to `ALGORITHM.md` (prose) and `MODEL_SPEC.md` (early design). This file states the algorithm in formal math.

## Problem

Given pool-seq reads from a population that is a mixture of $F$ known founder genomes (the panel), recover the **founder mixture vector** $\mathbf{h} \in \Delta^F$ and project it to **per-record alt-allele frequencies**.

For GrENE-Net specifically: $F = 231$ *A. thaliana* founders, $K \approx 8 \times 10^7$ unique panel k-mers genome-wide, $R \approx 3.4 \times 10^6$ atomized VCF records.

## Inputs (precomputed once per panel)

**Founder × k-mer indicator matrix:**
$$
\mathrm{cn} \in \{0,1\}^{F \times K},
\qquad
\mathrm{cn}_{f,k} =
\begin{cases}
1 & \text{founder } f \text{ carries k-mer } k \\
0 & \text{otherwise}
\end{cases}
$$
Built from PanGenie's `kmers.tsv.gz` + panel VCF + reference FASTA. Stored sparse (`scipy.sparse.csr`) at $\sim$30–50 GB for the full GrENE-Net panel.

**Founder × variant indicator matrix:**
$$
\mathrm{cn}_{\text{var}} \in \{0,1\}^{F \times R},
\qquad
\mathrm{cn}_{\text{var}, f, r} = \mathbb{1}[\text{founder } f \text{ carries the alt at variant } r]
$$
Built from the panel VCF GT field after `bcftools norm -m -` decomposition.

## Per-sample observation

For one pool-seq sample we count k-mer hits with jellyfish:
$$
c_k \in \mathbb{Z}_{\ge 0},
\qquad
c_k = \text{(number of read-hits to k-mer } k\text{)}, \quad k = 1, \ldots, K
$$
plus a coverage estimate $\lambda$ derived from total read length / genome size.

## Generative model

Each k-mer count is a Poisson draw whose rate sums founder contributions:
$$
\boxed{\;
c_k \;\sim\; \text{Poisson}\!\left( \lambda \cdot \sum_{f=1}^{F} h_f \cdot \mathrm{cn}_{f,k} \right)
\;}
\qquad k = 1, \ldots, K
$$
with the simplex constraint
$$
\mathbf{h} \in \Delta^{F-1}
\;\;\Longleftrightarrow\;\;
h_f \ge 0 \;\;\forall f, \quad \sum_{f=1}^{F} h_f = 1
$$

Define the per-k-mer expected (un-coverage-scaled) rate:
$$
\mu_k(\mathbf{h}) \;=\; \mathbf{h}^\top \mathrm{cn}_{:,k} \;=\; \sum_f h_f \, \mathrm{cn}_{f,k}
$$

The Poisson log-likelihood (dropping constants) is
$$
\ell(\mathbf{h}) \;=\; \sum_{k=1}^{K} \big[\, c_k \log \mu_k(\mathbf{h}) - \lambda \, \mu_k(\mathbf{h}) \,\big]
$$
This is concave in $\mathbf{h}$ (logarithm of a linear function plus a linear function), so the constrained MLE $\hat{\mathbf{h}} = \arg\max_{\mathbf{h} \in \Delta^{F-1}} \ell(\mathbf{h})$ exists and is unique up to over-determined founder pairs.

## EM derivation

**Latent variables.** Treat each unit of count at k-mer $k$ as a latent draw from one founder among those carrying $k$. Let $z_{k,f}$ = number of count units at $k$ attributed to founder $f$, with $\sum_f z_{k,f} = c_k$ and $z_{k,f} = 0$ if $\mathrm{cn}_{f,k} = 0$.

**E-step.** Given current $\mathbf{h}^{(t)}$, the posterior attribution probability for k-mer $k$, founder $f$ is
$$
P\!\left(f \mid k; \mathbf{h}^{(t)}\right)
\;=\;
\frac{h_f^{(t)} \, \mathrm{cn}_{f,k}}{\sum_{f'} h_{f'}^{(t)} \, \mathrm{cn}_{f',k}}
\;=\;
\frac{h_f^{(t)} \, \mathrm{cn}_{f,k}}{\mu_k(\mathbf{h}^{(t)})}
$$
so the expected attribution is
$$
\mathbb{E}\!\left[\, z_{k,f} \mid c_k, \mathbf{h}^{(t)} \,\right]
\;=\;
c_k \cdot \frac{h_f^{(t)} \, \mathrm{cn}_{f,k}}{\mu_k(\mathbf{h}^{(t)})}
$$

**M-step.** Maximize the expected complete log-likelihood w.r.t. $\mathbf{h}$ subject to $\sum_f h_f = 1$. The lagrangian gives
$$
\boxed{\;
h_f^{(t+1)}
\;\propto\;
\sum_{k=1}^{K}
\mathbb{E}\!\left[\, z_{k,f} \mid c_k, \mathbf{h}^{(t)} \,\right]
\;=\;
h_f^{(t)} \cdot \sum_{k=1}^{K} \mathrm{cn}_{f,k} \cdot \frac{c_k}{\mu_k(\mathbf{h}^{(t)})}
\;}
$$
followed by renormalization $h_f^{(t+1)} \leftarrow h_f^{(t+1)} / \sum_{f'} h_{f'}^{(t+1)}$.

This is the **standard multiplicative-update EM** (Lee–Seung NMF with simplex normalization, identical in form to RSEM/kallisto/salmon for RNA-seq). $\lambda$ cancels in the ratio.

## Vectorized implementation

Three matrix-vector products per iteration ($\mathrm{cn}$ stored sparse):
$$
\boldsymbol{\mu}^{(t)} \;=\; \mathrm{cn}^\top \mathbf{h}^{(t)}
\quad \text{(length } K\text{)}
$$
$$
\mathbf{w}^{(t)} \;=\; \mathbf{c} \,\oslash\, \boldsymbol{\mu}^{(t)}
\quad \text{(elementwise; length } K\text{)}
$$
$$
\widetilde{\mathbf{h}}^{(t+1)} \;=\; \mathbf{h}^{(t)} \,\odot\, \big(\mathrm{cn}\,\mathbf{w}^{(t)}\big)
\quad \text{(length } F\text{)}
$$
$$
\mathbf{h}^{(t+1)} \;=\; \widetilde{\mathbf{h}}^{(t+1)} \,/\, \big\Vert \widetilde{\mathbf{h}}^{(t+1)} \big\Vert_1
$$

In code (`poolfreq/src/em_solver.py:86-107`):

```python
denom   = h @ cn               # μ_k for all k (length K)
cw      = counts / denom       # c_k / μ_k     (length K)
em_term = h * (cn @ cw)        # h_f · Σ_k cn[f,k] · c_k/μ_k  (length F)
h_new   = em_term / em_term.sum()
```

Iterate until $\Vert \mathbf{h}^{(t+1)} - \mathbf{h}^{(t)} \Vert_2 < 10^{-6}$. Typically 30–80 iterations on full Chr1 from a uniform start. Float32 throughout; total wall time 5–15 min/sample.

## Optional priors

For per-window EM in recombinant pools, two regularizers are bolted onto the M-step (`em_solver.py:74-99`).

**Symmetric Dirichlet** with concentration $\alpha$:
$$
h_f^{(t+1)}
\;\propto\;
h_f^{(t)} \sum_k \mathrm{cn}_{f,k} \frac{c_k}{\mu_k} \;+\; (\alpha - 1)
$$
$\alpha = 1$ gives MLE; $\alpha > 1$ pulls toward uniform.

**Anchor toward a prior $\mathbf{h}_{\text{prior}}$** (e.g. the chrom-wide global $\mathbf{h}$), weighted by $\beta$:
$$
h_f^{(t+1)}
\;\propto\;
h_f^{(t)} \sum_k \mathrm{cn}_{f,k} \frac{c_k}{\mu_k} \;+\; \beta \, N \, h_{\text{prior},f}
$$
where $N = \sum_k c_k$ is the total k-mer count. $\beta = 0$ → pure MLE; $\beta = 1$ → prior worth as much as the data; sweet spot $\beta \in [0.05, 0.5]$.

## Projection step

Once $\hat{\mathbf{h}}$ is recovered, the per-record alt frequency is exact:
$$
\boxed{\;
\widehat{\mathrm{AF}}_r
\;=\;
\hat{\mathbf{h}}^\top \mathrm{cn}_{\text{var},:,r}
\;=\;
\sum_{f=1}^{F} \hat{h}_f \, \mathrm{cn}_{\text{var}, f, r}
\;}
$$
Linear in $\hat{\mathbf{h}}$, no estimation noise added. Same form as hapFIRE's $\text{AF} = \mathbf{h}_{\text{hap}}^\top \mathrm{H2S}$ projection (haplotype-to-SNP matrix), only the matrix is $\mathrm{cn}_{\text{var}}$ (founder × variant) and the inputs are k-mer counts instead of phased reads.

## Per-window mode

For evolved samples with mosaic ancestry, partition $[1, L]$ into windows $\{W_w\}$ and run the EM independently per window:
$$
\hat{\mathbf{h}}_w \;=\; \arg\max_{\mathbf{h} \in \Delta^{F-1}} \;\sum_{k \in W_w} \big[\, c_k \log \mu_k(\mathbf{h}) - \lambda \mu_k(\mathbf{h}) \,\big]
$$
Then for variant $r$ in window $w$,
$$
\widehat{\mathrm{AF}}_r \;=\; \hat{\mathbf{h}}_w^\top \mathrm{cn}_{\text{var},:,r}
$$

Produced TSVs: `cactus_em_recomb_window_{Nkb}.tsv`. The window size sweep (10/100/200 kb) is in `results/recomb_window_size_compare.tsv`; on EQUAL n50_g3 R² monotonically improves with window size (0.883 → 0.942 → 0.951), on SKEWED every per-window mode collapses below its global baseline.

## Properties

| Property | Why |
|---|---|
| Monotone-improving | Standard EM result for the Poisson likelihood |
| Stays on the simplex | $h_f^{(t)} = 0 \Rightarrow h_f^{(t+1)} = 0$ (multiplicative); explicit renormalization $\Vert\mathbf{h}\Vert_1 = 1$ |
| $\lambda$-free | Cancels in the M-step ratio; coverage estimate is sanity-only |
| Over-determined at $K \gg F$ | $80 \text{M} \gg 231$, so $\hat{\mathbf{h}}$ is identified up to founder pairs that share every k-mer (rare) |
| Coverage-saturated below 1× | Joint signal across $K$ k-mers is already huge at low per-k-mer count |

## Relationship to HARP and related methods

This is **not HARP code** — same shape, different evidence:

| | cactus_em | HARP (hapFIRE) |
|---|---|---|
| Latent | founder mixture $\mathbf{h}$ | haplotype mixture $\mathbf{h}_{\text{hap}}$ per LD block |
| Evidence | k-mer counts $c_k$ | per-base read pileup $P(\text{base} \mid \text{hap}, q)$ |
| Inference | Poisson EM, multiplicative update | per-base likelihood EM |
| Projection | $\hat{\mathbf{h}}^\top \mathrm{cn}_{\text{var}}$ | $\hat{\mathbf{h}}_{\text{hap}}^\top \mathrm{H2S}$ + CVXPY founder solve |
| Sees SVs? | Yes (k-mers cover bubbles) | No (per-base SNP-only) |

The *form* of the M-step
$$
h_f^{(t+1)} \;\propto\; h_f^{(t)} \sum_k \mathrm{cn}_{f,k} \frac{c_k}{\mu_k(\mathbf{h}^{(t)})}
$$
also recurs in:
- **Lee & Seung (1999)** — multiplicative-update NMF.
- **Dempster, Laird & Rubin (1977)** — original EM for finite mixtures.
- **Bray et al. 2016 (kallisto)**, **Patro et al. 2017 (salmon)**, **Li & Dewey 2011 (RSEM)** — RNA-seq quantification: each read maps to multiple transcripts; EM splits its mass by current transcript abundance estimates. Substitute "transcripts → founders" and "reads → k-mers" and the iteration is essentially identical.

That is why **`kallisto_em` slots into the methods registry** (`build_session_dataset.py:77-78`) as a pseudo-baseline: it solves the same multinomial mixture problem with the same multiplicative-update EM, just under different naming.

## References

- Kessner D et al. (2013) *Mol. Biol. Evol.* 30(5):1145-1158 — HARP.
- Lee D & Seung S (1999) *Nature* 401:788-791 — NMF multiplicative update.
- Dempster A, Laird N & Rubin D (1977) *J. R. Stat. Soc. B* 39:1-38 — EM algorithm.
- Bray N et al. (2016) *Nat. Biotechnol.* 34:525-527 — kallisto pseudoalignment + EM.
- Wu X et al. (2024) — hapFIRE (GrENE-Net Text S4).
- Ebler J et al. (2022) *Nat. Genet.* 54:518-525 — PanGenie.
