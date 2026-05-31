# Session handoff — kMate

> **Merged into [`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) (2026-05-30).**
> There is no separate "handoff" doc anymore. `docs/PIPELINE_STATE.md` is the
> single source of truth for production inputs (§0), the per-sample run recipe
> (§0.1), reproducibility/environment (§0.2), and what's deprecated. This stub
> stays only so old links don't break.

**Production in one line:** exactly ONE pipeline — **arch3** — decomposes the
panel into `panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz`; **both** K_pa
(`data/kmer_pa_231_arch3_filt2inv`) and V_pa
(`panel/arch3/chr{N}/var_pa_231_arch3_chr{N}`) build from that VCF and nothing
else. Old `founders_231_v3qc*` VCFs and `data/*_231_v3qc_v3*` matrices are
archived. Env: `kmate`. → full detail in `docs/PIPELINE_STATE.md` §0.
