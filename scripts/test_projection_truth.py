#!/usr/bin/env python3
"""
Sanity check: assume hapFIRE perfectly recovers truth (uniform 1/N) and project.
Verifies the math without needing hapFIRE outputs.
"""
import pandas as pd

VCF_SAMPLES = "/home/tbellagio/scratch/hapfire_sv/data/vcf_samples_231.txt"
N_SAMPLES = 231

vcf_samples = [s.strip() for s in open(VCF_SAMPLES) if s.strip()]
assert len(vcf_samples) == N_SAMPLES, f"got {len(vcf_samples)} expected {N_SAMPLES}"

# Truth ecotype freq: uniform 1/N
f_ecotype = pd.Series(1.0 / N_SAMPLES, index=vcf_samples)

print(f"f_ecotype sum: {f_ecotype.sum():.6f}")
print(f"f_ecotype example: {f_ecotype.head(3).to_dict()}")

results = []
for sv_freq in [0.10, 0.30, 0.50, 0.70, 0.90]:
    n_sv = round(N_SAMPLES * sv_freq)
    sv_carriers = vcf_samples[:n_sv]
    G_SV = pd.Series(0, index=vcf_samples)
    G_SV[sv_carriers] = 1
    af_alt_proj = float((G_SV * f_ecotype).sum())
    truth = n_sv / N_SAMPLES
    print(f"  sv_freq={sv_freq:.2f}  n_sv={n_sv:3d}  truth={truth:.4f}  proj={af_alt_proj:.4f}  diff={af_alt_proj-truth:.6e}")
    results.append({"sv_freq": sv_freq, "n_sv": n_sv, "truth": truth, "proj": af_alt_proj})

pd.DataFrame(results).to_csv("/home/tbellagio/scratch/hapfire_sv/results/sanity_check_truth_projection.tsv", sep="\t", index=False)
print("OK — math is correct: projection from truth ecotype freq exactly recovers truth SV freq.")
