#!/usr/bin/env python
"""cis-eQTL for the candidate genes, to overlay on the locus figures.

Ported from the MOI-LAB pipeline (`scripts/run_eqtl.py` + `pipelines/eqtl/run_eqtl_one_gene.R`,
see METHODS_FUNCTIONAL_TRACKS.md). Per gene: cut the 1001g biallelic SNP panel to the gene
body +-10 kb, put the expression phenotype in column 6 of the .fam, run `gemma -lm 4`.

This is the SNP version, chosen deliberately. Testing our own indel/SV against expression
would be the mechanistically informative design, but it is not available where it matters:
of the 107 founders carrying both a panel genotype and non-MU expression, GPX6's 1.16 kb
deletion is CALLED IN ONLY 29 (27%) with 1-4 ALT carriers, and FUS3 is the same
(sv-panel-support-asymmetry, landing on the two best candidates). Only CRK14, AT2G30000,
SCPL34 and AT4G13200 have both a call rate above 88% and enough carriers. So the SNP scan
is what can be run across all of them, at n = 471 rather than 107.

WHAT THIS CAN AND CANNOT SAY. Every candidate here is SNP-tagged at r^2 = 1.0, so a cis-eQTL
peak over our variant is expected whether or not the indel is causal. It establishes that
the gene is under cis regulation and that our variant sits on the regulatory haplotype. It
is NOT evidence that the indel itself changes expression.

Expression: TG_data_20180606.Rdata, `TG.genes$d_log2_batch` (batch-corrected log2), batch
`MU` dropped as upstream does, replicate samples averaged per accession (874 samples ->
507 accessions, 471 of them genotyped).

`gemma -lm 4` is a plain linear model with no kinship term, matching upstream. Population
structure is therefore uncorrected; within a 20 kb cis window that mostly rescales the
whole window rather than moving the peak, but the p-values are not calibrated and should
be read as a profile, not a test.

env: kmate for this script; plink2 + gemma from the gwas_tools env. Compute node, ~1 min/gene.
Writes results/eqtl/<AGI>.assoc.txt and results/eqtl/cis_eqtl_all.csv.
"""
from __future__ import annotations
import os, sys, subprocess, shutil
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, GEA_DIR)
sys.path.insert(0, f"{GEA_DIR}/r2_gea_nonsnp/phase1_replication")
import lib                                                       # noqa: E402

PROJ = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
EQ = f"{PROJ}/data/eqtl"
BFILE = f"{EQ}/1001gbi"
OUT = f"{HERE}/results/eqtl"
TMP = f"{OUT}/_tmp"
ENVBIN = "/global/home/users/tbellg/miniforge3/envs/gwas_tools/bin"
PLINK2, GEMMA = f"{ENVBIN}/plink2", f"{ENVBIN}/gemma"
WINDOW = 10_000
MIN_N = 50

CANDIDATES = {"AT4G11600": "GPX6", "AT3G26790": "FUS3", "AT2G30000": "AT2G30000",
              "AT4G13200": "AT4G13200", "AT4G23220": "CRK14", "AT3G51450": "SSL7",
              "AT5G23210": "SCPL34", "AT5G48410": "GLR1.3", "AT2G01180": "LPP1"}


def phenotypes() -> pd.DataFrame:
    """accession -> mean d_log2_batch per candidate gene, MU batch dropped."""
    E = pd.read_csv(f"{EQ}/cand_expression_d_log2_batch.csv", index_col=0)
    M = pd.read_csv(f"{EQ}/TG_meta.csv")
    M["index"] = M["index"].astype(str)
    keep = (M.batch_comb != "MU").values
    E = E.loc[:, keep]                       # columns are TG.meta$index, in meta order
    acc = M.loc[keep, "index"].values
    P = E.T.groupby(acc).mean()              # replicate accessions averaged
    P.index.name = "accession"
    return P


def run_gene(agi: str, sym: str, genes: pd.DataFrame, P: pd.DataFrame, fam: pd.DataFrame):
    g = genes[genes.gene == agi]
    if not len(g):
        print(f"  {sym}: no TAIR10 coordinates"); return None
    ch = int(str(g.chrom.iloc[0]).replace("Chr", ""))
    lo, hi = int(g.start.min()) - WINDOW, int(g.end.max()) + WINDOW
    if agi not in P.columns:
        print(f"  {sym}: not in the expression matrix"); return None

    pheno = fam.iid.map(P[agi]).astype(float)
    n = int(pheno.notna().sum())
    if n < MIN_N:
        print(f"  {sym}: only {n} genotyped samples with expression"); return None

    os.makedirs(TMP, exist_ok=True)
    stem = f"{TMP}/{agi}"
    r = subprocess.run([PLINK2, "--bfile", BFILE, "--chr", str(ch), "--from-bp", str(lo),
                        "--to-bp", str(hi), "--make-bed", "--out", stem],
                       capture_output=True, text=True)
    if r.returncode or not os.path.exists(f"{stem}.bim"):
        print(f"  {sym}: plink2 failed\n{r.stderr[-400:]}"); return None
    nsnp = sum(1 for _ in open(f"{stem}.bim"))

    # phenotype into column 6 of the subset .fam, in ITS order (plink2 may reorder)
    F = pd.read_csv(f"{stem}.fam", sep=r"\s+", header=None)
    ph = F[1].astype(str).map(P[agi]).astype(float)
    F[5] = ph.fillna(-9).values
    F.iloc[:, :6].to_csv(f"{stem}.fam", sep=" ", header=False, index=False)

    # GEMMA defaults (-miss 0.05, -maf 0.01) drop most of the window: the 1001g calls are
    # patchy and at 471 phenotyped samples a 5% missingness ceiling left only ~120-330 of
    # ~2,000-3,000 SNPs, including the tagging SNPs of our own variants -- which makes the
    # "is the eQTL on our haplotype" readout unanswerable. Relaxed to 30%.
    r = subprocess.run([GEMMA, "-bfile", stem, "-lm", "4", "-miss", "0.3", "-maf", "0.01",
                        "-o", agi],
                       capture_output=True, text=True, cwd=OUT)
    f = f"{OUT}/output/{agi}.assoc.txt"
    if not os.path.exists(f):
        print(f"  {sym}: gemma failed\n{(r.stderr or r.stdout)[-400:]}"); return None
    A = pd.read_csv(f, sep="\t")
    A["gene"], A["symbol"], A["chrom"] = agi, sym, f"Chr{ch}"
    A["nlp"] = -np.log10(A.p_wald.clip(lower=1e-300))
    print(f"  {sym:<10} {agi}  Chr{ch}:{lo:,}-{hi:,}  {nsnp:,} SNPs, n={int(ph.notna().sum())}"
          f"  peak -log10p = {A.nlp.max():.2f} at {int(A.loc[A.nlp.idxmax(),'ps']):,}")
    return A


def main():
    os.makedirs(OUT, exist_ok=True)
    P = phenotypes()
    print(f"expression: {P.shape[0]} accessions x {P.shape[1]} candidate genes")
    fam = pd.read_csv(f"{BFILE}.fam", sep=r"\s+", header=None, usecols=[0, 1],
                      names=["fid", "iid"], dtype=str)
    print(f"genotypes : {len(fam)} accessions; overlap = "
          f"{fam.iid.isin(P.index).sum()}\n")
    genes = lib.load_genes()
    out = [a for agi, sym in CANDIDATES.items()
           if (a := run_gene(agi, sym, genes, P, fam)) is not None]
    if out:
        A = pd.concat(out, ignore_index=True)
        A.to_csv(f"{OUT}/cis_eqtl_all.csv", index=False)
        print(f"\nwrote {OUT}/cis_eqtl_all.csv  ({len(A):,} SNP tests)")
    shutil.rmtree(TMP, ignore_errors=True)


if __name__ == "__main__":
    main()
