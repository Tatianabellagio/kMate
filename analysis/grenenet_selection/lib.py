"""Shared helpers for the GrENE-Net SV-GEA analysis (kMate outputs).

The goal (mirroring the phase-1 SNP paper, but for SVs): find structural
variants whose frequency rises in hot sites and falls in cold ones — climate-
dependent directional selection — and ask whether SVs add adaptive signal SNPs
miss. See README.md.

Data:
  - kMate per-sample AF TSVs: analysis/grenenet_selection/common/rerun_kfw_hb/evolved/<MLFH...>.tsv
      (--unit chrom + full-panel Kf_w; old-panel grenenet_kmate_arch3 deleted 2026-07-08)
      cols: chrom pos ref_len alt_len alt_freq info n_called se
  - founding p0 (gen 0): the SEEDMIX kMate outputs (mean over 8 reps).
  - sample metadata (site/plot/generation/coverage): Table_S5.
"""
from __future__ import annotations
import os, glob, json
import numpy as np
import pandas as pd
from scipy import stats

# PROJ derived from this file's own location (analysis/grenenet_selection/lib.py), not
# hardcoded, so the repo tolerates future scratch moves without a repeat of the
# 2026-07-22 break (user-scratch kmate tree migrated to project/fc_moilab scratch;
# the old hardcoded /global/scratch/users/tbellg/kmate silently stopped existing).
PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Repointed 2026-07-07 to the full-panel-Kf_w + haploblock (--unit chrom) rerun.
# Prev: rerun_perfounder (pre-Kf_w-fullpanel, 2026-07-06 — STALE). See
# analysis/grenenet_selection/common/rerun_kfw_hb/README.md.
GEA = f"{PROJ}/analysis/grenenet_selection"

# ---- panel version -----------------------------------------------------------------------------------------------
# KMATE_PANEL (env) = arch3v2 (default since 2026-10-22: the production panel, panel/arch3v2/README.md) or arch3 (the
# previous panel, panel/arch3/; only to reproduce arch3-era outputs). Every panel-dependent input and every output
# directory derived from it is defined here; with arch3v2 each output directory gets the suffix "_arch3v2", so the two
# versions never write to the same place (the arch3 outputs stay in the unsuffixed folders). Panel-independent files (TAIR10-derived
# tables, the multi-allelic catalogue both panels decompose) stay in their single location.
PANEL = os.environ.get("KMATE_PANEL", "arch3v2")
if PANEL not in ("arch3", "arch3v2"):
    raise ValueError(f"KMATE_PANEL must be arch3 or arch3v2, not {PANEL!r}")
VSUF = "" if PANEL == "arch3" else f"_{PANEL}"
# Scripts whose every panel-dependent path comes from this block (burden chain, audited 2026-10-20), as paths relative
# to GEA (not bare names: archived copies share some names). Under arch3v2 (the default) any other script refuses to
# start: scripts not yet ported still hard-code arch3 panel files and would mix arch3v2 allele frequencies (lib.AF_STORE, lib.OUT)
# with arch3 genotypes or write into arch3 folders. Make a script panel-aware (paths from here only, outputs under
# versioned dirs), have it audited, then add its path here.
PANEL_AWARE = {
    "genes/convergence/gea_go_enrichment.py", "genes/convergence/candidate_patterns.py", "genes/convergence/sv_defence_load.py",
    "genes/candidates/676_WDR5A/evidence_WDR5A.py", "genes/candidates/337_CPK32/evidence_CPK32.py", "genes/convergence/lineup_founder_control.py", "genes/convergence/review_sheets/fix_indel_carriers.py",
    "genes/convergence/review_sheets/shared_tolerant.py",
    "genes/convergence/review_sheets/fill_caches.py", "genes/convergence/review_sheets/dossier_noncoding.py",
    "genes/top_candidates/top_top_candidates/001_ALA6/selcoef_ALA6.py",
    *(f"genes/candidates/{d}/evidence_{d.split('_', 1)[1]}.py" for d in ("260_GORK", "237_BCA1", "405_EPF2", "269_PHOT1", "375_BLUS1", "373_CRK36", "314_ABA4", "368_CYP707A4", "406_DPBF2", "255_VIP6", "295_EMF1", "282_EAF1B", "287_AHL22", "273_RVE6", "390_MYB102", "396_BGLU27", "645_APS1", "313_NHX7", "293_HTT3", "261_BBX18", "327_ALDH12A1", "292_RGA", "014_PLDGAMMA2", "055_CYP89A9", "243_GSTU26", "332_BHLH23", "329_ACS6", "538_MFB16.4", "050_PLDZETA1", "017_ROSY1", "027_AT1G32030", "184_WOX2", "185_ARGAH1", "325_RBOHF")),   # functional review sheets, 2026-10-31
    *(f"genes/candidates/{d}/evidence_{d.split('_', 1)[1]}.py" for d in ("004_SKS16", "011_AT5G27880", "097_AT4G11350", "1078_HMGS", "120_AT3G18952", "121_CAF2", "128_DTX40", "133_AT3G43829", "173_EXO70E1", "235_LIG6", "239_RKF1", "244_CIPK9", "246_DRG", "248_AT4G13490", "249_BBX24", "250_AHL20", "254_SUB", "257_SBT5.4", "259_FUS3", "277_AGL14", "278_AIPP3", "279_FZR1", "280_FLXL2", "283_AGL5", "307_HBP5", "322_LFNR2", "335_ASG1", "340_WRKY57", "352_HHO1", "360_WRKY66", "365_CBP60G", "371_MKK9", "376_MYOB13", "380_UBC33", "394_FLZ18", "397_HEMB2", "407_AT3G11325", "445_AT1G10280", "449_TIR1", "451_RPL21", "458_PAT08", "476_D6PKL1", "484_TPT", "551_AT3G24610", "552_SUS4", "556_AHL26", "561_AT2G41390", "562_MORC4", "661_SPL13A", "829_MMG4.18", "832_PEX19-2", "848_CC1", "858_AT1G32920", "920_AT3G46650", "922_LTPG7", "926_CASP3", "952_CML38")),   # the rest of genes/unranked, 2026-11-01
    *(f"genes/top_candidates/{d}/evidence_{d.split('_', 1)[1]}.py" for d in ("017_SMXL6", "020_HRS1", "044_SAUR74", "050_AT1G32430", "051_AT2G17160", "052_MNF13.9", "053_RPP13", "063_TPS05", "070_AGL100", "077_PGIP1", "078_2ODDC23.1", "081_ANP2", "083_PNC2", "097_WRKY17", "104_AL3", "107_CYP71B21", "111_AT4G07868", "202_TRM21", "231_OHP1", "294_PUB8", "295_ARI12", "314_EAP3", "439_CML35")),   # sheets in place, 2026-11-01
    *(f"genes/ranked/{d}/evidence_{d.split('_', 1)[1]}.py" for d in ("005_CRK45", "006_TAF14B", "007_AT4G07390", "008_BLH9", "009_EFL3", "010_PER50", "015_NQR", "018_CSN2", "023_ADT4", "027_MEBL", "029_AT2G13623", "030_GBSS1", "032_AT2G36540", "049_BGLU43", "056_AT5G27820", "059_AT2G23210", "106_ERF094", "145_AT4G13450", "169_MRN17.10", "170_GPT2", "181_NAC061", "209_WRKY62", "215_GSH1", "239_AZI1", "240_MKD10.60", "245_SVP", "253_TCP15", "267_APR2", "343_FBL8", "349_AT1G68862", "719_AT2G17695")),   # sheets in place, 2026-11-01
    "genes/flc_callrate/flc_cr30_report.py", "qc/big_bubble_impact.py", "genes/flc_callrate/flc_advisor_tests.py", "genes/flc_callrate/flc_hapfire.py", "genes/flc_callrate/flc_short_read.py", "genes/top_candidates/presentation_examples.py",
    "genes/flc_callrate/flc_split_gea.py", "genes/flc_callrate/n_bubbles_scan.py",
    *(f"burden/{f}.py" for f in ("founder_snps", "build_tables", "allele_change", "allele_matching", "tests", "diag_allele_change_startfreq", "allele_change_variants")),
    *(f"r1_sv_negative_selection/{f}.py" for f in (      # the SV event layer and the burden's gene-level tables
        "_build_bubble_map", "_build_event_map", "_sv_origin_partition", "_sv_event_catalog", "_sv_event_carriers",
        "_build_event_units", "_sv_event_classes", "_sv_garden_af", "_make_graffite_vcf", "_te_mask_prep",
        "_te_classify", "_insertion_asm_genes", "_helixer_only_protein_check", "_sv_gcargo_classes",
        "_sv_allele_change_table", "_sv_classic_sfs", "_sv_sfs_context", "_sv_cds_genes", "_cds_gene_evolution",
        "_cds_gene_evolution_checks", "_cds_sv_carrier_expression", "_cds_sv_expression_lmm", "_cds_genes_matched", "_cds_genes_classes", "_cds_genes_classes_followup")),
    # defense slides (Chapter 3): read the paper notebooks' data cells / lib paths only
    "defense_slides/make_load_slides.py", "defense_slides/make_load_paper_slides.py", "defense_slides/make_intro_slides.py",
    "defense_slides/make_gea_manhattan_slides.py", "defense_slides/audit_example_candidates.py", "defense_slides/make_example_slides.py", "defense_slides/example_axes.py", "defense_slides/make_gea_atlas_slides.py", "defense_slides/make_gea_tagging_slides.py", "defense_slides/make_gea_pair_slides.py", "defense_slides/make_crk_slides.py", "defense_slides/make_pkp3_slide.py", "defense_slides/make_pkp3_cam5_style.py", "defense_slides/make_aba_slide.py", "defense_slides/make_aop_sites_scatter.py",                                            # PC1 SV Manhattan, data and names from notebooks/figure_gea_pc1_zoom.py
    # kMate-derived roots
    "r3_persite_gwas/ecotype_fitness.py", "common/build_af_store.py", "common/build_garden_climate.py", "common/build_p0.py", "common/build_pool_matrix.py",
    "r2_gea_nonsnp/phase1_replication/assemble_lastgen.py",
    # the GEA, steps 1-5 of r2_gea_nonsnp/phase1_replication/_submit_gea.sh (audited 2026-10-21)
    "r2_gea_nonsnp/phase1_replication/build_class_matrices.py", "r2_gea_nonsnp/phase1_replication/build_lfmm_input.py",
    "r2_gea_nonsnp/phase1_replication/multiaxis/reblock_cr50.py", "genes/convergence/event_gea/build_event_matrix.py",
    # the GEA, step 6 (audited 2026-10-21)
    "genes/convergence/class_signal_counts.py", "genes/convergence/annotate_gea_universe.py",
    "genes/convergence/gea_hit_landscape_data.py", "genes/convergence/gea_sv_cargo_data.py",
    "genes/convergence/gea_results_numbers.py", "r2_gea_nonsnp/phase1_replication/multiaxis/raw_block_significance.py",
    # the two LD paper figures, SNP tagging and LD decay (audited 2026-10-21)
    *(f"r3_persite_gwas/{f}.py" for f in ("build_tagging_masked", "build_tagging_events", "_build_tagging_masked_nb", "_audit_tagging")),
    "qc/_compute_ld_decay_panels.py", "qc/_build_ld_decay_panels_nb.py",
    # the panel-composition and Truvari-collapse paper figures (audited 2026-10-21)
    "qc/_build_panel_stats_nb.py", *(f"r1_sv_negative_selection/{f}.py" for f in (
        "truvari_sweep", "truvari_pos_range", "sv_sfs_alleles_vs_events", "_build_truvari_sweep_nb")),
    "r1_sv_negative_selection/_build_insertion_cargo_rm_nb.py",   # paper figures sv_classes/, sv_insertion_genes/ (audited 2026-10-21)
    # the v2 candidate list, steps 1-3 and the coding layer (audited 2026-10-22)
    "genes/convergence/v2/step1_gea_pool.py", "genes/convergence/v2/layer_coding.py",
    "genes/convergence/v2/step2_candidate_genes.py", "genes/convergence/v2/step3_functional_ranking.py",
    # port A, burden figures (2026-10-22):
    # port B, sfs_per_region chain (2026-10-22):
    "r1_sv_negative_selection/_sv_event_polarity.py", "r1_sv_negative_selection/_build_main_fig1_nb.py",
    # port C, GEA figure and hit landscape (2026-10-22):
    *(f"genes/convergence/{f}.py" for f in ("_build_gea_hit_landscape_nb", "gea_region_baseline", "gea_promoter_regulatory")),
    # port D, v2 candidate list steps 4-5 and layers (2026-10-22):
    *(f"genes/convergence/v2/{f}.py" for f in ("layer_atac", "layer_snp_visibility", "layer_expression", "layer_literature",
                                                "layer_eqtl", "layer_tfbs", "step4_merge_evidence", "step5_manifest")),
    # port F, Drive AF export (2026-10-22):
    "common/build_gen_matrices.py", "common/build_gen_matrix_sv.py", "common/append_n_called_to_indexes.py",
    "r2_gea_nonsnp/phase1_replication/build_class_matrices_seedmix.py",
    # port E, dossier layer (2026-10-22):
    *(f"genes/convergence/{f}.py" for f in (
        "candidate_dossier", "audit_candidate", "snp_dossier", "neighbour_expression", "carrier_expression_test",
        "run_eqtl_cis", "plot_expression", "plot_atac", "build_sv_hit_dossier", "build_snp_tagging", "atac_overlap",
        "tfbs_turnover", "coding_evidence", "dapseq_binding", "copy_evidence_to_paper_figures")),
    *(f"genes/convergence/v2/{f}.py" for f in (
        "locus_data", "locus_plots", "_build_locus_notebook", "build_catalog", "hidden_pool", "hidden_finalize")),
    *(f"genes/dissection/{f}.py" for f in (
        "plot_locus_combined", "plot_variant_garden_grid", "plot_variant_garden_ribbon", "plot_variant_garden_change")),
    "common/build_group_means.py",
    # port E, evidence-sheet helpers (2026-10-23): the per-folder panels the dossiers and evidence sheets use
    *(f"genes/convergence/{f}.py" for f in (
        "plot_sv_cargo", "plot_axes_manhattan", "flowering_time_test", "dubin_cold_expression", "tfbs_event_union",
        "drought_forest", "atlas_expression")),
    "notebooks/figure_gea_pkp3.py",                     # draws with 048_PKP3/evidence_PKP3.py (arch3v2)
    "notebooks/figure_gea_pc1_zoom.py", "notebooks/figure_gea_pc3_zoom.py",                 # the same evidence module, climate PC3, one zoom per chromosome
    "genes/te_expression/te_cis_expression.py", "genes/te_expression/te_cis_summary.py",         # SV TE content vs cis expression, genome-wide
    "genes/te_expression/te_region_vs_content.py",      # TE region vs TE content enrichment among SV hits
    "genes/te_expression/te_enrichment_confounders.py", # the same enrichment against linkage, genotyping and carrier origin
    # the dossier extras (2026-10-23): SNP-peak LD, neighbour expression, selection coefficients, the AOP configurations
    *(f"genes/convergence/{f}.py" for f in ("top_snp_sv_neighbours", "top_snp_ld", "neighbour_expression_figure")),
    # defense ch. 3 (2026-11-12): do SVs carry more of the positive (climate) signal than SNPs? tagging vs matched null,
    # single-causal fine-mapping share
    *(f"genes/convergence/{f}.py" for f in ("top_snp_sv_tagging_null", "sv_finemap_share", "snp_sv_ld_by_rank", "selcoef_sv_vs_snp",
                                            "direct_selection_finemap", "snp_threshold_overlap",
                                            "gea_common_threshold_conditional")),
    "genes/top_candidates/top_top_candidates/selection_coef/run_selcoef.py",
    "genes/top_candidates/top_top_candidates/selection_coef/run_selcoef_text_genes.py", "genes/top_candidates/top_top_candidates/selection_coef/allgen_s.py",   # the same s, every gene the GEA text names
    *(f"genes/top_candidates/top_top_candidates/177_AOP1/scripts/aop_{f}.py" for f in ("unit_map", "genotypes", "configurations", "configurations_simple")),
    *(f"genes/top_candidates/top_top_candidates/013_CRK_cluster/scripts/crk_cluster_{f}.py" for f in (   # the cluster's own
        "evidence", "region", "content", "sheet", "svlocus", "tubemap", "diversity")),                             # scripts (arch3v2)
    # the top_top genes' own evidence scripts, arch3v2 versions since 2026-10-23 (one per gene; the arch3 ones are backed up)
    *(f"genes/top_candidates/top_top_candidates/{p}" for p in (
        "001_ALA6/evidence_ALA6.py",
        "002_SPL9/evidence_SPL9.py",
        "016_CIPK3/evidence_CIPK3.py",
        "037_VTI13/evidence_VTI13.py",
        "047_AT1G52180/evidence_AT1G52180.py",
        "048_PKP3/evidence_PKP3.py",
        "066_AT3G27700/evidence_AT3G27700.py",
        "068_HIP1/evidence_HIP1.py",
        "073_MLE8.5/SIG5/evidence_SIG5.py",
        "073_MLE8.5/evidence_MLE8.5.py",
        "084_DTX5/evidence_DTX5.py",
        "086_AT1G49330/evidence_AT1G49330.py",
        "177_AOP1/evidence_AOP1.py",
        "216_CKB3/evidence_CKB3.py",
        "223_SDC/evidence_SDC.py",
        "227_PAN/evidence_PAN.py",
        "297_AFP2/evidence_AFP2.py",
        "303_ABI5/evidence_ABI5.py",
        "330_SMXL7/evidence_SMXL7.py",
        "437_MES3/evidence_MES3.py",
        "secondary/004_EARLI1/evidence_EARLI1.py",
        "secondary/011_TIC/evidence_TIC.py",
    )),
}


def _panel_guard():
    if PANEL == "arch3":
        return
    import sys, __main__
    if "ipykernel" in sys.modules:   # a notebook kernel: the guard cannot see which notebook, so it warns instead
        print(f"lib: KMATE_PANEL={PANEL} (paths from lib only; a notebook that hard-codes arch3 files mixes panels)",
              file=sys.stderr)
        return
    script = getattr(__main__, "__file__", "") or (sys.argv[0] if sys.argv else "")
    rel = os.path.relpath(os.path.realpath(script), os.path.realpath(GEA)) if script else ""
    if rel not in PANEL_AWARE:
        raise RuntimeError(f"KMATE_PANEL={PANEL}, but {rel or 'this session'} is not panel-aware (lib.PANEL_AWARE): "
                           "it may read arch3 panel files or write to arch3 folders. Make it panel-aware first.")


def check_version_dir(path: str, what: str) -> None:
    """Refuse an output/input directory that belongs to the other panel version (name must end in VSUF)."""
    name = os.path.basename(os.path.realpath(path).rstrip("/"))
    ok = name.endswith(VSUF) if VSUF else not name.endswith("_arch3v2")
    if not ok:
        raise SystemExit(f"{what} {path} does not belong to panel {PANEL} (KMATE_PANEL); refusing")


_panel_guard()


def var_pa_prefix(chrom: str) -> str:
    """{prefix}.{var_pa,var_called,meta}.npz of the panel version, chrom 'Chr1'..'Chr5' or 'chr1'..'chr5'."""
    c = chrom.lower()
    return f"{PROJ}/panel/{PANEL}/{c}/var_pa_231_{PANEL}_{c}"


def panel_vcf(chrom: str) -> str:
    """The panel version's final VCF of one chromosome."""
    c = chrom.lower()
    return f"{PROJ}/panel/{PANEL}/{c}/merged_231_{c}_final.vcf.gz"


def catalogue_vcf(chrom: str) -> str:
    """The annotated multi-allelic catalogue (bubbles with per-ALT atomic IDs); shared by every panel version."""
    c = chrom.lower()
    return f"{PROJ}/panel/arch3/{c}/{c}_135_annotated.sorted.vcf.gz"


KMER_PA_PREFIX = f"{PROJ}/data/kmer_pa_231_{PANEL}_filt2inv/kmer_pa"   # + _ChrN.{kmer_pa,meta}.npz
# kMate per-sample outputs of the GrENE-Net cohort (evolved samples and the 8 SEEDMIX replicates).
# arch3v2: common/rerun_arch3v2/ (submit_arch3v2.sh, same recipe); arch3: common/rerun_kfw_hb/ (its README.md).
_KMATE_RUN = "rerun_kfw_hb" if PANEL == "arch3" else f"rerun{VSUF}"
OUT = f"{GEA}/common/{_KMATE_RUN}/evolved"
SEEDMIX = f"{GEA}/common/{_KMATE_RUN}/seedmix"
POOL_MATRICES = f"{GEA}/common/results/pool_matrices{VSUF}"             # build_pool_matrix.py, assemble_lastgen.py
SV_ADAPTIVE = f"{GEA}/r1_sv_negative_selection/results/sv_adaptive{VSUF}"    # the SV event layer and its tables
SV_ADAPTIVE_SHARED = f"{GEA}/r1_sv_negative_selection/results/sv_adaptive"   # panel-independent files kept there
FITNESS_DIR = f"{GEA}/r3_persite_gwas/results/ecotype_fitness{VSUF}"         # ecotype_fitness.py
BURDEN_RESULTS = f"{GEA}/burden/results{VSUF}"
# the GEA and the candidate-gene layer (port E, 2026-10-22); arch3 = the unsuffixed folders
GEA_TILE = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/tile_cr50{VSUF}"   # LFMM tables, lfmm_<cls>_gen9_<axis>.csv
GEA_CLASS_MATRICES = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices_cr50{VSUF}"   # gen9.pools.csv, <cls>_gen9.records.csv
EVENT_GEA = f"{GEA}/genes/convergence/event_gea/cr50{VSUF}"            # nonsnp/sv_gen9_{units.csv,dp.npy}
CONV_RESULTS = f"{GEA}/genes/convergence/results{VSUF}"                # convergence/dossier outputs of the panel version
CONV_RESULTS_SHARED = f"{GEA}/genes/convergence/results"              # panel-independent: eqtl/ (1001G SNPs), symbol/protein caches
V2_LIST = f"{CONV_RESULTS}/candidate_list_v2_events_cr50"              # the v2 candidate list (v2/step1-5)
GENES_DIR = f"{GEA}/genes"
# one dossier tree for every panel (user, 2026-10-22: no per-panel duplicate folders); arch3v2 is the list it follows
DOSSIER_TOP = f"{GENES_DIR}/top_candidates"                           # dossier folders (candidate_dossier.py)
DOSSIER_RANKED = f"{GENES_DIR}/ranked"                                # second-tier dossier folders
DOSSIER_CATALOG = f"{GENES_DIR}/candidates"                           # catalog of the whole v2 list + INDEX.csv (v2/build_catalog.py)
DISSECTION_PLOTS = f"{GENES_DIR}/dissection/results/plots"       # raw renders of grid / ribbon / garden_change
DISSECTION_LOCI = f"{GENES_DIR}/dissection/loci"                 # raw renders of plot_locus_combined.py
DISSECTION_REGIONS = f"{GENES_DIR}/dissection/regions"           # plot_locus_combined.py's panel-VCF window extracts


def save_unlinked(fig, path: str, **kw) -> None:
    """fig.savefig(path) after removing any existing file at path. Renders under genes/dissection/ are hard-linked into
    the dossier folders (genes/top_candidates/, ranked/, ...), so writing over the same name in place would silently
    change the curated dossier copy too; unlinking first leaves that copy as it was."""
    if os.path.lexists(path):
        os.remove(path)
    fig.savefig(path, **kw)
GEN_MATRICES = f"{GEA}/common/results/gen_matrices{VSUF}"              # per-sample gen1-3 AF (common/build_gen_matrices.py)
GROUP_MEANS = f"{GEA}/common/results/group_means{VSUF}.npz"            # build_group_means() cache (hot/cold pilot dp, p0)
T5 = ("/global/scratch/users/tbellg/pang/grenenet_reads/"
      "Table_S5_sample_collection_sequencing_library.csv")
# Authoritative GrENE-Net sample table (has usesample, flowerscollected,
# fix_57_generation) — the phase-1 pooling keys. site_gen_plot pools built from
# its usesample==True rows reproduce the phase-1 merged_hapFIRE columns exactly.
SAMPLES_DATA = ("/global/scratch/projects/fc_moilab/projects/grenenet-phase1/"
                "frequency/hapFIRE_frequencies/samples_data_fix57.csv")
# per-site ERA5 bioclim (bio1-19), all 31 cohort sites (+ more); col `site`.
BIOCLIM = ("/global/scratch/projects/fc_moilab/projects/grenenet-phase1/"
           "drive_zenodo/data-intermediate/bioclimvars_experimental_sites_era5.csv")
AF_STORE = f"{GEA}/common/results/af_store{VSUF}"      # compact per-sample NPY store (build_af_store.py)
AF_SCALE, AF_NAN = 10000, 65535          # uint16 AF encoding (4-decimal + NaN sentinel)
# TAIR10 gene annotation (Chr1..Chr5, matches our SV `chrom` exactly): one row per
# gene. Sibling *_genes_transposons.gff adds TEs (the adaptive-SV class). Canonical
# copy lives in ~/ara_key_files (stable home location).
ARA_KEYS = "/global/home/users/tbellg/ara_key_files"
TAIR10_GENES = f"{ARA_KEYS}/TAIR10_GFF3_genes.gff"
TAIR10_GENES_TE = f"{ARA_KEYS}/TAIR10_GFF3_genes_transposons.gff"
# phase-1 hapFIRE LD haploblocks: 1.05M SNPs -> 16,674 blocks (id `chrom_idx`,
# e.g. '1_0'), TAIR10 Chr coords. The 12 kendall_* files share one SNP->block map.
LD_BLOCKS = ("/global/scratch/users/tbellg/gea_grene-net/ARCHIVE/"
             "linages_wza_picmin/kendall_0_w_id_n_blocks.csv")
# Finer LD blocks from BigLD on the 231-founder panel (kMate mcf90 recompute):
# per-chrom TSVs (chrom, start_pos, end_pos, n_variants), one r2 threshold each.
# clq0.9 = 58,376 blocks (vs 16,674 hapFIRE) -> the finer window definition for the
# phase-1 GEA re-run. Interval-based (unlike hapFIRE's nearest-SNP map), so they do
# NOT tile the genome: variants in inter-block gaps are unassigned (return '').
CLQ_BLOCKS_DIR = f"{GEA}/blocks/results/blocks_mcf90"

# ==================================================================================================================
# ⛔ EXCLUDED GARDENS -- site 33 (Lausanne) is excluded from EVERY analysis (decided 2026-10-21, following GrENE-net
# phase 1, which removed it by hand from its last-generation GEA, genomic offset and leave-one-out analyses).
# Why: no usable generation-1 data (plots 7 and 10 failed lab work, plot 1 = 7 flowers, QC-excluded) and no plot
# continuity: the six plots sequenced in generation 2 (4, 5, 6, 9, 11, 12) had zero year-1 survivors, while the
# year-1 plots (1, 7, 10) were dead in year 2 (phase-1 survival.csv). Its "generation 2" is therefore not two
# generations of selection in a plot. Applied here, in pool_table(), load_climate() and garden_climate(), so every
# pool, garden list, climate table and fitness table built through lib drops it. A script that reads the sample
# tables (SAMPLES_DATA, final_gen.csv, survival.csv, ...) without lib must drop EXCLUDED_SITES itself.
EXCLUDED_SITES = frozenset({33})
# ==================================================================================================================

# Pilot sites (extend as the cohort grows). site -> (label, role)
SITE_CLIMATE = {4: "hot", 54: "cold"}   # 4=Cadiz/Madrid region Spain, 54=Cologne DE
BIO_COLS = [f"bio{i}" for i in range(1, 20)]


def load_climate() -> pd.DataFrame:
    """Per-site bio1-19 (ERA5), indexed by int site. The full-GEA climate axis."""
    c = pd.read_csv(BIOCLIM)
    c["site"] = c["site"].astype(int)
    c = c[~c.site.isin(EXCLUDED_SITES)]               # ⛔ EXCLUDED_SITES (site 33)
    return c.set_index("site")[BIO_COLS]


GARDEN_CLIMATE = f"{GEA}/common/results/garden_climate.csv"     # common/build_garden_climate.py
def garden_climate() -> pd.DataFrame:
    """The 30 gardens with GEA data x 22 axes (bio1-19, pc1-3), indexed by int site. The one garden-climate table.
    The file holds the 30 (rebuilt 2026-10-22 from the arch3v2 pools); EXCLUDED_SITES are dropped here as well."""
    g = pd.read_csv(GARDEN_CLIMATE, index_col="site")
    return g[~g.index.isin(EXCLUDED_SITES)]


def climate_range(axis: str) -> tuple[float, float]:
    """(min, max) of a climate axis over the 30 gardens: the colour range every garden plot uses for that axis."""
    v = garden_climate()[axis]
    return float(v.min()), float(v.max())


FOUNDER_CLIMATE = f"{GEA}/common/results/founder_climate.csv"   # common/build_founder_climate.py
def founder_climate() -> pd.DataFrame:
    """The 231 panel founders x bio1-19 (home climate), indexed by ecotype id. The one founder-climate table."""
    return pd.read_csv(FOUNDER_CLIMATE, index_col="ecotypeid")


def founder_range(axis: str) -> tuple[float, float]:
    """(min, max) of a home-climate variable over the 231 founders: the colour range every founder plot uses."""
    v = founder_climate()[axis]
    return float(v.min()), float(v.max())


def decode_af(u: np.ndarray) -> np.ndarray:
    """uint16 AF store -> float32 alt_freq in [0,1]; AF_NAN sentinel -> NaN."""
    f = u.astype(np.float32) / AF_SCALE
    f[u == AF_NAN] = np.nan
    return f


# QC exclusion: samples with too little USABLE panel data (Chr1 nonzero-k-mer
# fraction < 0.10) — dead/contaminated libraries whose h/AF are garbage. NOT a
# sequencing-depth cut (depth doesn't isolate them; corr(depth,nzfrac)~0.57). See
# analysis/grenenet_selection/qc/results/qc_coverage_audit.csv + analysis/grenenet_selection/notebooks/qc_coverage_audit.ipynb. Applied globally so no pool/
# site/analysis sees them. Set 2026-07-07.
QC_EXCLUDE_FILE = f"{os.path.dirname(os.path.abspath(__file__))}/../../data/qc_lowcov_exclude.txt"
def qc_excluded() -> set[str]:
    """Sample IDs excluded by the low-usable-data QC (see QC_EXCLUDE_FILE)."""
    try:
        with open(QC_EXCLUDE_FILE) as fh:
            return {ln.split("\t")[0].strip() for ln in fh
                    if ln.strip() and not ln.startswith("#")}
    except FileNotFoundError:
        return set()


CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def genome_h(samp: str, base: str) -> np.ndarray | None:
    """genome-wide founder h for one sample = mean over chroms of the per-chrom h.

    Reads the --unit chrom cohort format `{base}/{samp}_{ch}.h_per_chrom.npz`, key
    `{ch}` (a single 231-vector per chrom).
    """
    gs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_per_chrom.npz"
        if not os.path.exists(f):
            return None
        gs.append(np.load(f, allow_pickle=True)[ch].astype(np.float64))
    return np.mean(gs, 0)


def list_samples(base: str = OUT) -> list[str]:
    """Genome-wide per-sample TSVs present in `base` (excludes per-chrom files and
    QC-excluded low-usable-data samples)."""
    excl = qc_excluded()
    return sorted(s for s in (os.path.basename(p)[:-4] for p in glob.glob(f"{base}/*.tsv")
                              if "_Chr" not in os.path.basename(p))
                  if s not in excl)


def sample_map(samples: list[str] | None = None) -> pd.DataFrame:
    """Map sample_id -> site, plot, date, year, generation, coverage, climate.

    Indexed by sample_id. `coverage` is Table_S5 weighted_mean_coverage.
    """
    t5 = pd.read_csv(T5).set_index("sampleid")
    if samples is None:
        samples = list_samples()
    keep = [s for s in samples if s in t5.index]
    # `coverage` = per-SAMPLE; `weighted_mean_coverage` = per-PLOT aggregate
    # (constant across a plot's generations) — keep both, default to per-sample.
    m = t5.loc[keep, ["site", "plot", "date", "year", "generation",
                      "coverage", "weighted_mean_coverage"]].copy()
    m = m.rename(columns={"weighted_mean_coverage": "wmean_cov_plot"})
    m["climate"] = m["site"].map(SITE_CLIMATE)
    return m


def cohort_meta(samples: list[str]) -> pd.DataFrame:
    """Full-GEA per-sample metadata: site/plot/generation/coverage + bio1-19.

    Indexed by sample_id, in the given `samples` order. Joins Table_S5 to the
    per-site ERA5 bioclim — the continuous climate axis for the 31-site GEA
    (supersedes the binary SITE_CLIMATE used by the 2-site pilot).
    """
    t5 = pd.read_csv(T5).set_index("sampleid")
    keep = [s for s in samples if s in t5.index]
    m = t5.loc[keep, ["site", "plot", "date", "year", "generation",
                      "coverage", "weighted_mean_coverage"]].copy()
    m = m.rename(columns={"weighted_mean_coverage": "wmean_cov_plot"})
    m["site"] = m["site"].astype(int)
    return m.join(load_climate(), on="site")


def pool_table(store: str = AF_STORE) -> pd.DataFrame:
    """Timepoint -> site_gen_plot pool map with flower weights, for samples in the
    store. One row per timepoint sample; columns: sample_id, pool, site, plot,
    generation (=fix_57_generation), flowerscollected, coverage. `pool` =
    "{site}_{generation}_{plot}" — the phase-1 merge unit (745 pools)."""
    sd = pd.read_csv(SAMPLES_DATA)
    sd = sd[sd["usesample"]].copy()
    sd = sd[~sd.site.astype(int).isin(EXCLUDED_SITES)]   # ⛔ EXCLUDED_SITES (site 33)
    have = {str(s) for s in np.load(f"{store}/samples.npy", allow_pickle=True)}
    have -= qc_excluded()                       # drop low-usable-data / contaminated libraries
    sd = sd[sd.sampleid.isin(have)]
    sd["generation"] = sd["fix_57_generation"].astype(int)
    sd["pool"] = (sd.site.astype(str) + "_" + sd.generation.astype(str)
                  + "_" + sd["plot"].astype(str))
    return sd[["sampleid", "pool", "site", "plot", "generation",
               "flowerscollected", "coverage"]].reset_index(drop=True)


def load_af(sample: str, base: str = OUT) -> pd.DataFrame:
    """Load one sample's per-record AF TSV."""
    return pd.read_csv(f"{base}/{sample}.tsv", sep="\t")


def rec_key(df: pd.DataFrame) -> pd.Series:
    """Stable per-record key (chrom:pos:ref_len:alt_len) for joining samples.

    kMate TSVs carry allele *lengths*, not sequences, so this keys on
    (locus, ref_len, alt_len). NOT unique: 8.70% of arch3v2 records share it with another record
    (49.5% of SVs). Join on record order (af_store index / `col`) instead."""
    return (df.chrom.astype(str) + ":" + df.pos.astype(str) + ":"
            + df.ref_len.astype(str) + ":" + df.alt_len.astype(str))


SV_MIN_BP = 50   # SV threshold: length change must exceed this many bp


def is_snp(df: pd.DataFrame) -> pd.Series:
    return (df.ref_len == 1) & (df.alt_len == 1)


def is_indel(df: pd.DataFrame) -> pd.Series:
    """Any non-SNP record (ref or alt not a single base) — includes small indels."""
    return (df.ref_len != 1) | (df.alt_len != 1)


def is_sv(df: pd.DataFrame, min_bp: int = SV_MIN_BP) -> pd.Series:
    """SV = large indel: length change |alt_len - ref_len| > min_bp (default 50 bp)."""
    return (df.alt_len - df.ref_len).abs() > min_bp


def sv_size(df: pd.DataFrame) -> pd.Series:
    """SV length proxy = |alt_len - ref_len| (0 for SNPs/MNPs of equal length)."""
    return (df.alt_len - df.ref_len).abs()


def build_p0(sv_only: bool = False, cache: bool = True) -> pd.Series:
    """Founding (gen-0) per-record AF = mean alt_freq over the 8 SEEDMIX reps.

    Returns a Series indexed by rec_key. Cached to analysis/grenenet_selection/.
    """
    tag = "sv" if sv_only else "all"
    cache_path = f"{GEA}/p0_seedmix_{tag}.pkl"   # pickle: no pyarrow dependency
    if cache and os.path.exists(cache_path):
        return pd.read_pickle(cache_path)
    fs = [f for f in glob.glob(f"{SEEDMIX}/SEEDMIX_S*.tsv")
          if "_Chr" not in os.path.basename(f)]
    if not fs:
        raise FileNotFoundError(f"no SEEDMIX TSVs in {SEEDMIX}")
    tot = None
    n = np.zeros(0)
    for f in fs:
        d = pd.read_csv(f, sep="\t")
        if sv_only:
            d = d[is_sv(d)]
        s = pd.Series(d.alt_freq.to_numpy(), index=rec_key(d))
        if tot is None:
            tot = s.fillna(0.0)
            n = s.notna().astype(float)
        else:
            tot = tot.add(s.fillna(0.0), fill_value=0.0)
            n = n.add(s.notna().astype(float), fill_value=0.0)
    p0 = (tot / n.replace(0, np.nan)).rename("p0")
    if cache:
        os.makedirs(GEA, exist_ok=True)
        p0.to_pickle(cache_path)
    return p0


def build_qc_table(samples: list[str] | None = None, base: str = OUT,
                   out_csv: str | None = None) -> pd.DataFrame:
    """One-time per-sample QC summary (reads each TSV's lean columns once).

    Per sample: record/SV/SNP counts, finite-AF fraction, AF means (all + SV),
    called-mask stats, + joined coverage/generation/climate + eff_n_founders.
    Cache to analysis/grenenet_selection/pilot_qc.csv so notebooks just load it.
    """
    if samples is None:
        samples = list_samples(base)
    m = sample_map(samples)
    rows = []
    for s in samples:
        d = pd.read_csv(f"{base}/{s}.tsv", sep="\t",
                        usecols=["ref_len", "alt_len", "alt_freq", "n_called"])
        sv = (d.ref_len != 1) | (d.alt_len != 1)
        fin = d.alt_freq.notna()
        enf = eff_n_founders(s, base)
        rows.append(dict(
            sample=s,
            n_records=len(d), n_sv=int(sv.sum()), n_snp=int((~sv).sum()),
            finite_af_frac=float(fin.mean()),
            af_mean=float(d.alt_freq.mean()),
            sv_af_mean=float(d.loc[sv, "alt_freq"].mean()),
            median_n_called=float(d.n_called.median()),
            frac_ncalled_lt50=float((d.n_called < 50).mean()),
            eff_n_founders_mean=float(np.nanmean(list(enf.values()))) if enf else np.nan,
        ))
    qc = pd.DataFrame(rows).set_index("sample")
    qc = qc.join(m[["site", "plot", "generation", "coverage", "climate"]])
    if out_csv is None:
        out_csv = f"{GEA}/pilot_qc.csv"
    os.makedirs(GEA, exist_ok=True)
    qc.to_csv(out_csv)
    return qc


def _af_array(path: str) -> np.ndarray:
    """Just the alt_freq column as a float array (records are panel-ordered)."""
    return pd.read_csv(path, sep="\t", usecols=["alt_freq"]).alt_freq.to_numpy(dtype=float)


def build_group_means(samples: list[str] | None = None, base: str = OUT,
                      cache: bool = True) -> pd.DataFrame:
    """Per-record mean alt_freq for each (climate, generation) group + founding p0.

    All TSVs (evolved + SEEDMIX) share the SAME record order, so we aggregate by
    POSITION (fast numpy), not string keys. Returns record meta
    (chrom,pos,ref_len,alt_len) + p0 (mean over SEEDMIX reps) + one column per
    group `<climate>_g<gen>` (NaN-aware mean) + `<...>_n`. This is the Δp
    foundation: Δp_group = group_col - p0. Cached to analysis/grenenet_selection/.
    """
    # numpy .npz cache: readable across pandas versions (the pickle cache broke
    # when written by pandas 3.x and read by 2.x — StringDtype pickle mismatch).
    cache_path = GROUP_MEANS        # the panel version's (record order = its af_store / kMate TSVs)
    if cache and os.path.exists(cache_path):
        z = np.load(cache_path, allow_pickle=False)
        names = [str(c) for c in z["_columns"]]
        df = pd.DataFrame({c: z[c] for c in names})
        df["chrom"] = df["chrom"].astype(str)
        return df
    if samples is None:
        samples = list_samples(base)
    m = sample_map(samples)
    # group_means is the hot/cold 2-site pilot aggregation (the Δp foundation): only
    # climate-classified samples (SITE_CLIMATE: hot/cold sites) contribute. Every other
    # site has climate=NaN and would KeyError the (climate,gen) accumulator below — and
    # groupby() silently drops NaN groups, so they were never meant to be aggregated.
    # (Old cache built clean only because list_samples then returned the pilot alone;
    # it now returns all 2168.) Drop NaN-climate samples up front.
    samples = [s for s in samples if s in m.index and pd.notna(m.loc[s, "climate"])]
    if not samples:
        raise ValueError("group_means: no climate-classified (hot/cold) samples found")
    m = m.loc[samples]
    meta = pd.read_csv(f"{base}/{samples[0]}.tsv", sep="\t",
                       usecols=["chrom", "pos", "ref_len", "alt_len"])
    N = len(meta)

    # founding p0 = NaN-aware mean over the 8 SEEDMIX reps
    sm = [f for f in glob.glob(f"{SEEDMIX}/SEEDMIX_S*.tsv")
          if "_Chr" not in os.path.basename(f)]
    s_sum = np.zeros(N); s_cnt = np.zeros(N)
    for f in sm:
        a = _af_array(f)
        # GUARD: by-position aggregation REQUIRES every TSV on the SAME segregating
        # panel. The stale 06-02 cache was built on the old 10.33M panel (a superset
        # carrying AC=0 monomorphic records) — a record-count mismatch here means an
        # old-panel TSV slipped in and would silently misalign the sum. Fail loudly.
        if len(a) != N:
            raise ValueError(
                f"panel mismatch: {f} has {len(a):,} records, expected {N:,} "
                f"(all TSVs must share the segregating panel order)")
        ok = np.isfinite(a)
        s_sum[ok] += a[ok]; s_cnt += ok
    out = meta.copy()
    out["p0"] = s_sum / np.where(s_cnt > 0, s_cnt, np.nan)

    # one accumulator per (climate, generation); each sample read once
    keys = list(m.groupby(["climate", "generation"]).groups.keys())
    acc = {k: [np.zeros(N), np.zeros(N), 0] for k in keys}
    for s in samples:
        k = (m.loc[s, "climate"], m.loc[s, "generation"])
        a = _af_array(f"{base}/{s}.tsv")
        if len(a) != N:                                    # same panel-alignment guard
            raise ValueError(
                f"panel mismatch: {s} has {len(a):,} records, expected {N:,} "
                f"(segregating panel); refusing to misalign by-position aggregation")
        ok = np.isfinite(a)
        acc[k][0][ok] += a[ok]; acc[k][1] += ok; acc[k][2] += 1
    for (clim, gen), (gsum, gcnt, n) in acc.items():
        out[f"{clim}_g{gen}"] = gsum / np.where(gcnt > 0, gcnt, np.nan)
        out[f"{clim}_g{gen}_n"] = n

    # NOTE on monomorphic sites: the AC=0 panel scaffolding that contaminated the old
    # cache came from the PRE-segregating 10.33M panel. It is excluded simply by reading
    # the segregating-panel TSVs (the input record set has no AC=0 records) — guaranteed
    # upstream, not here. We deliberately do NOT drop "all-absent" rows: on the segregating
    # panel those are founder-segregating variants merely unobserved in the 175-sample
    # hot/cold pilot (p0=0, dp=0) — legitimate panel records, kept so group_means stays
    # one-row-per-segregating-record and any key-join resolves. The alignment assertions
    # above are the real guard against a stale/mixed-panel input.

    if cache:
        os.makedirs(GEA, exist_ok=True)
        arrs = {}
        for c in out.columns:
            a = out[c].to_numpy()
            arrs[c] = a.astype("U10") if a.dtype == object else a
        np.savez(cache_path, _columns=np.array(list(out.columns)), **arrs)
    return out


def load_genes() -> pd.DataFrame:
    """TAIR10 protein-coding genes: chrom, start, end, gene (AT-ID), name.

    chrom is 'Chr1'..'Chr5' (matches the SV records). ChrC/ChrM dropped.
    """
    g = pd.read_csv(TAIR10_GENES, sep="\t", header=None,
                    names=["chrom", "src", "feat", "start", "end",
                           "score", "strand", "frame", "attr"])
    g = g[g.chrom.isin([f"Chr{i}" for i in range(1, 6)])].copy()
    g["gene"] = g.attr.str.extract(r"ID=([^;]+)")
    g["name"] = g.attr.str.extract(r"Name=([^;]+)")
    return g[["chrom", "start", "end", "strand", "gene", "name"]].reset_index(drop=True)


# ------------------------------------------------------------------ genomic location: the project's ONE rule
# Every analysis that says where a variant sits (burden, SFS, catalogue, GEA screen, candidate list) calls `locate`.
# Location: CDS > UTR (incl. exons of non-coding transcripts) > intron > reference TE > promoter (<= PROM_BP upstream,
# strand-aware) > intergenic; a record with REF or ALT longer than SV_MIN_BP (every SV, |ALT - REF| > SV_MIN_BP as
# `is_sv`, and long replacements) by its two flanking bases (`flanks`), anything else by its REF span. Gene models: the
# TAIR10 GFF without the TE remnants TAIR10 annotates as protein-coding genes (`tair10_te_derived_genes.csv`, `te_gene`)
# and without TAIR10's transposable_element_gene transcripts, for the location and for the gene named alike.
# PROM_BP: 1 kb, not 2 kb: the median Arabidopsis gene sits only 863 bp from its nearest upstream neighbour, so a 2 kb
# window lies inside another gene's body for 73.4% of genes against 54.2% at 1 kb (also TAIR's standard upstream set).
PROM_BP = 1000
LOCATIONS = ["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]
TE_DERIVED_GENES = f"{SV_ADAPTIVE_SHARED}/tair10_te_derived_genes.csv"   # _te_derived_genes.py (panel-independent)
_LAYERS = {}


def merged_intervals(df: pd.DataFrame) -> dict:
    """Non-overlapping sorted intervals per chrom: {chrom: (starts, ends)}."""
    out = {}
    for c, d in df.groupby("chrom"):
        d = d.sort_values("start"); s, e = [], []
        for a, b in zip(d.start.to_numpy(), d.end.to_numpy()):
            if s and a <= e[-1]: e[-1] = max(e[-1], b)
            else: s.append(a); e.append(b)
        out[c] = (np.array(s), np.array(e))
    return out


def hit_intervals(iv: dict, chrom, lo, hi) -> np.ndarray:
    """bool per query: does [lo, hi] overlap any merged interval of iv on its chrom."""
    chrom, lo, hi = np.asarray(chrom, dtype=object), np.asarray(lo, np.int64), np.asarray(hi, np.int64)
    res = np.zeros(lo.size, bool)
    for c in np.unique(chrom):
        if c not in iv: continue
        m = chrom == c; s, e = iv[c]
        j = np.searchsorted(s, hi[m], side="right") - 1
        ok = j >= 0; r = np.zeros(m.sum(), bool)
        r[ok] = e[j[ok]] >= lo[m][ok]
        res[m] = r
    return res


def gene_models() -> pd.DataFrame:
    """The TAIR10 GFF (genes + TEs, Chr1-5) without every row of the TE-derived "genes": the gene models of `locate`."""
    g = pd.read_csv(TAIR10_GENES_TE, sep="\t", header=None, comment="#",
                    names=["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"])
    g = g[g.chrom.isin(CHROMS)]
    tes = set(pd.read_csv(TE_DERIVED_GENES).query("te_gene").gene)
    g = g[~g.attr.str.extract(r"(?:ID|Parent)=([^;,]+)")[0].str.replace(r"\.\d+$", "", regex=True).isin(tes)].copy()
    g["gene"] = np.where(g.feat == "gene", g.attr.str.extract(r"ID=([^;]+)")[0], "")
    return g.reset_index(drop=True)


def _location_tables():
    if not _LAYERS:                                # read once per process
        g = gene_models()
        genes = g[g.feat == "gene"]
        prom = genes.copy()
        prom["start"] = np.where(genes.strand == "+", genes.start - PROM_BP, genes.end + 1)
        prom["end"] = np.where(genes.strand == "+", genes.start - 1, genes.end + PROM_BP)
        # transcript features count only for the rule's genes (not TAIR10's transposable_element_gene transcripts), and
        # the exons of a non-coding transcript (ncRNA, miRNA, tRNA, ... : no CDS) count with the UTR: transcribed, not coding
        tx = g.attr.str.extract(r"Parent=([^;,]+)")[0]
        host = tx.str.replace(r"\.\d+$", "", regex=True).isin(set(genes.gene))
        exon = g[(g.feat == "exon") & host]
        noncoding = exon[~tx[exon.index].isin(set(tx[g.feat == "CDS"]))]
        _LAYERS["IV"] = {"CDS": merged_intervals(g[(g.feat == "CDS") & host]),
                         "UTR": merged_intervals(pd.concat([g[g.feat.isin(["five_prime_UTR", "three_prime_UTR"]) & host], noncoding])),
                         "exon": merged_intervals(exon), "gene": merged_intervals(genes),
                         "TE": merged_intervals(g[g.feat == "transposable_element"]), "promoter": merged_intervals(prom)}
        _LAYERS["genes"] = genes[["chrom", "start", "end", "strand", "gene"]].reset_index(drop=True)
    return _LAYERS


def location_layers() -> dict:
    """Merged intervals of every layer of the rule: {CDS, UTR, exon, gene, TE, promoter: {chrom: (starts, ends)}}."""
    return _location_tables()["IV"]


def location_genes() -> pd.DataFrame:
    """The genes of the rule (chrom, start, end, strand, gene): the `gene` features of `gene_models`."""
    return _location_tables()["genes"]


def location_context(IV: dict, chrom, lo, hi, promoter_first: bool = False) -> np.ndarray:
    """Index into LOCATIONS of [lo, hi], highest-priority layer it touches. Gene body outranks TE and promoter; TE
    outranks promoter (26% of events in a TE also sit in a promoter window); promoter_first flips that one choice."""
    H = {k: hit_intervals(v, chrom, lo, hi) for k, v in IV.items()}
    cat = np.full(np.asarray(lo).size, 5, np.int8)
    if promoter_first:
        cat[H["TE"]] = 3; cat[H["promoter"]] = 4
    else:
        cat[H["promoter"]] = 4; cat[H["TE"]] = 3
    cat[H["gene"] & ~H["exon"]] = 2; cat[H["UTR"]] = 1; cat[H["CDS"]] = 0
    return cat


def flanks(pos, ref, alt):
    """The two reference bases next to an SV's variable segment, which both alleles share: REF and ALT are trimmed of the
    bases they share at their start (the VCF anchor) and then at their end, and the flanks are the reference bases just
    before and just after what is left of REF. A left-anchored record (REF and ALT start with the same base) has its left
    flank at pos; a right-anchored one (they end with the same base) its right flank at pos + ref_len - 1; a record
    sharing neither has them at pos - 1 and pos + ref_len."""
    pos = np.asarray(pos, np.int64); L, R = np.empty(pos.size, np.int64), np.empty(pos.size, np.int64)
    for i, (r, a) in enumerate(zip(ref, alt)):
        n = min(len(r), len(a)); k = 0
        while k < n and r[k] == a[k]: k += 1
        s = 0
        while s < n - k and r[-1 - s] == a[-1 - s]: s += 1
        L[i], R[i] = pos[i] + k - 1, pos[i] + len(r) - s
    return L, R


def _genes_at(chrom, x):
    """Per query point: the gene whose body contains it (the last-starting one), and every gene whose promoter window
    contains it, nearest TSS first (';'-joined)."""
    G = location_genes(); chrom = np.asarray(chrom, dtype=object); x = np.asarray(x, np.int64)
    body = np.full(x.size, "", dtype=object); prom = np.full(x.size, "", dtype=object)
    for c in np.unique(chrom):
        g = G[G.chrom == c].sort_values("start", kind="stable")
        gs, ge, gid = g.start.to_numpy(), g.end.to_numpy(), g.gene.to_numpy()
        plus = g.strand.to_numpy() == "+"
        ws = np.where(plus, gs - PROM_BP, ge + 1); tss = np.where(plus, gs, ge)    # promoter windows, all PROM_BP long
        o = np.argsort(ws, kind="stable"); ws, wid, wtss = ws[o], gid[o], tss[o]
        for i in np.flatnonzero(chrom == c):
            j = np.searchsorted(gs, x[i], side="right")
            cand = np.flatnonzero(ge[max(0, j - 50):j] >= x[i])
            if cand.size:
                body[i] = gid[max(0, j - 50) + cand[-1]]
            a, b = np.searchsorted(ws, x[i] - PROM_BP + 1), np.searchsorted(ws, x[i], side="right")
            if b > a:
                prom[i] = ";".join(wid[a:b][np.argsort(np.abs(wtss[a:b] - x[i]), kind="stable")])
    return body, prom


def located_by_flanks(ref_len, alt_len) -> np.ndarray:
    """bool per record: located by its two flanking bases, because REF or ALT is longer than SV_MIN_BP (an SV, |ALT - REF|
    > SV_MIN_BP, or a long replacement whose ALT is about as long as REF); the SV / short-indel CLASS stays |ALT - REF|."""
    return np.maximum(np.asarray(ref_len, np.int64), np.asarray(alt_len, np.int64)) > SV_MIN_BP


def _located(chrom, pos, ref_len, alt_len, ref, alt, flank_left, flank_right, promoter_first, IV):
    """(location index, lo, hi, SV rows, left- and right-flank priority rank) of `locate`."""
    ch = np.asarray(chrom, dtype=object); pos = np.asarray(pos, np.int64)
    rl = np.asarray(ref_len, np.int64); al = np.asarray(alt_len, np.int64)
    IV = location_layers() if IV is None else IV
    lo, hi = pos.copy(), pos + np.maximum(rl - 1, 0)
    sv = np.flatnonzero(located_by_flanks(rl, al))
    loc = location_context(IV, ch, lo, hi, promoter_first)
    lc = rc = np.zeros(0, np.int8)
    if sv.size:
        if flank_left is not None:
            L, R = np.asarray(flank_left)[sv].astype(np.int64), np.asarray(flank_right)[sv].astype(np.int64)
        elif ref is not None:
            L, R = flanks(pos[sv], [ref[i] for i in sv], [alt[i] for i in sv])
        else:
            raise ValueError(f"{sv.size} records with REF or ALT > {SV_MIN_BP} bp need flank_left/flank_right or ref/alt to be located")
        lo[sv], hi[sv] = L, R
        lc = location_context(IV, ch[sv], L, L, promoter_first); rc = location_context(IV, ch[sv], R, R, promoter_first)
        # each flank alone, not the span between them (that is the segment); the flank of higher PRIORITY wins, which is
        # the lower index only while TE outranks promoter
        rank = np.array([0, 1, 2, 4, 3, 5]) if promoter_first else np.arange(6)
        loc[sv] = np.where(rank[lc] <= rank[rc], lc, rc)
        lc, rc = rank[lc], rank[rc]                # `locate` compares the flanks by priority
    return loc, lo, hi, sv, lc, rc


def location_index(chrom, pos, ref_len, alt_len, ref=None, alt=None, flank_left=None, flank_right=None,
                   promoter_first: bool = False, IV: dict | None = None) -> np.ndarray:
    """`locate`'s location only, as an index into LOCATIONS (fast: no gene lookup). IV: other layers (e.g. shifted)."""
    return _located(chrom, pos, ref_len, alt_len, ref, alt, flank_left, flank_right, promoter_first, IV)[0]


def locate(chrom, pos, ref_len, alt_len, ref=None, alt=None, flank_left=None, flank_right=None,
           promoter_first: bool = False) -> pd.DataFrame:
    """Where each record sits, by the project's one rule (see above). Returns, in input order:
      location    LOCATIONS entry
      gene        the gene the location refers to: for promoter the gene whose promoter window holds the record (nearest
                  TSS); otherwise the gene whose body holds it, else that promoter gene; '' when neither. For an SV,
                  the flank that gave the location is read first
      prom_genes  every gene whose promoter window holds it, nearest TSS first (divergent promoters share upstream DNA)
      lo, hi      the positions located: the two flanks of an SV (each a point), the REF span of anything else
    A record located by its flanks (`located_by_flanks`) needs them: pass flank_left / flank_right (e.g. from
    `sv_event_catalog.csv`) or its ref / alt strings.
    """
    loc, lo, hi, sv, lc, rc = _located(chrom, pos, ref_len, alt_len, ref, alt, flank_left, flank_right, promoter_first, None)
    ch = np.asarray(chrom, dtype=object)
    gb = np.full(lo.size, "", dtype=object); gp = np.full(lo.size, "", dtype=object)
    nsv = np.setdiff1d(np.arange(lo.size), sv)
    if nsv.size:                                   # REF span (<= SV_MIN_BP bp, so no gene fits inside): its start, else its end
        for x in (lo[nsv], hi[nsv]):
            b, p = _genes_at(ch[nsv], x)
            gb[nsv] = np.where(gb[nsv] == "", b, gb[nsv]); gp[nsv] = np.where(gp[nsv] == "", p, gp[nsv])
    if sv.size:
        (bl, pl), (br, pr) = _genes_at(ch[sv], lo[sv]), _genes_at(ch[sv], hi[sv])
        first = lc <= rc                           # the gene of the flank that gave the location, else of the other one
        gb[sv] = np.where(first, np.where(bl != "", bl, br), np.where(br != "", br, bl))
        gp[sv] = np.where(first, np.where(pl != "", pl, pr), np.where(pr != "", pr, pl))
    gp0 = np.array([p.split(";")[0] for p in gp], dtype=object)
    prom = loc == LOCATIONS.index("promoter")      # a promoter call names its promoter's gene, never a body gene
    gene = np.where(prom, gp0, np.where(gb != "", gb, gp0))   # (a body without exon/intron/UTR/CDS call: an ncRNA exon)
    return pd.DataFrame({"location": np.array(LOCATIONS, dtype=object)[loc], "gene": gene, "prom_genes": gp,
                         "lo": lo, "hi": hi})


def location_layers_touched(chrom, pos, ref_len, alt_len, ref=None, alt=None, flank_left=None, flank_right=None,
                            IV: dict | None = None) -> dict:
    """Every layer a record touches, not only the winner: {CDS, UTR, intron, TE, promoter: bool array}, by the same
    positions as `locate` (an SV: either flank)."""
    IV = location_layers() if IV is None else IV
    _, lo, hi, sv, _, _ = _located(chrom, pos, ref_len, alt_len, ref, alt, flank_left, flank_right, False, IV)
    ch = np.asarray(chrom, dtype=object)
    def hits(c, a, b):
        H = {k: hit_intervals(v, c, a, b) for k, v in IV.items()}
        return {"CDS": H["CDS"], "UTR": H["UTR"], "intron": H["gene"] & ~H["exon"], "TE": H["TE"], "promoter": H["promoter"]}
    out = hits(ch, lo, hi)
    if sv.size:
        hl, hr = hits(ch[sv], lo[sv], lo[sv]), hits(ch[sv], hi[sv], hi[sv])
        for k in out:
            out[k][sv] = hl[k] | hr[k]
    return out


EXPR_DIR = f"{PROJ}/data/eqtl"
_EXPRESSION = {}                                            # cache of `expression()`s matrix


def expression(genes=None, accessions=None) -> pd.DataFrame:
    """Gene expression in the 1001 Genomes rosettes: the project's ONE expression rule, used by every analysis.

    CPM x the lab's gene-length factor x the lab's batch factor (SHARE/TUTORIALS/example_load_1001expression.R uses
    the same factors in d_log2_batch, batch MU dropped), with nothing set to NA: 0 CPM stays 0 (a 0 is non-expression,
    not missing data) and values below 1 CPM are kept. Unlike d_log2_batch, each gene keeps its own level (d_log2_batch
    is centred on the gene's mean); within a gene, log2(value) = d_log2_batch + a constant.
    Replicate libraries of an accession are averaged on the CPM scale. Made by `r1_sv_negative_selection/
    _export_1001t_cpm.R`. Returns genes x accessions (accession ids as str) in corrected CPM; for means or display take
    log2(CPM + 1), for tests between accessions use the CPM (ranks). genes / accessions: optional subsets (missing ones
    are left out).
    """
    if "X" not in _EXPRESSION:                                   # read the matrix once per process
        S = pd.read_csv(f"{EXPR_DIR}/TG_samples.csv")
        acc = pd.Series(S.accession.astype(str).to_numpy(), index=S["sample"])
        M = pd.read_csv(f"{EXPR_DIR}/TG_genes_cpm_noMU.tsv.gz", sep="\t", index_col=0)
        _EXPRESSION["X"] = M.T.groupby(acc.reindex(M.columns).to_numpy()).mean().T
    X = _EXPRESSION["X"]
    if genes is not None:
        X = X.loc[[g for g in dict.fromkeys(genes) if g in X.index]]
    X = X.copy()
    if accessions is not None:
        X = X[[a for a in dict.fromkeys(map(str, accessions)) if a in X.columns]]
    return X


def annotate_svs(df: pd.DataFrame, flank: int = 0, genes: pd.DataFrame | None = None
                 ) -> pd.DataFrame:
    """Tag each SV row (needs `chrom`,`pos`) with the gene(s) it overlaps.

    SV span = [pos, pos + ref_len) when ref_len present, else the point `pos`.
    `flank` bp widens the SV span on each side (e.g. 2000 for nearby/promoter
    hits). Adds columns: gene (first hit AT-ID), gene_name, n_genes (overlaps),
    genes_all (';'-joined). Vectorized per chrom via sorted-interval search.
    """
    if genes is None:
        genes = load_genes()
    out_gene = np.full(len(df), "", dtype=object)
    out_name = np.full(len(df), "", dtype=object)
    out_n = np.zeros(len(df), dtype=int)
    out_all = np.full(len(df), "", dtype=object)
    pos = df["pos"].to_numpy()
    rl = df["ref_len"].to_numpy() if "ref_len" in df else np.ones(len(df), int)
    lo = pos - flank
    hi = pos + np.where(rl > 1, rl, 1) + flank
    idx = df.index.to_numpy()
    for c, gc in genes.groupby("chrom"):
        gs = gc.start.to_numpy(); ge = gc.end.to_numpy()
        gid = gc.gene.to_numpy(); gnm = gc.name.fillna("").to_numpy()
        sel = np.where(df["chrom"].to_numpy() == c)[0]
        for j in sel:
            ov = np.where((gs <= hi[j]) & (ge >= lo[j]))[0]
            if ov.size:
                out_gene[j] = gid[ov[0]]; out_name[j] = gnm[ov[0]]
                out_n[j] = ov.size; out_all[j] = ";".join(gid[ov])
    res = df.copy()
    res["gene"] = out_gene; res["gene_name"] = out_name
    res["n_genes"] = out_n; res["genes_all"] = out_all
    return res


def assign_ld_blocks(chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Assign each record (chrom='Chr1'..'Chr5', pos) to its phase-1 hapFIRE LD
    block, by inheriting the block of the NEAREST genotyped SNP (per chrom).

    Returns a string array of block ids ('chrom_idx'); '' where the chrom has no
    SNP map. SNP density ~1/114 bp, so the nearest-SNP block is reliable.
    """
    b = pd.read_csv(LD_BLOCKS, usecols=["pos", "chrom", "block"])
    out = np.full(len(pos), "", dtype=object)
    pos = np.asarray(pos, dtype=np.int64)
    for ci, g in b.groupby("chrom"):
        m = chrom == f"Chr{ci}"
        if not m.any():
            continue
        g = g.sort_values("pos")
        sp = g.pos.to_numpy(); bl = g.block.to_numpy()
        p = pos[m]
        j = np.clip(np.searchsorted(sp, p), 0, len(sp) - 1)
        jm1 = np.clip(j - 1, 0, len(sp) - 1)
        # pick whichever flanking SNP is closer in bp
        use_prev = (j > 0) & (np.abs(p - sp[jm1]) <= np.abs(sp[j] - p))
        nearest = np.where(use_prev, jm1, j)
        out[np.where(m)[0]] = bl[nearest]
    return out


def assign_clq_blocks(chrom: np.ndarray, pos: np.ndarray, r2: float = 0.9) -> np.ndarray:
    """Assign each record (chrom='Chr1'..'Chr5', pos) to its BigLD clq{r2} block.

    Reads the per-chrom interval TSVs in CLQ_BLOCKS_DIR (chrom, start_pos,
    end_pos, n_variants). Unlike the hapFIRE nearest-SNP map, these are LD islands
    that do NOT tile the genome, so a variant landing in an inter-block gap gets ''
    (drop it upstream of WZA). Block ids are 'Chr{n}_{idx}', idx = the block's rank
    by start_pos within its chrom (stable, unique). Returns a string object array.
    """
    out = np.full(len(pos), "", dtype=object)
    pos = np.asarray(pos, dtype=np.int64)
    chrom = np.asarray(chrom, dtype=str)
    tag = f"clq{r2}"
    for ci in range(1, 6):
        c = f"Chr{ci}"
        f = f"{CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv"
        if not os.path.exists(f):
            continue
        g = pd.read_csv(f, sep="\t").sort_values("start_pos").reset_index(drop=True)
        st = g["start_pos"].to_numpy(np.int64)
        en = g["end_pos"].to_numpy(np.int64)
        m = np.where(chrom == c)[0]
        if len(m) == 0:
            continue
        p = pos[m]
        j = np.searchsorted(st, p, "right") - 1          # candidate block (last start <= p)
        inb = (j >= 0) & (p <= en[np.clip(j, 0, len(en) - 1)])
        idx = np.where(inb)[0]
        out[m[idx]] = [f"{c}_{int(j[k])}" for k in idx]
    return out


def collapse_to_blocks(df: pd.DataFrame, stat: str = "z_emp",
                       block_col: str = "block") -> pd.DataFrame:
    """Collapse a per-SV GEA frame to one LEAD SV per LD block (max |stat|).

    Adds n_sv_block (SVs in the block). Use for candidate dedup / Manhattan so a
    single low-recomb block (up to ~900 correlated SVs) can't dominate the list.
    NOTE: lead-by-|stat| is selection-biased — for unbiased GIF/QQ use a
    stat-independent representative instead (e.g. first by pos, or .sample()).
    """
    d = df.copy()
    d["_abs"] = d[stat].abs()
    d["n_sv_block"] = d.groupby(block_col)[stat].transform("size")
    lead = d.sort_values("_abs", ascending=False).groupby(block_col, as_index=False).first()
    return lead.drop(columns="_abs")


def bh(pv: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    m = len(pv); o = pv.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((pv[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def lamgc(pv: np.ndarray) -> float:
    """Genomic-control lambda from a p-value vector (median chi2_1 / expected median)."""
    return float(np.median(stats.chi2.isf(np.clip(pv, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1))


def multisite_gwas_raw(tag: str = "clq90_pc1") -> dict:
    """Load the multisite founder-GWAS Z/C (+recomputed JOINT/GLOBAL/CLIMATE) for further
    C-whitened contrasts (e.g. climate-cluster, locality-index) without a new GWAS run.

    Recomputes chi2_joint/p_joint/z_global/z_clim from the raw Z (M markers x S sites,
    genomic-controlled) and cross-site null covariance C in
    analysis/grenenet_selection/archive/window_hapfreq_retired/hapfreq/multisite_founder_gwas_{tag}.npz, using the identical
    formulas as founder_gwas_multisite.py. Self-checks the recomputed JOINT lambda_GC
    against the saved _meta.json so a silent Z/C misread can't drift unnoticed (the npz
    doesn't store the derived per-marker stats, and the sibling .csv is p_joint-sorted so
    it can't be joined back to Z positionally).
    """
    h = f"{GEA}/hapfreq"
    z = np.load(f"{h}/multisite_founder_gwas_{tag}.npz", allow_pickle=True)
    meta = json.load(open(f"{h}/multisite_founder_gwas_{tag}_meta.json"))
    Z, C, sites, bio1 = z["Z"], z["C"], z["sites"], z["bio1"]
    chrom, start, end, mac = z["chrom"], z["start"], z["end"], z["mac"]
    M, S = Z.shape
    Cinv = np.linalg.pinv(C)
    one = np.ones(S)
    dg = float(one @ Cinv @ one)
    z_global = (Z @ Cinv @ one) / np.sqrt(dg)
    c0 = (bio1 - bio1.mean()) / bio1.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one
    z_clim = (Z @ Cinv @ c) / np.sqrt(float(c @ Cinv @ c))
    chi2_joint = np.einsum("mi,ij,mj->m", Z, Cinv, Z)
    p_joint = stats.chi2.sf(chi2_joint, S)
    lam = lamgc(p_joint)
    assert abs(lam - meta["tests"]["JOINT"]["lam"]) < 1e-3, (
        f"recomputed JOINT lambda {lam:.4f} != saved meta {meta['tests']['JOINT']['lam']:.4f} "
        "-- Z/C recompute drifted from the original GWAS run")
    return dict(Z=Z, C=C, Cinv=Cinv, sites=sites, bio1=bio1, one=one, dg=dg,
                unit=np.array([f"{c_}:{s}-{e}" for c_, s, e in zip(chrom, start, end)]),
                chrom=chrom, mac=mac, chi2_joint=chi2_joint, p_joint=p_joint, q_joint=bh(p_joint),
                z_global=z_global, z_clim=z_clim, M=M, S=S, meta=meta)


# founder MAC count matching build_sv_landscape.py's MAF>=0.05 floor (~5.2% of 231
# founders) -- the bar that already defines the LD blocks and the SV landscape. Default for
# founder_panel_keep below; independent of maf_filter's own min_maf (used by the
# climate-cluster / locality-index marker ranking, currently 1% per project decision).
COMMON_MAC = 12


def maf_filter(raw: dict, min_maf: float = 0.01) -> dict:
    """Restrict a lib.multisite_gwas_raw() dict to markers with founder MAF >= min_maf.

    Rare-founder-allele markers have noisier kinship-corrected z-scores, so any NEW
    per-marker contrast (not already lambda-validated the way JOINT/GLOBAL/CLIMATE were)
    should be ranked only among markers clearing a standard MAF floor (default 1%, the
    common GWAS QC convention). min_maf is converted to a founder MAC via the cohort's
    n_founders (in raw["meta"]) so it's independent of the founder count. Per-site fields
    (S-length: sites/bio1/one/C/Cinv) are untouched; per-marker fields (M-length) are
    masked and M is updated.
    """
    min_mac = int(np.ceil(min_maf * raw["meta"]["n_founders"]))
    keep = raw["mac"] >= min_mac
    out = dict(raw)
    for k in ("Z", "unit", "chrom", "mac", "chi2_joint", "p_joint", "q_joint", "z_global", "z_clim"):
        out[k] = raw[k][keep]
    out["M"] = int(keep.sum())
    return out


def founder_panel_keep(chrom: np.ndarray, pos: np.ndarray, min_mac: int = COMMON_MAC,
                       called_min: float = 0.9) -> np.ndarray:
    """Founder-panel MAC/call-rate mask for arbitrary (chrom,pos) records.

    Same floor as COMMON_MAC / build_sv_landscape.py (MAF>=0.05, called-frac>=0.9 among
    the 231 founders) -- the bar that already defines the LD blocks and the SV landscape.
    Use to floor any OTHER per-variant catalog (e.g. the raw per-sample/pool kMate variant
    calls behind site_variant_temporal_scoef.py) before counting SV/indel/SNP composition
    or scoring selection on it: unlike the founder panel, that catalog carries no MAF floor
    of its own (only a noisy pool-level p0 reachability check), so >50% of its SV-class
    calls are literal founder singletons -- see KINSHIP_TEMPORAL_METHODS_RESEARCH context.
    Positions absent from the founder panel entirely do not pass (conservative: no MAC to
    vouch for them). Returns a bool array aligned to (chrom,pos), True where at least one
    panel record at that exact position passes.

    WARNING -- POSITION-LEVEL, NOT RECORD-LEVEL (bug found 2026-09-18). At a multi-allelic
    position every record passes if ANY record there clears the floor, so rare alleles ride in
    on a common sibling's frequency. Example: Chr2:7158959 holds 15 alleles of one ~2 kb
    insertion; one has 21 carriers, so an allele with 3 carriers passed a MAC>=12 filter.
    This mixes units (a record is admitted on its bubble's frequency but then analysed as an
    allele) and was concentrated in lineage-restricted SVs: 78% of "cold-specific" insertions
    in the MAC>=12 temporal set had < 12 carriers of their own. New code should use
    `founder_record_keep()`. Kept unchanged so existing results stay reproducible.
    """
    import warnings
    warnings.warn("founder_panel_keep() is position-level: at multi-allelic sites rare alleles pass "
                  "on a sibling's frequency. Use founder_record_keep() for a per-record filter.",
                  stacklevel=2)
    import scipy.sparse as sp
    out = np.zeros(len(pos), bool)
    chrom = np.asarray(chrom); pos = np.asarray(pos, dtype=np.int64)
    for ch in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        m = chrom == ch
        if not m.any():
            continue
        base = var_pa_prefix(ch)
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        lpos = meta["pos"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz"); vc = sp.load_npz(f"{base}.var_called.npz")
        n_alt = np.asarray(vp.sum(0)).ravel(); n_cal = np.asarray(vc.sum(0)).ravel()
        passmaf = (n_alt >= min_mac) & (n_alt <= vp.shape[0] - min_mac) & (n_cal >= called_min * vp.shape[0])
        order = np.argsort(lpos); lpos_s = lpos[order]; pass_s = passmaf[order]
        p = pos[m]
        j = np.searchsorted(lpos_s, p)
        found = (j < len(lpos_s)) & (lpos_s[np.clip(j, 0, len(lpos_s) - 1)] == p)
        res = np.zeros(len(p), bool)
        res[found] = pass_s[np.clip(j[found], 0, len(lpos_s) - 1)]
        out[m] = res
    return out



def founder_record_keep(kind: str, rows: np.ndarray, min_mac: int = COMMON_MAC,
                        called_min: float = 0.9, return_counts: bool = False):
    """Founder-panel MAC / call-rate filter PER RECORD (allele), for af_store rows.

    The record-level replacement for `founder_panel_keep()`. Each record is judged on its OWN
    carrier count and call rate, so a rare allele cannot pass on a common sibling at the same
    position.

    kind : "nonsnp" or "snp" -- which af_store index `rows` refer to.
    rows : int array of af_store row indices (into index_{kind}.npz).
    Relies on af_store's per-chromosome rows being the panel version's records of that kind in the
    same order (asserted here on ref/alt lengths for every chromosome).
    Pass: min_mac <= n_carriers <= N - min_mac  and  n_called >= called_min * N  (N = 231).
    Returns bool mask aligned to `rows`; with return_counts=True also (n_carriers, n_called).
    """
    import scipy.sparse as sp
    rows = np.asarray(rows, dtype=np.int64)
    idx = np.load(f"{AF_STORE}/index_{kind}.npz")
    ch = idx["chrom"].astype("U5")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    n_alt = np.zeros(rows.size, np.int64); n_cal = np.zeros(rows.size, np.int64)
    off = 0; N = None
    for c in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        n = int((ch == c).sum())
        base = var_pa_prefix(c)
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        issnp = (prl == 1) & (pal == 1)
        panel = np.where(issnp if kind == "snp" else ~issnp)[0]
        assert panel.size == n, f"{c}/{kind}: af_store {n} rows vs panel {panel.size} records"
        assert np.array_equal(prl[panel], rl[off:off + n]) and np.array_equal(pal[panel], al[off:off + n]), \
            f"{c}/{kind}: af_store/panel record order differs"
        sel = (rows >= off) & (rows < off + n)
        if sel.any():
            vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
            N = vp.shape[0]
            cols = panel[rows[sel] - off]
            n_alt[sel] = np.asarray((vp[:, cols] > 0).sum(0)).ravel()
            n_cal[sel] = np.asarray((vc[:, cols] > 0).sum(0)).ravel()
        off += n
    N = N or 231
    keep = (n_alt >= min_mac) & (n_alt <= N - min_mac) & (n_cal >= called_min * N)
    return (keep, n_alt, n_cal) if return_counts else keep

def matched_perm_test(values: np.ndarray, bins: np.ndarray, idx: np.ndarray,
                      nperm: int = 10000, seed: int = 0) -> tuple[float, float, float, float]:
    """Size-matched permutation test (the design used throughout sv_adaptive/): is
    mean(values[idx]) higher than expected from len(idx) draws matched on `bins` (e.g.
    block-size bin)?

    Each idx item's null draws come from its OWN bin, `nperm` times, so bigger/smaller
    items are compared to same-size controls. Returns (observed, null_mean, obs/null,
    one-sided p = P(null_mean_over_draws >= observed)).
    """
    rng = np.random.default_rng(seed)
    bin_members = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    obs = values[idx].mean()
    null = np.empty((len(idx), nperm))
    for k, i in enumerate(idx):
        null[k] = rng.choice(values[bin_members[bins[i]]], nperm)
    nm = null.mean(0)
    p = (1 + (nm >= obs).sum()) / (nperm + 1)
    ratio = obs / nm.mean() if nm.mean() != 0 else np.nan
    return float(obs), float(nm.mean()), float(ratio), float(p)


def eff_n_founders(sample: str, base: str = OUT) -> dict:
    """Per-chrom effective #founders = 1/sum(h^2) from the saved h vectors.

    The per-chrom driver writes one file per chrom: <sample>_Chr{N}.h_per_chrom.npz
    (global mode), each holding that chrom's h under key 'Chr{N}'.
    """
    out = {}
    for p in sorted(glob.glob(f"{base}/{sample}_Chr*.h_per_chrom.npz")):
        z = np.load(p, allow_pickle=True)
        for k in z.files:
            if k.startswith("Chr"):
                h = z[k]
                out[k] = float(1.0 / np.sum(h ** 2)) if h.size else np.nan
    return out
