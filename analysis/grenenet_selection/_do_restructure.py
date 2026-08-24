#!/usr/bin/env python
"""GEA restructure migration. Disposable — delete after the move lands.

  python _do_restructure.py            # dry run
  python _do_restructure.py --apply    # git mv + patch + verify
  python _do_restructure.py --verify   # re-run the import verifier only

Handles the thing that makes this non-trivial: every script reaches the shared
`lib.py` via a chain of os.path.dirname() around __file__. Moving a script one
level down (whether on its own, or because its whole directory moved) makes
that chain resolve one level too shallow and `import lib` breaks at runtime,
not at move time. So each moved .py gets exactly one more dirname(), and the
verifier then evaluates every chain for real and checks it lands on lib.py.
"""
import os, re, subprocess, sys, collections

GEA = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(GEA))
APPLY = "--apply" in sys.argv
VERIFY_ONLY = "--verify" in sys.argv

MAP = {}

MAP["r1_sv_negative_selection"] = """
_compute_s_dist_by_stratum.py _compute_s_climate_slope.py
_temporal_selection_snp_vs_nonsnp.py _temporal_s_snp_vs_nonsnp.py
_temporal_s_plots_snp_vs_nonsnp.py _temporal_s_enrich_initqty.py
_temporal_sel_drift_maf.py _nonsnp_temporal_category.py
site_variant_temporal_scoef.py _sv_temporal_direct.py
_build_temporal_s_consolidated_nb.py _build_temporal_s_nofilter_nb.py
_compute_sfs_shift_by_site.py _build_sfs_shift_nb.py
_compute_sfs_time_site4.py _build_sfs_time_site4_nb.py
_compute_parallelism.py _picmin.py _site_parallelism.py
_plot_site_parallelism.py _build_parallelism_picmin_nb.py
site_sv_enrichment.py aggregate_sites_enrichment.py
_render_site4_enrichment_fig.py _enrich_threshold_sweep.py
_audit_s_classes.py noise_check_sv_s.py _sv_callqual_artifact.py
_founder_load_test.py _build_sv_selection_audit_nb.py
_build_sv_polarity_manhattan_nb.py
run_site_enrichment.sbatch _extract_vcf_callqual.sh
SV_TEMPORAL_PURGING_SUMMARY.md BAYPASS_TEMPORAL_PLAN.md
""".split()

MAP["r2_gea_nonsnp"] = """
build_kendall.py _build_kendall_nb.py _kendall_site_collapsed_test.py
build_lfmm_input.py _build_lfmm_nb.py build_lfmm_scree.py
_build_kselection_nb.py build_struct_pca.py
build_mixedmodel.py _build_mixedmodel_nb.py
build_two_stage_gea.py _build_twostage_nb.py build_two_stage_pooled.py
build_lmm_scoef.py
build_wza.py wza_script.py _build_wza_nb.py build_significant_blocks.py
nonsnp_only_genes.py nonsnp_only_genes_describe.py
_build_nonsnp_only_genes_nb.py go_enrichment_nonsnp.py genes_from_regions.py
check_lea_elip_persite.py _build_snp_vs_nonsnp_viz_nb.py
_build_persite_new_peaks_nb.py _build_persite_new_peaks_sv_nb.py
_build_candidate_nb.py plot_candidates.py plot_phase1style.py
plot_gi_locus.py plot_micropangenome.py
run_lfmm_kscan.sh run_lfmm_ksweep.sh push_gea_to_drive.sh
phase1_replication wza_investigation cam5_replication
""".split()

MAP["r3_persite_gwas"] = """
build_selection_trait.py founder_genotype.py extract_shortread_geno.py
build_class_grms.py build_untagged_grm.py
class_split_gwas.py _build_class_gwas_persite_nb.py
_build_class_gwas_multitrait_nb.py plot_class_gwas_pngs.py
ecotype_fitness.py ecotype_fitness_gwas.py ecotype_fitness_enrichment.py
founder_persite_sv_enrichment.py
varexp_selection.py varexp_untagged.py varexp_bioclim.py lasso_bioclim.py
build_nonsnp_tagging.py build_tagging_masked.py _build_tagging_masked_nb.py
compute_tagging_null.py verify_tagging_permutation_null.py
run_class_grms.sbatch run_class_split_gwas.sbatch run_ecotype_gwas.sbatch
lasso_bioclim.sbatch run_tagging_filter_sweep.sh tagging_sweep.sbatch
VAREXP_SELECTION_HANDOFF.md
""".split()

# not the main analysis: hit interrogation + the "are SVs drivers or
# passengers" investigation, at both the fine-mapping and haplotype units
MAP["extras"] = """
driver_passenger
ecotype_block_sv_enrichment.py ecotype_hap_sv_enrichment.py
ecotype_hap_sv_regional.py ecotype_hap_sv_rotation.py
ecotype_hap_sv_rotation2.py
_sv_hap_context.py _sv_hap_freqrobust.py _sv_haplotype_axes_sweep.py
_sv_hap_rotationnull.py _sv_founder_direction.py _sv_founder_mechanism.py
_sv_passenger_test.py _sv_winning_genetics.py _audit_sv_fitness.py
_winners_ratio.py _winners_sv_depletion.py _plot_winner_ratio_allsites.py
_plot_winner_ratio_vs_purging.py _plot_winners_sv_depletion.py
sv_hap_context.sbatch sv_hap_freqrobust.sbatch
sv_haplotype_axes_sweep.sbatch sv_hap_rotationnull.sbatch
audit_sv_fitness.sbatch
""".split()

MAP["common"] = """
build_af_store.py build_gen_matrices.py build_gen_matrix.py
build_pool_matrix.py build_p0.py build_sample_h_cache.py
build_af_store_array.sh build_gen_matrices_array.sh
rebuild_group_means.sbatch run_sample_h_cache.sbatch
push_af_matrices_to_drive.sh
rerun_kfw_hb
""".split()

MAP["blocks"] = """
blocks_tiling.py recompute_blocks.py dynamic_ld_blocks.py
make_coarse_blocks.py build_hap_membership.py
block_cluster_pc1ve.py blockcoherence_data.py block_founder_vs_evolved.py
block_haplotype_counts.py block_kmer_coverage.py block_missing_sensitivity.py
block_panel_support_tag.py block_unit_frontier.py
coherence_vs_floor.py delta_p_coherence.py eval_block_coherence.py
verify_founder_ve.py _check_clq90_blocks.py
unit_distributions.py plot_unit_distributions.py
build_sv_landscape.py plot_sv_landscape.py
_build_blockcoherence_nb.py _build_blocks_units_nb.py
blocks_recompute.sbatch blocks_recompute_mcf90.sbatch
blocks_recompute_p80.sbatch blocks_clq05_mcf90.sbatch
clq_sweep_mcf90.sbatch cluster_mcf90.sbatch coherence_vs_floor.sbatch
diag_genomewide.sbatch dynld.sbatch unit_dist.sbatch
hap_membership.sbatch hap_membership_clq50.sbatch
hap_membership_clq90.sbatch hap_membership_clq90nosv.sbatch
install_bigld.sbatch
WINDOW_UNIT_VALIDATION.md
hap_blocks bigld_env
""".split()

MAP["qc"] = """
seedmix_identifiability.py seedmix_kmer_identifiability.py _build_9977_nb.py
compare_window_vs_global_af.py plot_coverage_distribution.py
_build_qc_notebook.py _compute_panel_overlap_grenenet.py
_build_panel_overlap_grenenet_nb.py _build_panel_stats_nb.py
_build_support_nb.py
seedmix_validation
""".split()

MAP["archive"] = ["PIPELINE_B_POOLED_MODEL.md"]

EXTERNAL = {"_build_floor_derivation_nb.py":
            os.path.join(REPO, "benchmarks", "localonly_p231")}

STAYS = {"lib.py", "README.md", "GLOBAL_MODE_DECISION.md", "EXPORT_MANIFEST.md",
         "_RESTRUCTURE_PLAN.md", "_do_restructure.py"}


def tracked(rel):
    return subprocess.run(
        ["git", "-C", REPO, "ls-files", "--error-unmatch",
         os.path.join("analysis/grenenet_selection", rel)],
        capture_output=True).returncode == 0


# --- the dirname-chain rewriter -------------------------------------------
def _chain(n, abspath):
    inner = "os.path.abspath(__file__)" if abspath else "__file__"
    return "os.path.dirname(" * n + inner + ")" * n


def deepen(text, old_path):
    """Add one os.path.dirname() ONLY to chains that currently resolve to GEA.

    Scripts routinely carry two inserts: one for their own directory (so they
    can import sibling modules) and one for GEA (so they can import lib).
    Deepening both would fix `import lib` and break the sibling import. So
    each chain is evaluated at the file's CURRENT location and rewritten only
    if it points at GEA.
    """
    out, hits, shift = text, 0, 0
    for m in CHAIN_RE.finditer(text):
        try:
            got = eval(m.group(1), {"os": os, "__file__": old_path})
        except Exception:
            continue
        # the lib-reaching chain is precisely the one landing where lib.py is
        if not os.path.isfile(os.path.join(got, "lib.py")):
            continue
        new = f"sys.path.insert(0, os.path.dirname({m.group(1)}))"
        out = out[:m.start() + shift] + new + out[m.end() + shift:]
        shift += len(new) - (m.end() - m.start())
        hits += 1
    return out, hits


CHAIN_RE = re.compile(
    r"sys\.path\.insert\(\s*0\s*,\s*((?:os\.path\.dirname\()+"
    r"os\.path\.abspath\(__file__\)\)+|(?:os\.path\.dirname\()+__file__\)+)\s*\)")


def verify():
    """Evaluate every dirname-chain for real and check it reaches lib.py.

    A file passes if ANY of its chains lands on lib.py — scripts legitimately
    carry a second insert pointing at their own directory.
    """
    bad, ok, skipped = [], 0, 0
    for root, dirs, files in os.walk(GEA):
        dirs[:] = [d for d in dirs if d not in ("__pycache__", "archive")]
        for f in files:
            if not f.endswith(".py"):
                continue
            p = os.path.join(root, f)
            try:
                t = open(p).read()
            except (UnicodeDecodeError, OSError):
                continue
            if not re.search(r"^\s*(import lib\b|from lib import)", t, re.M):
                continue
            chains = CHAIN_RE.findall(t)
            if not chains:
                skipped += 1          # uses an absolute/cwd path instead
                continue
            resolved = []
            for c in chains:
                try:
                    resolved.append(eval(c, {"os": os, "__file__": p}))
                except Exception:
                    pass
            if any(os.path.isfile(os.path.join(r, "lib.py")) for r in resolved):
                ok += 1
            else:
                bad.append(f"{os.path.relpath(p, GEA)} -> {resolved}")
    print(f"\n=== verify: {ok} scripts resolve `import lib` correctly, "
          f"{skipped} use an absolute/cwd path, {len(bad)} BROKEN")
    for b in bad:
        print("    BROKEN", b)
    return 1 if bad else 0


def main():
    if VERIFY_ONLY:
        return verify()

    seen, problems, untracked = collections.Counter(), [], []
    for dest, items in MAP.items():
        for it in items:
            seen[it] += 1
            if not os.path.exists(os.path.join(GEA, it)):
                problems.append(f"MISSING  {it} (-> {dest})")
            elif os.path.isfile(os.path.join(GEA, it)) and not tracked(it):
                untracked.append(f"{it} -> {dest}")
    for it, n in seen.items():
        if n > 1:
            problems.append(f"DUPLICATE {it} appears {n}x")

    assigned = set(seen) | STAYS | set(EXTERNAL)
    orphans = [f for f in sorted(os.listdir(GEA))
               if os.path.isfile(os.path.join(GEA, f))
               and f.rsplit(".", 1)[-1] in ("py", "sh", "sbatch", "md")
               and f not in assigned]

    print(f"=== plan: {sum(len(v) for v in MAP.values())} items -> "
          f"{len(MAP)} dirs (+{len(EXTERNAL)} external) ===")
    for dest, items in MAP.items():
        print(f"  {dest:28s} {len(items)}")
    if problems:
        print("\n!!! PROBLEMS — nothing will move:")
        for p in problems: print("   ", p)
        return 1
    if untracked:
        print(f"\n--- untracked, SKIPPED (stay at top level): {len(untracked)}")
        for u in untracked: print("   ", u)
    if orphans:
        print(f"\n--- unassigned, stay at top level: {len(orphans)}")
        for o in orphans: print("   ", o)
    if not APPLY:
        print("\n(dry run — pass --apply to execute)")
        return 0

    # ---- collect what will move, and patch BEFORE moving, so each chain can
    # ---- be evaluated at the location it was written for
    plan = []                        # (src_abs, dest_abs, is_dir)
    for dest, items in MAP.items():
        for it in items:
            src = os.path.join(GEA, it)
            if os.path.isfile(src) and not tracked(it):
                continue
            plan.append((src, os.path.join(GEA, dest, it), os.path.isdir(src)))
    for it, dpath in EXTERNAL.items():
        if tracked(it):
            plan.append((os.path.join(GEA, it),
                         os.path.join(dpath, it), False))

    targets = []
    for src, _, is_dir in plan:
        if is_dir:
            for root, dirs, files in os.walk(src):
                dirs[:] = [d for d in dirs if d != "__pycache__"]
                targets += [os.path.join(root, f)
                            for f in files if f.endswith(".py")]
        elif src.endswith(".py"):
            targets.append(src)

    npatch = ntot = 0
    for p in targets:
        try:
            t = open(p).read()
        except (UnicodeDecodeError, OSError):
            continue
        new, hits = deepen(t, p)
        if hits:
            open(p, "w").write(new)
            npatch += 1
            ntot += hits
    print(f"deepened {ntot} GEA-pointing sys.path chains across {npatch} of "
          f"{len(targets)} .py files (own-dir chains left alone)")

    moved_files = {}
    for src, dst, is_dir in plan:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        subprocess.run(["git", "-C", REPO, "mv", src, dst], check=True)
        if not is_dir:
            moved_files[os.path.basename(src)] = os.path.relpath(dst, GEA)
    print(f"moved {len(moved_files)} files, "
          f"{sum(1 for _, _, d in plan if d)} directories")

    refs = 0
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in
                   (".git", "__pycache__", "old_docs")]
        for f in files:
            if not f.endswith((".py", ".sh", ".sbatch", ".md", ".ipynb")):
                continue
            p = os.path.join(root, f)
            try:
                t = open(p).read()
            except (UnicodeDecodeError, OSError):
                continue
            o = t
            for it, newrel in moved_files.items():
                if it.endswith(".md"):
                    continue
                t = re.sub(r"(analysis/grenenet_selection/)" + re.escape(it) + r"\b",
                           r"\1" + newrel, t)
            if t != o:
                open(p, "w").write(t)
                refs += 1
    print(f"rewrote script-path references in {refs} files")
    return verify()


sys.exit(main())
