import pandas as pd
import scipy.stats
import numpy as np
import sys, argparse
from scipy.stats import norm
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype

def WZA(gea, statistic, MAF_filter=0.0):
    ## gea - the name of the pandas dataFrame with the gea results
    ## statistic - the name of the column with your p-values
    ## MAF_filter - the lowest MAF you will tolerate
    ## NOTE, this function assumes that the DataFrame has a column named pbar_qbar
    
    ## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)
    gea[statistic] = gea[statistic].clip(lower=1e-15)
    gea[statistic] = gea[statistic].replace(1, 1-1e-3)
    
    # Convert the p-values into 1-sided Z scores (hence the 1 - p-values)
    gea["z_score"] = scipy.stats.norm.ppf(1 - np.array(gea[statistic], dtype=float))
    gea["pbar_qbar"] = gea["MAF"] * (1 - gea["MAF"])

    ## Apply the MAF filter
    gea_filt = gea[gea["MAF"] >= MAF_filter].copy()
    
    if gea_filt.shape[0] != gea.shape[0]:
        print(f"Window filtered out due to MAF. Initial SNPs: {gea.shape[0]}, After MAF filter: {gea_filt.shape[0]}")
    
    if gea_filt.shape[0] == 0:
        return np.nan
    
    ## Calculate the numerator and the denominator for the WZA
    gea_filt["weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z_score"]
    gea_filt["weiZ_den"] = gea_filt["pbar_qbar"] ** 2
    numerator = gea_filt["weiZ_num"].sum()
    denominator = np.sqrt(gea_filt["weiZ_den"].sum())

    if denominator == 0:
        print(f"Denominator is zero for window. Returning NaN for window.")
        return np.nan

    ## We've calculated the num. and the den., let's take the ratio
    weiZ = numerator / denominator

    ## Return the final dataframe
    return weiZ

def top_candidate(gea, thresh):
    hits = (gea["pVal"] < thresh).sum()
    snps = gea.shape[0]
    top_candidate_p = scipy.stats.binomtest(hits, snps, thresh, alternative="greater").pvalue
    return top_candidate_p, hits

def adjust_WZA_with_spline(wza_df, roller=50, minEntries=40, deg=2, sd_fit="poly",
                           mean_fit="interp", mean_poly_deg=5):
    wza_df.to_csv('before_filtering_wza_df.csv')
    # remove null Z values - they won't help us
    wza_t = wza_df[~wza_df.Z.isnull()].reset_index()
    wza_s = wza_t.sort_values('SNPs')

    problematic_windows = wza_df[wza_df["SNPs"] < minEntries]
    print("Problematic windows with fewer SNPs than minEntries:")
    problematic_windows.to_csv('problematic_windows.csv')
    
    rolled_Z_vars = wza_s.Z.rolling(window=roller, min_periods=minEntries).var()
    masking_array = ~rolled_Z_vars.isnull()
    rolled_Z_sd = np.sqrt(rolled_Z_vars)[masking_array]
    rolled_Z_means = wza_s.Z.rolling(window=roller, min_periods=minEntries).mean()[masking_array]
    rolled_mean_SNP_number = wza_s.SNPs.rolling(window=roller, min_periods=minEntries).mean()[masking_array]

    if rolled_Z_vars.isnull().any():
        print(f"WARNING: Rolling window calculation resulted in NaN for some windows.")

    Xs = np.asarray(rolled_mean_SNP_number, dtype=float)
    target = np.asarray(wza_df["SNPs"], dtype=float)
    if sd_fit == "isotonic":
        # kMate fix (2026-07-21; validated in wza_investigation/wza_sd_fix_test.ipynb
        # via a permutation null + an independent audit): on clq0.9/mcf90 blocks the
        # SNP-count tail is so sparse that BOTH the canonical deg-2 and the author's
        # deg-7 SD polynomials EXTRAPOLATE pathologically near the cap — deg-2's
        # parabola turns over and UNDER-predicts SD (positive-but-tiny -> fabricated
        # Z_pVal==0 for unremarkable blocks; negative -> NaN), deg-7 explodes. The null
        # SD-vs-SNP-count is genuinely monotone-up-then-plateau (permutation-null
        # isotonic-R2~=0.98), and deg-2 fabricates p==0 even with climate SHUFFLED
        # (signal-free) while isotonic fabricates none. So fit SD with a
        # monotone-non-decreasing isotonic regression (flat beyond support). See
        # STATUS_clq90 §0.
        #
        # The MEAN fit is a separate choice (`mean_fit`), and it matters just as much:
        # p = 1 - Phi((Z - mean)/sd). The original isotonic branch predicted the mean with
        # np.interp -- piecewise-linear through every rolling point, i.e. NO smoothing.
        # Because run_wza.py feeds raw p and main() rank-transforms it genome-wide
        # (`csv["pVal"] = csv[stat].rank()/n`), E[z]=0 by construction and mean(Z) has
        # almost no real dependence on block size -- so interp just tracks rolling noise.
        # Measured out-of-sample (notebooks/cap_poly_decision.ipynb §4) it is the WORST of
        # six candidates in all 8 blockdef x class cases, ~25-40% worse than any smoother;
        # isotonic_auto / deg2 / const agree within a couple of percent of each other.
        #
        # PRODUCTION (2026-07-28, FINAL): mean_fit="poly_clamp", mean_poly_deg=5.
        # UPSTREAM FITS THE TREND HERE TOO (an unclamped deg-2, general_WZA_script.py:86), so
        # fitting it is the faithful choice; `const` was a kMate invention and is retired.
        # Judged on accuracy at the LARGE blocks -- the only place the candidates differ, and
        # a place no aggregate RMSE can see (the 8-15 blocks past the rolling-support edge are
        # ~0.02% of rolling points) -- across all 24 blockdef x model x class cells:
        #        mean |predicted - empirical|   worst cell
        #   deg5_clamp        1.60                4.41     <- PRODUCTION
        #   isotonic_auto     1.90                5.05     (also: direction flips, see below)
        #   deg2_clamp        2.15                7.26     (= upstream's degree, clamped)
        #   const             2.81               10.35     <- WORST, retired
        # Half-split OOS RMSE had favoured `const`, but that is a bias-variance artifact: each
        # half holds only ~5 of the ~10 largest blocks, so the validation TARGET is very noisy
        # there and a zero-variance predictor wins even when the trend is real. Production
        # fits on ALL blocks, where the trend is well estimated -- and it replicates (sign of
        # the large-block slope agrees in >=90% of half-splits for 6/24 cells, and is visibly
        # present in more).
        # `isotonic_auto` is rejected outright: its Spearman direction FLIPS between adjacent
        # climate axes in 20 of 24 cells (mean_fit_direction_audit.py, 480 scans; |rho| as low
        # as 0.0005 -- "significant" only because n~40k rolling points), which would apply
        # opposite-trending mean corrections across the 60 scans the gene unions are built on.
        from sklearn.isotonic import IsotonicRegression
        ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
        ir.fit(Xs, np.asarray(rolled_Z_sd, dtype=float))
        sd_predictions = ir.predict(target)
        means = np.asarray(rolled_Z_means, dtype=float)
        if mean_fit == "isotonic_auto":
            # direction chosen by the data; collapses to ~constant when there is no trend
            im = IsotonicRegression(increasing="auto", out_of_bounds="clip")
            im.fit(Xs, means)
            mean_predictions = im.predict(target)
        elif mean_fit == "const":
            mean_predictions = np.full_like(target, means.mean())
        elif mean_fit == "poly_clamp":
            # PRODUCTION (2026-07-28). Upstream fits the mean with a deg-2 polynomial
            # (general_WZA_script.py:86) -- i.e. it FITS THE TREND; `const` was a kMate
            # invention and measured WORST of four at the large blocks (mean |error| 2.81,
            # worst cell 10.35, vs 1.60/4.41 for this). Degree 5 + CLAMPING are the only
            # departures, both forced by our block-size range:
            #   * clamping (fit inside rolling support, hold the boundary value outside)
            #     is what prevents the unbounded blow-up: plain deg-2 predicts mean(Z) =
            #     -289 at the largest snp block where the empirical value is -10.5.
            #   * degree 5 beats 2/7/10 on large-block accuracy; degrees >=15 Runge-
            #     oscillate near the sparse upper end (top-band RMSE 11-216) and are unusable.
            # x is standardized so the fit is numerically well conditioned.
            mu = Xs.mean(); sg = Xs.std() or 1.0
            pm = np.poly1d(np.polyfit((Xs - mu) / sg, means, deg=mean_poly_deg))
            mean_predictions = pm((np.clip(target, Xs.min(), Xs.max()) - mu) / sg)
        elif mean_fit == "poly":
            mean_predictions = np.poly1d(np.polyfit(Xs, means, deg=deg))(target)
        else:  # "interp" -- pre-2026-07-28 behaviour, kept for reproducing old outputs
            mean_predictions = np.interp(target, Xs, means)
    else:
        # canonical polynomial interpolation of the rolling mean/SD (deg=2 Booker, or 7).
        # NO SAFETY FLOOR (removed 2026-07-03): a negative-SD tail yields NaN (excluded
        # downstream), never fabricated significance. For sparse-tail block definitions
        # (e.g. clq0.9/mcf90) prefer sd_fit="isotonic" — the polynomial extrapolates
        # badly there (see the isotonic branch above / STATUS_clq90 §0).
        sd_predictions = np.poly1d(np.polyfit(Xs, np.asarray(rolled_Z_sd, dtype=float), deg=deg))(target)
        mean_predictions = np.poly1d(np.polyfit(Xs, np.asarray(rolled_Z_means, dtype=float), deg=deg))(target)

    with np.errstate(invalid="ignore"):
        wza_p_values = [1 - norm.cdf(wza_df["Z"][i], loc=mean_predictions[i], scale=sd_predictions[i]) for i in range(wza_df.shape[0])]
    wza_df["Z_pVal"] = wza_p_values

    return wza_df

def main():
    ## Define command line args
    parser = argparse.ArgumentParser(description="A script that implements the WZA, a method for combining evidence across closely linked SNPs in GEA studies.")

    parser.add_argument("--correlations", "-c", required=True, dest="correlations", type=str, help="The file containing the correlations")
    parser.add_argument("--summary_stat", "-s", required=True, dest="summary_stat", type=str, help="The name of the column you are analysing")
    parser.add_argument("--window", "-w", required=True, dest="window", type=str, help="The name of column containing the windows you want to analyse")
    parser.add_argument("--output", required=True, dest="output", type=str, help="The name of the output file")
    parser.add_argument("--sample_snps", required=False, dest="sample_snps", type=int, default=0, help="[OPTIONAL] Give the number of SNPs you want to downsample to.")
    parser.add_argument("--resamples", required=False, dest="resamples", type=int, default=100, help="[OPTIONAL] Number of times to resample WZA scores")
    parser.add_argument("--min_snps", required=False, dest="min_snps", type=int, default=2, help="[OPTIONAL] Minimum number of SNPs per window")
    parser.add_argument("--large_i_small_p", required=False, action="store_true", help="[OPTIONAL] Extreme values of the summary stat you're using are large values.")
    parser.add_argument("--top_candidate_threshold", required=False, dest="top_candidate_threshold", type=float, default=99, help="[OPTIONAL] Percentile threshold for top-candidate test")
    parser.add_argument("--verbose", "-v", required=False, action="store_true", help="[OPTIONAL] Verbose mode")
    parser.add_argument("--MAF", required=False, dest="MAF", type=str, help="[OPTIONAL] MAF column name.")
    parser.add_argument("--sep", required=False, dest="sep", type=str, default="\t", help="What separator do you use in your file?")
    parser.add_argument("--retain", required=False, dest="retain", nargs="+", type=str, help="Columns to add to output file")
    parser.add_argument("--no_SNP_number_correction", required=False, action="store_true", help="Provide this flag if you just want the raw WZA scores")
    parser.add_argument("--empiricalP", required=False, action="store_true", help="Flag for empirical p-values")
    parser.add_argument("--maf_filter", required=False, dest="maf_filter", type=float, default=0.0, help="[OPTIONAL] MAF cutoff applied in main (LOCAL COPY: default 0 = keep all; orig hardcoded 0.05)")
    parser.add_argument("--poly_deg", required=False, dest="poly_deg", type=int, default=2, help="[OPTIONAL] SNP-number-correction polynomial degree (LOCAL: 2=canonical Booker; phase-1 used 7)")
    parser.add_argument("--roller", required=False, dest="roller", type=int, default=50, help="[OPTIONAL] Rolling-window size for SNP-number correction (canonical 50)")
    parser.add_argument("--min_entries", required=False, dest="min_entries", type=int, default=40, help="[OPTIONAL] Min entries per rolling window (canonical 40; phase-1 used 10)")
    parser.add_argument("--sd_fit", required=False, dest="sd_fit", type=str, default="poly", choices=["poly", "isotonic"], help="[LOCAL] SNP-number-correction fit: 'poly' (deg via --poly_deg; canonical Booker) or 'isotonic' (monotone-non-decreasing SD; robust for sparse-tail block defs like clq0.9/mcf90 where poly extrapolates to fabricated p==0/NaN)")
    parser.add_argument("--mean_fit", required=False, dest="mean_fit", type=str, default="interp", choices=["interp", "isotonic_auto", "const", "poly", "poly_clamp"], help="[LOCAL, only used when --sd_fit isotonic] how to predict the MEAN of Z vs SNP count. 'poly_clamp' = PRODUCTION 2026-07-28 (degree via --mean_poly_deg, default 5; fitted inside rolling support and held flat outside). Upstream uses an UNCLAMPED deg-2 here (general_WZA_script.py:86); clamping prevents mean(Z)=-289 at our largest blocks. 'const' scored WORST of four at large blocks (mean |err| 2.81 vs 1.60). 'isotonic_auto' rejected: direction flips across axes in 20/24 cells. 'interp' = pre-2026-07-28, unsmoothed.")

    parser.add_argument("--mean_poly_deg", required=False, dest="mean_poly_deg", type=int, default=5, help="[LOCAL] polynomial degree for --mean_fit poly_clamp. 5 is production: it beats 2/7/10 on large-block accuracy; >=15 Runge-oscillates near the sparse tail and is unusable.")

    args = parser.parse_args()

    ## Print all arguments passed at the beginning
    print("\nArguments passed:")
    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")


    csv = pd.read_csv(args.correlations, sep=args.sep, engine="python")
    print(list(csv))
    if args.window not in list(csv):
        print("The window variable you provided is not in the dataframe you gave")
        return
    if args.summary_stat not in list(csv):
        print("The summary statistic variable you provided is not in the dataframe you gave")
        return

    if "MAF" in list(csv):
        pass
    else:
        csv["MAF"] = csv[args.MAF].copy()

    maf_filter = args.maf_filter   # LOCAL COPY: parameterized (orig hardcoded 0.05)
    ## Apply MAF filter here...
    csv = csv[csv["MAF"] > maf_filter]

    if args.empiricalP:
        csv["pVal"] = csv[args.summary_stat].copy()
    else:
        if args.large_i_small_p:
            if args.summary_stat == "RDA":
                csv["pVal"] = 1 - (csv[args.summary_stat] ** 2).rank() / csv.shape[0]
            else:
                csv["pVal"] = 1 - csv[args.summary_stat].rank() / csv.shape[0]
        else:
            csv["pVal"] = csv[args.summary_stat].rank() / csv.shape[0]

    csv_genes = csv[csv[args.window] != "None"]

    if args.verbose:
        print("here's a peek at the input data")
        print(csv_genes.head())

    csv_gb_gene = csv_genes.groupby(args.window)
    print(args.sample_snps)
    if args.sample_snps == -1:
        csv_gb_gene_SNP_count = csv_genes.groupby(args.window)
        num_SNP_list = np.array([s[1].shape[0] for s in csv_gb_gene_SNP_count if s[1].shape[0] >= args.min_snps])
        max_SNP_count = int(np.percentile(num_SNP_list[num_SNP_list != 0], 75))

        if args.verbose:
            print("Using the 75th percentile number of SNPs as the maximum in each gene:", max_SNP_count)
    elif args.sample_snps == 0:
        max_SNP_count = 1e6
        if args.verbose:
            print("The maximum number of SNPs in each gene:", max_SNP_count)
    else:
        max_SNP_count = args.sample_snps  # Use the value passed in the argument
        if args.verbose:
            print("The maximum number of SNPs in each gene:", max_SNP_count)
    

    all_genes = []
    count = 0
    for g in csv_gb_gene:
        count += 1
        gene = g[0]
        gene_df = g[1].copy()

        original_snp_count = gene_df.shape[0]  # Track the original number of SNPs
        ## Perform the WZA on the annotations in the contig using parametric p-values
        if original_snp_count <= max_SNP_count:
            wza = WZA(gene_df, "pVal")
            snp_count_used = original_snp_count  # Use the original SNP count since no downsampling occurred
        else:
            snp_count_used = max_SNP_count  # This reflects the reduced number of SNPs after downsampling
            wza = np.array([WZA(gene_df.sample(max_SNP_count), "pVal") for i in range(args.resamples)]).mean()

        if pd.isna(wza):
            print(f"WARNING: WZA calculation resulted in NaN for gene: {gene}")

        top_candidate_p, hits = top_candidate(gene_df, 1 - (args.top_candidate_threshold / 100))

        if args.verbose:
            print(f"\nGene #: {count}\tgene: {gene}\tWZA: {wza}\tTC: {top_candidate_p}")

        output = {
            "gene": gene,
            "SNPs": snp_count_used,
            "hits": hits,
            "Z": wza,
            "top_candidate_p": top_candidate_p
        }

        all_genes.append(output)

    if args.retain is not None:
        WZA_DF_temp = pd.DataFrame(all_genes)
        print(WZA_DF_temp.SNPs.var())

        if args.verbose:
            print("\nAdding retained columns to the final dataframe")

        retained_df_list = []
        for r in args.retain:
            if is_string_dtype(csv[r]):
                retained_df_list.append(csv.groupby(args.window)[r].apply(lambda x: x.iloc[0]))
            elif is_numeric_dtype(csv[r]):
                retained_df_list.append(csv.groupby(args.window)[r].mean())

        retained_df = pd.concat(retained_df_list, axis=1)
        WZA_DF_tmp = pd.concat([WZA_DF_temp.set_index("gene"), retained_df], axis=1).reset_index()

        if "index" in list(WZA_DF_tmp):
            WZA_DF_tmp = WZA_DF_tmp[WZA_DF_tmp["index"] != "None"]
        WZA_DF_tmp.rename(index={"index": "gene"}, inplace=True)

        if args.no_SNP_number_correction:
            WZA_DF_tmp.to_csv(args.output, index=False)
            return

        if WZA_DF_tmp.SNPs.var() == 0:
            print("There is no variation in SNP number among your windows")
            print("SNP number correction will achieve nothing, outputting raw WZA scores")
            WZA_DF_tmp.to_csv(args.output, index=False)
            return

        WZA_DF = adjust_WZA_with_spline(WZA_DF_tmp, roller=args.roller, minEntries=args.min_entries, deg=args.poly_deg, sd_fit=args.sd_fit, mean_fit=args.mean_fit, mean_poly_deg=args.mean_poly_deg)
    else:
        WZA_DF_tmp = pd.DataFrame(all_genes)
        if WZA_DF_tmp.SNPs.var() == 0:
            print("NOTE!\nThere is no variation in SNP number among your windows")
            print("SNP number correction will achieve nothing, outputting raw WZA scores")
            WZA_DF_tmp.to_csv(args.output, index=False)
            return

        if args.no_SNP_number_correction:
            WZA_DF_tmp.to_csv(args.output, index=False)
            return

        WZA_DF = adjust_WZA_with_spline(WZA_DF_tmp, roller=args.roller, minEntries=args.min_entries, deg=args.poly_deg, sd_fit=args.sd_fit, mean_fit=args.mean_fit, mean_poly_deg=args.mean_poly_deg)

    WZA_DF.to_csv(args.output, index=False)

main()