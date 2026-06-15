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

def adjust_WZA_with_spline(wza_df, roller=50, minEntries=40, deg=2):
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

    # Generating weights for polynomial function (deg, default 2 = canonical) - standard deviation
    sd_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_sd, deg=deg)
    sd_polynomial_model = np.poly1d(sd_weights)

    # Generating weights for polynomial function (deg, default 2 = canonical) - mean
    mean_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_means, deg=deg)
    mean_polynomial_model = np.poly1d(mean_weights)


    sd_predictions = sd_polynomial_model(wza_df["SNPs"])
    # SAFETY FLOOR (kMate local fix; see Booker email thread): in sparse large-window
    # tails the degree-2 polynomial can predict negative/zero SD -> norm.cdf returns NaN.
    # Floor the predicted SD at the smallest strictly-positive prediction (fallback: the
    # empirical SD of Z). Distribution-agnostic; avoids switching the polynomial degree.
    _pos = sd_predictions[sd_predictions > 0]
    _floor = _pos.min() if _pos.size else np.nanstd(wza_df["Z"].to_numpy())
    sd_predictions = np.clip(sd_predictions, a_min=_floor, a_max=None)

    mean_predictions = mean_polynomial_model(wza_df["SNPs"])

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

        WZA_DF = adjust_WZA_with_spline(WZA_DF_tmp, roller=args.roller, minEntries=args.min_entries, deg=args.poly_deg)
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

        WZA_DF = adjust_WZA_with_spline(WZA_DF_tmp, roller=args.roller, minEntries=args.min_entries, deg=args.poly_deg)

    WZA_DF.to_csv(args.output, index=False)

main()