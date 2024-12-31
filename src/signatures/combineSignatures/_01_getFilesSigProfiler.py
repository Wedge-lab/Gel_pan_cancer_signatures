import os
import sys

import numpy as np
import pandas as pd
from dotenv import load_dotenv

load_dotenv()
RESULT_DIR = os.getenv("RESULT_DIR")

if __name__ == "__main__":
    results_dir, cohort_list, output_dir, sig_type, max_sigs = sys.argv[1:]

    cohorts = pd.read_csv(cohort_list, sep="\t")
    full_stats_df = pd.DataFrame()

    for index, row in cohorts.iterrows():
        # Load stats data from runs and reformat columns
        stats_df = pd.read_csv(f"{row.cohort_dir}/All_solutions_stat.csv")
        stats_df["best_solution"] = stats_df.Signatures.map(lambda x: "*" in x).astype(
            int
        )
        stats_df["Signatures"] = np.array(
            [str(sig).rstrip("*") for sig in list(stats_df.Signatures)]
        ).astype(int)
        stats_df["Matrix Frobenius%"] = (
            stats_df["Matrix Frobenius%"].map(lambda x: x.rstrip("%")).astype(float)
        )
        stats_df.rename(
            {
                "Signatures": "signatures",
                "Matrix Frobenius%": "matrix_frobenius",
                "Minimum Stability": "min_stability",
                "Stability (Avg Silhouette)": "avg_stability",
            },
            axis=1,
            inplace=True,
        )
        # Add file paths for Signatures and Standard Errors
        stats_df["signatures_file"] = [
            f"{row.cohort_dir}/All_Solutions/{sig_type}_{sig}_Signatures/Signatures/{sig_type}_S{sig}_Signatures.txt"
            for sig in stats_df.signatures
        ]
        stats_df["signatures_SE_file"] = [
            f"{row.cohort_dir}/All_Solutions/{sig_type}_{sig}_Signatures/Signatures/{sig_type}_S{sig}_Signatures_SEM_Error.txt"
            for sig in stats_df.signatures
        ]
        stats_df["activities_file"] = [
            f"{row.cohort_dir}/All_Solutions/{sig_type}_{sig}_Signatures/Activities/{sig_type}_S{sig}_NMF_Activities.txt"
            for sig in stats_df.signatures
        ]
        stats_df["samples_file"] = [
            f"{row.cohort_dir}/Samples.txt" for sig in stats_df.signatures
        ]
        stats_df["cohort"] = row.cohort

        # Decomposed Signatures
        if sig_type == "SBSSV":
            stats_df["signatures_decomposed_file"] = [
                f"{row.cohort_dir}/All_Solutions/{sig_type}_{sig}_Signatures/Decompose_Solution/Signatures/Decompose_Solution_Signatures.txt"
                for sig in stats_df.signatures
            ]
            stats_df["activities_decomposed_file"] = [
                f"{row.cohort_dir}/All_Solutions/{sig_type}_{sig}_Signatures/Decompose_Solution/Activities/Decompose_Solution_Activities.txt"
                for sig in stats_df.signatures
            ]
        elif sig_type == "SBSCNV":
            # CAN'T RUN AIC WITH CNV BECAUSE SIGNATURE DECOMPOSITION AND ASSIGNMENT HASN'T HAPPENED FOR ALL SOLUTIONS :
            stats_df["signatures_decomposed_file"] = [
                f"{RESULT_DIR}/results/signatures/decomposedSignatures/decomposed_CNV48_COSMICref/{row.cohort}/{sig}/Decompose_Solution/Signatures/Decompose_Solution_Signatures.txt"
                for sig in stats_df.signatures
            ]
            stats_df["activities_decomposed_file"] = [
                f"{RESULT_DIR}/results/signatures/decomposedSignatures/decomposed_CNV48_COSMICref/{row.cohort}/{sig}/Decompose_Solution/Activities/Decompose_Solution_Activities.txt"
                for sig in stats_df.signatures
            ]
        else:
            stats_df["signatures_decomposed_file"] = [
                f"{row.cohort_dir}/{sig}_Signatures/{sig_type}/Suggested_Solution/COSMIC_{sig_type}_Decomposed_Solution/Signatures/COSMIC_{sig_type}_Signatures.txt"
                for sig in stats_df.signatures
            ]
            stats_df["activities_decomposed_file"] = [
                f"{row.cohort_dir}/{sig}_Signatures/{sig_type}/Suggested_Solution/COSMIC_{sig_type}_Decomposed_Solution/Activities/COSMIC_{sig_type}_Activities_refit.txt"
                for sig in stats_df.signatures
            ]

        stats_df["max_sigs"] = max_sigs
        if (sig_type == "SBSCNV") & (row.cohort == "Prostate"):
            stats_df["max_sigs"] = 8
        full_stats_df = pd.concat((full_stats_df, stats_df))

    # Save cohort level filenames
    full_stats_df.to_csv(
        f"{output_dir}/all_cohort_stats.tsv", index=False, header=True, sep="\t"
    )
