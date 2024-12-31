import os
import sys

import numpy as np
import pandas as pd
import scipy.stats
from dotenv import load_dotenv

from signatures.plotting.combinedSignatures import sig_dirs

load_dotenv()
DATA_DIR = os.getenv("DATA_DIR")

if __name__ == "__main__":
    samples_file, activities_file = sys.argv[1:]

    # Load samples file
    samples_df = pd.read_csv(samples_file, sep="\t")
    print("Samples: ", len(samples_df))

    print("Load signature activity files")
    activities_df = pd.DataFrame(index=samples_df["sample_id"])
    signatures_files = {
        key: f"{sig_dir}/Combined_Solution_Activities.tsv"
        for key, sig_dir in sig_dirs.items()
    }
    for sig_set, sig_file in signatures_files.items():
        sig_df = pd.read_csv(sig_file, sep="\t").set_index("Samples")
        sig_df[f"T{sig_set}"] = np.sum(sig_df, axis=1)
        activities_df = (
            pd.merge(
                activities_df, sig_df, how="left", left_index=True, right_index=True
            )
            .fillna(0)
            .astype(int)
        )
    activities_df.index.rename("sample_id", inplace=True)
    print(activities_df.head())
    print("Activities: ", activities_df.shape)

    # Get activities for sample set
    activities_df = pd.merge(
        samples_df[["sample_id"]],
        activities_df,
        left_on="sample_id",
        right_index=True,
        how="inner",
    ).set_index("sample_id")
    print("Combined: ", activities_df.shape)

    # Generate mock signature
    samples_df.set_index(["sample_id", "group"], inplace=True)
    samples_df = (samples_df - np.mean(samples_df)) / np.std(samples_df)
    samples_df["intercept"] = np.ones(len(samples_df))
    # Negative Binomial scipy parameterisation - n=theta, p=theta/(mu+theta)
    theta = 1
    betaNB = np.ones(samples_df.shape[1]) * 1
    mu_NB = np.exp(betaNB @ np.array(samples_df).T)
    k_NB = scipy.stats.nbinom.rvs(theta, theta / (mu_NB + theta))
    # Binomial for zero inflation
    betaZI = np.ones(samples_df.shape[1]) / 2
    expit = lambda x: 1 / (1 + np.exp(-x))
    mu_ZI = expit(betaNB @ np.array(samples_df).T)
    k_ZI = scipy.stats.binom.rvs(1, mu_ZI)
    activities_df["MUTmock"] = k_NB * k_ZI
    print(
        f"Mock range: {np.min(activities_df['MUTmock'])} - {np.max(activities_df['MUTmock'])}"
    )

    # Save activities file
    activities_df.to_csv(activities_file, sep="\t")
    # activities_df[["MUTmock"]].to_csv(activities_file, sep="\t")
