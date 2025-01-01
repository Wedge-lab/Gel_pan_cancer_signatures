import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

from signatures.config import load_environment
from signatures.sampleCuration.pan_cancer_sample import getSamples

load_environment()
DATA_DIR = os.getenv("DATA_DIR")

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    samples_file, targets_file = sys.argv[1:]

    # Get sample dataframe
    print("Load samples")
    sample_df = getSamples(hyper=True)
    sample_df["sample_id"] = (
        sample_df.participant_id.astype(str)
        + "_"
        + sample_df.tumour_sample_platekey
        + "_"
        + sample_df.germline_sample_platekey
    )

    # Remove samples where sex doesn't match organ
    sample_df = sample_df[
        (
            ~(
                sample_df.tumour_tissue.map(lambda x: x in ["Ovary", "Uterus"])
                & (sample_df.is_female == -1)
            )
        )
        & (
            ~(
                sample_df.tumour_tissue.map(lambda x: x in ["Testis", "Prostate"])
                & (sample_df.is_female == 1)
            )
        )
    ]
    print("Removed male samples for Ovary, Uterus and female for Testis, Prostate")
    print(f"Number of samples: {len(sample_df)}")

    # Get genotype grid
    genotype_df = pd.read_csv(
        f"{DATA_DIR}/pancancer_mutationTimeR_updated_CNACalls_2023_09_16.tsv",
        usecols=[
            "participant_id",
            "tumour_sample_platekey",
            "germline_sample_platekey",
            "is_WGD",
        ],
        sep="\t",
    )  # 'ploidy',
    # genotype_df['ploidy'] = (genotype_df.ploidy>2.5).astype(int)
    genotype_df["is_WGD"] = genotype_df.is_WGD.astype(int)
    genotype_df["sample_id"] = (
        genotype_df.participant_id.astype(str)
        + "_"
        + genotype_df.tumour_sample_platekey
        + "_"
        + genotype_df.germline_sample_platekey
    )
    genotype_df = genotype_df.drop(
        np.intersect1d(
            genotype_df.keys(),
            [
                "participant_id",
                "tumour_sample_platekey",
                "germline_sample_platekey",
                "day_sampling",
                "sample_platekey",
                "tumour_pseudo_id",
            ],
        ),
        axis=1,
    )

    # Get the number of chromothripsis and chromoplexy mutations
    sv_dir = "/re_gecip/cancer_pan/rhoulston/analysisResults/SVClassification/"
    sv_counts = pd.DataFrame()
    for cohort in os.listdir(sv_dir):
        if cohort != "README.md":
            sv_counts = pd.concat(
                (
                    sv_counts,
                    pd.read_csv(
                        f"{sv_dir}/{cohort}/summaryLevel/sampleClusterClassCounts.tsv",
                        sep="\t",
                    ),
                )
            )
    sv_counts["sample_id"] = (
        sv_counts.participant_id.astype(str)
        + "_"
        + sv_counts.tumour_sample_platekey
        + "_"
        + sv_counts.germline_sample_platekey
    )

    # sv_counts = sv_counts[~pd.isnull(sv_counts['n_chromothripsis'])]
    # cluster_types = ['n_chromoplexy', 'n_chromothripsis', 'n_tandem_duplication']
    cluster_types = ["n_tandem_duplication"]

    for cluster_type in cluster_types:
        sv_counts[re.sub("^n_", "is_", cluster_type)] = (
            sv_counts[cluster_type] > 0
        ).astype(int)
    genotype_df = pd.merge(
        genotype_df,
        sv_counts[
            ["sample_id"]
            + cluster_types
            + [re.sub("^n_", "is_", cluster_type) for cluster_type in cluster_types]
        ],
        on="sample_id",
        how="inner",
    )

    # Get the number of Kataegis events
    ktg_dir = (
        f"{DATA_DIR}/kataegis_analysis/results/SigProfilerClusters_persample/summary/"
    )
    ktg_counts = pd.read_csv(
        f"{ktg_dir}/All_10983_clustered_events_table.tsv",
        sep="\t",
        usecols=["sample_id", "kataegis"],
    )
    ktg_counts.rename({"kataegis": "n_kataegis"}, axis=1, inplace=True)
    ktg_counts = ktg_counts[~pd.isnull(ktg_counts["n_kataegis"])]
    cluster_types = ["n_kataegis"]
    for cluster_type in cluster_types:
        ktg_counts[re.sub("^n_", "is_", cluster_type)] = (
            ktg_counts[cluster_type] > 0
        ).astype(int)
    genotype_df = pd.merge(
        genotype_df,
        ktg_counts[
            ["sample_id"]
            + cluster_types
            + [re.sub("^n_", "is_", cluster_type) for cluster_type in cluster_types]
        ],
        on="sample_id",
        how="inner",
    )

    # Sample list which is also in NCRAS
    sample_list = np.intersect1d(genotype_df.sample_id, sample_df.sample_id)
    sample_df = sample_df.set_index("sample_id").loc[sample_list]
    genotype_df = genotype_df.set_index("sample_id").loc[sample_list]

    # Confounding variables
    variables = ["log_age", "is_female"] + [f"pc{i}" for i in range(1, 4)]

    # Construct X matrix
    print("Construct X")
    X = pd.DataFrame(
        {
            variable: sample_df[variable]
            for variable in variables
            + [
                "tumour_tissue",
                "tumour_group",
                "participant_id",
                "tumour_sample_platekey",
                "germline_sample_platekey",
            ]
        }
    )

    # Get non-nan values
    [
        print(key, np.sum(~np.isnan(X[key])))
        for key in X.columns
        if X[key].dtype in [float, int]
    ]
    cuts = np.prod(
        np.array(
            [~np.isnan(X[key]) for key in X.columns if X[key].dtype in [float, int]]
        ),
        axis=0,
    ).astype(bool)
    print(f"Non-nan: {np.sum(cuts)}")

    # Apply cuts
    X = X[cuts]

    # Save
    X["sample_id"] = (
        X.participant_id.astype(str)
        + "_"
        + X.tumour_sample_platekey
        + "_"
        + X.germline_sample_platekey
    )
    X.rename(columns={"tumour_group": "group"}, inplace=True)
    print(X.columns)
    X[["sample_id", "group"] + variables].to_csv(samples_file, index=False, sep="\t")

    # Filter genotype_df for cuts
    print(genotype_df.keys())
    print(genotype_df.head())
    genotype_df = (genotype_df.loc[X.sample_id]).astype(int)
    print(len(sample_df), len(genotype_df))

    # Save targets file
    genotype_df.to_csv(targets_file, sep="\t")
