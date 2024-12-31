import sys
import warnings

import numpy as np
import pandas as pd
import scipy.stats
import treatments

from signatures.associations.sampleCuration.pan_cancer_sample import getSamples

warnings.filterwarnings("ignore")


if __name__ == "__main__":
    samples_file, targets_file, n = sys.argv[1:]
    n = int(n)

    # Get sample dataframe
    print("Load samples")
    sample_df = getSamples(hyper=True)
    print("Primary only!")
    sample_df = sample_df[sample_df.tumour_type == "PRIMARY"]
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
    print("Removed male samples for Ovary,Uterus and female for Testis,Prostate")
    print(f"Number of samples: {len(sample_df)}")

    # Get and save treatment grid
    treatment_df = treatments.getTreatment2(
        sample_df[
            [
                "participant_id",
                "tumour_sample_platekey",
                "germline_sample_platekey",
                "sample_id",
                "day_sampling",
            ]
        ],
        base="",
    ).reset_index()
    treatment_df.rename(
        columns=dict(
            zip(
                treatment_df.keys()[3:],
                [x.replace(" ", "_") for x in treatment_df.keys()[3:]],
            )
        ),
        inplace=True,
    )
    treatment_df = treatment_df.drop(
        np.intersect1d(
            treatment_df.keys(),
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

    # Sample list which is also in NCRAS
    sample_list = np.intersect1d(treatment_df.sample_id, sample_df.sample_id)
    sample_df = sample_df.set_index("sample_id").loc[sample_list]
    treatment_df = treatment_df.set_index("sample_id").loc[sample_list]

    # Confounding variables
    variables = ["log_age", "is_female"] + [
        f"pc{i}" for i in range(1, 4)
    ]  # logit_purity',

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

    # Filter treatment_df for cuts
    print(treatment_df.keys())
    treatment_df = (treatment_df.loc[X.sample_id] > 0).astype(int)
    print(len(sample_df), len(treatment_df))

    # Generate mock targets - for 1%,5%,20%,50% hit rate
    for hit_rate in [10, 5, 1]:
        treatment_df[f"MOCK{hit_rate:03d}"] = scipy.stats.binom.rvs(
            n, hit_rate / 1000, size=len(treatment_df)
        )
    # Generate mock targets - perturbation of genes within cohorts
    for treat in [
        "CARBOPLATIN",
        "FLUOROURACIL",
        "RADIOTHERAPY",
    ]:  # "rtprescribeddose"]:
        indices = np.arange(len(treatment_df))
        MOCKtreat = np.zeros(len(treatment_df), dtype=int)
        for group in np.unique(X.group):
            MOCKtreat[X.group == group] = treatment_df[treat][
                np.random.choice(
                    indices[X.group == group],
                    size=np.sum(X.group == group),
                    replace=False,
                )
            ]
        treatment_df[f"MOCK{treat}"] = MOCKtreat

    # Save targets file
    treatment_df.to_csv(targets_file, sep="\t")
