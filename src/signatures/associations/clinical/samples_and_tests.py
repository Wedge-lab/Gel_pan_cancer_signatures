import os
import sys

import pandas as pd
from dotenv import load_dotenv

load_dotenv()
DATA_DIR = os.getenv("DATA_DIR")
SAMPLE_LIST = os.getenv("SAMPLE_LIST")

if __name__ == "__main__":
    samples_file, activities_file = sys.argv[1:]

    # Initialise test dataframe
    test_df = pd.DataFrame(columns=["group", "target", "type"])

    # Load samples file
    sample_df = pd.read_csv(samples_file, sep="\t")

    # Get tumour group
    sample_list_df = pd.read_csv(
        SAMPLE_LIST,
        usecols=[
            "participant_id",
            "tumour_sample_platekey",
            "germline_sample_platekey",
            "signature_extraction_group",
        ],
        delim_whitespace=True,
    ).rename(columns={"age_sampling": "age"})
    sample_list_df["sample_id"] = (
        sample_list_df.participant_id.astype(str)
        + "_"
        + sample_list_df.tumour_sample_platekey
        + "_"
        + sample_list_df.germline_sample_platekey
    )
    sample_df = pd.merge(
        sample_df,
        sample_list_df[
            ["sample_id", "tumour_sample_platekey", "signature_extraction_group"]
        ],
        on="sample_id",
        how="left",
    )

    # Get tumour groups
    subgroup_size = sample_df.groupby(["signature_extraction_group", "group"]).size()
    for subgroup in subgroup_size.subgroup:
        sample_df[subgroup] = (sample_df.subgroup == subgroup).astype(int)
    subgroup_size.reset_index().rename(
        {"signature_extraction_group": "group", "group": "target"}
    )
    subgroup_size["type"] = "binomial"

    # Get: grade, stage, ...
    tumour_xm = pd.read_csv(f"{DATA_DIR}/GEL_NCRAS_XM.tsv", sep="\t")
    VOI = [
        "grade",
        "tumoursize",
        "stage_best",
        "er_status",
        "pr_status",
        "her2_status",
        "nodesexcised",
        "nodesinvolved",
        "laterality",
        "behaviour_coded_desc",
        "npi",  # Notting Prognostic Index - size, stage and lymph node involvement
        "dukes",
        "figo",  # Stage
        "t_best",  # Size and stage of primary
        "n_best",  # Lymph node metastasis
        "m_best",  # Distant metastases
        "gleason_combined",  # Grade
        "basisofdiagnosis",
        "clarks",
        "breslow",  # Melanoma thickness and depth
        "multifocal",  # Breast with multiple tumours from one primary
    ]
    av_tumour = pd.read_csv(
        f"{DATA_DIR}/V11_reheadered/av_tumour.tsv",
        sep="\t",
        usecols=[
            "tumour_pseudo_id",
        ]
        + VOI,
    )

    sample_df = pd.merge(
        sample_df[["tumour_sample_platekey", "extraction_group"]],
        pd.merge(
            tumour_xm[["tumour_sample_platekey", "tumour_pseudo_id"]],
            av_tumour,
            on="tumour_pseudo_id",
            how="left",
        ),
        on="tumour_sample_platekey",
        how="left",
    )
