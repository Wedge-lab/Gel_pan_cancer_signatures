import os
import re
import sys

import numpy as np
import pandas as pd

from signatures.config import load_environment
from signatures.utils import regexSearch

load_environment()
DATA_DIR = os.getenv("DATA_DIR")

if __name__ == "__main__":
    samples_file, targets_file = sys.argv[1:]

    # Load samples file
    sample_df = pd.read_csv(samples_file, sep="\t")
    sample_df["tumour_sample_platekey"] = sample_df.sample_id.map(
        lambda x: "_".join(x.split("_")[1:3])
    )

    # Crossmatch between GMCs and NCRAS
    tumour_xm = pd.read_csv(f"{DATA_DIR}/GEL_NCRAS_XM.tsv", sep="\t")
    # Drop duplicate platekeys
    tumour_xm["grade_val"] = np.where(
        (tumour_xm.grade.astype(str) == "nan") | (tumour_xm.grade.astype(str) == "GX"),
        1,
        0,
    )
    tumour_xm["stage_val"] = np.where(tumour_xm.stage_best.astype(str) == "nan", 1, 0)
    tumour_xm["random"] = np.random.rand(len(tumour_xm))
    tumour_xm = tumour_xm.sort_values(
        ["days_to_sampling", "grade_val", "stage_val", "random"]
    ).drop_duplicates("tumour_sample_platekey")
    tumour_xm = tumour_xm[["tumour_sample_platekey", "tumour_pseudo_id"]]

    # Variables of interest
    VOI = [
        "grade",
        "stage_best",  #'tumoursize',
        "er_status",
        "pr_status",
        "her2_status",
        #'nodesexcised', 'nodesinvolved', 'laterality',
        #'behaviour_coded_desc',
        "npi",  # Notting Prognostic Index - size, stage and lymph node involvement
        "dukes",
        "figo",  # Stage
        #'t_best', # Size and stage of primary
        #'n_best', # Lymph node metastasis
        #'m_best', # Distant metastases
        "gleason_combined",  # Grade
        #'basisofdiagnosis',
        "clarks",
        "breslow",  # Melanoma thickness and depth
        #'multifocal' # Breast with multiple tumours from one primary
    ]

    # Load NCRAS data
    av_tumour = pd.read_csv(
        f"{DATA_DIR}/V11_reheadered/av_tumour.tsv",
        sep="\t",
        usecols=[
            "tumour_pseudo_id",
        ]
        + VOI
        + [
            "site_coded_desc",
        ],
    )
    sample_av_df = pd.merge(
        sample_df[["sample_id", "tumour_sample_platekey", "group"]],
        pd.merge(
            tumour_xm[["tumour_sample_platekey", "tumour_pseudo_id"]],
            av_tumour,
            on="tumour_pseudo_id",
            how="left",
        ),
        on="tumour_sample_platekey",
        how="left",
    )

    # Categorical variables
    sample_av_df["stage_best"] = sample_av_df.stage_best.astype(str).map(
        lambda x: regexSearch("([0-9]+)[A-Z]*", x, 1)
        if re.match("([0-9]+)[A-Z]*", x)
        else regexSearch("([A-Z])", x, 1)
        if re.match("([A-Z])", x)
        else "nan"
    )
    sample_av_df["figo"] = sample_av_df.figo.astype(str).map(
        lambda x: regexSearch("([0-9]+)[a-z]*", x, 1)
        if re.match("([0-9]+)[a-z]*", x)
        else regexSearch("([IV]+)[A-Z]+", x, 1)
        if re.match("([IV]+)[A-Z]+", x)
        else "nan"
    )
    VOI_mapping = {
        "grade": {"G1": 0, "G2": 0, "G3": 1, "G4": 1, "GX": np.nan},
        "stage_best": {
            "0": 0,
            "1": 0,
            "2": 0,
            "3": 1,
            "4": 1,
            "6": np.nan,
            "A": np.nan,
            "B": np.nan,
            "C": np.nan,
            "U": np.nan,
            "X": np.nan,
        },
        "dukes": {"A": 0, "B": 0, "C": 1, "C1": 1, "C2": 1, "D": 1, "X": np.nan},
        "figo": {
            "1": 0,
            "2": 0,
            "3": 1,
            "4": 1,
            "I": 0,
            "III": 1,
            "IV": 1,
            "nan": np.nan,
        },
        "t_best": {"0": 0, "1": 0, "2": 0, "3": 1, "4": 1, "nan": np.nan},
        "m_best": {"0": 0, "1": 0, "2": 0, "3": 1, "4": 1, "nan": np.nan},
        "n_best": {"0": 0, "1": 0, "2": 0, "3": 1, "4": 1, "nan": np.nan},
        "er_status": {"B": np.nan, "N": 0, "P": 1, "Pm": 1, "U": np.nan, "X": np.nan},
        "pr_status": {"B": np.nan, "N": 0, "P": 1, "Pm": 1, "U": np.nan, "X": np.nan},
        "her2_status": {"B": np.nan, "N": 0, "P": 1, "Pm": 1, "U": np.nan, "X": np.nan},
    }
    sample_av_df.replace(VOI_mapping, inplace=True)

    # Triple negative breast cancer
    triple_negative = np.sum(
        np.array(sample_av_df[["er_status", "pr_status", "her2_status"]]), axis=1
    )
    triple_negative[triple_negative > 0] = 1
    triple_negative = 1 - triple_negative
    sample_av_df["triple_negative"] = triple_negative

    # Continuous variables
    sample_av_df["breslow"] = sample_av_df["breslow"].map(
        lambda x: float(x)
        if not "-" in str(x)
        else np.sum(np.array(str(x).split("-")).astype(float)) / 2
    )
    binary_split = {
        "tumoursize": 30.0,
        "npi": 4.4,
        "gleason_combined": 7.0,
        "clarks": 4.0,
        "breslow": 3.125,
    }
    for key in ["npi", "gleason_combined", "clarks", "breslow"]:  # 'tumoursize',
        sample_av_df[key] = np.where(
            np.isnan(sample_av_df[key]), np.nan, sample_av_df[key] > binary_split[key]
        )

    print(sample_av_df.head())

    # CRC and Testis sites
    crc_colon_rectum = sample_av_df.site_coded_desc.astype(str).map(
        lambda x: 1 if "colon" in x.lower() else 0 if "rectum" in x.lower() else np.nan
    )
    crc_colon_rectum[sample_av_df.group != "ColoRect-AdenoCA"] = np.nan
    sample_av_df["crc_colon_rectum"] = crc_colon_rectum

    # Seminoma, teritoma?
    testis_desc = sample_av_df.site_coded_desc.astype(str).map(
        lambda x: 1
        if "descended" in x.lower()
        else 0
        if "undescended" in x.lower()
        else np.nan
    )
    # testis_desc[sample_av_df.group!="Testis-GCT"] = np.nan
    sample_av_df["testis_desc"] = testis_desc

    # Save sample df
    sample_av_df[
        [
            "sample_id",
        ]
        + VOI
        + ["triple_negative", "crc_colon_rectum", "testis_desc"]
    ].to_csv(targets_file, sep="\t", index=False)
