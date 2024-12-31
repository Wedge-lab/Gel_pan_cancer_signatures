"Crossmatch GEL with NCRAS using curated sample table from Alex Cornish and Dan Chubb"

import os

import numpy as np
import pandas as pd
from dotenv import load_dotenv

load_dotenv()
DATA_DIR = os.getenv("DATA_DIR")
RESULT_DIR = os.getenv("RESULT_DIR")

group_mapping = {
    "BileDuct-AdenoCA": ["adenocarcinoma"],  # cholangiocarcinoma
    "Bladder-TCC": ["transitional cell"],
    "Breast-DuctalCA": ["duct"],
    "Breast-LobularCA": ["lobular"],
    "CNS-Astro": ["astrocytoma"],
    "CNS-GBM-IDHmut": ["idh-mutant"],
    "CNS-GBM-IDHwt": ["idh-wildtype"],
    "CNS-Menin": ["meningioma"],
    "CNS-Oligo": ["oligodendroglioma"],
    "ColoRect-AdenoCA": ["adenocarcinoma"],
    "Connective-Chondro": ["chondrosarcoma"],
    "Connective-Leiomyo": ["leiomyosarcoma"],
    "Connective-Liposarc": ["liposarcoma"],
    "Connective-Myxofibro": ["myxo", "fibro"],
    "Connective-Osteosarc": ["osteosarcoma"],
    "Connective-SCS": ["spindle cell"],
    "Connective-SS": ["synovial"],
    "Eso-AdenoCA": ["adenocarcinoma"],
    "Haem-ALL": ["lymphoblastic leuk"],
    "Haem-AML": [
        "acute myeloid leuk",
        "acute monocytic leuk",
        "acute myelomonocytic leuk",
        "acute promyelocytic leuk",
    ],
    "Haem-CLL": ["chronic"],
    "Haem-CML": ["chronic myelogenous leuk", "chronic myelogenous leuk"],
    "Haem-MM": ["multiple myeloma"],
    "Haem-MPN": ["myeloproliferative neoplasm"],
    "HeadNeck-SCC": ["squamous cell"],
    "Kidney-CCRCC": ["adenocarcinoma"],
    "Kidney-ChRCC": ["chromophobe"],
    "Kidney-PRCC": ["adenocarcinoma", "papillary"],
    "Liver-HCC": ["hepatocellular"],
    "Lung-AdenoCA": ["adenocarcinoma"],
    "Lung-LargeCell": ["large cell"],
    "Lung-SCC": ["squamous cell"],
    "Lung-SmallCell": ["small cell"],
    "Mes-Mesothelioma": [],
    "Ovary-AdenoCA": ["adenocarcinoma"],
    "Panc-AdenoCA": ["adenocarcinoma"],
    "Prost-AdenoCA": ["adenocarcinoma"],
    "Skin-Melanoma": ["melanoma"],
    "Stomach-AdenoCA": ["adenocarcinoma"],
    "Testis-GCT": ["germ cell", "seminoma"],
    "Uterus-AdenoCA": ["adenocarcinoma"],
}

tissue_mapping = {
    "BileDuct": ["bile duct"],
    "Bladder": ["bladder", "ureter", "urethra"],
    "Breast": ["breast", "nipple"],
    "CNS": [
        "brain",
        "cerebellum",
        "cerebrum",
        "frontal lobe",
        "occipital lobe",
        "parietal lobe",
        "temporal lobe",
    ],
    "ColoRect": ["colon", "caecum", "appendix", "rectum"],
    "Connective": [
        "bone",
        "connective",
        "femur",
        "humerus",
        "vertebra",
        "mandible",
        "finger",
        "myometrium",
        "pelvis",
        "peritoneum",
        "sternum",
        "scapula",
        "subcutaneous tissue",
        "thigh",
    ],
    "Eso": ["oesophagus", "oesophageal"],
    "Haem": ["bone marrow", "leukemia", "leukaemia"],
    "HeadNeck": [
        "tongue",
        "mouth",
        "cheek",
        "glottis",
        "gum",
        "palate",
        "pharynx",
        "tonsil",
        "sinus",
        "lip",
        "nasal cavity",
        "nasopharynx",
        "oropharynx",
        "molar",
    ],
    "Kidney": ["kidney"],
    "Liver": ["liver"],
    "Lung": ["lung", "bronchus"],
    "Mes": ["pleura"],
    "Ovary": ["ovary", "fallopian tube"],
    "Panc": ["pancreas", "pancreatic duct"],
    "Prost": ["prostate"],
    "Skin": ["skin", "melanoma", "ear", "eyelid"],
    "Stomach": ["stomach", "gastric antrum", "cardia", "pylorus"],
    "Testis": ["testis"],
    "Uterus": ["uterus", "endometrium"],
}


def crossmatchGelNcras():
    # Import curated sample list from Alex and Dan
    sample_df = pd.read_csv(
        f"{RESULT_DIR}/sample_lists/sample_list_2021_06_29.tsv", delim_whitespace=True
    ).rename(columns={"age_sampling": "age"})

    # Import NCRAS data and crossmatch on tumour_pseudo_id
    sact = pd.read_csv(
        f"{DATA_DIR}/V11/V11_reheadered/sact.tsv",
        sep="\t",
        usecols=["participant_id", "tumour_pseudo_id", "start_date_of_cycle"],
    )
    av_treatment = pd.read_csv(
        f"{DATA_DIR}/V11/V11_reheadered/av_treatment.tsv",
        sep="\t",
        usecols=["participant_id", "tumour_pseudo_id", "eventdate", "eventdesc"],
    )
    av_tumour = pd.read_csv(
        f"{DATA_DIR}/V11/V11_reheadered/av_tumour.tsv",
        sep="\t",
        usecols=[
            "participant_id",
            "tumour_pseudo_id",
            "diagnosisdatebest",
            "date_first_surgery",
        ],
    )
    # Outer join all NCRAS tables on tumour_pseudo_id
    av_xm = pd.merge(
        av_tumour,
        pd.merge(
            sact[["tumour_pseudo_id", "start_date_of_cycle"]],
            av_treatment[["tumour_pseudo_id", "eventdate", "eventdesc"]],
            on="tumour_pseudo_id",
            how="outer",
        ),
        on="tumour_pseudo_id",
        how="outer",
    )

    # Crossmatch GEL and NCRAS on participant_id and day of diagnosis
    av_sample_xm = pd.merge(
        sample_df,
        av_xm,
        left_on=["participant_id", "day_diagnosis"],
        right_on=["participant_id", "diagnosisdatebest"],
        how="left",
    )

    # Get key NCRAS dates
    tumour_xm = pd.DataFrame()
    av_dates = [
        "date_first_surgery",
        "diagnosisdatebest",
        "start_date_of_cycle",
        "eventdate",
    ]
    for av_date in av_dates:
        tmp = av_sample_xm[
            [
                "participant_id",
                "tumour_sample_platekey",
                "tumour_pseudo_id",
                "day_sampling",
                "eventdesc",
                av_date,
            ]
        ].rename({av_date: "av_date"}, axis=1)
        tmp["date_source"] = av_date
        tumour_xm = pd.concat((tumour_xm, tmp.copy()))

    # Keep crossmatches with date difference <28 days
    tumour_xm["days_to_sampling"] = (
        pd.to_datetime(tumour_xm.day_sampling) - pd.to_datetime(tumour_xm["av_date"])
    ).dt.days
    tumour_xm.drop_duplicates(
        ["tumour_sample_platekey", "tumour_pseudo_id", "date_source", "av_date"],
        inplace=True,
    )
    tumour_xm["abs_days_to_sampling"] = np.abs(tumour_xm.days_to_sampling)
    # Drop any crossmatches with date difference greater than 28d
    # Drop duplicates of platekey and pseudo_id
    tumour_xm = (
        tumour_xm[tumour_xm.abs_days_to_sampling < 28]
        .sort_values("abs_days_to_sampling")
        .drop_duplicates(["tumour_sample_platekey", "tumour_pseudo_id"])
    )
    len(tumour_xm), len(np.unique(tumour_xm.tumour_sample_platekey))

    # Get histology, grade and stage from NCRAS
    av_tumour = pd.read_csv(
        f"{DATA_DIR}/V11/V11_reheadered/av_tumour.tsv",
        sep="\t",
        usecols=[
            "tumour_pseudo_id",
            "histology_coded_desc",
            "site_coded_desc",
            "grade",
            "stage_best",
        ],
    )
    sample_av_df = pd.merge(
        sample_df[["tumour_sample_platekey", "tumour_group"]],
        pd.merge(
            tumour_xm[
                [
                    "tumour_sample_platekey",
                    "tumour_pseudo_id",
                    "date_source",
                    "days_to_sampling",
                    "eventdesc",
                ]
            ],
            av_tumour,
            on="tumour_pseudo_id",
            how="left",
        ),
        on="tumour_sample_platekey",
        how="left",
    )
    sample_av_df["tumour_tissue"] = sample_av_df["tumour_group"].map(
        lambda x: x.split("-")[0]
    )

    # Match histology to tumour group
    sample_av_df["histology_clean"] = sample_av_df["histology_coded_desc"].map(
        lambda x: str(x).lower()
    )
    sample_av_df["site_clean"] = sample_av_df["site_coded_desc"].map(
        lambda x: str(x).lower()
    )
    sample_av_df["group_match"] = sample_av_df[
        ["tumour_group", "histology_clean"]
    ].apply(
        lambda x: np.sum(
            [
                (y in x["histology_clean"]) | (y == "any")
                for y in group_mapping[x["tumour_group"]] + ["nan"]
            ]
        )
        > 0,
        axis=1,
    )
    # Match site to tumour tissue
    sample_av_df["tissue_match"] = sample_av_df[["tumour_tissue", "site_clean"]].apply(
        lambda x: np.sum(
            [y in x["site_clean"] for y in tissue_mapping[x["tumour_tissue"]]]
        )
        > 0,
        axis=1,
    )

    # Subset using tissue and group matches
    subset = (sample_av_df.group_match) & (sample_av_df.tissue_match)
    sample_curated = sample_av_df[subset]

    return sample_curated  # [['tumour_sample_platekey', 'tumour_pseudo_id']]


def unique_crossmatch(sample_df):
    # GEL-NCRAS crossmatch
    GEL_NCRAS_XM = pd.read_csv(
        f"{RESULT_DIR}/sampleCuration/GEL_NCRAS_XM.tsv", sep="\t"
    )
    # NCRAS data
    av_tumour = pd.read_csv(f"{DATA_DIR}/V11/V11_reheadered/av_tumour.tsv", sep="\t")

    sample_xm_df = pd.merge(
        sample_df,
        pd.merge(
            GEL_NCRAS_XM[["tumour_pseudo_id", "tumour_sample_platekey"]],
            av_tumour,
            how="inner",
            on="tumour_pseudo_id",
        ),
        how="left",
        on=["participant_id", "tumour_sample_platekey"],
    )

    av_keys = [
        "sex",
        "age_sampling",
        "stage_best",
        "grade",
        "purity",
        "er_status",
        "pr_status",
        "her2_status",
        "npi",
        "dukes",
        "figo",
        "gleason_combined",
        "clarks",
        "breslow",
    ]

    # Clean duplicated samples
    duplicate_keys = np.unique(sample_xm_df.tumour_sample_platekey)[
        np.unique(sample_xm_df.tumour_sample_platekey, return_counts=True)[1] > 1
    ]

    sample_xm_df.set_index("tumour_sample_platekey", inplace=True)
    for sample in duplicate_keys:
        for key in av_keys:
            # Conditional if entries are different, otherwise we don't need to worry
            if len(np.unique(sample_xm_df.loc[sample][key].astype(str))) > 1:
                non_nan = sample_xm_df.loc[sample][key][
                    sample_xm_df.loc[sample][key].astype(str) != "nan"
                ]

                # If only one non-nan: use that value
                if len(non_nan) == 1:
                    sample_xm_df.loc[sample, key] = non_nan[0]
                else:  # Set as nan due to conflict
                    sample_xm_df.loc[sample, key] = np.nan
    sample_xm_df.reset_index(inplace=True)

    # Drop the duplicated samples
    sample_xm_df.drop_duplicates(
        [
            "tumour_sample_platekey",
        ]
        + av_keys,
        inplace=True,
    )

    return sample_xm_df


if __name__ == "__main__":
    output_dir = f"{RESULT_DIR}/sampleCuration"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    tumour_xm = crossmatchGelNcras()

    tumour_xm.to_csv(f"{output_dir}/GEL_NCRAS_XM.tsv", sep="\t", index=False)
