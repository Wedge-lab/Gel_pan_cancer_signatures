# Import standard modules
import os
import sys

import numpy as np
import pandas as pd
from dotenv import load_dotenv

from signatures.associations.sampleCuration.pan_cancer_sample import getSamples

load_dotenv()
DATA_DIR = os.getenv("DATA_DIR")


def crossmatchGELNCRAS(sample_platekeys=None):
    # Cancer analysis table
    cancer_analysis_df = (
        pd.read_csv(f"{DATA_DIR}/V11_reheadered/cancer_analysis.tsv", sep="\t")
        .reset_index()
        .rename({"index": "cancer_analysis_index"}, axis=1)
    )
    cancer_analysis_df["sample_platekey"] = (
        "tumo"
        + cancer_analysis_df.tumour_sample_platekey
        + "_norm"
        + cancer_analysis_df.germline_sample_platekey
    )
    cancer_analysis_df.set_index("sample_platekey", inplace=True)

    # Subset data with sample platekeys
    if sample_platekeys is not None:
        cancer_analysis_df = cancer_analysis_df.loc[sample_platekeys]
    cancer_analysis_df.reset_index(inplace=True)

    # AV tumour table
    av_tumour_df = (
        pd.read_csv(f"{DATA_DIR}/V11_reheadered/av_tumour.tsv", sep="\t")
        .reset_index()
        .rename({"index": "av_tumour_index"}, axis=1)
    )

    # Crossmatch on participant_id
    av_cancer_xm = pd.merge(
        av_tumour_df[
            [
                "av_tumour_index",
                "participant_id",
                "diagnosisdatebest",
                "tumour_pseudo_id",
            ]
        ],
        cancer_analysis_df[
            [
                "cancer_analysis_index",
                "participant_id",
                "sample_platekey",
                "tumour_clinical_sample_time",
            ]
        ],
        on="participant_id",
        how="inner",
    )

    # Get diagnosis - sampling date difference
    for key in ["diagnosisdatebest", "tumour_clinical_sample_time"]:
        av_cancer_xm[key] = pd.to_datetime(av_cancer_xm[key], format="%Y-%m-%d")
    av_cancer_xm["sample_delay"] = (
        av_cancer_xm["tumour_clinical_sample_time"] - av_cancer_xm["diagnosisdatebest"]
    ).dt.total_seconds() / (3600 * 24)

    # Remove matches with longest delays
    # av_cancer_xm = av_cancer_xm.sort_values(['tumour_pseudo_id','sample_delay']).drop_duplicates('participant_id')
    av_cancer_xm = av_cancer_xm.sort_values("sample_delay").drop_duplicates(
        "sample_platekey"
    )

    return av_cancer_xm[
        [
            "participant_id",
            "sample_platekey",
            "tumour_pseudo_id",
            "diagnosisdatebest",
            "tumour_clinical_sample_time",
            "sample_delay",
            "cancer_analysis_index",
            "av_tumour_index",
        ]
    ]


def getTreatment(sample_df, base=""):
    # Get NCRAS - GEL crossmatch
    print("NCRAS")
    sample_df["sample_platekey"] = (
        "tumo"
        + sample_df.tumour_sample_platekey
        + "_norm"
        + sample_df.germline_sample_platekey
    )
    sample_df = pd.merge(
        sample_df,
        crossmatchGELNCRAS(sample_platekeys=sample_df.sample_platekey)[
            ["sample_platekey", "tumour_pseudo_id"]
        ],
        how="inner",
        on="sample_platekey",
    )
    print(len(sample_df))

    # ICD codes
    print("ICD codes")
    icd10 = pd.read_csv(
        f"{DATA_DIR}/ICD10/icd102016syst_codes.txt",
        sep=";",
        header=None,
        usecols=[4, 7, 8],
    )
    # icd10 = icd10[icd10[8].map(lambda x: 'breast' in x.lower())]
    icd10 = pd.concat(
        [
            icd10[[4, 8]].rename(columns={4: "code", 8: "desc"}),
            icd10[[7, 8]].rename(columns={7: "code", 8: "desc"}),
        ]
    ).drop_duplicates(("code", "desc"))

    # SACT table
    print("SACT")
    sact_df = pd.read_csv(f"{DATA_DIR}/V11_reheadered/sact.tsv", sep="\t")
    sact_df = pd.merge(
        sact_df,
        icd10,
        how="left",
        left_on="primary_diagnosis",
        right_on="code",
        suffixes=("", "_icd"),
    )
    sact_df = sact_df[sact_df.desc.astype(str) != "nan"]

    # Merge sample with treatment dataframe on participant_id
    treatment_df = pd.merge(
        sample_df, sact_df, on="tumour_pseudo_id", how="left", suffixes=("", "_sact")
    )
    for key in [
        "date_decision_to_treat",
        "start_date_of_regimen",
        "start_date_of_cycle",
        "administration_date",
        "date_of_final_treatment",
        "day_sampling",
    ]:
        treatment_df[key] = pd.to_datetime(treatment_df[key], format="%Y-%m-%d")

    # Only consider treatments after sampling
    treatment_df = treatment_df[
        treatment_df.day_sampling > treatment_df.start_date_of_cycle
    ]

    # Get duration of treatment (if final date not given we use sampling date)
    treatment_df["treatment_duration"] = (
        (
            treatment_df[["day_sampling", "date_of_final_treatment"]].min(axis=1)
            - treatment_df["start_date_of_regimen"]
        ).dt.total_seconds()
        / (3600 * 24)
    ).astype(int)
    print(len(treatment_df))

    # Get unique treatments
    treatment_map = {
        "AC": "DOXORUBICIN + CYCLOPHOSPHAMIDE",
        "BEP": "BLEOMYCIN + ETOPOSIDE + CISPLATIN",
        "CHOP": "CYCLOPHOSPHAMIDE + DOXORUBICIN + VINCRISTINE + PREDNISOLON",
        "CHOP_R": "CYCLOPHOSPHAMIDE + DOXORUBICIN + VINCRISTINE + PREDNISOLON + RITUXIMAB",
        "CTD": "THALIDOMIDE + CYCLOPHOSPHAMIDE + DEXAMETHASONE",
        "CYCLOPHOSPHAMIDE_HIGH_DOSE": "CYCLOPHOSPHAMIDE",
        "CYTARABINE_HIGH_DOSE": "CYTARABINE",
        "CARBO": "CARBOPLATIN",
        "GEMCARBO": "GEMCITABINE + CARBOPLATIN",
        "TRASTUZUMAB EMTANSINE": "TRASTUZUMAB + EMTANSINE",
        "NAB-PACLITAXEL": "ALBUMIN + PACLITAXEL",
        "FEC": "FLUOROURACIL + EPIRUBICIN + CYCLOPHOSPHAMIDE",
        "FEC 100": "FLUOROURACIL + EPIRUBICIN + CYCLOPHOSPHAMIDE",
        "FEC 60 OR 75": "FLUOROURACIL + EPIRUBICIN + CYCLOPHOSPHAMIDE",
        "EC": "EPIRUBICIN + CYCLOPHOSPHAMIDE",
        "EOX": "EPIRUBICIN + OXALIPLATIN + CAPECITABINE",
        "ECX": "EPIRUBICIN + CISPLATIN + CAPECITABINE",
        "ECF": "EPIRUBICIN + CISPLATIN + FLUOROURACIL",
        "MDG": "FLUOROURACIL + FOLINIC ACID",
        "TCH": "TAXOTERE + CARBOPLATIN + HERCEPTIN",
        "PCV": "PROCARBAZINE + LOMUZTINE + VINCRISTINE",
        "MITOMYCIN_INTRAVESICULAR": "MITOMYCIN",
        "METHOTREXATE_HIGH_DOSE": "METHOTREXATE",
        "METHOTREXATE_INTRATHECAL": "METHOTREXATE",
        "HCX": "TRASTUZUMAB + CISPLATIN + CAPECITABINE",
        "FU": "FLUOROURACIL",
    }
    treatment_df.analysis_group = treatment_df.analysis_group.map(
        lambda x: " + ".join(
            [
                treatment_map[t.strip()]
                if t.strip() in treatment_map.keys()
                else t.strip()
                for t in x.split("+")
            ]
        )
    )

    # Create grid of treatment durations
    treatments = np.unique(
        [
            treatment
            for treatments in np.unique(treatment_df["analysis_group"])
            for treatment in treatments.split(" + ")
        ]
    )
    treatment_grid = sample_df.set_index("sample_platekey")
    for treatment in treatments:
        subset = treatment_df.analysis_group.map(lambda x: treatment in x.split(" + "))
        sample_platekeys = treatment_df["sample_platekey"][subset]

        treatment_grid[treatment] = np.zeros(len(treatment_grid), dtype=int)
        treatment_grid[treatment].loc[sample_platekeys] = np.array(
            treatment_df["treatment_duration"][subset]
        )
    print(treatment_grid.shape)

    print("RTDS")
    # Radio therapy data
    rtds_df = pd.read_csv(f"{DATA_DIR}/V11_reheadered/rtds.tsv", sep="\t")
    rtds_df = pd.merge(
        rtds_df,
        icd10,
        how="left",
        left_on="radiotherapydiagnosisicd",
        right_on="code",
        suffixes=("", "_icd"),
    )
    rtds_df = rtds_df[rtds_df.desc.astype(str) != "nan"]

    # Merge with sample_df
    radio_df = pd.merge(sample_df, rtds_df, on="participant_id", how="left")

    for key in [
        "day_sampling",
        "apptdate",
        "decisiontotreatdate",
        "earliestclinappropriatedate",
        "treatmentstartdate",
        "proceduredate",
    ]:
        radio_df[key] = pd.to_datetime(radio_df[key], format="%Y-%m-%d")

    # Only consider treatments after sampling
    radio_df = radio_df[radio_df.day_sampling > radio_df.treatmentstartdate]
    print(len(radio_df))

    participant_dosage = (
        radio_df[["participant_id", "rtprescribeddose"]].groupby("participant_id").sum()
    )
    print(participant_dosage.head())
    print(treatment_grid.head())

    treatment_grid = pd.merge(
        treatment_grid,
        participant_dosage,
        how="left",
        left_on="participant_id",
        right_index=True,
    )
    treatment_grid.rtprescribeddose[np.isnan(treatment_grid.rtprescribeddose)] = 0.0

    return treatment_grid


def getTreatment2(sample_df, base=""):
    # File of treatments for samples
    treatment_prior_file = f"{DATA_DIR}/processedClinicalData/treatmentPriorRecords.tsv"
    treatment_prior = pd.read_csv(treatment_prior_file, sep="\t")
    treatment_df = pd.merge(
        sample_df,
        treatment_prior.drop(
            [
                "participant_id",
                "germline_sample_platekey",
                "tumour_normal_id",
                "date_sampling",
            ],
            axis=1,
        ),
        on="tumour_sample_platekey",
        how="inner",
    )
    treatment_df["treatment"] = np.where(
        treatment_df.treatment_type == "RADIOTHERAPY", "RADIOTHERAPY", treatment_df.drug
    )

    # Generate treatment grid
    treatment_df["exposed"] = 1
    treatment_grid = (
        pd.pivot_table(
            treatment_df,
            index="tumour_sample_platekey",
            columns="treatment",
            values="exposed",
        )
        .fillna(0)
        .astype(int)
    )

    # Merge with sample_df and retrieve all samples
    dtypes = dict(
        zip(
            treatment_grid.keys(),
            [
                np.int32,
            ]
            * len(treatment_grid.keys()),
        )
    )
    treatment_grid = (
        pd.merge(
            sample_df,
            treatment_grid,
            left_on="tumour_sample_platekey",
            right_index=True,
            how="left",
        )
        .set_index("tumour_sample_platekey")
        .fillna(0)
        .astype(dtypes)
    )

    return treatment_grid


if __name__ == "__main__":
    # filename to save covariates into
    filename_output = sys.argv[1]

    # GET DATA - include hypermutated samples
    sample_df = getSamples(
        PC=True,
        signatures=True,
        primary=True,
        ethnicity=False,
        hyper=False,
        ER=False,
        MSI=False,
    )

    # Confounding variables
    variables = (
        ["log_age", "is_female"] + [f"pc{i}" for i in range(1, 4)]
    )  # 'log_age', ['age','log_age','logit_purity','is_female'] + [f"pc{i}" for i in range(1,6)]

    # Get non-nan values
    [
        print(key, np.sum(~np.isnan(X[key])))
        for key in variables
        if X[key].dtype in [float, int]
    ]
    cuts = np.prod(
        np.array(
            [~np.isnan(X[key]) for key in variables if X[key].dtype in [float, int]]
        ),
        axis=0,
    ).astype(bool)
    print(f"Non-nan: {np.sum(cuts)}")

    # Apply cuts
    sample_df = sample_df[cuts]

    # Save
    sample_df.rename(
        columns={"hyper_status": "hypermutation", "tumour_tissue": "category"}
    )
    sample_df[
        ["tumour_sample_platekey", "hypermutation", "category"] + variables
    ].to_csv(filename_output, index=False, sep="\t")

    # Get and save treatment grid
    treatment_grid = getTreatment(
        sample_df[
            [
                "participant_id",
                "tumour_sample_platekey",
                "germline_sample_platekey",
                "day_sampling",
            ]
        ]
    )
    treatment_grid.to_csv(f"./treatment_duration.tsv", sep="\t")
