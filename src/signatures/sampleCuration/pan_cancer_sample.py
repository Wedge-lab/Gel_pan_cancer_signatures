# Import standard modules
import os
import sys

import numpy as np
import pandas as pd
from dotenv import load_dotenv

load_dotenv()
RESULT_DIR = os.getenv("RESULT_DIR")
DATA_DIR = os.getenv("DATA_DIR")
SAMPLE_LIST = os.getenv("SAMPLE_LIST")


def getSamples(
    PC=True,
    ethnicity=False,
    hyper=False,
    signatures=True,
    ER=False,
    MSI=False,
    primary=True,
):
    # Import sample data
    sample_keys = [
        "tumour_sample_platekey",
        "germline_sample_platekey",
        "age_sampling",
        "tumour_group",
        "tumour_type",
        "sex",
        "purity",
        "participant_id",
        "pseudonymous_id",
        "day_diagnosis",
        "day_sampling",
        "day_death",
        "day_last_followup",
    ]
    sample_df = pd.read_csv(
        f"{SAMPLE_LIST}",
        usecols=sample_keys,
        delim_whitespace=True,
    ).rename(columns={"age_sampling": "age"})
    sample_df["sample_platekey"] = (
        "tumo"
        + sample_df.tumour_sample_platekey
        + "_norm"
        + sample_df.germline_sample_platekey
    )
    print(f"pan-cancer sample list {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    # Aggregate GVCF sample stats
    if PC:
        agg_keys = ["platekey", "pc1", "pc2", "pc3", "pc4", "pc5"]
        agg_df = pd.read_csv(
            f"{DATA_DIR}/aggV2/additional_data/aggregate_gvcf_sample_stats/aggregate_gvcf_sample_stats_v10_78195.tsv",
            usecols=agg_keys,
            sep="\t",
        )
        sample_df = pd.merge(
            sample_df,
            agg_df,
            how="inner",
            left_on="germline_sample_platekey",
            right_on="platekey",
        )
        print(f"PC data {len(sample_df)}")
        print(
            f"Total: {len(sample_df)}, \
        Unique participant: {len(np.unique(sample_df['participant_id']))}, \
        Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
        )

    # Participant information
    if ethnicity:
        participant_df = pd.read_csv(
            f"{DATA_DIR}/prs_signatures/data/participant_2022-03-21_15-52-34.tsv",
            sep="\t",
        )
        sample_df = pd.merge(
            sample_df,
            participant_df,
            how="inner",
            left_on="participant_id",
            right_on="Participant Id",
        )
        sample_df["white_european"] = sample_df["Participant Ethnic Category"].map(
            lambda x: x
            in ["White: Any other White background", "White: British", "White: Irish"]
        )
        sample_df["white_european"] = (
            sample_df["white_european"]
            & (
                np.abs(
                    (
                        sample_df["pc1"]
                        - np.mean(sample_df["pc1"][sample_df["white_european"]])
                    )
                    / np.std(sample_df["pc1"][sample_df["white_european"]])
                )
                < 6
            )
            & (
                np.abs(
                    (
                        sample_df["pc2"]
                        - np.mean(sample_df["pc2"][sample_df["white_european"]])
                    )
                    / np.std(sample_df["pc2"][sample_df["white_european"]])
                )
                < 6
            )
        )
        print(f"Participant (ethnicity) data: {len(sample_df)}")
        print(
            f"White european subset: {np.sum(sample_df.white_european)}/{len(sample_df)} - not applied"
        )
        print(
            f"Total: {len(sample_df)}, \
        Unique participant: {len(np.unique(sample_df['participant_id']))}, \
        Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
        )
        sample_df = sample_df[sample_df.white_european]

    # Signature data
    if signatures:
        sig_df = pd.read_csv(f"{RESULT_DIR}/SIGmats/v2_draft/SBS288_mat.tsv", sep="\t")
        sig_df.set_index("Samples", inplace=True)
        signatures = sig_df.columns
        sig_df["TSMC"] = sig_df[signatures].agg("sum", axis=1)
        sig_df["tumour_sample_platekey"] = sig_df.index.map(
            lambda x: "_".join(x.split("_")[1:3])
        )
        sample_df = pd.merge(
            sample_df,
            sig_df[["TSMC", "tumour_sample_platekey"]],
            how="inner",
            on="tumour_sample_platekey",
        )
        print(f"Total somatic mutation signature: {len(sample_df)}")
        print(
            f"Total: {len(sample_df)}, \
        Unique participant: {len(np.unique(sample_df['participant_id']))}, \
        Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
        )

    # Hypermutations
    if hyper:
        sample_df["tumour_tissue"] = sample_df.tumour_group.map(
            lambda x: x.split("-")[0]
        )
        if not hyper:
            sample_df = sample_df[
                (
                    sample_df.tumour_tissue.map(
                        lambda x: x
                        not in [
                            "Uterus",
                            "CNS",
                            "Connective",
                            "Kidney",
                            "Haem",
                            "Prost",
                        ]
                    )
                )
                | (sample_df.TSMC < 14900)
            ]
            sample_df = sample_df[
                (
                    sample_df.tumour_tissue.map(
                        lambda x: x
                        not in [
                            "ColoRect",
                            "Stomach",
                            "Ovary",
                            "Breast",
                            "HeadNeck",
                            "Eso",
                        ]
                    )
                )
                | (sample_df.TSMC < 40800)
            ]
            print(f"Hypermutations: {len(sample_df)}")

    # Further sample data
    if ER:
        additional_keys = ["participant_id", "er_status"]
        additional_df = pd.read_csv(
            f"{DATA_DIR}/v14/v14_reheadered/av_tumour.tsv",
            usecols=additional_keys,
            sep="\t",
        )
        sample_df = pd.merge(
            sample_df, additional_df, how="left", on="participant_id"
        ).drop_duplicates("tumour_sample_platekey")
        print(f"Additional (er status) sample data: {len(sample_df)}")
        print(
            f"Total: {len(sample_df)}, \
        Unique participant: {len(np.unique(sample_df['participant_id']))}, \
        Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
        )

    # CRC MSI status
    if MSI:
        crc_df = pd.read_csv(
            f"{DATA_DIR}/0.sampleLists/projectTable.somaticLandscape.allUnique.2020-10-06.tsv",
            sep="\t",
            usecols=[
                "participant_id",
                "tumour_sample_platekey",
                "germline_sample_platekey",
                "subtype",
            ],
        )
        sample_df = pd.merge(
            sample_df,
            crc_df[["tumour_sample_platekey", "subtype"]],
            on="tumour_sample_platekey",
            how="left",
        )
        # All MSI status
        mmr_df = pd.read_csv(
            "/re_gecip/cancer_pan/dchubb/combined_results.txt",
            sep="\t",
            usecols=["Position", "msi_status"],
        )
        sample_df = pd.merge(
            sample_df,
            mmr_df,
            left_on="tumour_sample_platekey",
            right_on="Position",
            how="left",
        )
        # Uterus MSI status
        print(f"MSI status: {len(sample_df)}")
        print(
            f"Total: {len(sample_df)}, \
        Unique participant: {len(np.unique(sample_df['participant_id']))}, \
        Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
        )

    # Primary tumour
    if primary:
        sample_df = sample_df[
            (sample_df.tumour_group.map(lambda x: x.split("-")[0]) == "Skin")
            | (sample_df.tumour_type != "METASTASES")
        ]
        print(f"Primary tumour: {len(sample_df)}")
        print(
            f"Total: {len(sample_df)}, \
        Unique participant: {len(np.unique(sample_df['participant_id']))}, \
        Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
        )

        # Melanoma sample site
        # Melanoma sample site
        melanoma_df = pd.read_csv(
            f"{DATA_DIR}/MALIGNANT_MELANOMAv10_clinical.tsv",
            sep="\t",
            usecols=[
                "tumour_pseudo_id",
                "platekeys",
                "gel_topographies_names",
                "treatment_lowest",
            ],
        )
        melanoma_df.sort_values(["platekeys", "treatment_lowest"], inplace=True)
        melanoma_df.drop_duplicates("platekeys", keep="first", inplace=True)
        melanoma_df = melanoma_df[
            (melanoma_df.treatment_lowest < 30) & (melanoma_df.treatment_lowest > -1)
        ]
        melanoma_df["tumour_sample_platekey"] = melanoma_df["platekeys"].map(
            lambda x: "_".join(x.split("_")[1:3])
        )

        sample_df = pd.merge(
            sample_df, melanoma_df, how="left", on="tumour_sample_platekey"
        )
        sample_df["skin_site"] = (
            sample_df["gel_topographies_names"]
            .astype(str)
            .map(lambda x: "skin" in x.lower())
        )
        sample_df = sample_df[
            (sample_df.tumour_group.map(lambda x: x.split("-")[0]) != "Skin")
            | (sample_df.skin_site)
        ]
        print(f"Skin site tumour: {len(sample_df)}")
        print(
            f"Total: {len(sample_df)}, \
        Unique participant: {len(np.unique(sample_df['participant_id']))}, \
        Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
        )

        print(
            f"Total: {len(melanoma_df)}, \
        Unique platekey: {len(np.unique(melanoma_df['tumour_sample_platekey']))}"
        )

    # Variable transformations
    # Sex
    sample_df["is_female"] = np.where(sample_df["sex"] == "FEMALE", 1, -1)
    # # Tissue type - onehot encode
    sample_df["tumour_tissue"] = np.array(
        [t.split("-")[0] for t in sample_df["tumour_group"]]
    )
    sample_df["tumour_subtissue"] = np.array(
        [t.split("-")[1] for t in sample_df["tumour_group"]]
    )
    # Purity
    sample_df["logit_purity"] = np.log(
        (sample_df["purity"] + 1e-10) / (1 - sample_df["purity"] + 1e-10)
    )
    print(
        f"Nan purity {np.sum(np.isnan(sample_df.purity))} Nan logit_purity {np.sum(np.isnan(sample_df.logit_purity))}"
    )
    # Log Age
    sample_df["log_age"] = np.log(sample_df["age"])

    return sample_df


def mergeSamples(file_loc="", hyper=False):
    # Import sample data
    print("Sample data")
    sample_keys = [
        "tumour_sample_platekey",
        "germline_sample_platekey",
        "age_sampling",
        "tumour_group",
        "tumour_type",
        "sex",
        "purity",
        "participant_id",
        "pseudonymous_id",
        "day_sampling",
    ]
    sample_df = pd.read_csv(
        f"{RESULT_DIR}/sample_lists_incl_SEGs/sample_list.tsv",
        usecols=sample_keys,
        delim_whitespace=True,
    ).rename(columns={"age_sampling": "age"})
    print(f"pan-cancer sample list {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    # Aggregate GVCF sample stats
    agg_keys = ["platekey", "pc1", "pc2", "pc3", "pc4", "pc5"]
    agg_df = pd.read_csv(
        f"{DATA_DIR}/aggV2/additional_data/aggregate_gvcf_sample_stats/aggregate_gvcf_sample_stats_v10_78195.tsv",
        usecols=agg_keys,
        sep="\t",
    )
    sample_df = pd.merge(
        sample_df,
        agg_df,
        how="inner",
        left_on="germline_sample_platekey",
        right_on="platekey",
    )
    print(f"PC data {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    # Signature data
    sig_df = pd.read_csv(f"{RESULT_DIR}/SIGmats/v2_draft/SBS288_mat.tsv", sep="\t")
    sig_df.set_index("Samples", inplace=True)
    signatures = sig_df.columns
    sig_df["TSMC"] = sig_df[signatures].agg("sum", axis=1)
    sig_df["tumour_sample_platekey"] = sig_df.index.map(
        lambda x: "_".join(x.split("_")[1:3])
    )
    sample_df = pd.merge(
        sample_df,
        sig_df[["TSMC", "tumour_sample_platekey"]],
        how="inner",
        on="tumour_sample_platekey",
    )
    print(f"Total somatic mutation signature: {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    # Hypermutations
    sample_df["tumour_tissue"] = sample_df.tumour_group.map(lambda x: x.split("-")[0])
    sample_df["tumour_subtissue"] = sample_df.tumour_group.map(
        lambda x: x.split("-")[1]
    )
    gruber_groups = {
        "BileDuct": "Hepatopancreatobiliary",
        "Panc": "Hepatopancreatobiliary",
        "Liver": "Hepatopancreatobiliary",
        "Mes": "Lung",
        "Lung": "Lung",
        "Stomach": "UpperGastro",
        "Eso": "UpperGastro",
    }
    sample_df.tumour_tissue.replace(gruber_groups, inplace=True)

    # Further sample data
    additional_keys = ["participant_id", "er_status"]
    additional_df = pd.read_csv(
        f"{DATA_DIR}/v14/v14_reheadered/av_tumour.tsv",
        usecols=additional_keys,
        sep="\t",
    )
    sample_df = pd.merge(
        sample_df, additional_df, how="left", on="participant_id"
    ).drop_duplicates("tumour_sample_platekey")
    print(f"Additional (er status) sample data: {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    # All MSI status
    mmr_df = pd.read_csv(
        "/re_gecip/cancer_pan/dchubb/combined_results.txt",
        sep="\t",
        usecols=["Position", "msi_status"],
    )
    sample_df = pd.merge(
        sample_df,
        mmr_df,
        left_on="tumour_sample_platekey",
        right_on="Position",
        how="left",
    )

    # CRC MSI status
    crc_df = pd.read_csv(
        f"{DATA_DIR}/0.sampleLists/projectTable.somaticLandscape.allUnique.2020-10-06.tsv",
        sep="\t",
        usecols=[
            "participant_id",
            "tumour_sample_platekey",
            "germline_sample_platekey",
            "subtype",
        ],
    )
    sample_df = pd.merge(
        sample_df,
        crc_df[["tumour_sample_platekey", "subtype"]],
        on="tumour_sample_platekey",
        how="left",
    )
    crc_msi_status = sample_df.subtype.copy()
    crc_msi_status[
        (sample_df.subtype.astype(str) == "nan")
        & (sample_df.tumour_tissue == "ColoRect")
    ] = sample_df["msi_status"][
        (sample_df.subtype.astype(str) == "nan")
        & (sample_df.tumour_tissue == "ColoRect")
    ]
    sample_df["crc_msi_status"] = crc_msi_status
    sample_df["crc_msi_status"].replace(["NEG", "POS"], ["MSS", "MSI"], inplace=True)
    # Endometrial MSI status
    endo_df = sample_df[
        sample_df.tumour_group.map(lambda x: x.split("-")[0] == "Uterus")
    ][["tumour_sample_platekey", "msi_status"]]
    endo_msi_status = pd.read_csv(
        f"{RESULT_DIR}/Endo_MSI/combined_results.txt", sep="\t"
    )
    endo_df = pd.merge(
        endo_df,
        endo_msi_status,
        left_on="tumour_sample_platekey",
        right_on="Position",
        how="left",
        suffixes=("_pan", "_endo"),
    ).astype({"msi_status_pan": str, "msi_status_endo": str})
    endo_tcga_status = pd.read_csv(
        f"{RESULT_DIR}/Endo_MSI/subtyping_TCGA.tsv", sep="\t"
    )
    endo_tcga_status["tumour_sample_platekey"] = endo_tcga_status.sample_id.map(
        lambda x: x.split(".")[2]
    )
    endo_df = pd.merge(
        endo_df, endo_tcga_status, on="tumour_sample_platekey", how="left"
    ).astype({"TCGA_subtype": str})
    endo_msi_status = endo_df.TCGA_subtype.copy()
    endo_msi_status[endo_df.TCGA_subtype == "nan"] = endo_df["msi_status_endo"][
        endo_df.TCGA_subtype == "nan"
    ]
    endo_msi_status[
        (endo_df.TCGA_subtype == "nan") & (endo_df.msi_status_endo == "nan")
    ] = endo_df["msi_status_pan"][
        (endo_df.TCGA_subtype == "nan") & (endo_df.msi_status_endo == "nan")
    ]
    endo_df["endo_msi_status"] = endo_msi_status
    endo_df["endo_msi_status"].replace(
        ["CNA-high", "CNA-low", "NEG", "POS"],
        ["MSS", "MSS", "MSS", "MSI"],
        inplace=True,
    )
    sample_df = pd.merge(
        sample_df,
        endo_df[["tumour_sample_platekey", "endo_msi_status"]],
        how="left",
        on="tumour_sample_platekey",
    )
    sample_df.msi_status[sample_df.endo_msi_status == "MSS"] = "NEG"
    sample_df.msi_status[
        (sample_df.endo_msi_status == "MSI") | (sample_df.endo_msi_status == "POL")
    ] = "POS"
    print(f"MSI status: {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    # Hypermutation status
    hyper_status = np.zeros(len(sample_df), dtype=int)
    Endo_set = ["Uterus", "CNS", "Connective", "Kidney", "Haem", "Prost"]
    CRC_set = ["ColoRect", "Ovary", "Breast", "HeadNeck", "UpperGastro"]
    # MSI
    hyper_status[
        (sample_df.tumour_tissue.map(lambda x: x in Endo_set))
        & (
            (np.log10(sample_df.TSMC) > 4.20)
            | (sample_df.msi_status == "POS")
            | (sample_df.endo_msi_status == "MSI")
        )
    ] = 1
    hyper_status[
        (sample_df.tumour_tissue.map(lambda x: x in CRC_set))
        & (
            (np.log10(sample_df.TSMC) > 4.61)
            | (sample_df.msi_status == "POS")
            | (sample_df.crc_msi_status == "MSI")
        )
    ] = 1
    hyper_status[
        (sample_df.tumour_tissue.map(lambda x: x not in Endo_set + CRC_set))
        & (sample_df.msi_status == "POS")
    ] = 1
    # POLE
    hyper_status[
        (sample_df.tumour_tissue.map(lambda x: x in Endo_set))
        & ((np.log10(sample_df.TSMC) > 5.53) | (sample_df.endo_msi_status == "POL"))
    ] = 2
    hyper_status[
        (sample_df.tumour_tissue.map(lambda x: x in CRC_set))
        & (
            (np.log10(sample_df.TSMC) > 5.59)
            | (sample_df.crc_msi_status.replace("MSI&POL", "POL") == "POL")
        )
    ] = 2
    sample_df["hyper_status"] = hyper_status
    if not hyper:
        sample_df = sample_df[sample_df.hyper_status == 0]

    # Variable transformations
    # Sex
    sample_df["is_female"] = np.where(sample_df["sex"] == "FEMALE", 1, 0)
    # Purity
    sample_df["logit_purity"] = np.log(
        (sample_df["purity"] + 1e-10) / (1 - sample_df["purity"] + 1e-10)
    )
    # Log Age
    sample_df["log_age"] = np.log(sample_df["age"])

    # Primary tumoue
    sample_df = sample_df[
        (sample_df.tumour_tissue == "Skin") | (sample_df.tumour_type != "METASTASES")
    ]
    print(f"Primary tumour: {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    # Melanoma sample site
    melanoma_df = pd.read_csv(
        f"{DATA_DIR}/MALIGNANT_MELANOMAv10_clinical.tsv",
        sep="\t",
        usecols=[
            "tumour_pseudo_id",
            "platekeys",
            "gel_topographies_names",
            "treatment_lowest",
        ],
    )
    melanoma_df.sort_values(["platekeys", "treatment_lowest"], inplace=True)
    melanoma_df.drop_duplicates("platekeys", keep="first", inplace=True)
    melanoma_df = melanoma_df[
        (melanoma_df.treatment_lowest < 30) & (melanoma_df.treatment_lowest > -1)
    ]
    melanoma_df["tumour_sample_platekey"] = melanoma_df["platekeys"].map(
        lambda x: "_".join(x.split("_")[1:3])
    )

    sample_df = pd.merge(
        sample_df, melanoma_df, how="left", on="tumour_sample_platekey"
    )
    sample_df["skin_site"] = (
        sample_df["gel_topographies_names"]
        .astype(str)
        .map(lambda x: "skin" in x.lower())
    )
    sample_df = sample_df[(sample_df.tumour_tissue != "Skin") | (sample_df.skin_site)]
    print(f"Skin site tumour: {len(sample_df)}")
    print(
        f"Total: {len(sample_df)}, \
    Unique participant: {len(np.unique(sample_df['participant_id']))}, \
    Unique platekey: {len(np.unique(sample_df['tumour_sample_platekey']))}"
    )

    print(
        f"Total: {len(melanoma_df)}, \
    Unique platekey: {len(np.unique(melanoma_df['tumour_sample_platekey']))}"
    )

    return sample_df


if __name__ == "__main__":
    # filename to save covariates into
    filename_output = sys.argv[1]

    # GET DATA - include hypermutated samples
    print("Load samples")
    sample_df = getSamples(hyper=True)

    # Remove samples where sex doesn't match organ
    sample_df = sample_df[
        (
            ~(
                sample_df.tumour_tissue.map(
                    lambda x: x in ["Breast", "Ovary", "Uterus"]
                )
                & (sample_df.is_female == 0)
            )
        )
        & (
            ~(
                sample_df.tumour_tissue.map(lambda x: x in ["Testis", "Prostate"])
                & (sample_df.is_female == 1)
            )
        )
    ]
    print("Removed male samples for Breast,Ovary,Uterus and female for Testis,Prostate")
    print(f"Number of samples: {len(sample_df)}")

    # Get mutation rates
    # Germline
    germline_df = (
        pd.read_csv(
            f"{DATA_DIR}/cancGeneHits/germline/output/aggv2-germline_cancer-gene_hits.tsv",
            sep="\t",
        )
        .set_index("gene_id")
        .T
    )
    # OncoKB
    somatic_df = pd.read_csv(
        f"{DATA_DIR}/cancGeneHits/somatic/output/OncoKB_somatic_cancer-gene_hits.tsv",
        sep="\t",
        index_col=0,
    )
    # LOH
    loh_df = pd.read_csv(
        f"{DATA_DIR}/cancGeneHits/somatic/output/lohfrac_gene_sample_matrix.tsv",
        sep="\t",
    ).set_index("tumour_sample_platekey")
    # Add in platekeys which failed battenberg and assume no LoH
    failed_battenberg_platekeys = np.setxor1d(
        np.intersect1d(loh_df.index, somatic_df.index), somatic_df.index
    )
    loh_df = pd.concat(
        (
            loh_df,
            pd.DataFrame(
                np.zeros((len(failed_battenberg_platekeys), len(loh_df.keys()))),
                index=failed_battenberg_platekeys,
                columns=loh_df.keys(),
            ),
        )
    )
    # Only tumour suppressor genes
    cancer_genes = pd.read_csv(f"{DATA_DIR}/cancer_gene_census.csv")
    genes = cancer_genes["Gene Symbol"][
        cancer_genes["Role in Cancer"].map(lambda x: "TSG" in str(x))
    ]
    # Get binary two hit model data
    tumour_sample_platekeys = list(sample_df.tumour_sample_platekey)
    germline_sample_platekeys = list(sample_df.germline_sample_platekey)
    germline_df = germline_df.loc[germline_sample_platekeys][genes].set_index(
        sample_df.tumour_sample_platekey
    )
    loh_df = loh_df.loc[tumour_sample_platekeys][genes].set_index(
        sample_df.tumour_sample_platekey
    )
    somatic_df = somatic_df.loc[tumour_sample_platekeys][genes].set_index(
        sample_df.tumour_sample_platekey
    )
    twohit_df = pd.DataFrame(
        np.array(germline_df).astype(int)
        + np.array(loh_df > 0.5).astype(int)
        + np.array(somatic_df),
        columns=genes,
    )
    sample_df["hit_rate"] = np.sum(np.array(twohit_df), axis=1)

    # Confounding variables
    variables = (
        ["log_age", "is_female"] + [f"pc{i}" for i in range(1, 4)]
    )  # + ['hit_rate'] # 'log_age', ['age','log_age','logit_purity','is_female'] + [f"pc{i}" for i in range(1,6)]

    # Construct X matrix
    print("Construct X")

    # Get survival times and whether patient is still alive
    sample_df["last_followup_time"] = (
        pd.to_datetime(sample_df["day_last_followup"])
        - pd.to_datetime(sample_df.day_sampling)
    ).dt.days
    sample_df["death_time"] = (
        pd.to_datetime(sample_df["day_death"]) - pd.to_datetime(sample_df.day_sampling)
    ).dt.days
    sample_df["survival_time"] = np.where(
        ~np.isnan(sample_df.death_time),
        sample_df.death_time,
        sample_df.last_followup_time,
    )
    sample_df["vital_status"] = np.where(~np.isnan(sample_df.death_time), 1, 0)
    dates_of_interest = ["survival_time", "vital_status"]

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
            + dates_of_interest
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
    X[["sample_id", "group"] + dates_of_interest + variables].to_csv(
        filename_output, index=False, sep="\t"
    )
