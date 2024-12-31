import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dotenv import load_dotenv

from signatures.plotting.combinedSignatures import loadSignatures, signatureRenamer

load_dotenv()
DATA_DIR = os.getenv("DATA_DIR")
RESULT_DIR = os.getenv("RESULT_DIR")

mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["font.family"] = "STIXGeneral"
plt.rc("axes", labelsize=16)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)


if __name__ == "__main__":
    # Get renaming dictionary
    combined_acts, combined_sigs, relabel_map = loadSignatures()
    rename_dict = signatureRenamer(relabel_map)

    ### Load sample list
    samples_file = f"{DATA_DIR}/sample_lists_incl_SEGs/sample_list.tsv"
    sample_df = pd.read_csv(
        samples_file,
        usecols=[
            "participant_id",
            "tumour_sample_platekey",
            "germline_sample_platekey",
            "tumour_group",
            "signature_extraction_group",
        ],
        sep="\t",
    )
    sample_df["tumour_tissue"] = sample_df.tumour_group.map(lambda x: x.split("-")[0])
    sample_df["sample_id"] = (
        sample_df.participant_id.astype(str)
        + "_"
        + sample_df.tumour_sample_platekey
        + "_"
        + sample_df.germline_sample_platekey
    )

    # Activities
    combined_act_df = pd.DataFrame(index=list(sample_df.sample_id))
    for sig_type in ["SBS288", "DBS78", "ID83", "CNV48", "SV32"]:
        combined_act_df = pd.merge(
            combined_act_df,
            combined_acts[sig_type],
            left_index=True,
            right_index=True,
            how="left",
        ).fillna(0)

    inactivation_genes = {
        "HRD": ["BRCA1", "BRCA2", "PALB2", "BRIP1", "RAD51B"],
        "APOBEC": ["APOBEC3B"],
        "MMR": ["MSH6", "MSH2", "MLH1", "PMS2"],
    }
    repair_signatures = {
        "HRD": ["SBS3", "ID6", "CN17"],
        "APOBEC": ["SBS2", "SBS13"],
        "MMR": ["SBS15", "SBS26", "SBS44"],
    }

    # Germline
    germline_df = (
        pd.read_csv(
            f"{RESULT_DIR}/cancGeneHits/germline_30/output/aggv2-germline_cancer-gene_hits.tsv",
            sep="\t",
        )
        .set_index("gene_id")
        .T
    )
    # OncoKB
    somatic_df = pd.read_csv(
        f"{RESULT_DIR}/cancGeneHitssomatic/output/OncoKB_somatic_cancer-gene_hits.tsv",
        sep="\t",
        index_col=0,
    )
    # LOH
    loh_df = pd.read_csv(
        f"{RESULT_DIR}/cancGeneHits/somatic/output/lohfrac_gene_sample_matrix.tsv",
        sep="\t",
    ).set_index("tumour_sample_platekey")
    # Match up the indices
    loh_df = (
        pd.merge(
            sample_df[["tumour_sample_platekey", "sample_id"]],
            loh_df,
            how="left",
            left_on="tumour_sample_platekey",
            right_index=True,
        )
        .set_index("sample_id")
        .drop("tumour_sample_platekey", axis=1)
        .fillna(0)
    )
    somatic_df = (
        pd.merge(
            sample_df[["tumour_sample_platekey", "sample_id"]],
            somatic_df,
            how="left",
            left_on="tumour_sample_platekey",
            right_index=True,
        )
        .set_index("sample_id")
        .drop("tumour_sample_platekey", axis=1)
        .fillna(0)
    )
    germline_df = (
        pd.merge(
            sample_df[
                ["germline_sample_platekey", "tumour_sample_platekey", "sample_id"]
            ],
            germline_df,
            how="left",
            left_on="germline_sample_platekey",
            right_index=True,
        )
        .set_index("sample_id")
        .drop(["tumour_sample_platekey", "germline_sample_platekey"], axis=1)
        .fillna(0)
    )
    index = np.intersect1d(loh_df.index, somatic_df.index)
    index = np.intersect1d(index, np.array(germline_df.index))

    # List of all genes to consider
    gene_list = [
        gene for gene_list in inactivation_genes.values() for gene in gene_list
    ]

    two_hit = (
        (loh_df.loc[index][gene_list] > 0).astype(int)
        + (germline_df.loc[index][gene_list]).astype(int)
        + somatic_df.loc[index][gene_list].astype(int)
    )

    two_hit = pd.merge(
        two_hit, combined_act_df, left_index=True, right_index=True, how="inner"
    )

    for mut_type in ["HRD", "MMR", "APOBEC"]:
        print(f"Fraction of samples with {mut_type} gene inactivation")
        two_hit[mut_type + "_gene"] = (two_hit[inactivation_genes[mut_type]] > 1).sum(
            axis=1
        ) > 0
        two_hit[mut_type + "_sig"] = (two_hit[repair_signatures[mut_type]] > 0).sum(
            axis=1
        ) > (1 if mut_type == "HRD" else 0)
        two_hit[mut_type + "_both"] = (
            two_hit[mut_type + "_gene"] & two_hit[mut_type + "_sig"]
        )

        inactivated = (
            pd.merge(
                two_hit[[mut_type + "_gene", mut_type + "_sig", mut_type + "_both"]],
                pd.concat(
                    (
                        sample_df[
                            ["sample_id", "signature_extraction_group"]
                        ].set_index("sample_id"),
                        pd.DataFrame(
                            {"signature_extraction_group": "All"},
                            index=sample_df.sample_id,
                        ),
                    )
                ),
                left_index=True,
                right_index=True,
            )
            .groupby("signature_extraction_group")
            .sum()
        )
        inactivated["total"] = (
            pd.concat(
                (
                    sample_df,
                    pd.DataFrame(
                        {"signature_extraction_group": "All"}, index=sample_df.sample_id
                    ),
                )
            )
            .groupby("signature_extraction_group")
            .size()
        )
        inactivated = pd.merge(
            inactivated,
            inactivated / np.array(inactivated["total"])[:, None],
            left_index=True,
            right_index=True,
            suffixes=("", "_f"),
        )
        inactivated["frac_sig_gene"] = (
            inactivated[mut_type + "_both"] / inactivated[mut_type + "_sig"]
        )
        inactivated["frac_gene_sig"] = (
            inactivated[mut_type + "_both"] / inactivated[mut_type + "_gene"]
        )
        inactivated["total"] = sample_df.groupby("signature_extraction_group").size()

        inactivated.to_csv(
            f"{RESULT_DIR}/repairPathways/{mut_type}_counts.tsv", sep="\t"
        )
