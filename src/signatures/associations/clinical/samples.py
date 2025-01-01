# Import standard modules
import os
import sys

import numpy as np
import pandas as pd

from signatures.config import load_environment
from signatures.sampleCuration.pan_cancer_sample import getSamples

load_environment()
DATA_DIR = os.getenv("DATA_DIR")

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

    # Get mutation rates
    # Germline
    germline_df = (
        pd.read_csv(f"{DATA_DIR}/aggv2-germline_cancer-gene_hits.tsv", sep="\t")
        .set_index("gene_id")
        .T
    )
    # OncoKB
    somatic_df = pd.read_csv(
        f"{DATA_DIR}/OncoKB_somatic_cancer-gene_hits.tsv", sep="\t", index_col=0
    )
    # LOH
    loh_df = pd.read_csv(
        f"{DATA_DIR}/lohfrac_gene_sample_matrix.tsv", sep="\t"
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
    X[["sample_id", "group"] + variables].to_csv(filename_output, index=False, sep="\t")
