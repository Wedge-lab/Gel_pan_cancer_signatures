import os
import sys

import numpy as np
import pandas as pd
import scipy.stats
from dotenv import load_dotenv

load_dotenv()
DATA_DIR = os.getenv("DATA_DIR")

if __name__ == "__main__":
    samples_file, targets_file = sys.argv[1:]

    # Sample list
    sample_df = pd.read_csv(samples_file, sep="\t")
    ids = sample_df.sample_id.map(
        lambda x: [
            x.split("_")[0],
            "_".join(x.split("_")[1:3]),
            "_".join(x.split("_")[3:]),
        ]
    )
    sample_df = pd.concat(
        (
            sample_df,
            pd.DataFrame(
                ids.tolist(),
                index=ids.index,
                columns=[
                    "participant_id",
                    "tumour_sample_platekey",
                    "germline_sample_platekey",
                ],
            ),
        ),
        axis=1,
    )

    # OncoKB
    somatic_df = pd.read_csv(
        f"{DATA_DIR}/OncoKB_somatic_cancer-gene_hits.tsv", sep="\t", index_col=0
    )

    # Only tumour suppressor genes
    cancer_genes = pd.read_csv("{DATA_DIR}/cancer_gene_census.csv")
    genes = cancer_genes["Gene Symbol"][
        cancer_genes["Role in Cancer"].map(lambda x: "TSG" in str(x))
    ]

    # Get binary germline hit model data
    tumour_sample_platekeys = list(sample_df.tumour_sample_platekey)
    somatic_df = somatic_df.loc[tumour_sample_platekeys][genes].set_index(
        sample_df.sample_id
    )
    somatic_df = pd.DataFrame(
        np.array(somatic_df > 0).astype(int), index=sample_df.sample_id, columns=genes
    )
    print(np.unique(somatic_df, return_counts=True))
    somatic_df.index.name = "sample_id"

    # DNA types from database
    DNA_repair = pd.merge(
        pd.read_csv(f"{DATA_DIR}/human-dna-repair-genes.tsv", sep="\t"),
        pd.DataFrame(index=np.array(genes)),
        left_on="Gene",
        right_index=True,
        how="inner",
    )
    DNA_repair_types = {
        "Chromatin Structure and Modification": "CSM",
        "DNA polymerases ": "POL",
        "Direct reversal of damage": "DRD",
        "Editing and processing nucleases": "EPN",
        "Fanconi anemia": "FA",
        "Genes defective in diseases associated with sensitivity to DNA damaging agents": "DDA",
        "Homologous recombination": "HR",
        "Other conserved DNA damage response genes": "OTHER",
    }
    DNA_repair.replace({"Type": DNA_repair_types}, inplace=True)
    for type in np.unique(DNA_repair["Type"]):
        type_hit = np.zeros(len(somatic_df), dtype=int)
        type_hit[
            np.sum(
                somatic_df[[gene for gene in DNA_repair.Gene[DNA_repair.Type == type]]]
                > 0,
                axis=1,
            )
            > 0
        ] = 1
        somatic_df[f"{type}"] = type_hit

    # Generate mock targets - for 1%,5%,20%,50% hit rate
    for hit_rate in [5, 15, 50]:
        somatic_df[f"MOCK{hit_rate:02d}"] = scipy.stats.binom.rvs(
            1, hit_rate / 100, size=len(somatic_df)
        )

    # Generate mock targets - perturbation of genes within cohorts
    for gene in ["NTHL1", "BRCA2", "MGMT", "POL", "FA", "HR", "MMR"]:
        indices = np.arange(len(somatic_df))
        MOCKgene = np.zeros(len(somatic_df), dtype=int)
        for group in np.unique(sample_df.group):
            sample = np.random.choice(
                indices[sample_df.group == group],
                size=np.sum(sample_df.group == group),
                replace=False,
            )
            MOCKgene[sample_df.group == group] = somatic_df[gene][sample]
        somatic_df[f"MOCK{gene}"] = MOCKgene

    # Save targets file
    somatic_df.to_csv(targets_file, sep="\t")
