"combine signatures between cohorts - iteratively add cohort extracted signatures to pan-cancer COSMIC list"

import inspect
import os
import string
import sys

import numpy as np
import pandas as pd
import SigProfilerAssignment
import SigProfilerExtractor
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerAssignment import decompose_subroutines

from signatures.config import load_environment
from signatures.utils import regexSearch

print(f"SigProfilerAssignment file location: {inspect.getfile(Analyze.decompose_fit)}")
print(f"Extractor version {SigProfilerExtractor.__version__}")
print(f"Assignment version {SigProfilerAssignment.__version__}")


load_environment()
DATA_DIR = os.getenv("DATA_DIR")


def posteriorSignatures(sig_values, sig_se):
    # Get posterior covariance
    posterior_cov = 1 / np.sum(1 / sig_se**2, axis=1)

    # Inverse weighted mean
    posterior_mean = posterior_cov * np.sum(sig_values / sig_se**2, axis=1)
    posterior_mean /= np.sum(posterior_mean)

    return posterior_mean, np.sqrt(posterior_cov)


def cosineSimilarity(sigsX, sigsY):
    X = np.array(sigsX)
    Y = np.array(sigsY)
    # Reshape signatures if only one provided
    if len(X.shape) == 1:
        X = X[:, None]
    if len(Y.shape) == 1:
        Y = Y[:, None]

    # Generate dataframe of cosine similarities
    X = X / np.sqrt(np.sum(X**2, axis=0))[None, :]
    Y = Y / np.sqrt(np.sum(Y**2, axis=0))[None, :]

    return X.T @ Y


def cosineSimilarityReduce(signatures, threshold=0.8):
    """cosineSimilarityReduce
    Map signatures together where cosine similarity>0.8
    """

    # Generate dataframe of cosine similarities
    sigmat = (
        np.array(signatures)
        / np.sqrt(np.sum(np.array(signatures) ** 2, axis=0))[None, :]
    )
    cosine_similarity_mtx = sigmat.T @ sigmat

    # Dataframe of cosine similarities for signature pairs
    XX, YY = np.meshgrid(signatures.keys(), signatures.keys())
    cossim_df = pd.DataFrame(
        {
            "cos_sim": np.triu(cosine_similarity_mtx).flatten(),
            "sig_x": XX.flatten(),
            "sig_y": YY.flatten(),
        }
    )
    # Drop same signatures and lower triangle
    cossim_df = cossim_df[
        (cossim_df["sig_x"] != cossim_df["sig_y"]) & (cossim_df["cos_sim"] > threshold)
    ]
    # Drop signatures in same cohort
    cossim_df = cossim_df[
        cossim_df.sig_x.map(lambda x: x.split("-")[1])
        != cossim_df.sig_y.map(lambda x: x.split("-")[1])
    ]
    # Sort by descending cosine similarity
    cossim_df = cossim_df.sort_values("cos_sim", ascending=False)

    print(cossim_df)

    # Create a mapping of new signatures across tissue types
    mapping = {}
    i = 0
    while len(cossim_df) > 0:
        mapping[i] = [cossim_df.sig_x.iloc[0], cossim_df.sig_y.iloc[0]]
        cossim_df = cossim_df[1:]

        cohort_set = [x.split("-")[1] for x in mapping[i]]
        drop_rows = []

        for row in cossim_df.iterrows():
            if (row[1].sig_x in mapping[i]) & (
                row[1].sig_y.split("-")[1] not in cohort_set
            ):
                if (
                    len(
                        set(mapping[i])
                        - set(cossim_df.sig_x[cossim_df.sig_y == row[1].sig_y])
                    )
                    == 0
                ):
                    mapping[i].append(row[1].sig_y)
                    cohort_set.append(row[1].sig_y.split("-")[1])
                drop_rows.append(row[0])
            elif (row[1].sig_y in mapping[i]) & (
                row[1].sig_x.split("-")[1] not in cohort_set
            ):
                if (
                    len(
                        set(mapping[i])
                        - set(cossim_df.sig_y[cossim_df.sig_x == row[1].sig_x])
                    )
                    == 0
                ):
                    mapping[i].append(row[1].sig_x)
                    cohort_set.append(row[1].sig_x.split("-")[1])
                drop_rows.append(row[0])

        for row in cossim_df.iterrows():
            if (row[1].sig_x in mapping[i]) & (row[1].sig_y in mapping[i]):
                drop_rows.append(row[0])

        cossim_df = cossim_df.drop(np.unique(drop_rows), axis=0)

        i += 1

    mapping_sigs = [x for key in mapping for x in mapping[key]]
    for sig in signatures.columns:
        if not sig in mapping_sigs:
            mapping[i] = [sig]
            i += 1

    return mapping


def orderSignatures(signatures):
    """orderSignatures - sort signatures into a sensible ordering
    1) mutation type - SBS, DBS, ID
    2) signature number - 1,2,3,...
    3) novel signature cohort - Bladder, Breast...
    4) signature letter - 7a,7b,7c
    """

    mut_type = np.array(
        [regexSearch("T*([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig, 1) for sig in signatures]
    )  # type: ignore
    mut_index = np.array(
        [regexSearch("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig, 2) for sig in signatures]
    ).astype(int)  # type: ignore
    mut_exten = np.array(
        [regexSearch("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig, 3) for sig in signatures]
    )  # type: ignore
    mut_tissue = np.array(
        [
            regexSearch("([A-Za-z]+)([0-9]+)([a-zA-Z]*)-([a-zA-Z_]+)", sig, 4)
            if "-" in sig
            else "Z"
            for sig in signatures
        ]
    )

    order_exten = np.zeros(len(mut_exten))
    order_exten[np.argsort(mut_exten)] = np.arange(len(mut_exten))

    order_tissue = np.zeros(len(mut_tissue))
    order_tissue = np.arange(len(np.unique(mut_tissue)))[
        np.unique(mut_tissue, return_inverse=True)[1]
    ]

    # print(mut_type)

    order = np.argsort(
        pd.Series(mut_type)
        .replace(
            ["SBS", "DBS", "ID", "CN", "CNV", "RefSigR", "RS", "SV"],
            [1, 2, 3, 4, 5, 6, 7, 8],
        )
        .astype(int)
        + mut_index / 1e3
        + order_tissue / 1e5
        + order_exten / 1e7
        + np.array([sig[0] == "T" for sig in signatures]).astype(int) * 100
    )
    return order


def deconvolve_sigs(
    signatures,
    output,
    signature_database,
    nnls_add_penalty=0.05,
    nnls_remove_penalty=0.01,
    initial_remove_penalty=0.05,
    genome_build="GRCh37",
    cosmic_version=3.2,
    make_plots=True,
    collapse_to_SBS96=True,
    connected_sigs=True,
    verbose=False,
    new_signature_thresh_hold=0.8,
):
    layer_directory2 = output + "/Decompose_Solution"
    if not os.path.exists(layer_directory2):
        os.makedirs(layer_directory2)

    try:
        if not os.path.exists(output):
            os.makedirs(output)
    except:
        print("The {} folder could not be created".format("output"))

        #################
        # Decomposition
        #################
    try:
        processAvg = pd.read_csv(signatures, sep="\t", index_col=0)
    except:
        try:
            processAvg = signatures
        except:
            sys.exit(
                "Error in formatting of input signatures, Pass a text file of signatures in the format of denovo signatures"
            )

    listOfSignatures = list(processAvg.columns)
    index = processAvg.index

    mutation_type = str(processAvg.shape[0])
    if mutation_type == "78":
        mutation_context = "DBS78"
    elif mutation_type == "83":
        mutation_context = "ID83"
    elif mutation_type == "48":
        mutation_context = "CNV48"
    else:
        mutation_context = "SBS" + mutation_type

    # creating list of mutational type to sync with the vcf type input
    signature_names = decompose_subroutines.make_letter_ids(
        idlenth=processAvg.shape[1], mtype=mutation_context
    )

    originalProcessAvg = processAvg.copy()
    if processAvg.shape[0] == 1536 and collapse_to_SBS96 == True:
        # collapse the 1596 context into 96 only for the deocmposition
        processAvg = pd.DataFrame(processAvg, index=index)
        processAvg = processAvg.groupby(processAvg.index.str[1:8]).sum()
    if processAvg.shape[0] == 288 and collapse_to_SBS96 == True:
        # collapse the 288 context into 96 only for the deocmposition
        processAvg = pd.DataFrame(processAvg, index=index)
        processAvg = processAvg.groupby(processAvg.index.str[2:9]).sum()
    processAvg = np.array(processAvg)

    if len(pd.read_csv(signature_database, nrows=1, index_col=0, sep="\t").keys()) == 0:
        print("No signatures in reference")
        # Start from scratch
        final_signatures = {}
        final_signatures["globalsigids"] = []
        final_signatures["newsigids"] = signature_names  # decomposed_solution.columns
        final_signatures["dictionary"] = {
            sigid: [sigid] for sigid in final_signatures["newsigids"]
        }
        final_signatures["activity_percentages"] = [
            [100] for sigid in final_signatures["newsigids"]
        ]
        final_signatures["globalsigs"] = processAvg[[]]
        final_signatures["newsigs"] = np.array(processAvg)
        return final_signatures

    final_signatures = decompose_subroutines.signature_decomposition(
        processAvg,
        mutation_type,
        layer_directory2,
        genome_build=genome_build,
        signature_database=signature_database,
        mutation_context=mutation_context,
        add_penalty=0.05,
        connected_sigs=connected_sigs,
        remove_penalty=0.01,
        make_decomposition_plots=make_plots,
        originalProcessAvg=originalProcessAvg,
        new_signature_thresh_hold=new_signature_thresh_hold,
        sig_exclusion_list=[],
    )

    return final_signatures


if __name__ == "__main__":
    best_solutions_file, cosmic_file, output_dir, sig_type, sample_id_type = sys.argv[
        1:
    ]

    # Load file of best solutions
    best_solution = pd.read_csv(best_solutions_file, sep="\t")

    # Load COSMIC reference signatures
    cosmic = pd.read_csv(cosmic_file, sep="\t", index_col=0)

    # Initialise dataframes
    decomposed_sigs = pd.DataFrame()
    decomposed_sigs_se = pd.DataFrame()
    mutations = pd.DataFrame()
    cohort_signatures = {}

    # Get best solutions signatures and activities files
    for index, row in best_solution.iterrows():
        # Get mutations
        cohort_mutations = pd.read_csv(row.samples_file, sep="\t", index_col=0)

        # Get signatures
        cohort_sigs = pd.read_csv(row.signatures_file, sep="\t", index_col=0)
        # cohort_sigs.rename(columns={key:f"{key}-{row.cohort}" for key in cohort_sigs.columns if key[:len(sig_type)]==sig_type}, inplace=True)
        cohort_sigs.rename(
            columns={key: f"{key}-{row.cohort}" for key in cohort_sigs.columns},
            inplace=True,
        )

        # Get signature uncertainties
        cohort_se = pd.read_csv(row.signatures_SE_file, sep="\t", index_col=0)
        # cohort_se.rename(columns={key:f"{key}-{row.cohort}" for key in cohort_se.columns if key[:len(sig_type)]==sig_type}, inplace=True)
        cohort_se.rename(
            columns={key: f"{key}-{row.cohort}" for key in cohort_se.columns},
            inplace=True,
        )

        # Add cohort mutations and signatures to allSignatures table
        columns = [
            col for col in cohort_sigs.columns if not col in decomposed_sigs.columns
        ]
        decomposed_sigs = pd.merge(
            decomposed_sigs,
            cohort_sigs[columns],
            how="right",
            left_index=True,
            right_index=True,
        )
        decomposed_sigs_se = pd.merge(
            decomposed_sigs_se,
            cohort_se,
            how="right",
            left_index=True,
            right_index=True,
        )
        mutations = pd.merge(
            mutations, cohort_mutations, how="right", left_index=True, right_index=True
        )

        # Signatures extracted de-novo in the cohort
        cohort_signatures[row.cohort] = list(cohort_sigs.columns)

    # Save de novo signatures
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    mutations.to_csv(f"{output_dir}/Sample_Mutations.tsv", sep="\t")
    decomposed_sigs.to_csv(f"{output_dir}/Decomposed_Solution_Signatures.tsv", sep="\t")
    decomposed_sigs_se.to_csv(
        f"{output_dir}/Decomposed_Solution_Signature_se.tsv", sep="\t"
    )

    # Standard error of combined signatures
    combined_signature_se = pd.DataFrame(index=decomposed_sigs_se.index)
    grouped_signatures_cosmic = {}

    iteration = 0
    while iteration < 50:
        output_iter = f"{output_dir}/iter_{iteration}/"
        if not os.path.exists(output_iter):
            os.makedirs(output_iter)

        cosmic.to_csv(f"{output_iter}/cosmic_extended.tsv", sep="\t")

        # Run decomposition if reference includes signatures
        print(f"Iter {iteration}. Decomposing signatures...")
        decomposition = deconvolve_sigs(
            f"{output_dir}/Decomposed_Solution_Signatures.tsv",
            output_iter,
            f"{output_iter}/cosmic_extended.tsv",
            genome_build="GRCh38",
            verbose=False,
            new_signature_thresh_hold=0.8,
            make_plots=False,
            collapse_to_SBS96=sig_type == "SBS288",
        )

        # Relabel signatures with original labels
        relabel = dict(zip(decomposition["dictionary"].keys(), decomposed_sigs.columns))
        novel_signatures = [relabel[sig] for sig in decomposition["newsigids"]]
        print("Novel: ", novel_signatures)

        # Continue until no more novel signatures
        if len(novel_signatures) == 0:
            break

        # Group novel signatures by cosine similarity
        grouped_signatures = cosineSimilarityReduce(decomposed_sigs[novel_signatures])
        for i in range(len(grouped_signatures)):
            if i < 26:
                grouped_signatures[f"{string.ascii_uppercase[i]}"] = (
                    grouped_signatures.pop(i)
                )
            else:
                grouped_signatures[
                    f"{string.ascii_uppercase[i//26+1]}{string.ascii_uppercase[i%26]}"
                ] = grouped_signatures.pop(i)
        print(grouped_signatures)

        # Get mean of signature combinations and add to reference table
        grouped_signature_values = pd.DataFrame(index=cosmic.index)
        grouped_signature_se = pd.DataFrame(index=decomposed_sigs_se.index)
        novel_candidates = pd.DataFrame()
        for new_sig, old_sigs in grouped_signatures.items():
            if len(old_sigs) == 1:
                if sig_type == "SBS288":
                    grouped_signature_values[new_sig] = (
                        decomposed_sigs[old_sigs[0]]
                        .groupby(decomposed_sigs.index.str[2:9])
                        .sum()[grouped_signature_values.index]
                    )
                else:
                    grouped_signature_values[new_sig] = decomposed_sigs[old_sigs[0]][
                        grouped_signature_values.index
                    ]
                grouped_signature_se[new_sig] = decomposed_sigs_se[old_sigs[0]]
            else:
                signature, signature_se = posteriorSignatures(
                    decomposed_sigs[old_sigs], decomposed_sigs_se[old_sigs]
                )
                if sig_type == "SBS288":
                    grouped_signature_values[new_sig] = signature.groupby(
                        decomposed_sigs.index.str[2:9]
                    ).sum()[grouped_signature_values.index]
                else:
                    grouped_signature_values[new_sig] = signature[
                        grouped_signature_values.index
                    ]
                grouped_signature_se[new_sig] = signature_se

            if cosmic.shape[1] > 0:
                max_cossim = np.max(
                    cosineSimilarity(grouped_signature_values[new_sig], cosmic)
                )
            else:
                max_cossim = 0
            novel_candidates = pd.concat(
                [
                    novel_candidates,
                    pd.DataFrame(
                        [[new_sig, old_sigs, len(old_sigs), max_cossim]],
                        columns=[
                            "new_sig",
                            "old_sigs",
                            "cohort_count",
                            "max_cosmic_similarity",
                        ],
                    ),
                ]
            )
        # Order groups by max cosine similarity with a reference signature -> take signature with lowest similarity
        novel_sig = (
            novel_candidates.sort_values("max_cosmic_similarity").iloc[0].new_sig
        )
        print(novel_sig)
        cosmic[f"{sig_type}{string.ascii_uppercase[iteration]}-Pan"] = (
            grouped_signature_values[novel_sig]
        )
        combined_signature_se[f"{sig_type}{string.ascii_uppercase[iteration]}-Pan"] = (
            grouped_signature_se[novel_sig]
        )
        grouped_signatures_cosmic[
            f"{sig_type}{string.ascii_uppercase[iteration]}-Pan"
        ] = grouped_signatures[novel_sig]

        iteration += 1

    # Save combined signature uncertainties
    combined_signature_se.to_csv(
        f"{output_dir}/Combined_Solution_Signatures_SE.tsv", sep="\t"
    )
    # Save extended cosmic file
    cosmic.to_csv(f"{output_dir}/Combined_Solution_Signatures.tsv", sep="\t")

    # Initialise dataframes
    combined_acts = pd.DataFrame()

    for index, row in best_solution.iterrows():
        print(f"Decomposing {row.cohort} to new reference...")

        # Run decomposition
        output = f"{output_dir}/{row.cohort}/"
        if not os.path.exists(output):
            os.makedirs(output)

        # Run decompose fit
        Analyze.decompose_fit(
            row.samples_file,
            output,
            signatures=row.signatures_file,
            signature_database=f"{output_dir}/Combined_Solution_Signatures.tsv",
            genome_build="GRCh38",
            verbose=False,
            new_signature_thresh_hold=0.8,
            make_plots=False,
            collapse_to_SBS96=sig_type == "SBS288",
        )

        cohort_acts = pd.read_csv(
            f"{output}/Decompose_Solution/Activities/Decompose_Solution_Activities.txt",
            sep="\t",
        ).set_index("Samples")  # .rename(columns=relabel)
        combined_acts = pd.concat((combined_acts, cohort_acts)).fillna(0)

        print(f"...done")

    # Change sample IDs to match other tumour types
    if sample_id_type == "tumour_sample_platekey":
        cancer_analysis = cancer_analysis = pd.read_csv(
            f"{DATA_DIR}/V11/V11_reheadered/cancer_analysis.tsv",
            sep="\t",
            usecols=[
                "participant_id",
                "tumour_sample_platekey",
                "germline_sample_platekey",
            ],
        )
        cancer_analysis["sample_id"] = (
            cancer_analysis["participant_id"].astype(str)
            + "_"
            + cancer_analysis["tumour_sample_platekey"]
            + "_"
            + cancer_analysis["germline_sample_platekey"]
        )
        combined_acts = (
            pd.merge(
                combined_acts.reset_index(),
                cancer_analysis[["sample_id", "tumour_sample_platekey"]],
                left_on="Samples",
                right_on="tumour_sample_platekey",
                how="left",
            )
            .drop(["Samples", "tumour_sample_platekey"], axis=1)
            .rename({"sample_id": "Samples"}, axis=1)
            .set_index("Samples")
        )
    elif sample_id_type == "pid_germline_tumour":
        combined_acts.set_index(
            combined_acts.index.map(
                lambda x: x.split(".")[0]
                + "_"
                + "_".join(x.split(".tumo")[1].split("_norm"))
                if "tumo" in x
                else "_".join(np.array(x.split("_"))[[0, 3, 4, 1, 2]])
            ),
            inplace=True,
        )
        print(combined_acts.head())

    # Save results
    signature_order = combined_acts.keys()[orderSignatures(combined_acts.keys())]
    combined_acts[signature_order].to_csv(
        f"{output_dir}/Combined_Solution_Activities.tsv", sep="\t", index=True
    )

    # Save signature group mappings
    pd.DataFrame(
        grouped_signatures_cosmic.items(), columns=["new_signatures", "old_signatures"]
    ).to_csv(f"{output_dir}/Combined_Solution_Mapped.tsv", sep="\t", index=False)
