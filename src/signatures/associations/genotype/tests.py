import sys

import pandas as pd

if __name__ == "__main__":
    samples_file, signatures_file, targets_file, tests_file, tests_file_binary = (
        sys.argv[1:]
    )
    min_sig_count = 10
    min_tgt_count = 5

    # Load in signatures, targets and groups
    signatures = pd.read_csv(signatures_file, sep="\t").set_index("sample_id")
    targets = pd.read_csv(targets_file, sep="\t").set_index("sample_id")
    groups = pd.read_csv(samples_file, usecols=["sample_id", "group"], sep="\t")

    # Convert to binary
    signatures = (signatures > 0).astype(int).T
    targets = (targets > 0).astype(int).T

    # Number of non-zero activities for each signature-group pair
    signature_counts = (
        pd.merge(
            signatures.T,
            groups.set_index("sample_id"),
            left_index=True,
            right_index=True,
            how="inner",
        )
        .groupby("group")
        .sum()
        .reset_index()
        .melt(id_vars=["group"], value_name="count_sig", var_name="signature")
    )

    # Number of mutated samples in each gene-group pair
    target_group = pd.merge(
        targets.T,
        groups.set_index("sample_id"),
        left_index=True,
        right_index=True,
        how="inner",
    ).groupby("group")
    target_counts = pd.merge(
        target_group.sum()
        .reset_index()
        .melt(id_vars=["group"], value_name="count_tgt", var_name="target"),
        target_group.count()
        .reset_index()
        .melt(id_vars=["group"], value_name="total_tgt", var_name="target"),
    )

    # Only perform tests with at least min_sig_count non-zero signatures and min_tgt_count mutated and non-mutated samples
    set_counts = pd.merge(signature_counts, target_counts, on="group", how="inner")
    set_counts = set_counts[
        (set_counts["count_sig"] > min_sig_count)
        & (set_counts["count_tgt"] > min_tgt_count)
        & (set_counts["total_tgt"] - set_counts["count_tgt"] > min_tgt_count)
    ]

    # Target is binomial with n total possible hits
    types = pd.DataFrame.from_dict(
        {
            key: "binomial"
            if key[:3] == "is_"
            else "poisson"
            if key[:2] == "n_"
            else ""
            for key in set_counts.target.unique()
        },
        orient="index",
        columns=["type"],
    )
    set_counts = pd.merge(
        set_counts, types, left_on="target", right_index=True, how="inner"
    )

    # Save test file
    print(f"{len(set_counts)} tests")
    header = ["signature", "group", "target", "type"]
    set_counts[header].sort_values(header).to_csv(tests_file, sep="\t", index=False)

    # Test file for non-zero logistic regression
    set_counts = set_counts[
        (set_counts["total_tgt"] - set_counts["count_sig"] > min_sig_count)
    ]
    print(f"{len(set_counts)} binary tests")
    header = ["signature", "group", "target"]
    set_counts[header].sort_values(header).to_csv(
        tests_file_binary, sep="\t", index=False
    )
