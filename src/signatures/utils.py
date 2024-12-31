import re

import numpy as np
import pandas as pd


def regexSearch(pattern, string, index=1):
    try:
        return re.search(pattern, string).group(index)
    except AttributeError:
        raise ValueError(f"Pattern {pattern} not found in {string}")


def orderSignatures(signatures):
    """orderSignatures - sort signatures into a sensible ordering
    1) mutation type - SBS, DBS, ID
    2) signature number - 1,2,3,...
    3) novel signature cohort - Bladder, Breast...
    4) signature letter - 7a,7b,7c
    """

    mut_type = np.array(
        [regexSearch("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig, 1) for sig in signatures]
    )
    mut_index = np.array(
        [regexSearch("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig, 2) for sig in signatures]
    ).astype(int)
    mut_exten = np.array(
        [regexSearch("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig, 3) for sig in signatures]
    )
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

    print(mut_type)

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
    )
    return order


def BenjiminiHochberg(pvalues, alpha=0.05):
    m = len(pvalues)
    k = np.arange(m)

    order = np.argsort(pvalues)
    passed = pvalues[order] <= alpha * k / m

    return order[passed]


def BH_threshold(pvalues, alpha=0.05):
    passed_indices = BenjiminiHochberg(np.array(pvalues), alpha=alpha)
    if len(passed_indices) > 0:
        return np.max(np.array(pvalues)[passed_indices] + 1e-10)
    else:
        return 0


def resultsQuery(
    results_df,
    include=[],
    drop=[],
    sorter="resample_pvalue",
    keys=[
        "target",
        "signature",
        "group",
        "resample_pvalue",
        "wilks_pvalue_zinb",
        "wilks_pvalue_log0",
        "model_alt_zinb",
        "target_means_alt_zinb",
        "target_covs_alt_zinb",
        "power80_zinb",
        "resample_power80",
    ],
):
    subset = np.zeros(len(results_df), dtype=bool)

    for triarchy in include:
        subset[
            ((results_df.group == triarchy[0]) | (triarchy[0] == ""))
            & ((results_df.signature == triarchy[1]) | (triarchy[1] == ""))
            & ((results_df.target == triarchy[2]) | (triarchy[2] == ""))
        ] = True

    for triarchy in drop:
        subset[
            ((results_df.group == triarchy[0]) | (triarchy[0] == ""))
            & ((results_df.signature == triarchy[1]) | (triarchy[1] == ""))
            & ((results_df.target == triarchy[2]) | (triarchy[2] == ""))
        ] = False

    return results_df[keys][subset].sort_values(sorter)
