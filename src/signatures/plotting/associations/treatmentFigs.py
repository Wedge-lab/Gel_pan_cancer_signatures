import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats
from dotenv import load_dotenv
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable

from signatures.plotting.combinedSignatures import (
    loadSignatures,
    publish_fig,
    signatureRenamer,
)
from signatures.utils import BenjiminiHochberg

load_dotenv()
RESULT_DIR = os.getenv("RESULT_DIR")
FIGURE_DIR = os.getenv("FIGURE_DIR")

mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["font.family"] = "STIXGeneral"
plt.rc("axes", labelsize=16)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)

default_colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def correlationGrid(tgt_df, target_subset, p_threshold):
    # Generate target grid
    unique_targets, inverse_target = np.unique(tgt_df.target, return_inverse=True)
    unique_dept, inverse_dept = np.unique(tgt_df.dependent, return_inverse=True)

    target_grid = np.zeros((len(unique_targets), len(unique_dept)))

    target_grid[inverse_target, inverse_dept] = np.where(
        tgt_df.wilks_pvalue < p_threshold, tgt_df.target_means_alt, 0
    )
    # target_grid[inverse_target, inverse_dept] = np.sign(tgt_df[tgt_subset].target_means_alt)\
    #                                                    *(-np.log10(tgt_df[tgt_subset].wilks_pvalue))
    target_grid = pd.DataFrame(target_grid, index=unique_targets, columns=unique_dept)

    # Insert missing columns
    target_grid[np.setdiff1d(target_subset, target_grid.keys())] = np.nan
    # Insert missing rows
    target_grid = target_grid.T
    target_grid[np.setdiff1d(target_subset, target_grid.keys())] = np.nan
    target_grid = target_grid.T

    return target_grid


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


if __name__ == "__main__":
    run_name = "treatment_dCRT"
    results_dir = f"{RESULT_DIR}/results/associations/treatmentSigs/"

    # Rename list of treatments
    treatment_list = {
        "rtprescribeddose": "Radiotherapy",
        "RADIOTHERAPY": "Radiotherapy",
        "OXALIPLATIN": "OXALIPLATIN",
        "CARBOPLATIN": "CARBOPLATIN",
        "CISPLATIN": "CISPLATIN",
        "FLUOROURACIL": "FLUOROURACIL",
        "CYCLOPHOSPHAMIDE": "CYCLOPHOSPHAMIDE",
        "EPIRUBICIN": "EPIRUBICIN",
        "BEVACIZUMAB": "BEVACIZUMAB",
        "BOSUTINIB": "BOSUTINIB",
        "CAPECITABINE": "CAPECITABINE",
        "CETUXIMAB": "CETUXIMAB",
        "DASATINIB": "DASATINIB",
        "DOCETAXEL": "DOCETAXEL",
        "DOXORUBICIN": "DOXORUBICIN",
        "FOLINIC_ACID": "FOLINIC ACID",
        "GEMCITABINE": "GEMCITABINE",
        "HYDROXYCARBAMIDE": "HYDROXYCARBAMIDE",
        "IFOSFAMIDE": "IFOSFAMIDE",
        "IMATINIB": "IMATINIB",
        "IRINOTECAN": "IRINOTECAN",
        "PACLITAXEL": "PACLITAXEL",
        "PEMETREXED": "PEMETREXED",
        "PERTUZUMAB": "PERTUZUMAB",
        "TEMOZOLOMIDE": "TEMOZOLOMIDE",
        "TRASTUZUMAB": "TRASTUZUMAB",
    }
    treatment_list = {
        key: treatment_list[key][0].upper() + treatment_list[key][1:].lower()
        for key in treatment_list
    }

    # Input data
    sample_df = pd.read_csv(
        f"{results_dir}/{run_name}/input/samples.tsv", sep="\t"
    ).set_index("sample_id")
    target_df = pd.read_csv(
        f"{results_dir}/{run_name}/input/targets.tsv", sep="\t"
    ).set_index("sample_id")
    test_df = pd.read_csv(f"{results_dir}/{run_name}/input/tests.tsv", sep="\t")
    signature_df = pd.read_csv(
        f"{results_dir}/{run_name}/input/signatures.tsv", sep="\t"
    ).set_index("sample_id")

    results_zinb = pd.read_csv(
        f"{results_dir}/{run_name}/output/signature_target_assoc_nb.csv"
    )
    results_log0 = pd.read_csv(
        f"{results_dir}/{run_name}/output/signature_target_assoc_logistic.csv"
    )
    results_target = pd.read_csv(
        f"{results_dir}/{run_name}/output/target_target_assoc.csv"
    )

    mock_match = np.vectorize(lambda x: bool(re.match(".*mock.*", x, re.IGNORECASE)))
    results_combined = pd.merge(
        results_zinb,
        results_log0,
        on=("target", "signature", "group"),
        suffixes=("_zinb", "_log0"),
    )  # .replace({'group':group_labels})

    results_mock = results_combined[
        (mock_match(results_combined.target)) | (mock_match(results_combined.signature))
    ].copy()

    results_combined["zscore"] = results_combined.target_means_alt_zinb / np.sqrt(
        results_combined.target_covs_alt_zinb
    )
    results_combined["log10_pvalue_zinb"] = -np.log10(results_combined.resample_pvalue)
    results_combined["log10_pvalue_log0"] = -np.log10(
        results_combined.wilks_pvalue_log0
    )
    results_combined["exp_means_alt_zinb"] = np.exp(
        results_combined["target_means_alt_zinb"]
    )
    results_combined = pd.merge(
        results_combined,
        pd.DataFrame(treatment_list.keys(), columns=["target"]),
        on="target",
        how="inner",
    ).replace({"target": treatment_list})

    # Get renaming dictionary
    combined_acts, combined_sigs, relabel_map = loadSignatures()
    rename_dict = signatureRenamer(relabel_map)
    results_combined.replace({"signature": rename_dict}, inplace=True)

    # Rename connective as sarcoma
    results_combined.group = results_combined.group.str.replace("Connective", "Sarcoma")

    ####### -------- Z-score - Q plot ------- #########
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    xlim = [-7, 7]

    subset = (
        (results_combined.resample_pvalue > 0)
        & (np.abs(results_combined.target_means_alt_zinb) < 1000)
        & (results_combined.signature.map(lambda x: x[0] != "T"))
        & ~mock_match(results_combined.target)
    )
    # (results_combined.model_alt_zinb.isin(['negbin', 'glm.binom', 'binomial', 'poisson']))

    passed_indices = BenjiminiHochberg(
        np.array(results_combined["resample_pvalue"]), 0.001
    )
    p_threshold = np.max(np.array(results_combined["resample_pvalue"])[passed_indices])

    # zscore = results_combined.target_means_alt/np.sqrt(results_combined.target_covs_alt)
    plt.scatter(
        results_combined.target_means_alt_zinb[subset],
        -np.log10(results_combined.resample_pvalue[subset]),
        c=-np.log10(results_combined.wilks_pvalue_log0[subset]),
        s=4,
        cmap="viridis_r",
        vmax=12,
    )

    ylim = (0, ax.axes.get_ylim()[1])

    plt.plot(xlim, [-np.log10(p_threshold), -np.log10(p_threshold)], "--r")
    plt.plot([0, 0], [-np.log10(p_threshold), ylim[1]], "--r")

    plt.xlim(xlim)
    plt.gca().set_ylim(bottom=0, top=ylim[1])

    cbar = plt.colorbar()
    cbar.set_label(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{logistic})$")

    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{ZINB-dCRT})$")

    pos_sets = {
        "target,group": [],
        "signature,target": [
            ("DBS2", "Radiotherapy", p_threshold, 0, -1.2),
            ("ID5", "Radiotherapy", p_threshold, 0, 4),
            ("ID8", "Radiotherapy", p_threshold, 2, -1),
        ],
        "signature,group": [("DBS5", "ColoRect-AdenoCA", 1e-10, 3, -0.5)],
    }
    for key, combos in pos_sets.items():
        for i, set_row in enumerate(combos):
            print(set_row)
            rset = results_combined[
                subset
                & (results_combined[key.split(",")[0]] == set_row[0])
                & (results_combined[key.split(",")[1]] == set_row[1])
                & (results_combined.resample_pvalue < set_row[2])
            ]
            posy = np.mean(-np.log10(rset.resample_pvalue)) + set_row[3]
            if set_row[4] < 0:
                ha = "right"
            else:
                ha = "left"
            plt.text(
                set_row[4], posy, f"{set_row[0]}, {set_row[1]}", ha=ha, va="center"
            )
            for index, row in rset.iterrows():
                plt.plot(
                    [set_row[4], row.target_means_alt_zinb],
                    [posy, -np.log10(row.resample_pvalue)],
                    ":k",
                    alpha=0.3,
                )

    # Highlight hits
    select = [
        ("Oxaliplatin", "DBS5", "ColoRect-AdenoCA", 0, -0.5),
        ("Cisplatin", "DBS5", "Lung-AdenoCA", 3, -2.8),
    ]
    for i, set_row in enumerate(select):
        row = results_combined[
            (results_combined.target == set_row[0])
            & (results_combined.signature == set_row[1])
            & (results_combined.group == set_row[2])
        ].iloc[0]
        pos = (row.target_means_alt_zinb, -np.log10(row.resample_pvalue))
        text_pos = (pos[0] + set_row[4], pos[1] + set_row[3])
        plt.text(
            text_pos[0],
            text_pos[1],
            f"{row.signature} ({row.target}) \n{row.group}",
            ha="right" if set_row[4] < 0 else "left",
            va="center",
        )
        plt.plot([pos[0], text_pos[0]], [pos[1], text_pos[1]], ":k", alpha=0.3)

    ellipse = Ellipse(
        (5.2, 38),
        width=3.6,
        height=42,
        facecolor="none",
        edgecolor="grey",
        linestyle="--",
    )
    ax.add_patch(ellipse)
    plt.text(5.2, 45, "Breast-DuctalCA,\nSBS26/SBS44", ha="center", va="center")

    publish_fig("sig_treatment_bqplot", publish=publish)

    ###### ------ Mock associations ------ ######
    subset = (results_mock.resample_pvalue > 0) & (
        results_mock.model_alt_zinb.isin(["negbin", "glm.binom", "binomial"])
    )

    mock_size = np.vectorize(lambda x: bool(re.match("MOCK[0-9]+", x)))
    subsets = [
        subset & (results_mock.signature == "SBSmock"),
        subset & mock_size(results_mock.target),
    ]
    titles = ["Signature Mock", "Target Resample", "Target Perturbation"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

    ymax = 10

    for i, subset in enumerate(subsets):
        plt.sca(axes[i])

        n_test = np.sum(subset)
        k = np.arange(1, n_test + 1)

        pv_expected = np.linspace(0, 1 - 1 / n_test, n_test) + 1 / (2 * n_test)
        pv_expected = np.arange(1 / n_test, 1 + 1e-10, 1 / n_test)
        percentiles = scipy.stats.beta.ppf(
            np.array([0.05, 0.5, 0.95])[:, None], k[None, :], n_test + 1 - k[None, :]
        )

        Y = np.sort(-np.log10(results_mock[subset]["wilks_pvalue_zinb"]))
        Y[Y > ymax] = ymax + 0.01
        plt.scatter(
            -np.log10(percentiles[1][::-1])[Y < ymax],
            Y[Y < ymax],
            s=20,
            label=f"Wilks",
            c=default_colours[3],
            marker=".",
        )  # , marker=markers[i])
        plt.scatter(
            -np.log10(percentiles[1][::-1])[Y > ymax],
            Y[Y > ymax],
            s=20,
            c=default_colours[3],
            marker="^",
        )

        Y = np.sort(-np.log10(results_mock[subset]["resample_pvalue"]))
        plt.scatter(
            -np.log10(percentiles[1][::-1]),
            Y,
            s=20,
            label=f"dCRT",
            c=default_colours[0],
            marker=".",
        )  # , marker=markers[i])
        plt.scatter(
            -np.log10(percentiles[1][::-1])[Y > ymax],
            np.zeros(np.sum(Y > ymax)) + ymax - 0.1,
            s=20,
            c=default_colours[0],
            marker="^",
        )

        k = np.arange(1, n_test + 1)
        percentiles = scipy.stats.beta.ppf(
            np.array([0.05, 0.5, 0.95])[:, None], k[None, :], n_test + 1 - k[None, :]
        )
        plt.plot(
            -np.log10(percentiles[1][::-1]), -np.log10(percentiles[1][::-1]), c="k"
        )
        plt.fill_between(
            -np.log10(percentiles[1][::-1]),
            -np.log10(percentiles[0][::-1]),
            -np.log10(percentiles[2][::-1]),
            alpha=0.3,
            color="k",
        )

        plt.gca().set_xlim(left=0)
        plt.gca().set_ylim(bottom=0)
        plt.gca().set_ylim(top=ymax + 0.1)

        plt.title(titles[i], fontsize=20, pad=10)

        if i == 1:
            plt.legend(loc="upper left")
            plt.xlabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{expected})$")
        if i == 0:
            plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{observed})$")

    plt.subplots_adjust(wspace=0.01)

    publish_fig("treatment_mock_qqplots_x3", publish=publish)

    ###### -------- Therapy-therapy correlations ------- ########
    tgt_df = pd.read_csv(f"{results_dir}/{run_name}/output/target_target_assoc.csv")
    tgt_df.wilks_pvalue[tgt_df.wilks_pvalue == 0] = 1e-300
    tgt_df = tgt_df

    tgt_df.replace(
        {"dependent": treatment_list, "target": treatment_list}, inplace=True
    )

    # Select subset to work with
    subset = (
        (results_combined.resample_pvalue > 0)
        & (results_combined.model_alt_zinb == "negbin")
        & (np.abs(results_combined.target_means_alt_zinb) < 1000)
        & results_combined.signature.map(lambda x: x[0] != "T")
        & (results_combined.target != "TGFBR2")
    )

    # Subset of Genes to plot
    target_subset = np.unique(
        results_combined.target[
            subset & (results_combined.resample_pvalue < p_threshold)
        ]
    )
    target_subset = target_subset[~mock_match(target_subset)]

    order_treat = list(treatment_list.values())

    target_indices = np.intersect1d(
        target_subset, list(treatment_list.values()), return_indices=True
    )
    target_subset = target_subset[np.argsort(target_indices[2])]

    fig, axes = plt.subplots(1, 3, figsize=(50, 17), sharey=True)
    fig.subplots_adjust(right=0.95, wspace=0.03, hspace=0.12)

    fs = 30

    groups = ["ColoRect-AdenoCA", "Breast-DuctalCA", "CNS-GBM-IDHwt"]
    for i, group in enumerate(groups):
        plt.sca(axes[i])

        p_threshold = np.max(
            np.array(tgt_df.wilks_pvalue)[
                BenjiminiHochberg(np.array(tgt_df.wilks_pvalue), 0.01)
            ]
        )

        tgt_subset = tgt_df.group == group
        target_grid = correlationGrid(
            tgt_df[tgt_df.group == group], target_subset, 0.01
        )[target_subset].loc[target_subset]

        im = plt.pcolor(target_grid, vmin=-5, vmax=5, cmap="bwr")

        plt.title(f"Treatment-treatment associations in {group}", fontsize=fs * 1.4)

        axes[i].set_xticks(np.arange(len(target_subset)) + 0.5, minor=False)
        axes[i].set_xticks(np.arange(len(target_subset)), minor=True)
        axes[i].set_yticks(np.arange(len(target_subset)) + 0.5, minor=False)
        axes[i].set_yticks(np.arange(len(target_subset)), minor=True)

        axes[i].set_xticklabels(
            target_subset, rotation=90, fontsize=fs * 1.2, minor=False
        )
        axes[i].set_yticklabels(target_subset, fontsize=fs * 1.2, minor=False)
        axes[i].grid(True, which="minor", axis="both", linestyle="-", color="k")

    # add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
    cb_ax = fig.add_axes([0.96, 0.15, 0.01, 0.7])

    cbar = fig.colorbar(im, cax=cb_ax)
    cbar.set_label(r"$\beta$", fontsize=fs * 1.2)
    cbar.ax.tick_params(labelsize=fs)

    publish_fig(f'treatment_grid_beta_{"-".join(groups)}', publish=FIGURE_DIR)

    #### ------ Treatment data figure ------ #######
    # Load input data
    sample_df = pd.read_csv(
        f"{results_dir}/{run_name}/input/samples.tsv", sep="\t"
    ).set_index("sample_id")
    target_df = pd.read_csv(
        f"{results_dir}/{run_name}/input/targets.tsv", sep="\t"
    ).set_index("sample_id")
    test_df = pd.read_csv(f"{results_dir}/{run_name}/input/tests.tsv", sep="\t")
    target_df = (target_df > 0).astype(int)
    signature_df = pd.read_csv(
        f"{results_dir}/{run_name}/input/signatures.tsv", sep="\t"
    ).set_index("sample_id")
    len(sample_df), len(target_df), len(signature_df)
    sample_df["group"] = sample_df.group.str.replace("Connective", "Sarcoma")

    target_df.rename(treatment_list, axis=1, inplace=True)

    mock_match = np.vectorize(lambda x: bool(re.match(".*mock.*", x, re.IGNORECASE)))

    new_grid = pd.merge(
        sample_df[["group"]], target_df, left_index=True, right_index=True, how="left"
    )

    treatments_ = new_grid.keys()[7:]
    treatments_ = treatments_[~mock_match(treatments_)]
    treatment_xm = np.intersect1d(
        treatments_, list(treatment_list.values()), return_indices=True
    )
    order = treatment_xm[1][np.argsort(treatment_xm[2])]
    treatments_ = treatments_[order]

    groups = np.unique(new_grid["group"])
    group_xm = np.intersect1d(groups, np.unique(new_grid.group), return_indices=True)
    order = group_xm[1][np.argsort(group_xm[2])]
    groups = groups[order]

    count_grid = np.zeros((len(treatments_), len(groups)))

    for i, treat in enumerate(treatments_):
        for j, group in enumerate(groups):
            count_grid[i, j] = np.sum(new_grid[treat][new_grid.group == group] > 0)

    low_count = (count_grid <= 5) & (count_grid > 0)
    count_grid[(count_grid <= 5) & (count_grid > 0)] = 5

    ### ----- Treatment event counts ---- ###
    subset = np.sum(count_grid > 5, axis=1) > 0
    fig, ax = plt.subplots(
        1, 1, figsize=(12, 12 * np.sum(subset) / count_grid.shape[1])
    )

    im = plt.pcolor(
        count_grid[subset][::-1], cmap="inferno_r", norm=LogNorm(), vmin=5, vmax=100
    )

    low_count_idx = np.unravel_index(
        np.argwhere(np.array(low_count)[subset][::-1]),
        np.array(low_count)[subset][::-1].shape,
    )[1]
    plt.scatter(
        low_count_idx[:, 1] + 0.5,
        low_count_idx[:, 0] + 0.5,
        marker="o",
        color="k",
        s=10,
    )

    fs = 16
    ax.set_yticks(np.arange(len(treatments_[subset]) + 1), minor=True)
    ax.set_xticks(np.arange(len(groups) + 1), minor=True)
    ax.grid(True, which="minor", axis="both", linestyle="-", color="k")

    ax.set_xticks(np.arange(len(groups)) + 0.5, minor=False)
    ax.set_xticklabels(groups, rotation=90, fontsize=fs, minor=False)
    ax.set_yticks(np.arange(len(treatments_[subset])) + 0.5, minor=False)
    ax.set_yticklabels(treatments_[subset][::-1], rotation=0, fontsize=fs, minor=False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad="2%")
    # norm = mpl.colors.Normalize(vmin=0,vmax=10)
    cbar = fig.colorbar(im, cax=cax, orientation="vertical", ticks=[1, 10, 100, 1000])
    cbar.set_label("Number of samples", fontsize=fs)
    # cbar.ax.set_yticks()

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=f"$\mathrm{{count}} \leq 5$",
            markerfacecolor="k",
            markersize=10,
        )
    ]
    plt.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.0, 1.0),
        loc="lower right",
        fontsize=20,
        frameon=False,
    )
    publish_fig(f"treatment_counts", publish=FIGURE_DIR)
