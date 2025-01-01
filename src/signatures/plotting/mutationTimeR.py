import os

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats
from matplotlib import cm

from signatures.config import load_environment
from signatures.plotting.combinedSignatures import (
    loadSignatures,
    map_colors,
    publish_fig,
)
from signatures.utils import orderSignatures, regexSearch

load_environment()
DATA_DIR = os.getenv("DATA_DIR")
RESULT_DIR = os.getenv("RESULT_DIR")
FIGURE_DIR = os.getenv("FIGURE_DIR")

mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["font.family"] = "STIXGeneral"
plt.rc("axes", labelsize=16)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)

default_colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

sig_rename = {"SBS4608": "SBS288", "SBS288": "SBS288", "DBS78": "DBS78", "ID83": "ID83"}


def loadMatrices(
    sample_list,
    matrix_dir=f"{DATA_DIR}/MutationTimeRmatrices/output/matrix_dir/output/",
    mut_times=["clonal_[NA]", "clonal_[early]", "clonal_[late]", "subclonal"],
    sig_types=["SBS288", "DBS78", "ID83"],
    groups=[],
):
    full_matrix_df = {
        sig_type: {mut_time: pd.DataFrame() for mut_time in mut_times}
        for sig_type in sig_types
    }
    for sig_type in sig_types:
        print(sig_type)
        mut_type = regexSearch("([A-Z]+)[0-9]+", sig_type, 1)
        for mut_time in mut_times:
            print(mut_time)
            if len(groups) > 0:
                for group in groups:
                    # Load in matrix if exists
                    try:
                        matrix_df = pd.read_csv(
                            f"{matrix_dir}/{mut_type}/{group}_{mut_time}.{sig_type}.all",
                            sep="\t",
                            index_col=0,
                        )
                    except pd.errors.EmptyDataError:
                        continue
                    if sig_type == "SBS4608":
                        matrix_df.set_index(
                            matrix_df.index.map(
                                lambda x: regexSearch(
                                    "[A-Z]:[ACTG]([ACTG]\[[ACTG]>[ACTG]\][ACTG])[ACTG]",
                                    x,
                                    1,
                                )
                            ),
                            inplace=True,
                        )
                        matrix_df = matrix_df.groupby("MutationType").sum()

                    # Add matrix
                    full_matrix_df[sig_type][mut_time] = pd.concat(
                        (full_matrix_df[sig_type][mut_time], matrix_df.T)
                    )

                # Fill in missing samples
                full_matrix_df[sig_type][mut_time] = (
                    pd.merge(
                        sample_list[[]],
                        full_matrix_df[sig_type][mut_time],
                        left_index=True,
                        right_index=True,
                        how="left",
                    )
                    .fillna(0)
                    .astype(int)
                )
            else:
                if sig_type == "SBS288":
                    matrix_df = pd.read_csv(
                        f"{matrix_dir}/count_matrix_{mut_time}_{sig_type}.tsv",
                        sep="\t",
                        index_col=0,
                    )
                    matrix_df.set_index(
                        matrix_df.index.map(
                            lambda x: regexSearch(
                                "[A-Z]:([ACTG]\[[ACTG]>[ACTG]\][ACTG])", x, 1
                            )
                        ),
                        inplace=True,
                    )
                    matrix_df = matrix_df.groupby("MutationType").sum()
                    full_matrix_df[sig_type][mut_time] = matrix_df.T
                else:
                    full_matrix_df[sig_type][mut_time] = pd.read_csv(
                        f"{matrix_dir}/count_matrix_{mut_time}_{sig_type}.tsv",
                        sep="\t",
                        index_col=0,
                    ).T

        idx = full_matrix_df[sig_type]["clonal_[NA]"].index
        full_matrix_df[sig_type]["clonal"] = (
            full_matrix_df[sig_type]["clonal_[NA]"].loc[idx]
            + full_matrix_df[sig_type]["clonal_[early]"].loc[idx]
            + full_matrix_df[sig_type]["clonal_[late]"].loc[idx]
        )

    return full_matrix_df


def getSigMutationRates(
    sample_list,
    combined_acts,
    combined_sigs,
    full_matrix_df,
    mut_times=["clonal", "subclonal", "clonal_[early]", "clonal_[late]"],
    numerator="subclonal",
    denominator="clonal",
    sig_types=["SBS288", "DBS78", "ID83"],
):
    # Generate fraction matrix for each mutation type
    R = {key: pd.DataFrame() for key in mut_times}
    # Generate fraction matrix for each mutation type
    R_sample = {key: pd.DataFrame() for key in mut_times}
    for sig_type in sig_types:
        acts = combined_acts[sig_rename[sig_type]].copy()
        acts.set_index(
            acts.index.map(lambda x: "_".join(x.split("_")[1:3])), inplace=True
        )

        sigs = combined_sigs[sig_rename[sig_type]][acts.keys()].copy()

        # platekeys = np.intersect1d(sample_list[sample_list.purity>0.7].index, acts.index)
        platekeys = np.intersect1d(sample_list.index, acts.index)
        channels = full_matrix_df[sig_type][numerator].keys()
        signatures = sigs.keys()

        # Matrix to transform from platekeys to tumour groups by summing over samples
        groups, group_idx = np.unique(
            sample_list.loc[platekeys].group, return_inverse=True
        )
        group_matrix = np.zeros((len(groups), len(platekeys)), dtype=int)
        group_matrix[group_idx, np.arange(len(platekeys))] = 1

        # Get array of signatures and activities
        acts = np.array(acts.loc[platekeys])
        sigs = np.array(sigs.loc[channels])

        # Fraction of each mutation channel in each samples which is caused by the signature
        f = (acts[:, :, None] * sigs.T[None, :, :]) / (acts @ sigs.T)[:, None, :]

        for timer_type in R_sample.keys():
            # Stacked mutations for clonal and subclonal
            muts = np.array(
                full_matrix_df[sig_type][timer_type].loc[platekeys, channels]
            )
            R_sample[timer_type] = pd.merge(
                R_sample[timer_type],
                pd.DataFrame(
                    (f @ muts[:, :, None])[:, :, 0], columns=signatures, index=platekeys
                ),
                left_index=True,
                right_index=True,
                how="outer",
            )

        # Sum over groups
        acts = group_matrix @ acts
        # Fraction of each mutation channel in each samples which is caused by the signature
        f = (acts[:, :, None] * sigs.T[None, :, :]) / (acts @ sigs.T)[:, None, :]

        for timer_type in R.keys():
            # Stacked mutations for clonal and subclonal
            muts = group_matrix @ np.array(
                full_matrix_df[sig_type][timer_type].loc[platekeys, channels]
            )
            R[timer_type] = pd.merge(
                R[timer_type],
                pd.DataFrame(
                    (f @ muts[:, :, None])[:, :, 0], columns=signatures, index=groups
                ),
                left_index=True,
                right_index=True,
                how="outer",
            )

    for timer_type in R.keys():
        R[timer_type] = R[timer_type].T

    return R, R_sample


def getFractions(
    sample_list,
    combined_acts,
    R,
    R_sample,
    numerator="subclonal",
    denominator="clonal",
    sig_types=["SBS288", "DBS78", "ID83"],
):
    groups = np.unique(sample_list.group)

    nonzero_counts = (
        pd.merge(
            sample_list[["group"]],
            R_sample[numerator] + R_sample[denominator] > 0,
            left_index=True,
            right_index=True,
            how="inner",
        )
        .groupby("group")
        .sum()
    )

    rank_dfs = {
        sig_type: pd.DataFrame(index=combined_acts[sig_rename[sig_type]].keys())
        for sig_type in sig_types
    }
    stat_dfs = {
        sig_type: pd.DataFrame(
            np.zeros((len(combined_acts[sig_rename[sig_type]].keys()), len(groups))),
            index=combined_acts[sig_rename[sig_type]].keys(),
            columns=groups,
        )
        for sig_type in sig_types
    }
    pv_dfs = {
        sig_type: pd.DataFrame(
            np.zeros((len(combined_acts[sig_rename[sig_type]].keys()), len(groups))),
            index=combined_acts[sig_rename[sig_type]].keys(),
            columns=groups,
        )
        for sig_type in sig_types
    }

    for sig_type in sig_types:
        print(sig_type)

        signatures = combined_acts[sig_rename[sig_type]].keys()

        for group in groups:
            group_platekeys = np.array(
                sample_list.reset_index().set_index("group").loc[group].platekey
            )
            frac = R_sample[numerator].loc[group_platekeys] / (
                R_sample[numerator].loc[group_platekeys]
                + R_sample[denominator].loc[group_platekeys]
            )
            frac = frac[signatures]

            # Get non-nan fractions
            nonnan_frac = np.array(frac)[~np.isnan(np.array(frac))]

            # Get signature labels for fractions
            sig_labels = np.meshgrid(frac.keys(), frac.index)[0]
            sig_labels = sig_labels[~np.isnan(np.array(frac))]

            # Get sig_labels ordered by fraction
            sig_labels_ordered = sig_labels[np.argsort(nonnan_frac)]

            group_rank_df = (
                pd.DataFrame(
                    {
                        "sig": sig_labels_ordered,
                        group: np.arange(1, len(sig_labels) + 1)
                        / (len(sig_labels) + 1),
                    }
                )
                .groupby("sig")
                .mean()
            )

            rank_dfs[sig_type] = pd.merge(
                rank_dfs[sig_type],
                group_rank_df,
                left_index=True,
                right_index=True,
                how="left",
            )

            for sig in signatures:
                # This is a significance test
                if np.sum(sig_labels == sig) > 0:
                    mwu = scipy.stats.mannwhitneyu(
                        nonnan_frac[sig_labels == sig], nonnan_frac[sig_labels != sig]
                    )
                    stat_dfs[sig_type].loc[sig, group] = mwu.statistic
                    pv_dfs[sig_type].loc[sig, group] = mwu.pvalue
                else:
                    stat_dfs[sig_type].loc[sig, group] = np.nan
                    pv_dfs[sig_type].loc[sig, group] = np.nan

                # Also need some kind of z-score to plot!

    return nonzero_counts, rank_dfs, stat_dfs, pv_dfs


def plotFractions(
    sample_list,
    nonzero_counts,
    combined_acts,
    R,
    R_sample,
    rank_dfs,
    pv_dfs,
    stat_dfs,
    numerator="subclonal",
    denominator="clonal",
    relabel_map=None,
    label="Subclonal fraction",
):
    # High count signatures
    high_count_signatures = nonzero_counts.keys()[np.sum(nonzero_counts > 20) > 0]
    high_count_groups = np.unique(sample_list.group)[
        np.unique(sample_list.group, return_counts=True)[1] > 50
    ]

    print(
        [
            len(
                np.intersect1d(
                    high_count_signatures, combined_acts[sig_rename[sig_type]].keys()
                )
            )
            for sig_type in rank_dfs
        ]
    )

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(len(high_count_signatures) / 2, 2 * len(high_count_groups) / 2),
        gridspec_kw={
            "width_ratios": [
                len(
                    np.intersect1d(
                        high_count_signatures,
                        combined_acts[sig_rename[sig_type]].keys(),
                    )
                )
                for sig_type in rank_dfs
            ]
        },
        sharey=True,
    )
    fig.subplots_adjust(right=0.95, wspace=0.05, hspace=0.2)
    fs = 24

    size_factor = lambda x: np.log10(x + 1) ** 1.5 * 70

    # sample the colormaps that you want to use. Use 128 from each so we get 256
    # colors in total
    PVMIN = -np.log10(0.05)
    PVMAX = 10
    colors2 = cm.bwr(
        np.linspace(0.1, 0.9, 128)
    )  # plt.cm.bone_r(np.linspace(0., 1, 128))
    colors2[
        int(np.floor((1 - PVMIN / PVMAX) / 2 * 128)) : int(
            np.ceil((1 + PVMIN / PVMAX) / 2 * 128)
        )
    ] = [
        0.8,
    ] * 4
    pv_map = mcolors.LinearSegmentedColormap.from_list("my_colormap", colors2)

    for j in range(2):
        for i, sig_type in enumerate(list(rank_dfs.keys())):
            ax = axes[j, i]
            plt.sca(ax)

            signatures = np.intersect1d(
                high_count_signatures, combined_acts[sig_rename[sig_type]].keys()
            )
            signatures = signatures[orderSignatures(signatures)]

            R_num = R[numerator].loc[signatures, high_count_groups].copy()
            R_den = R[denominator].loc[signatures, high_count_groups].copy()
            rank_array = np.array(
                rank_dfs[sig_type].loc[signatures, high_count_groups].copy()
            )
            sizes = np.array(nonzero_counts.T.loc[signatures, high_count_groups])

            pv_array = np.array(
                -np.log10(pv_dfs[sig_type].loc[signatures, high_count_groups].copy())
            ) * np.sign(rank_array - 0.5)

            x_labels = signatures
            y_labels = high_count_groups

            xx, yy = np.meshgrid(np.arange(len(x_labels)), np.arange(len(y_labels)))
            im = plt.scatter(
                xx.flatten() + 0.5,
                yy.flatten() + 0.5,
                s=size_factor(sizes).T.flatten(),
                c=np.array(R_num / (R_num + R_den)).T.flatten()
                if j == 0
                else pv_array.T.flatten(),
                vmin=0.0 if j == 0 else -PVMAX,
                vmax=1.0 if j == 0 else PVMAX,
                cmap="viridis_r" if j == 0 else pv_map,
            )

            ax.set_yticks(np.arange(len(y_labels)) + 0.5, minor=False)
            ax.set_xticks(np.arange(len(x_labels)) + 0.5, minor=False)
            ax.set_yticklabels(y_labels, rotation=0, fontsize=fs, minor=False)
            ax.set_xticklabels(x_labels, rotation=90, fontsize=fs, minor=False)
            ax.set_yticks(np.arange(len(y_labels)), minor=True)
            ax.set_xticks(np.arange(len(x_labels)), minor=True)
            ax.grid(
                True, which="minor", axis="both", linestyle="-", color="k", alpha=0.5
            )

            ax.set_ylim(0, len(y_labels))
            ax.set_xlim(0, len(x_labels))

            # add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
            cb_ax = (
                fig.add_axes([0.96, 0.55, 0.015, 0.3])
                if j == 0
                else fig.add_axes([0.96, 0.14, 0.015, 0.3])
            )

            # cbar = plt.colorbar()
            cbar = fig.colorbar(im, cax=cb_ax)
            # cbar.set_label(r"$\mathrm{sgn}(\beta) \,\times\,-\log_{10}(\mathrm{pvalue})$", fontsize=36)
            cbar.set_label(
                label if j == 0 else r"$-\log_{10}(P\mathrm{value})$", fontsize=40
            )
            cbar.ax.tick_params(labelsize=fs * 1.2)

            if relabel_map is not None:
                for isig, sig in enumerate(signatures):
                    if sig in list(relabel_map["degasperi"][sig_rename[sig_type]].old):
                        ax.get_xticklabels()[isig].set_color(map_colors["degasperi"])
                    elif sig in list(relabel_map["novel"][sig_rename[sig_type]].new):
                        ax.get_xticklabels()[isig].set_color(map_colors["novel"])

            if (i == 0) & (j == 0):
                print("Adding legend")
                l1 = plt.scatter([], [], s=size_factor(1), c="k", alpha=0.7)
                l2 = plt.scatter([], [], s=size_factor(10), c="k", alpha=0.7)
                l3 = plt.scatter([], [], s=size_factor(100), c="k", alpha=0.7)
                l4 = plt.scatter([], [], s=size_factor(1000), c="k", alpha=0.7)
                legend = ax.legend(
                    [l1, l2, l3, l4],
                    [f"$1$", f"$10$", f"$100$", f"$1000$"],
                    loc="lower left",
                    bbox_to_anchor=(0, 1),
                    fontsize=fs,
                    ncol=4,
                    frameon=False,
                    title="Number of samples",
                )
                legend.get_title().set_fontsize("29")


if __name__ == "__main__":
    # Combined signatures as reference
    combined_acts, combined_sigs, relabel_map = loadSignatures()

    # Load in samples
    sample_list = pd.read_csv(
        f"{RESULT_DIR}/signatures/MutationTimeRmatricesVCFs/input/manifest.tsv",
        sep="\t",
        names=["platekey", "group", "vcf", "timer_vcf", "id"],
    )
    # Replace Connective with Sarcoma
    sample_list.group = sample_list.group.str.replace("Connective", "Sarcoma")
    groups = np.unique(sample_list.group)
    sample_list.set_index("platekey", inplace=True)

    # Load mutation matrices for MutationTimeR groups
    full_matrix_df = loadMatrices(
        sample_list,
        sig_types=["SBS288", "DBS78", "ID83"],
        matrix_dir=f"{RESULT_DIR}/signatures/MutationTimeRmatricesVCFs/output/",
    )

    # Transform mutation matrices from mutation classes to signatures
    R, R_sample = getSigMutationRates(
        sample_list, combined_acts, combined_sigs, full_matrix_df
    )

    # Get rank statistics of signatures in groups
    nonzero_counts = {}
    rank_dfs = {}
    stat_dfs = {}
    pv_dfs = {}
    (
        nonzero_counts["subclonal"],
        rank_dfs["subclonal"],
        stat_dfs["subclonal"],
        pv_dfs["subclonal"],
    ) = getFractions(
        sample_list,
        combined_acts,
        R,
        R_sample,
        numerator="subclonal",
        denominator="clonal",
    )
    # Get rank statistics of signatures in groups
    nonzero_counts["late"], rank_dfs["late"], stat_dfs["late"], pv_dfs["late"] = (
        getFractions(
            sample_list,
            combined_acts,
            R,
            R_sample,
            numerator="clonal_[late]",
            denominator="clonal_[early]",
        )
    )

    # Save a single table
    timer_df = pd.DataFrame()
    for sig_type in ["SBS288", "DBS78", "ID83"]:
        pv_df = (
            pv_dfs["subclonal"][sig_type]
            .stack()
            .reset_index()
            .rename(
                {"level_0": "signature", "level_1": "group", 0: "pv_subclonal"}, axis=1
            )
        )
        rank_df = (
            rank_dfs["subclonal"][sig_type]
            .stack()
            .reset_index()
            .rename(
                {"level_0": "signature", "level_1": "group", 0: "rank_subclonal"},
                axis=1,
            )
        )
        pv_df = pd.merge(
            pv_df,
            pv_dfs["late"][sig_type]
            .stack()
            .reset_index()
            .rename({"level_0": "signature", "level_1": "group", 0: "pv_late"}, axis=1),
            on=["signature", "group"],
            how="inner",
        )
        rank_df = pd.merge(
            rank_df,
            rank_dfs["late"][sig_type]
            .stack()
            .reset_index()
            .rename(
                {"level_0": "signature", "level_1": "group", 0: "rank_late"}, axis=1
            ),
            on=["signature", "group"],
            how="inner",
        )
        timer_df = pd.concat(
            (timer_df, pd.merge(pv_df, rank_df, on=["signature", "group"], how="inner"))
        )
    timer_df.to_csv(f"{RESULT_DIR}/exportData/mutationTimeR.tsv", sep="\t", index=False)

    # Plot clonal fraction and rank statistics
    plotFractions(
        sample_list,
        nonzero_counts["subclonal"],
        combined_acts,
        R,
        R_sample,
        rank_dfs["subclonal"],
        pv_dfs["subclonal"],
        stat_dfs["subclonal"],
        numerator="subclonal",
        denominator="clonal",
        relabel_map=relabel_map,
        label="Subclonal fraction",
    )
    publish_fig("subclonal_fraction_rank_mutationtimer", publish=FIGURE_DIR)

    # Plot clonal fraction and rank statistics
    plotFractions(
        sample_list,
        nonzero_counts["late"],
        combined_acts,
        R,
        R_sample,
        rank_dfs["late"],
        pv_dfs["late"],
        stat_dfs["late"],
        numerator="clonal_[late]",
        denominator="clonal_[early]",
        relabel_map=relabel_map,
        label="Late clonal fraction",
    )
    publish_fig("lateclonal_fraction_rank_mutationtimer", publish=FIGURE_DIR)
