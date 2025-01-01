import os
import re

import matplotlib as mpl
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
    signatureRenamer,
)
from signatures.utils import BH_threshold, orderSignatures

load_environment()
RESULT_DIR = os.getenv("RESULT_DIR")
FIGURE_DIR = os.getenv("FIGURE_DIR")

mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["font.family"] = "STIXGeneral"
plt.rc("axes", labelsize=16)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)
default_colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

if __name__ == "__main__":
    run_dir = f"{RESULT_DIR}/genotypeSigs/genotype_ktg_dCRT/"

    # Load results tables
    results_zinb = pd.read_csv(f"{run_dir}/output/signature_target_assoc_nb.csv")
    results_log0 = pd.read_csv(f"{run_dir}/output/signature_target_assoc_logistic.csv")
    results_target = pd.read_csv(f"{run_dir}/output/target_target_assoc.csv")
    mock_match = np.vectorize(lambda x: bool(re.match(".*mock.*", x, re.IGNORECASE)))
    results_combined = pd.merge(
        results_zinb,
        results_log0,
        on=("target", "signature", "group"),
        suffixes=("_zinb", "_log0"),
    )
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

    # Get renaming dictionary
    combined_acts, combined_sigs, relabel_map = loadSignatures()
    rename_dict = signatureRenamer(relabel_map)
    results_combined.replace({"signature": rename_dict}, inplace=True)

    # Rename connective as sarcoma
    results_combined.group = results_combined.group.str.replace("Connective", "Sarcoma")

    results_combined.resample_pvalue += 1e-200

    # Get non-mock individual signature results
    subset = (
        (~mock_match(results_combined.target))
        & (~mock_match(results_combined.signature))
        & (results_combined.signature.map(lambda x: x[0] != "T"))
        & (results_combined.resample_pvalue > 0)
    )
    results_combined = results_combined[subset]

    # Aggregate results over groups
    results_combined["pv_cov"] = (
        results_combined["target_means_alt_zinb"]
    ) ** 2 / scipy.stats.chi2.isf(results_combined["resample_pvalue"], df=1)

    results_combined["pv_IV"] = 1 / results_combined["pv_cov"]
    results_combined["pv_IVWmean"] = (
        results_combined["target_means_alt_zinb"] * results_combined["pv_IV"]
    )
    results_groupsum = (
        results_combined[["signature", "target", "pv_IV", "pv_IVWmean"]]
        .groupby(["target", "signature"])
        .sum()
    )
    results_groupsum["target_agg_means_alt_zinb"] = (
        results_groupsum["pv_IVWmean"] / results_groupsum["pv_IV"]
    )
    results_groupsum["target_agg_covs_alt_zinb"] = 1 / results_groupsum["pv_IV"]
    results_groupsum["target_agg_z_alt_zinb"] = results_groupsum[
        "target_agg_means_alt_zinb"
    ] * np.sqrt(results_groupsum["pv_IV"])

    # Get signatures with at least association
    pv_sws = BH_threshold(np.array(results_combined.resample_pvalue), alpha=0.01)
    unique_signatures = results_combined.signature.unique()
    unique_signatures = unique_signatures[orderSignatures(unique_signatures)]
    unique_signatures = {
        mut_type: [sig for sig in unique_signatures if sig[: len(mut_type)] == mut_type]
        for mut_type in ["SBS", "DBS", "ID", "CN", "SV"]
    }

    # Get the number of studywide significant results with the same sign as the aggregated result
    results_combined = pd.merge(
        results_combined,
        results_groupsum[["target_agg_z_alt_zinb"]].reset_index(),
        on=["target", "signature"],
        how="inner",
    )
    significant_results_count = (
        results_combined[
            (results_combined["resample_pvalue"] <= pv_sws)
            & (
                results_combined["target_means_alt_zinb"]
                * results_combined["target_agg_z_alt_zinb"]
            )
            > 0
        ]
        .groupby(["target", "signature"])
        .size()
    )

    fig, axes = plt.subplots(
        len(unique_signatures),
        5,
        figsize=(10, 18),
        gridspec_kw={
            "height_ratios": [
                len(unique_signatures[mut_type]) for mut_type in unique_signatures
            ][::-1]
        },
        sharex="col",
        sharey="row",
    )
    plt.subplots_adjust(hspace=0.02, wspace=0.02)
    genotype_labels = {
        "is_WGD": "WGD",
        "n_chromothripsis": "Chromothripsis",
        "n_chromoplexy": "Chromoplexy",
        "n_tandem_duplication": "Tandem duplication",
        "n_kataegis": "Kataegis",
    }

    max_cohort = 10
    norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.5 + max_cohort, clip=True)
    discrete_cmap = mpl.colors.ListedColormap(
        cm.viridis_r(np.linspace(0, 1, max_cohort + 1))
    )
    mapper = cm.ScalarMappable(norm=norm, cmap=discrete_cmap)

    for j, genotype in enumerate(genotype_labels):
        for i, mut_type in enumerate(unique_signatures.keys()):
            z_max = 12

            results_genotype = results_groupsum.loc[genotype]
            results_genotype = pd.merge(
                results_genotype,
                pd.DataFrame(index=unique_signatures[mut_type]),
                left_index=True,
                right_index=True,
                how="right",
            )

            ax = axes[4 - i, j]
            plt.sca(ax)
            colors = mapper.to_rgba(
                pd.merge(
                    pd.DataFrame(significant_results_count.loc[genotype]),
                    pd.DataFrame(index=unique_signatures[mut_type]),
                    left_index=True,
                    right_index=True,
                    how="right",
                ).fillna(0)[0]
            )
            z_scores = results_genotype.loc[
                unique_signatures[mut_type], "target_agg_z_alt_zinb"
            ]
            out_of_bounds = np.abs(z_scores) > z_max
            im = plt.scatter(
                z_scores[~out_of_bounds],
                np.arange(len(unique_signatures[mut_type]))[~out_of_bounds],
                c=colors[~out_of_bounds],
                s=15,
                zorder=0,
            )
            plt.scatter(
                np.sign(z_scores)[out_of_bounds & (z_scores > 0)] * (z_max * 1.03),
                np.arange(len(unique_signatures[mut_type]))[
                    out_of_bounds & (z_scores > 0)
                ],
                c=colors[out_of_bounds & (z_scores > 0)],
                s=15,
                zorder=0,
                marker=">",
            )
            plt.scatter(
                np.sign(z_scores)[out_of_bounds & (z_scores < 0)] * (z_max * 1.03),
                np.arange(len(unique_signatures[mut_type]))[
                    out_of_bounds & (z_scores < 0)
                ],
                c=colors[out_of_bounds & (z_scores < 0)],
                s=15,
                zorder=0,
                marker="<",
            )

            ylim = [-0.5, len(unique_signatures[mut_type]) - 0.5]
            plt.ylim(ylim)
            plt.plot([0, 0], ylim, "--k", alpha=0.5)

            xlim = [-z_max * 1.06, z_max * 1.06]
            plt.xlim(xlim)

            ax.set_yticks(np.arange(len(unique_signatures[mut_type])))
            ax.set_yticklabels(unique_signatures[mut_type], rotation=0, fontsize=10)

            ax.set_yticks(np.arange(len(unique_signatures[mut_type])) + 0.5, minor=True)
            ax.grid(True, which="minor", axis="y", linestyle="-", color="k", alpha=0.3)

            if j == 0:
                sig_type = {
                    "SBS": "SBS288",
                    "DBS": "DBS78",
                    "ID": "ID83",
                    "CN": "CNV48",
                    "SV": "SV32",
                }[mut_type]
                for isig, sig in enumerate(unique_signatures[mut_type]):
                    if sig in list(relabel_map["degasperi"][sig_type].old):
                        ax.get_yticklabels()[isig].set_color(map_colors["degasperi"])
                    elif sig in list(relabel_map["novel"][sig_type].new):
                        ax.get_yticklabels()[isig].set_color(map_colors["novel"])

            if i == 0:
                plt.xlabel(f"$Z$")

        plt.title(genotype_labels[genotype])

    # Add a single colorbar to the right of all subplots
    cbar_ax = fig.add_axes([0.93, 0.3, 0.02, 0.4])  # [left, bottom, width, height]
    cbar = plt.colorbar(mapper, cax=cbar_ax)
    cbar.set_label("Number of cohorts")
    cbar.set_ticks(np.arange(11))

    publish_fig("genomic_alteration_assoc", publish=FIGURE_DIR)
