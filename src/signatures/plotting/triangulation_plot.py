import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import Normalize

from signatures.config import load_environment

load_environment()
FIGURE_DIR = os.getenv("FIGURE_DIR")

mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["font.family"] = "STIXGeneral"
plt.rc("axes", labelsize=16)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)

default_colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def publish_fig(filename, publish="./"):
    plt.savefig(
        f"{publish}/{filename}.png",
        bbox_inches="tight",
        dpi=300,
        facecolor="w",
        transparent=False,
    )
    plt.savefig(
        f"{publish}/{filename}.pdf",
        bbox_inches="tight",
        format="pdf",
        facecolor="w",
        transparent=False,
    )


if __name__ == "__main__":
    # # Generate the full figure
    # 1) Need a colourbar
    # 2) What is the figure showing? How are results pancancer selected?
    # 3) Add columns for genomic features to come.
    # 4) Change x-labels
    relabel_x = {
        "Chromatin Structure and Modification": "Chromatin struct. & mod.",
        "DNA polymerases ": "DNA polymerases",
        "Homologous recombination": "HR",
        "grade": "Grade",
        "stage_best": "Stage",
        "er_status": "ER",
        "pr_status": "PR",
        "triple_negative": "Triple negative",
        "survival": "Survival",
        "subclonal": "Subclonal",
        "late": "Late clonal",
        "is_WGD": "WGD",
        "n_chromothripsis": "Chromothripsis",
        "n_chromoplexy": "Chromoplexy",
        "n_tandem_duplication": "Tandem duplications",
        "n_kataegis": "Kataegis",
    }

    group_labels = {
        "UV": ["SBS7a", "SBS7b", "SBS7c", "SBS7d", "DBS1"],
        "Smoking": ["SBS4", "DBS2", "ID3"],
        "MMR": ["SBS44", "SBS26", "SBS15", "SBS57"],  # ,"ID12","DBS12"],
        "Clock": ["SBS1", "SBS5", "ID1", "ID2"],
        r"$POLE$": ["SBS10a", "SBS10b", "SBS28", "DBS10", "DBS3"],
        "APOBEC": ["SBS2", "SBS13"],
        "HRD": ["SBS3", "ID6", "SV3"],
    }

    result_types = ["gene", "genotype", "treatment", "clinical", "survival", "timer"]
    titles = {
        "gene": "Gene inactivation",
        "genotype": "Genomic alterations",
        "treatment": "Treatment",
        "clinical": "Status",
        "survival": "Survival",
        "timer": "Timing",
    }
    result_colors = {
        "gene": 0,
        "genotype": 0,
        "treatment": 0,
        "clinical": 0,
        "survival": 0,
        "timer": 0,
    }
    cmaps = {0: ["bwr", -10, 10], 1: ["PRGn", -5, 5], 2: ["RdBu", -1, 1]}
    triang_dict = {}
    for j, result_type in enumerate(result_types + ["index", "colour_id"]):
        triang_dict[result_type] = pd.read_csv(
            f"/re_gecip/shared_allGeCIPs/pancancer_signatures/results/figs/triangulation/{result_type}.tsv",
            sep="\t",
            index_col=0,
        )
        not_null_count = (~pd.isnull(triang_dict[result_type])).sum()
        triang_dict[result_type] = triang_dict[result_type][
            not_null_count[not_null_count != 0].index
        ]
        print(triang_dict[result_type].keys())
    number_notnull = np.zeros(len(triang_dict[result_types[0]]), dtype=int)
    for j, result_type in enumerate(result_types):
        number_notnull += pd.isnull(triang_dict[result_type]).sum(axis=1)
    sig_subset = number_notnull > 0

    fig, ax = plt.subplots(1, 1, figsize=(12, 24))
    fs = 12

    xpos = 0
    xtick_labels = []
    for j, result_type in enumerate(result_types):
        cmap_variables = cmaps[result_colors[result_type]]

        x = np.arange(triang_dict[result_type][sig_subset].shape[1])
        y = np.arange(triang_dict[result_type][sig_subset].shape[0])
        xx, yy = np.meshgrid(x, y)
        plt.scatter(
            xpos + xx.flatten() + 0.5,
            yy.flatten() + 0.5,
            c=np.array(triang_dict[result_type][sig_subset]).flatten(),
            cmap=cmap_variables[0],
            vmin=cmap_variables[1],
            vmax=cmap_variables[2],
        )
        ylim = [0, triang_dict[result_type].shape[0]]
        print(result_type, ylim)

        xtick_labels += list(
            triang_dict[result_type][sig_subset].rename(relabel_x, axis=1).keys()
        )
        plt.text(
            xpos + triang_dict[result_type][sig_subset].shape[1] / 2 - 0.5,
            len(y) + 0.5,
            titles[result_type],
            ha="left",
            va="bottom",
            rotation=45,
            fontsize=12,
        )

        xpos += triang_dict[result_type][sig_subset].shape[1]

        plt.plot([xpos, xpos], ylim, "-k", linewidth=2)

    ylim = [0, triang_dict[result_types[-1]].shape[0]]
    plt.ylim(ylim)

    plt.sca(ax)
    cax = fig.add_axes(
        (
            ax.get_position().x1 + 0.01,
            ax.get_position().y0,
            0.02,
            ax.get_position().height,
        )
    )
    cbar = fig.colorbar(
        cm.ScalarMappable(
            norm=Normalize(vmin=cmaps[0][1], vmax=cmaps[0][2]), cmap=cmaps[0][0]
        ),
        cax=cax,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(rf"Z-score aggregated across cohorts", fontsize=fs * 1.5)

    y = np.arange(triang_dict[result_types[-1]][sig_subset].shape[0])

    plt.sca(ax)
    ax.set_xticks(np.arange(xpos), minor=True)
    ax.set_yticks(np.arange(len(y)), minor=True)
    ax.grid(True, which="minor", axis="both", linestyle="-", color="k", alpha=0.5)

    ax.set_xticks(np.arange(xpos) + 0.5, minor=False)
    ax.set_yticks(np.arange(len(y)) + 0.5, minor=False)

    ax.set_xticklabels(xtick_labels, fontsize=16, minor=False, rotation=90)
    ax.set_yticklabels(
        triang_dict["index"][sig_subset].index, fontsize=10, minor=False, rotation=0
    )
    plt.xlim(0, xpos)
    plt.ylim(ylim)

    for i, col_i in enumerate(triang_dict["colour_id"].colour_id[sig_subset]):
        ax.get_yticklabels()[i].set_color(default_colours[col_i + 1])

    for label, gene_group in group_labels.items():
        intersection = np.intersect1d(
            triang_dict["index"][sig_subset].index, gene_group, return_indices=True
        )
        plt.text(
            -4,
            float(np.mean(intersection[1])) + 0.5,
            label,
            ha="center",
            va="center",
            fontsize=fs * 1.2,
            color=default_colours[
                triang_dict["colour_id"].colour_id[sig_subset].loc[gene_group].iloc[0]
                + 1
            ],
            rotation=45,
        )

    publish_fig("triangulation", publish=FIGURE_DIR)
