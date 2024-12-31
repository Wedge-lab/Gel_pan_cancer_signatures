import os
import pandas as pd, numpy as np, re, tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from signatures.plotting.combinedSignatures import publish_fig, sig_dirs
from signatures.utils import regexSearch
from dotenv import load_dotenv

load_dotenv()
FIGURE_DIR = os.getenv("FIGURE_DIR")
DATA_DIR = os.getenv("DATA_DIR")

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)

# # AIC figures
def plotAIC(cohort_stats, sol_stats, title="",
            m=32, sig_type="SV32", min_stability=1.0, k_max=15, ax=None,
            legend=True, xlabel=True, label=None):

    # Plot results
    if ax is None:
        fig,ax=plt.subplots(1,1,figsize=(8,5))
    else:
        plt.sca(ax)

    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ln1=plt.plot(np.arange(1,k_max+1), cohort_stats["avg_stability"],
                "s-.", c=default_colors[0], label="Avg stability")
    plt.plot([0,k_max+1], [0.8,0.8], '--k', alpha=0.5)
    plt.ylabel("Stability / cosine similarity")
    if xlabel: plt.xlabel("Total signatures")

    ln_cos=plt.plot(np.arange(1,k_max+1), cohort_stats["Mean Cosine Distance"],
                    "o:", c=default_colors[3], label="Mean cosine distance")

    best = int(regexSearch("([0-9]+)\*",sol_stats.Signatures[sol_stats.Signatures.map(lambda x: "*" in x)].iloc[0], 1))
    ylim = ax.get_ylim()
    plt.fill_betweenx(ax.get_ylim(), [best-0.25,]*2, [best+0.25,]*2,
                    color=default_colors[3], alpha=0.25)

    ax.tick_params(axis='y', labelcolor=default_colors[0])
    plt.ylim(ylim)
    if label is not None: plt.text(k_max/2, 0.99, label, ha='center', va='top', fontsize=20)

    # Generate a new Axes instance, on the twin-X axes (same position)
    ax2 = ax.twinx()

    ln2=plt.plot(np.arange(1,k_max+1), cohort_stats['aic'],
                "*--", label="AIC", c=default_colors[9])
    ylim=ax2.get_ylim()
    plt.fill_betweenx(ylim, [np.argmin(cohort_stats['aic'])+0.7,]*2,
                                    [np.argmin(cohort_stats['aic'])+1.3,]*2,
                    color=default_colors[9], alpha=0.25)
    plt.ylabel("AIC")
    plt.ylim(ylim)

    ax2.tick_params(axis='y', labelcolor=default_colors[9])

    # added these three lines
    lns = ln1+ln2+ln_cos
    labs = [l.get_label() for l in lns]
    if legend:
        ax.legend(lns, labs, loc='lower left', bbox_to_anchor=(0.,1.), fontsize=12)
    plt.xlim(0.5,k_max+0.5)
    plt.title(title, fontsize=16)



if __name__=='__main__':

    sv_dir = f"{DATA_DIR}/Results_SigProfiler_v2"

    all_cohort_stats = pd.read_csv(f"{sig_dirs['SV32']}/all_cohort_stats_aic.tsv",
                            sep="\t", usecols=['signatures', 'cohort', 'avg_stability', 'Mean Cosine Distance', 'aic'])

    k_max=15
    cohorts = np.unique(all_cohort_stats.cohort)
    for cohort in tqdm.tqdm(cohorts, total=len(cohorts)):
        cohort_stats = all_cohort_stats[all_cohort_stats.cohort==cohort].sort_values("signatures")

        # Get solution statistics from SigProfiler
        sol_stats = pd.read_csv(f"{sv_dir}/{cohort.replace('_a','A').replace('_n','N').replace('_g','G')}_S1_S15_R500_PoisNoisRes_T/SV32/All_solutions_stat.csv")

        plotAIC(cohort_stats, sol_stats, m=32, sig_type="SV32", min_stability=1.0, k_max=k_max, ax=None,
                    title=cohort.replace("Connective", "Sarcoma"))

        # Save figure
        publish_fig(f"{cohort.replace('Connective', 'Sarcoma')}_SV32_aic", publish=FIGURE_DIR)


    cohort_set = ['Breast', 'Ovary', 'Uterus']
    fig, axes = plt.subplots(len(cohort_set), 1, figsize=(8,5*len(cohort_set)), sharex=True)
    plt.subplots_adjust(hspace=0.05)
    for i,cohort in enumerate(cohort_set):
        ax = axes[i]; plt.sca(ax)

        cohort_stats = all_cohort_stats[all_cohort_stats.cohort==cohort].sort_values("signatures")
        # Get solution statistics from SigProfiler
        sol_stats = pd.read_csv(f"{sv_dir}/{cohort.replace('_a','A').replace('_n','N').replace('_g','G')}_S1_S15_R500_PoisNoisRes_T/SV32/All_solutions_stat.csv")

        plotAIC(cohort_stats, sol_stats, m=32, sig_type="SV32", min_stability=1.0, k_max=k_max, ax=ax,
                    legend=i==0, xlabel=i==len(cohort_set)-1, label=cohort)

    # Save figure
    publish_fig(f"{'-'.join(cohort_set)}_SV32_aic", publish=FIGURE_DIR)