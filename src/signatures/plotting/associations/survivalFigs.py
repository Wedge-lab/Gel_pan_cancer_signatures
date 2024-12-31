import os
import scipy, scipy.stats
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from lifelines.statistics import logrank_test
from signatures.plotting.combinedSignatures import publish_fig, loadSignatures, signatureRenamer, map_colors
from signatures.utils import BH_threshold, orderSignatures
from dotenv import load_dotenv

load_dotenv()
RESULT_DIR = os.getenv("RESULT_DIR")
FIGURE_DIR = os.getenv("FIGURE_DIR")

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)

# ## Plot scatter survival figures
def plotScatter(results_df, pv_sws=0.05, sigs_=None, tissues_=None, pv_threshold=0.05,
                 size_factor=18, range_multiplier=3, fs=20, test_threshold=0.05, ax=None,
                 rename_sigs=None):

    results_pivoted = {}
    table_keys = ['wald_pvalue', 'target_mean_alt', 'target_phtest_alt', 'target_ljungtest_alt']
    for key in table_keys:
        results_pivoted[key] = pd.pivot_table(results_df, values=key,
                                              index=['group'], columns=['target'])
    results_pivoted['wald_pvalue'][results_pivoted['target_phtest_alt']<test_threshold] = np.nan
    #results_pivoted['wald_pvalue'][results_pivoted['target_ljungtest_alt']<test_threshold] = np.nan

    tissue_subset = np.nansum(np.array(results_pivoted['wald_pvalue'])<pv_threshold, axis=1)>0
    if tissues_ is not None:
        tissue_subset = tissue_subset&np.array([tissue in tissues_ \
                                                for tissue in results_pivoted['wald_pvalue'].index])
    tissues_ = results_pivoted['wald_pvalue'].index[tissue_subset][::-1]
    for tissue in tissues_:
        if not tissue in results_pivoted['wald_pvalue'].index:
            for key in table_keys:
                results_pivoted[key] = results_pivoted[key].T
                results_pivoted[key][tissue] = np.nan
                results_pivoted[key] = results_pivoted[key].T
    if sigs_ is None:
        sigs_subset = np.nansum(np.array(results_pivoted['wald_pvalue'])<pv_threshold, axis=0)>0
        sigs_ = results_pivoted['wald_pvalue'].columns[sigs_subset][
            orderSignatures(results_pivoted['wald_pvalue'].columns[sigs_subset])]
    else:
        for sig in sigs_:
            if not sig in results_pivoted['wald_pvalue'].columns:
                for key in table_keys:
                    results_pivoted[key][sig] = np.nan


    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=(20,20*len(tissues_)/len(sigs_)))#, sharex=True, sharey=True)
    plt.sca(ax)

    pv_table = np.array(results_pivoted['wald_pvalue'][sigs_].loc[tissues_])
    mu_table = np.array(results_pivoted['target_mean_alt'][sigs_].loc[tissues_])

    # Interquartile range
    dynamic_rng = np.diff(np.nanpercentile(np.array(mu_table).flatten(),
                                           np.array([25,75])))[0]

    # im = plt.pcolor(z_table, cmap=cm.coolwarm, vmin=-dynamic_rng*3, vmax=dynamic_rng*3)
    tried_grid = 1-np.isnan(pv_table).astype(int)
    tried_grid[np.abs(np.array(pv_table))<pv_sws] = 5
    plt.pcolor(tried_grid, cmap='bone_r', vmin=0, vmax=10)

    XX,YY = np.meshgrid(np.arange(pv_table.shape[0]), np.arange(pv_table.shape[1]))
    im = plt.scatter(YY.flatten()+0.5, XX.flatten()+0.5,
                     c=mu_table.T.flatten(),
                     s=-np.log10(pv_table.T.flatten())*size_factor,
                     cmap=cm.bwr, vmin=-dynamic_rng*range_multiplier, vmax=dynamic_rng*range_multiplier)

    ax.set_yticks(np.arange(len(tissues_)), minor=True)
    ax.set_xticks(np.arange(len(sigs_)), minor=True)
    ax.grid(True, which='minor', axis='both', linestyle='-', color='k')

    ax.set_xticks(np.arange(len(sigs_))+0.5, minor=False)
    ax.set_xticklabels(sigs_, rotation=90, fontsize=fs, minor=False);
    ax.set_yticks(np.arange(len(tissues_))+0.5, minor=False)

    tissue_labels = []
    test_keys = [key for key in results_df.keys() if ("phtest" in key)|("ljungtest" in key)]
    for tissue in tissues_:
        if scipy.stats.binom.sf(np.sum(np.array(results_df[results_df.group==tissue][test_keys])<0.01),
                                np.sum(results_df.group==tissue)*len(test_keys), 0.01)<0.05/len(sigs_):
            tissue_labels += [fr"*{tissue}*",]
        else:
            tissue_labels += [tissue,]
    ax.set_yticklabels(tissue_labels, rotation=0, fontsize=fs, minor=False);

    if rename_sigs is not None:
        for isig,sig in enumerate(sigs_):
            for sig_type in rename_sigs['degasperi']:
                if sig in list(rename_sigs['degasperi'][sig_type].old):
                    ax.get_xticklabels()[isig].set_color(map_colors['degasperi'])
                elif sig in list(rename_sigs['novel'][sig_type].new):
                    ax.get_xticklabels()[isig].set_color(map_colors['novel'])

    #plt.title(cov_labels[cov], fontsize=fs*2.5)
    plt.xlim(0,len(sigs_))
    plt.ylim(0,len(tissues_))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.05)
    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=fs*1.3)
    cbar.set_label(rf"$\beta$", fontsize=fs*1.3, rotation=0)

    legend_marks = np.array([0.05, 1e-3, 1e-5])
    legend_marks[np.argmin(np.abs(np.log10(legend_marks/pv_sws)))] = pv_sws
    legend_marks = np.sort(legend_marks)[::-1]
    legend_elements = [plt.scatter([],[],s=-size_factor*np.log10(mark), c='k', alpha=0.7) for mark in legend_marks]
    ax.legend(legend_elements,
              [fr'$P = {mark/10**(np.floor(np.log10(mark))):.1f}\times10^{{{int(np.floor(np.log10(mark)))}}}$'
                for mark in legend_marks],
              loc='lower right', bbox_to_anchor=(1,0.98),
              fontsize=fs*1.3, ncol=3, frameon=False)
    return sigs_, tissues_


# Make Kaplan-Meier lines
# Greenwood log-log confidence intervals
# https://www.math.wustl.edu/~sawyer/handouts/greenwood.pdf
def survival_line(deaths, times, alpha=0.05, bounds=True):

    za2 = scipy.stats.norm.isf(alpha/2)

    n_ = np.arange(len(times))[::-1]+1
    survival = np.cumprod(1-deaths/n_)

    if not bounds:
        return survival

    variance = 1/np.log(survival)**2 * np.cumsum(deaths/(n_*(n_-deaths)))
    cplus = np.log(-np.log(survival)) + za2*np.sqrt(variance)
    cminus = np.log(-np.log(survival)) - za2*np.sqrt(variance)

    lower = np.exp(-np.exp(cminus))
    upper = np.exp(-np.exp(cplus))

    return survival, lower, upper

def plot_KM(V, T, ax=None, label=None, c=None, alpha=0.05, count_times=[0]):
    order = np.argsort(T)
    survival, lower, upper = survival_line(V.copy()[order], T.copy()[order], alpha=alpha)

    y = np.insert(np.repeat(survival, 2), (0,0), (1,1))
    y_l = np.insert(np.repeat(lower, 2), (0,0), (np.nan,np.nan))
    y_u = np.insert(np.repeat(upper, 2), (0,0), (np.nan,np.nan))
    x = np.insert(np.repeat(T[order], 2), (0,2*len(T)), (0,np.inf))

    plt.plot(x, y, label=label, c=c)
    plt.fill_between(x, y_l, y_u, alpha=0.3, color=c)
    plt.scatter(T[order][V[order]==0], survival[V[order]==0], s=30, marker='|', c=c)

    if count_times is None:
        count_times = ax.get_xticks()
    counts = [len(T)-np.sum(T<time) for time in count_times]

    return count_times, counts

def plotSurvivalCurves(sample_df,
                        tgts = ["SBS17b","SBS2","CN17"],
                        groups = ["ColoRect-AdenoCA","Connective-Chondro", "Bladder-TCC"]):

    fig, axes = plt.subplots(1,len(groups), figsize=(20,5), sharey=True)
    plt.subplots_adjust(wspace=0.05)

    def_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fs=20
    for i, group in enumerate(groups):
        tgt = tgts[i]
        ax = axes[i]
        plt.sca(ax)

        # Merge targets with samples
        subset_df = sample_df[sample_df.group==group].copy()

        # Turn targets into binary parameters
        binary_tgt_arr = np.array(subset_df[targets])
        median = np.median(binary_tgt_arr, axis=0)
        binary_tgt_arr[binary_tgt_arr<=median] = 0
        binary_tgt_arr[binary_tgt_arr>median] = 1
        subset_df[targets] = binary_tgt_arr.astype(int)

        count_times = None # np.arange(0,np.max(subset_df.survival_time),500)

        ### KM curve split by target
        count_times, counts = plot_KM(np.array(subset_df[subset_df[tgt]==0].vital_status),
                np.array(subset_df[subset_df[tgt]==0].survival_time),
                ax=axes[i], label=f"No {tgt}", c=def_colors[0], alpha=0.05,
                        count_times=count_times)
        count_times = count_times[(count_times>=0)&\
                                (count_times<=np.max(subset_df.survival_time))]

        for it, time in enumerate(count_times):
            plt.text(time, -0.25, counts[it], ha='center', va='center', fontsize=fs*0.8, color=def_colors[0])

        count_times, counts = plot_KM(np.array(subset_df[subset_df[tgt]==1].vital_status),
                np.array(subset_df[subset_df[tgt]==1].survival_time),
                ax=axes[i], label=f"{tgt}", c=def_colors[1], alpha=0.05,
                        count_times=count_times)
        for it, time in enumerate(count_times):
            plt.text(time, -0.33, counts[it], ha='center', va='center', fontsize=fs*0.8, color=def_colors[1])

        results=logrank_test(np.array(subset_df[subset_df[tgt]==0].survival_time),
                                np.array(subset_df[subset_df[tgt]==1].survival_time),
                                event_observed_A=np.array(subset_df[subset_df[tgt]==0].vital_status),
                                event_observed_B=np.array(subset_df[subset_df[tgt]==1].vital_status))

        plt.text(np.max(subset_df.survival_time)*0.99, 0.1, fr"log-rank $P$ = {results.p_value:.2e}",
                ha='right', va='bottom', fontsize=fs*0.8)
        row = results_df[(results_df.group==group)&(results_df.target==tgt)].iloc[0]
        plt.text(np.max(subset_df.survival_time)*0.99, 0.02,
                fr"CPH $\beta = {row.target_mean_alt:.2f}_{{\,{row.target_95L_alt:.2f}}}^{{\,{row.target_95U_alt:.2f}}}$",
                ha='right', va='bottom', fontsize=fs*0.8)

        plt.title(group, fontsize=fs)

        if i==0:
            plt.text(-200, -0.25, "No sig:", color=def_colors[0], ha='right', va='center', fontsize=fs*0.8)
            plt.text(-200, -0.33, "Sig:", color=def_colors[1], ha='right', va='center', fontsize=fs*0.8)

        plt.xlabel("Time since sampling (days)", fontsize=fs)
        plt.ylim(0,1)
        plt.xlim(0,np.max(subset_df.survival_time))
        ax.set_xticks(count_times)

        plt.legend(loc='lower left')


def survivalFigure(sample_df, results_df, targets,
                    tgts = ["SBS17b","SBS2","CN17"],
                    groups = ["ColoRect-AdenoCA","Connective-Chondro", "Bladder-TCC"],
                    rename_sigs=None):

    _, axes = plt.subplots(2,1,figsize=(20,18), gridspec_kw={'height_ratios':[3,3]})
    plt.subplots_adjust(hspace=0.4)

    plt.text(0.07, 0.87, r"a)", fontsize=24, weight='bold', transform=plt.gcf().transFigure)
    plt.text(0.07, 0.44, r"b)", fontsize=24, weight='bold', transform=plt.gcf().transFigure)

    plt.sca(axes[0])

    divider = make_axes_locatable(axes[0])

    def_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fs=20
    for i, group in enumerate(groups):
        tgt = tgts[i]

        if i==0: ax = axes[0]
        else: ax = divider.append_axes("right", size=f"100%", pad=0.1)
        plt.sca(ax)

        # Merge targets with samples
        subset_df = sample_df[sample_df.group==group].copy()

        # Turn targets into binary parameters
        binary_tgt_arr = np.array(subset_df[targets])
        median = np.median(binary_tgt_arr, axis=0)
        binary_tgt_arr[binary_tgt_arr<=median] = 0
        binary_tgt_arr[binary_tgt_arr>median] = 1
        subset_df[targets] = binary_tgt_arr.astype(int)

        ### KM curve split by target
        count_times, counts = plot_KM(np.array(subset_df[subset_df[tgt]==0].vital_status),
                np.array(subset_df[subset_df[tgt]==0].survival_time),
                ax=ax, label=f"No {tgt}", c=def_colors[0], alpha=0.05,
                        count_times=None)
        count_times = count_times[(count_times>=0)&\
                                (count_times<=np.max(subset_df.survival_time))]

        for it, time in enumerate(count_times):
            plt.text(time, -0.20, counts[it], ha='center', va='center', fontsize=fs*0.8, color=def_colors[0])

        count_times, counts = plot_KM(np.array(subset_df[subset_df[tgt]==1].vital_status),
                np.array(subset_df[subset_df[tgt]==1].survival_time),
                ax=axes[i], label=f"{tgt}", c=def_colors[1], alpha=0.05,
                        count_times=count_times)
        for it, time in enumerate(count_times):
            plt.text(time, -0.25, counts[it], ha='center', va='center', fontsize=fs*0.8, color=def_colors[1])

        results=logrank_test(np.array(subset_df[subset_df[tgt]==0].survival_time),
                                np.array(subset_df[subset_df[tgt]==1].survival_time),
                                event_observed_A=np.array(subset_df[subset_df[tgt]==0].vital_status),
                                event_observed_B=np.array(subset_df[subset_df[tgt]==1].vital_status))

        plt.text(np.max(subset_df.survival_time)*0.99, 0.1, fr"log-rank $P$ = {results.p_value:.1e}",
                ha='right', va='bottom', fontsize=fs*0.8)
        row = results_df[(results_df.group==group)&(results_df.target==tgt)].iloc[0]
        plt.text(np.max(subset_df.survival_time)*0.99, 0.02,
                fr"CPH $\beta = {row.target_mean_alt:.2f}_{{\,{row.target_95L_alt:.2f}}}^{{\,{row.target_95U_alt:.2f}}}$",
                ha='right', va='bottom', fontsize=fs*0.8)

        plt.title(group, fontsize=fs)

        if i==0:
            plt.text(-200, -0.20, "No sig:", color=def_colors[0], ha='right', va='center', fontsize=fs*0.8)
            plt.text(-200, -0.25, "Sig:", color=def_colors[1], ha='right', va='center', fontsize=fs*0.8)

        plt.xlabel("Time since sampling (days)", fontsize=fs)
        if i>0:
            ax.set_yticks([])
        plt.ylim(0,1)
        plt.xlim(0,np.max(subset_df.survival_time))
        ax.set_xticks(count_times)

        plt.legend(loc='lower left')

    # Association scatter
    subset = (results_df['wald_pvalue']!=0)&(results_df['target_phtest_alt']>0.05)&(results_df['target_ljungtest_alt']>0.05)
    results_df = results_df[subset]
    pv_sws = BH_threshold(results_df.wald_pvalue, alpha=0.01)
    kwargs = {}
    _, _ = plotScatter(results_df, pv_sws=pv_sws, pv_threshold=0.1,
                                size_factor=20, range_multiplier=3,
                                ax=axes[1],
                                fs=18, test_threshold=0.01, rename_sigs=rename_sigs,
                                **kwargs)
    axes[1].set_title("CPH incl. grade", loc='left', fontsize=fs*1.5)


if __name__=='__main__':

    # Get renaming dictionary
    combined_acts, combined_sigs, relabel_map = loadSignatures()
    rename_dict = signatureRenamer(relabel_map)

    run_name="survival_log_grade-4_unique"

    # Import sample data
    sample_df = pd.read_csv(f"{RESULT_DIR}/survival/{run_name}/input/samples.tsv",
                            sep="\t", index_col=0)
    # Import tgtnature data
    target_df = pd.read_csv(f"{RESULT_DIR}/survival/{run_name}/input/signatures.tsv",
                            sep="\t", index_col=0)
    target_df.rename(rename_dict, axis=1, inplace=True)

    # Merge targets with samples
    targets = target_df.keys()
    sample_df = pd.merge(sample_df, target_df, how='inner', left_index=True, right_on='sample_id')

    # Load association results
    results_dir = f"{RESULT_DIR}/survival/{run_name}/output"
    results_df = pd.read_csv(f"{results_dir}/survival_results.tsv", sep="\t")

    results_df.replace({"target":rename_dict}, inplace=True)

    # Replace Connective with Sarcoma
    sample_df.group = sample_df.group.str.replace("Connective","Sarcoma")
    results_df.group = results_df.group.str.replace("Connective","Sarcoma")

    # Estimate pvalue
    results_df['wilks_pvalue'] = scipy.stats.chi2.sf(2*(results_df['CPHlikelihood_alt']-results_df['CPHlikelihood_null']), df=1)
    results_df['wald_pvalue'] = scipy.stats.chi2.sf((results_df.target_mean_alt/results_df.target_err_alt)**2, df=1)


    survivalFigure(sample_df, results_df, targets,
                            tgts=["SBS17b","CN17"],
                            groups=["ColoRect-AdenoCA", "Bladder-TCC"],
                            rename_sigs=relabel_map)

    publish_fig(f"survival_panel_{run_name}", publish=FIGURE_DIR)
