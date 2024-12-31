#!/usr/bin/env python
# coding: utf-8

import os
import scipy, scipy.stats
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from signatures.utils import BH_threshold
from signatures.plotting.combinedSignatures import loadSignatures, signatureRenamer, map_colors
from signatures.utils import orderSignatures
from dotenv import load_dotenv

load_dotenv()
FIGURE_DIR = os.getenv("FIGURE_DIR")
RESULT_DIR = os.getenv("RESULT_DIR")

default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)

def plot_scatter(results_df, pv_sws=0.05,
                 sigs_=None, tissues_=None,
                 pvalue='wald_pvalue', effect_size='target_mean_alt',
                 vmax=None, size_factor=13, fs=20, figsize=20, ax=None,
                 rename_sigs=None):

    results_pivoted = {}
    for key in [pvalue, effect_size]:
        results_pivoted[key] = pd.pivot_table(results_df, values=key,
                                              index=['group'], columns=['target'])

    if tissues_ is None:
        tissue_subset = np.nansum(np.array(results_pivoted[pvalue])<0.05, axis=1)>0
        tissues_ = results_pivoted[pvalue].index[tissue_subset][::-1]
    if sigs_ is None:
        sigs_subset = np.nansum(np.array(results_pivoted[pvalue])<0.05, axis=0)>0
        sigs_ = results_pivoted[pvalue].columns[sigs_subset][
            orderSignatures(results_pivoted[pvalue].columns[sigs_subset])]
    else:
        for sig in sigs_:
            if not sig in results_pivoted[pvalue].columns:
                results_pivoted[pvalue][sig] = np.nan
                results_pivoted[effect_size][sig] = np.nan

    if ax is None:
        fig, axes = plt.subplots(1,1,figsize=(figsize,figsize*len(tissues_)/len(sigs_)))#, sharex=True, sharey=True)
        ax = axes
    plt.sca(ax)

    pv_table = np.array(results_pivoted[pvalue][sigs_].loc[tissues_])
    mu_table = np.array(results_pivoted[effect_size][sigs_].loc[tissues_])

    # Interquartile range
    if vmax is None:
        dynamic_rng = np.diff(np.nanpercentile(np.array(mu_table).flatten(),
                                               np.array([25,75])))[0]
    else: dynamic_rng=vmax/3

    # im = plt.pcolor(z_table, cmap=cm.coolwarm, vmin=-dynamic_rng*3, vmax=dynamic_rng*3)
    tried_grid = 1-np.isnan(pv_table).astype(int)
    tried_grid[np.abs(np.array(pv_table))<pv_sws] = 5
    plt.pcolor(tried_grid, cmap='bone_r', vmin=0, vmax=10)

    XX,YY = np.meshgrid(np.arange(pv_table.shape[0]), np.arange(pv_table.shape[1]))
    im = plt.scatter(YY.flatten()+0.5, XX.flatten()+0.5,
                     c=mu_table.T.flatten(),
                     s=np.min(np.vstack((-np.log10(pv_table.T.flatten())*size_factor,
                                          np.ones(len(pv_table.flatten()))*150)), axis=0),
                     cmap=cm.bwr, vmin=-dynamic_rng*3, vmax=dynamic_rng*3)

    ax.set_yticks(np.arange(len(tissues_)), minor=True)
    ax.set_xticks(np.arange(len(sigs_)), minor=True)
    ax.grid(True, which='minor', axis='both', linestyle='-', color='k')

    ax.set_xticks(np.arange(len(sigs_))+0.5, minor=False)
    ax.set_xticklabels(sigs_, rotation=90, fontsize=fs, minor=False);
    ax.set_yticks(np.arange(len(tissues_))+0.5, minor=False)
    ax.set_yticklabels(tissues_, rotation=0, fontsize=fs, minor=False);

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

    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Scatter',
                              markerfacecolor='k', alpha=0.5, markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Scatter',
                              markerfacecolor='k', alpha=0.5, markersize=20),
                       Line2D([0], [0], marker='o', color='w', label='Scatter',
                              markerfacecolor='k', alpha=0.5, markersize=30)]
    sizes = np.array([0.05,0.00001,pv_sws])
    l1 = plt.scatter([],[],s=-size_factor*np.log10(sizes[0]), c='k', alpha=0.7)
    l2 = plt.scatter([],[],s=-size_factor*np.log10(sizes[1]), c='k', alpha=0.7)
    l3 = plt.scatter([],[],s=-size_factor*np.log10(sizes[2]), c='k', alpha=0.7)
    order = np.argsort(sizes)[::-1]
    ax.legend([[l1,l2,l3][idx] for idx in order],
              np.array([
               f'$P=0.05$',
               fr'$P=10^{{-5}}$',
               fr'$P={pv_sws/10**(np.floor(np.log10(pv_sws))):.1f}\times10^{{{int(np.floor(np.log10(pv_sws)))}}}$'
              ])[order],
              loc='lower right', bbox_to_anchor=(1,1),
              fontsize=fs*1.3, ncol=3, frameon=False)
    return sigs_, tissues_


def publish_fig(filename, publish="./"):
    plt.savefig(f'{publish}/{filename}.png', bbox_inches='tight', dpi=300, facecolor='w', transparent=False)
    plt.savefig(f'{publish}/{filename}.pdf', bbox_inches='tight', format='pdf', facecolor='w', transparent=False)

if __name__=='__main__':

    # ## Load results
    run_name = "clinical_tumour-group_InDelQC_nb"

    results_dir = f"{RESULT_DIR}/results/associations/clinicalSigs/{run_name}/output"
    results_zinb = pd.read_csv(f"{results_dir}/signature_target_assoc_nb.csv")
    results_log0 = pd.read_csv(f"{results_dir}/signature_target_assoc_logistic.csv")
    results_target = pd.read_csv(f"{results_dir}/target_target_assoc.csv")
    fig_dir = f"{results_dir}/figures"

    input_dir = f"{RESULT_DIR}/results/clinicalSigs/{run_name}/input/"
    samples_file = f"{input_dir}/samples.tsv"
    samples_cohort_file = f"{input_dir}/samples_cohort.tsv"
    signatures_file = f"{input_dir}/signatures.tsv"
    targets_file = f"{input_dir}/targets.tsv"
    targets_cohort_file = f"{input_dir}/targets_cohort.tsv"

    results_df = pd.merge(results_zinb, results_log0,
                                on=('target', 'signature', 'group'),
                                suffixes=("_zinb", '_log0'), how='outer')

    subset = (results_df.resample_pvalue>0)
    results_df = results_df[subset]

    results_df['zscore'] = results_df.target_means_alt_zinb/np.sqrt(results_df.target_covs_alt_zinb)
    results_df['log10_pvalue_zinb'] = -np.log10(results_df.resample_pvalue)
    results_df['log10_pvalue_log0'] = -np.log10(results_df.wilks_pvalue_log0)

    # Get renaming dictionary
    combined_acts, combined_sigs, relabel_map = loadSignatures()
    rename_dict = signatureRenamer(relabel_map)
    results_df.replace({'signature':rename_dict}, inplace=True)

    # Rename connective as sarcoma
    results_df.group = results_df.group.str.replace("Connective","Sarcoma")

    pvalue = "resample_pvalue"#"wilks_pvalue_log0"#wilks_pvalue_zinb'
    effect_size = 'target_means_alt_zinb'#"target_means_alt_log0"#

    transform = "_log"
    fig, axes = plt.subplots(3,1,figsize=(25,25), gridspec_kw={'height_ratios':[4,3,8]})
    plt.subplots_adjust(hspace=0.4)
    fs=20

    plt.text(0.07, 0.89, r"a)", fontsize=24, weight='bold', transform=plt.gcf().transFigure)
    plt.text(0.07, 0.65, r"b)", fontsize=24, weight='bold', transform=plt.gcf().transFigure)
    plt.text(0.07, 0.45, r"c)", fontsize=24, weight='bold', transform=plt.gcf().transFigure)
    plt.text(0.6, 0.45, r"d)", fontsize=24, weight='bold', transform=plt.gcf().transFigure)

    results_df["target_icovs_pvalue"] = 1/(results_df[effect_size]**2/\
                                        scipy.stats.chi2.isf(results_df[pvalue], df=1))
    # results_df["target_icovs_pvalue"] = 1/results_df["target_covs_alt_zinb"]
    results_df["target_muicovs_pvalue"] = results_df[effect_size]*results_df["target_icovs_pvalue"]
    signature_sets = {"HRD":["SBS3","DBS6","ID6","CN17","SV3"],
                    "APOBEC":["SBS2","SBS13"],
                    "MMR":["SBS1", "SBS5", "ID2","SBS15","SBS26","SBS44"],
                    r"$POLE$":["SBS10a", "SBS10b", "SBS28", "DBS3","DBS10"],
                    "Smoking":["SBS4","SBS92","DBS2", "ID3"]}
    signature_meanings = {sig:key for key in signature_sets.keys() for sig in signature_sets[key]}

    results_sigtypes = pd.DataFrame(results_df.replace({"signature":signature_meanings})[
        ["group","target","signature","target_muicovs_pvalue","target_icovs_pvalue"]
    ].groupby(["group","target","signature"]).sum()).reset_index()
    results_sigtypes["mean"] = results_sigtypes["target_muicovs_pvalue"]/results_sigtypes["target_icovs_pvalue"]
    results_sigtypes["std"] = np.sqrt(1/results_sigtypes["target_icovs_pvalue"])


    all_signatures = np.hstack([combined_acts[sig_set].keys()[
                np.sum(pd.merge(pd.read_csv(samples_file, sep="\t", usecols=['sample_id', 'group']),
                                combined_acts[sig_set]>0,
                                left_on='sample_id', right_index=True, how='inner')\
                    .drop('sample_id', axis=1).groupby('group').sum()>20)>0] for sig_set in combined_acts])
    all_signatures = all_signatures[orderSignatures(all_signatures)]


    plt.sca(axes[0])
    plt.title("Grade", fontsize=fs*1.2, loc='left')

    pv_threshold=0.05
    target = 'grade'

    tissues_=None; sigs_=None
    target_df = results_df[results_df.target==target].copy()

    p_threshold = BH_threshold(target_df[pvalue], alpha=0.01)

    _ = plot_scatter(target_df.drop('target', axis=1).rename({"signature":"target"}, axis=1),
                pvalue=pvalue, effect_size=effect_size,
                pv_sws=p_threshold, size_factor=20, fs=16, figsize=23,
                    ax=axes[0], rename_sigs=relabel_map,
                    sigs_=all_signatures)
    sigs1 = _[0]

    divider = make_axes_locatable(axes[2])
    groups = [group for group in  np.unique(results_sigtypes["group"])[::-1] if "-" in group]

    target_names = {"stage_best":"Stage", "grade":"Grade"}

    for it, target in enumerate(["stage_best", "grade"]):
        if it==0: ax = axes[2]
        else: ax = divider.append_axes("right", size=f"100%", pad=0.1)
        plt.sca(ax)

        for isig,sig in enumerate(signature_sets.keys()):
            subset = ((results_sigtypes.signature==sig)&(results_sigtypes.target==target))
            # ypos = np.intersect1d(np.array(results_sigtypes[subset]["group"]), groups, return_indices=True)[2][::-1]
            inter = np.intersect1d(groups, np.array(results_sigtypes[subset]["group"]), return_indices=True)
            ypos = inter[1][np.argsort(inter[2])]
            plt.errorbar(results_sigtypes[subset]["mean"],
                        ypos+0.3-isig*0.15,
                        yerr=None, xerr=results_sigtypes[subset]["std"]*2,
                        fmt="o", label=sig)

        ax.set_yticks(np.arange(len(groups))-0.5, minor=True)
        ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.5)
        ax.set_yticks(np.arange(len(groups)), minor=False)
        ax.set_yticklabels(groups, rotation=0, fontsize=fs, minor=False);

        plt.plot([0,0],[-0.5,len(groups)+0.5], '--k')
        plt.ylim(-0.5,len(groups)-0.5)
        if it==1: plt.legend(loc='upper left', frameon=True, bbox_to_anchor=(-0.25,1.1), ncol=1)

        plt.title(target_names[target], fontsize=fs*1.2)

        plt.xlim(-4,4)

        plt.xlabel(r"$\beta$", fontsize=fs)
        if it==1:
            ax.set_yticks([])


    ax = axes[1]
    ax.set_title("Tumour characteristics", fontsize=fs*1.2, loc='left')
    tumour_specific_targets = {'er_status':"ER", 'pr_status':"PR", 'her2_status':"HER2", 'triple_negative':"TNBC",
                            'breslow':"Breslow", 'dukes':"Dukes", 'figo':"FIGO",
                            'gleason_combined':"Gleason", 'npi':"NPI"}
    target_df = results_df[results_df.target.map(lambda x: x in tumour_specific_targets.keys())]

    target_df['group_target'] = target_df.group+" "+target_df.target.replace(tumour_specific_targets)

    p_threshold = BH_threshold(results_df[pvalue], alpha=0.01)

    # Reorder groups
    groups = np.unique(target_df.group_target)
    inverse_idx = np.unique(np.array([group.split(" ")[1] for group in groups]), return_inverse=True)[1]
    intersection = np.intersect1d(np.array([group.split(" ")[1] for group in np.unique(groups)]),
                np.array(list(tumour_specific_targets.values())), return_indices=True)
    groups = groups[np.argsort(intersection[2][inverse_idx])]
    groups = groups[np.argsort([group.split(" ")[0] for group in groups])]

    _ = plot_scatter(target_df.drop(['target', 'group'], axis=1)\
                        .rename({"signature":"target", "group_target":"group"}, axis=1),
                effect_size=effect_size,
                pvalue=pvalue,
                tissues_=groups[::-1], pv_sws=p_threshold, vmax=2, size_factor=20, figsize=24, fs=16,
                ax=axes[1], rename_sigs=relabel_map,
                sigs_=all_signatures)
    sigs2 = _[0]


    target_df = results_df[results_df.target.map(lambda x: x in tumour_specific_targets.keys())]
    target_df['group_target'] = target_df.group+" "+target_df.target.replace(tumour_specific_targets)

    target_sigtypes = pd.DataFrame(target_df.replace({"signature":signature_meanings})[
        ["group_target","signature","target_muicovs_pvalue","target_icovs_pvalue"]
    ].groupby(["group_target","signature"]).sum()).reset_index()

    target_sigtypes["mean"] = target_sigtypes["target_muicovs_pvalue"]/target_sigtypes["target_icovs_pvalue"]
    target_sigtypes["std"] = np.sqrt(1/target_sigtypes["target_icovs_pvalue"])

    ax = divider.append_axes("right", size=f"100%", pad=3.5)
    plt.sca(ax)
    plt.title("Tumour characteristics", fontsize=fs*1.2)

    for isig,sig in enumerate(list(signature_sets.keys())):
        subset = (target_sigtypes.signature==sig)
        inter = np.intersect1d(groups[::-1], np.array(target_sigtypes[subset]["group_target"]), return_indices=True)
        ypos = inter[1][np.argsort(inter[2])]
        plt.errorbar(target_sigtypes[subset]["mean"],
                    ypos+0.3-isig*0.15,
                    yerr=None, xerr=target_sigtypes[subset]["std"]*2,
                    fmt="o", label=sig)

    ax.set_yticks(np.arange(len(groups))-0.5, minor=True)
    ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.5)
    ax.set_yticks(np.arange(len(groups)), minor=False)
    ax.set_yticklabels(groups[::-1], rotation=0, fontsize=fs, minor=False);

    plt.plot([0,0],[-0.5,len(groups)+0.5], '--k')
    plt.ylim(-0.5,len(groups)-0.5)
    # plt.legend(loc='upper right')

    plt.xlim(-2.5,2.5)

    plt.xlabel(r"$\beta$", fontsize=fs)
    publish_fig("clinical_panel", publish=FIGURE_DIR)


    # ## Histology

    run_name = "clinical_tumour-group_InDelQC_nb"

    results_dir = f"{RESULT_DIR}/results/clinicalSigs/{run_name}/output"
    results_zinb = pd.read_csv(f"{results_dir}/signature_target_assoc_nb.csv")
    results_log0 = pd.read_csv(f"{results_dir}/signature_target_assoc_logistic.csv")
    results_target = pd.read_csv(f"{results_dir}/target_target_assoc.csv")
    fig_dir = f"{results_dir}/figures"

    results_dir_cohort = f"/re_gecip/shared_allGeCIPs/pancancer_signatures/results/associations/clinicalSigs/{run_name}/output_cohort"
    results_zinb = pd.concat((results_zinb, pd.read_csv(f"{results_dir_cohort}/signature_target_assoc_zinb.csv")))
    results_log0 = pd.concat((results_log0, pd.read_csv(f"{results_dir_cohort}/signature_target_assoc_logistic.csv")))
    results_target = pd.concat((results_target, pd.read_csv(f"{results_dir_cohort}/target_target_assoc.csv")))

    results_df = pd.merge(results_zinb, results_log0,
                                on=('target', 'signature', 'group'),
                                suffixes=("_zinb", '_log0'))#.replace({'group':group_labels})

    # Select subset of runs
    subset = (results_df.resample_pvalue>0)&\
            (results_df.model_alt_zinb=='negbin')&\
            (np.abs(results_df.target_means_alt_zinb)<1000)&\
            (results_df.resample_pvalue<0.01)#results_df.signature.map(lambda x: x[0]!='T')&\

    results_df['zscore'] = results_df.target_means_alt_zinb/np.sqrt(results_df.target_covs_alt_zinb)
    results_df['log10_pvalue_zinb'] = -np.log10(results_df.resample_pvalue)
    results_df['log10_pvalue_log0'] = -np.log10(results_df.wilks_pvalue_log0)

    # Rename connective as sarcoma
    results_df.group = results_df.group.str.replace("Connective","Sarcoma")
    results_df.target = results_df.target.str.replace("Connective","Sarcoma")
    results_df.replace({"signature":rename_dict}, inplace=True)
    results_df.target = results_df.target.map(lambda x: x.replace(".", "-"))

    # Select tumour groups
    histology_tests = ['Breast-DuctalCA',
                        'CNS-Astro', 'CNS-GBM-IDHmut', 'CNS-GBM-IDHwt', 'CNS-Menin', 'CNS-Oligo',
                        'Sarcoma-Chondro', 'Sarcoma-Leiomyo', 'Sarcoma-Liposarc',
                        'Sarcoma-Myxofibro', 'Sarcoma-Osteosarc', 'Sarcoma-SCS', 'Sarcoma-SS',
                        'Haem-ALL', 'Haem-AML', 'Haem-CML', 'Haem-MM', 'Haem-MPN',
                        'Kidney-CCRCC', 'Kidney-ChRCC', 'Kidney-PRCC',
                        'Lung-AdenoCA', 'Lung-LargeCell', 'Lung-SCC', 'Lung-SmallCell',
                        'crc_colon_rectum']
    results_df = results_df[results_df.target.map(lambda x: x in histology_tests)]

    pvalue = 'wilks_pvalue_log0'
    significant_signatures = pd.DataFrame(results_df[results_df[pvalue]<0.05/len(results_df)]\
                .groupby(['group','signature']).size()).reset_index().set_index('group')

    groups = np.unique(results_df.group)
    ncol = 3
    fig, axes = plt.subplots(int(np.ceil(len(groups)/ncol)),ncol,figsize=(10*ncol,5*np.ceil(len(groups)/ncol)))
    plt.subplots_adjust(hspace=0.4, wspace=0.5)

    fs=14

    for igrp, group in enumerate(groups):
        print(group)
        ax = axes[igrp//ncol, igrp%ncol]
        plt.sca(ax)
        significant_results = results_df[(results_df[pvalue]<0.05)&(results_df.group==group)]

        target_signatures = np.array([sig for sig in significant_results.signature])# if sig[0]!='T'])
        target_signatures = target_signatures[orderSignatures(target_signatures)]
        significant_results = significant_results.set_index('signature').loc[target_signatures]
        significant_results = significant_results[significant_results.target_covs_alt_log0<100]

        unique_targets = np.unique(significant_results['target'])
        unique_signatures = np.unique(significant_results.reset_index()['signature'])
        unique_signatures = unique_signatures[orderSignatures(unique_signatures)]
        # print(unique_signatures)

        for itgt,target in enumerate(unique_targets):
            subset = np.intersect1d(unique_signatures,
                                    np.array(significant_results[significant_results.target==target].index),
                                    return_indices=True)
            plt.errorbar(subset[1]+0.2*itgt/len(unique_targets)-0.1,
                        significant_results[significant_results.target==target]['target_means_alt_log0'].iloc[subset[2]],
                        2*np.sqrt(significant_results[significant_results.target==target]['target_covs_alt_log0']).iloc[subset[2]], fmt='o',
                        label="-".join(target.split("-")[1:]) if target[:3]!='crc' else "Colon")
        plt.title(target.split("-")[0] if target[:3]!='crc' else "Colon vs Rectum",
                fontsize=fs*1.2)
        xlim = [-0.5,len(unique_signatures)-0.5]
        plt.plot(xlim,[0,0], '--k')
        plt.xlim(xlim)

        plt.legend(bbox_to_anchor=(1.,1.), loc='upper left', frameon=False)

        ax.set_xticks(np.arange(len(unique_signatures))-0.5, minor=True)
        ax.grid(True, which='minor', axis='x', linestyle='-', color='k', alpha=0.3)

        ax.set_xticks(np.arange(len(unique_signatures)), minor=False)
        ax.set_xticklabels(unique_signatures, rotation=90, fontsize=fs, minor=False);

        if relabel_map is not None:
            for isig,sig in enumerate(unique_signatures):
                for sig_type in relabel_map['degasperi']:
                    if sig in list(relabel_map['degasperi'][sig_type].old):
                        ax.get_xticklabels()[isig].set_color(map_colors['degasperi'])
                    elif sig in list(relabel_map['novel'][sig_type].new):
                        ax.get_xticklabels()[isig].set_color(map_colors['novel'])

        plt.ylabel(r"$\beta$", fontsize=fs*1.2)

    for igrp in range(igrp+1, int(ncol*np.ceil(len(groups)/ncol))):
        ax = axes[igrp//ncol, igrp%ncol]
        plt.sca(ax); plt.axis('off')

    publish_fig("histology_rate_comparison", publish=FIGURE_DIR)

