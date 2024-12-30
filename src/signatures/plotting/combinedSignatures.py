import sys, os
import pandas as pd, numpy as np, re, tqdm
import scipy, scipy.stats

import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
import matplotlib.pyplot as plt
from pylab import cm
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from fisher import pvalue_npy
import scipy.cluster.hierarchy as shc
from matplotlib.patches import Rectangle
from signatures.combineSignatures._03_combineSignatures import orderSignatures, cosineSimilarity
from dotenv import load_dotenv
from signatures.utils import regexSearch

load_dotenv()
COMBINED_SIGS_DIR = os.getenv('COMBINED_SIGS_DIR')
REF_DIR = os.getenv('REF_DIR')
FIGURE_DIR = os.getenv('FIGURE_DIR')
DATA_DIR = os.getenv('DATA_DIR')

plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)

default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Signature directories
sig_dirs  = { "SBS288":f"{COMBINED_SIGS_DIR}/combinedSignatures_SBS288",
               "DBS78": f"{COMBINED_SIGS_DIR}/combinedSignatures_DBS78_V4",
               "ID83": f"{COMBINED_SIGS_DIR}/combinedSignatures_ID83",
               "CNV48": f"{COMBINED_SIGS_DIR}/combinedSignatures_CNV48",
               "SV32":  f"{COMBINED_SIGS_DIR}/combinedSignatures_SV32"}


# COSMIC and other reference files
references = {"SBS288":f"{REF_DIR}/COSMIC_v3.3.1_SBS_GRCh38.txt",
                "DBS78":f"{REF_DIR}/COSMIC_v3.3_DBS_GRCh38.txt",
                "ID83":f"{REF_DIR}/COSMIC_v3.3_ID_GRCh37.txt",
                "CNV48":f"{REF_DIR}/COSMIC_v3.3_CN_GRCh37.txt",
                "SV32":f"{REF_DIR}/Breast560_rearrangement.signatures.tsv"}

map_colors = {'degasperi':'darkorange',
              'novel':'firebrick'}

def publish_fig(filename, publish="./"):
    plt.savefig(f'{publish}/{filename}.png',
                   bbox_inches='tight', dpi=300, facecolor='w', transparent=False)
    plt.savefig(f'{publish}/{filename}.pdf',
                   bbox_inches='tight', format='pdf', facecolor='w', transparent=False)

def loadSignatures(sv_rename_deg=False):

    """
    Load signatures and activities from SBS to SVs
    Relabel signatures according to COSMIC / SV reference
    """

    # Load signatures and sample activities
    combined_sigs = {}
    combined_acts = {}
    for sig_type in sig_dirs:
        combined_sigs[sig_type] = pd.read_csv(f"{sig_dirs[sig_type]}/Combined_Solution_Signatures.tsv", sep="\t", index_col=0)
        combined_acts[sig_type] = pd.read_csv(f"{sig_dirs[sig_type]}/Combined_Solution_Activities.tsv", sep="\t", index_col=0).fillna(0)

    # Get reference signatures
    reference_sigs = {sig_type: pd.read_csv(references[sig_type], sep="\t", index_col=0) \
                        for sig_type in combined_sigs}

    # Maximum reference signature number for each type (e.g. 24 for CNV)
    max_sigs = {sig_type:np.max(np.array([regexSearch("([A-Z]+)([0-9]+)([a-z]*)", sig, 2) \
                                for sig in reference_sigs[sig_type].keys()]).astype(int))
                    for sig_type in combined_sigs}

    # Map signatures to COSMIC / Nik-Zainal reference
    reference_map = {}
    for sig_type in reference_sigs:
        reference_map[sig_type] = iterativeMapping(combined_sigs[sig_type], reference_sigs[sig_type],
                                            mut_type=regexSearch("([A-Z]+)[0-9]+", sig_type, 1).replace("CNV", "CN")) # type: ignore

    # Change SV6 to SV6a as we will also have an SV6b from Degasperi
    if sv_rename_deg: reference_map['SV32']['new'] = reference_map['SV32'].new.replace({'SV6':'SV6a'})

    # Get reference Degasperi 2022/2020 signatures
    degasperi_sigs = {"SBS288": pd.read_excel('/re_gecip/cancer_pan/aeverall/data/signatures/science.abl9283_tables_s1_to_s33.xlsx',
                                            f'Table S21').set_index('mutationClass'),
                        "DBS78": pd.read_excel('/re_gecip/cancer_pan/aeverall/data/signatures/science.abl9283_tables_s1_to_s33.xlsx',
                                                f'Table S22').set_index('mutationClass'),
                        "SV32": pd.read_csv('/re_gecip/shared_allGeCIPs/pancancer_signatures/data/RefSigv1_Rearr.tsv',
                                                sep="\t")}
    for sig_type in degasperi_sigs:
        degasperi_sigs[sig_type] = degasperi_sigs[sig_type].loc[combined_sigs[sig_type].index]

    # Replace signature labels from COSMIC
    mapped_sigs = {}
    for sig_type in combined_sigs:
        combined_acts[sig_type].rename(reference_map[sig_type].set_index('old')['new'].to_dict(), axis=1, inplace=True)
        combined_sigs[sig_type].rename(reference_map[sig_type].set_index('old')['new'].to_dict(), axis=1, inplace=True)
        mapped_sigs[sig_type] = reference_map[sig_type].new

    # Map signatures to Degasperi 2020/2022
    if sv_rename_deg:
        degasperi_sv_map = iterativeMapping(combined_sigs['SV32'].drop(reference_map['SV32'].old, axis=1),
                                                degasperi_sigs['SV32'],
                                                mut_type=regexSearch("([A-Z]+)[0-9]+", 'SV32', 1))

        combined_acts['SV32'].rename(degasperi_sv_map.set_index('old')['new'].to_dict(), axis=1, inplace=True)
        combined_sigs['SV32'].rename(degasperi_sv_map.set_index('old')['new'].to_dict(), axis=1, inplace=True)
        mapped_sigs['SV32'] = np.hstack((mapped_sigs['SV32'], degasperi_sv_map.new))
    else:
        degasperi_sv_map = None

    novel_map = {}
    for sig_type in combined_sigs:
        mut_type = regexSearch('([A-Z]+)([0-9]+)', sig_type, 1).replace("CNV", "CN") # type: ignore
        # if mut_type=='CNV': mut_type='CN'
        # Count up from max_sigs
        i = max_sigs[sig_type]+1
        # Initialise novel map dictionary
        novel_map[sig_type] = {"X":"X"}


        # Get new signature labels based on COSMIC max label and current labels
        for sig in combined_sigs[sig_type].drop(mapped_sigs[sig_type], axis=1).keys():
            new_sig = f"{mut_type}{i}"
            while new_sig in mapped_sigs[sig_type]:
                i+=1
                new_sig = f"{mut_type}{i}"
            novel_map[sig_type][sig] = new_sig
            i+=1

        novel_map[sig_type] = pd.DataFrame.from_dict(novel_map[sig_type], orient='index').drop("X").reset_index()\
                                        .rename({'index':'old', 0:'new'}, axis=1)

    # Rename signatures using novel labels
    for sig_type in novel_map:
        combined_acts[sig_type].rename(novel_map[sig_type].set_index('old')['new'].to_dict(), axis=1, inplace=True)
        combined_sigs[sig_type].rename(novel_map[sig_type].set_index('old')['new'].to_dict(), axis=1, inplace=True)

    # Map signatures to Degasperi 2020/2022
    degasperi_map = {}
    for sig_type in degasperi_sigs:
        degasperi_map[sig_type] = iterativeMapping(combined_sigs[sig_type].drop(reference_map[sig_type].new, axis=1),
                                                degasperi_sigs[sig_type],
                                                mut_type=regexSearch("([A-Z]+)[0-9]+", sig_type, 1).replace("CNV", "CN"))
    if sv_rename_deg: degasperi_map['SV32'] = pd.concat((degasperi_map['SV32'], degasperi_sv_map))
    for sig_type in ['ID83', 'CNV48']: degasperi_map[sig_type] = pd.DataFrame.from_dict({"old":[],"new":[]})

    return combined_acts, combined_sigs, {'reference':reference_map, 'degasperi':degasperi_map, 'novel': novel_map}

def signatureRenamer(relabel_map, sv_rename_deg=False):

    ### Rename signatures
    rename_df = pd.DataFrame()
    for sig_type in ['SBS288', 'DBS78', 'ID83', 'CNV48', 'SV32']:
        rename_df = pd.concat((rename_df, relabel_map['reference'][sig_type]))
        rename_df = pd.concat((rename_df, relabel_map['novel'][sig_type]))
    if sv_rename_deg: rename_df = pd.concat((rename_df, relabel_map['degasperi']['SV32']))

    return dict(zip(rename_df.old, rename_df.new))

def iterativeMapping(signatures, reference, threshold=0.8, mut_type='SV'):

    """ Iterative mapping
    Get matching of novel signatures to reference signatures
    """

    # Get matching of novel signatures to reference signatures
    nov_sig_list = list(np.sort(signatures.keys()))
    ref_sig_list = list(np.sort(reference.keys()))
    cossim_mtx = cosineSimilarity(np.array(signatures[np.sort(signatures.keys())]),
                                              np.array(reference[np.sort(reference.keys())]))
    mapping={}
    while (len(nov_sig_list)>0)&(len(ref_sig_list)>0):
        best_match = np.unravel_index(np.argmax(cossim_mtx), (len(nov_sig_list), len(ref_sig_list)))

        if cossim_mtx[best_match]>threshold:
            mapping[nov_sig_list[best_match[0]]] = ref_sig_list[best_match[1]]

            nov_sig_list.remove(nov_sig_list[best_match[0]])
            ref_sig_list.remove(ref_sig_list[best_match[1]])

            cossim_mtx = np.delete(cossim_mtx, best_match[0], axis=0)
            cossim_mtx = np.delete(cossim_mtx, best_match[1], axis=1)
        else:
            break

    # Save index mapping
    # mapping = {key: val.replace("RefSigR","SV") for key,val in mapping.items()}
    mapping = {key: mut_type+regexSearch("[A-Za-z]+([0-9]+[a-z]*)",val, 1) for key,val in mapping.items()}
    return pd.DataFrame.from_dict({'old':mapping.keys(), 'new':mapping.values()})


def plotSignature(ax, signature, mutation_labels, mutation_type="SBS", xticks=True, text="", fs=12,
                 pad=0.1):

    """
    Generate bar plot of the signature
    - Shape of figure, x labels and bar colours depend on the mutation type
    """

    plt.sca(ax)

    if mutation_type=='SBS':
        context_1 = signature.index.map(lambda x: regexSearch("([A-Z])\[([A-Z>]+)\]([A-Z])",x, 1))
        snv = signature.index.map(lambda x: regexSearch("([A-Z])\[([A-Z>]+)\]([A-Z])",x, 2))
        context_2 = signature.index.map(lambda x: regexSearch("([A-Z])\[([A-Z>]+)\]([A-Z])",x, 3))
        order = np.argsort(np.unique(snv, return_inverse=True)[1]*100 +                           np.unique(context_1, return_inverse=True)[1]*10 +                           np.unique(context_2, return_inverse=True)[1])
        snv = snv[order]; context_1 = context_1[order]; context_2 = context_2[order]
        signature = signature.iloc[order]

        super_type = [f"{context_1[i]}[{x}]" for i,x in enumerate(snv)]
        sub_type = context_2
        mutation_labels = np.array(sub_type)

        type_counts = np.repeat(16, 6)
        colors = np.repeat(np.array(['deepskyblue','k','red','silver','limegreen','pink']), type_counts)
        plt.xlim(-1,96)
        label_freq=1#4

        pad=pad

    elif mutation_type=='DBS':
        super_type = signature.index.map(lambda x: regexSearch("([A-Z]+)>([A-Z]+)",x, 1)+">")
        sub_type = signature.index.map(lambda x: regexSearch("([A-Z]+)>([A-Z]+)",x, 2))
        mutation_labels = np.array(sub_type)

        type_counts = np.array([9,6,9,6,9,6,6,9,9,9])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        colors = np.repeat(np.array(['deepskyblue', 'steelblue','limegreen','forestgreen','salmon','firebrick',
                                     'orange','darkorange', 'violet', 'indigo']), type_counts)
        plt.xlim(-1,78)
        label_freq=1#4

        pad=pad*1.5

    elif mutation_type=='ID':
        super_type = [":".join(label.split(":")[:3])+"" for label in mutation_labels]
        sub_type = [":".join(label.split(":")[3:]) for label in mutation_labels]
        mutation_labels = np.array(sub_type)

        type_counts = np.array([12,12,24,24,11])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        colors = np.repeat(np.array(['darkorange','forestgreen','firebrick','deepskyblue','darkviolet']), type_counts)
        label_freq=1#4
        plt.xlim(-1,83)

        pad=pad*1

    elif mutation_type=='CNV':
        super_type = [":".join(label.split(":")[:2])+"" for label in mutation_labels]
        sub_type = [":".join(label.split(":")[2:]) for label in mutation_labels]
        mutation_labels = np.array(sub_type)

        label_freq=1#2
        type_counts = np.array([3,5,5,5,5,5,5,5,5,5])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        plt.xlim(-1,48)

        pad=pad*0.8

    elif mutation_type=='SV':
        super_type = ["_".join(label.split("_")[:2])+"" for label in mutation_labels]
        sub_type = ["_".join(label.split("_")[2:]) for label in mutation_labels]
        mutation_labels = np.array(sub_type)

        label_freq=1
        type_counts = np.array([5,5,5,1,5,5,5,1])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        plt.xlim(-1,32)

        pad=pad*2.2

    else: raise ValueError("mutation_type unknown")

    plt.bar(np.arange(len(signature)), signature, color=colors)
    plt.xticks(np.arange(len(signature)), mutation_labels)

    label_unique, label_counts = np.unique(super_type, return_counts=True)
    label_counts = label_counts[np.argsort(np.intersect1d(super_type, label_unique, return_indices=True)[1])]
    if xticks:
        ax.xaxis.set_tick_params(which='both', labelbottom=True,
                                 labelsize=fs*0.6*(100/len(mutation_labels))**0.15, rotation=90)
        for j,label in enumerate(ax.xaxis.get_ticklabels()):
            if (j+2)%label_freq!=0: label.set_visible(False)
        for i1 in np.cumsum(label_counts)-label_counts/2:
            plt.text(i1-0.5, -pad*np.diff(ax.get_ylim()),
                     fr"{super_type[int(i1)]}", ha='center', va='top',
                     rotation=90, fontsize=fs*0.8)#, transform=plt.gcf().transFigure)

    ylim = ax.get_ylim()
    for posx in np.cumsum(type_counts)[:-1]: plt.plot([posx-0.5,]*2,ylim,'-k', alpha=0.3)
    plt.ylim(ylim)

    plt.text(1, ylim[1]*0.9, text, ha='left', va='center', fontsize=fs)



def plotNovelSignatures(combined_sigs, relabel_map, sv_rename_deg=False):

    """
    Plot all novel signatures in subplots
    """
    hpixels=300
    hpad=13

    plot_types = ['SBS288','DBS78', 'ID83','CNV48', 'SV32']

    # width_ratios = [[len(relabel_map['degasperi'][sig_type])+len(relabel_map['novel'][sig_type])+1 for sig_type in plot_types[:3]],
    #                 [len(relabel_map['degasperi'][sig_type])+len(relabel_map['novel'][sig_type])+1 for sig_type in plot_types[3:]]]
    width_ratios = [[len(relabel_map['novel'][sig_type])+1 for sig_type in plot_types[:3]],
                    [len(relabel_map['novel']['CNV48'])+1,len(combined_sigs.keys())+1]]
    for j in range(2):
        width_ratios[j] = np.insert((np.cumsum(width_ratios[j])/np.sum(width_ratios[j])*hpixels).astype(int), 0,0)
        width_ratios[j][0]-=hpad
        width_ratios[j][-1]+=hpad
    print(width_ratios)

    fig = plt.figure(constrained_layout=False, figsize=(70,40))
    gs = fig.add_gridspec(nrows=hpixels, ncols=100, left=0.05, right=0.48,
                        wspace=0.1, hspace=0.1)

    # all_axes = []
    for i, sig_type in enumerate(plot_types):

        j=0 if i<3 else 1

        # if sig_type=='SV32': signatures = np.hstack((relabel_map['reference'][sig_type].new, relabel_map['degasperi'][sig_type].old, relabel_map['novel'][sig_type].new))
        # else: signatures = np.hstack((relabel_map['degasperi'][sig_type].old, relabel_map['novel'][sig_type].new))
        if (sig_type=='SV32')&sv_rename_deg: signatures = np.hstack((relabel_map['reference'][sig_type].new, relabel_map['degasperi'][sig_type].old, relabel_map['novel'][sig_type].new))
        elif (sig_type=='SV32'): signatures = np.hstack((relabel_map['reference'][sig_type].new, relabel_map['novel'][sig_type].new))
        else: signatures = np.hstack((relabel_map['novel'][sig_type].new))
        signatures = np.array(signatures)[orderSignatures(signatures)]
        max_val = np.max(np.array(combined_sigs[sig_type]))

        # type_axes = []
        for isig, sig in enumerate(signatures):
            hmin=width_ratios[j][i%3]+hpad
            hmax=width_ratios[j][i%3+1]-hpad
            if j==0:
                ax = fig.add_subplot(gs[
                    int(hmin+isig*(hmax-hmin)/len(signatures)):hmin+int((isig+1)*(hmax-hmin)/len(signatures))-1,0:55
                ])
            else:
                ax = fig.add_subplot(gs[
                    int(hmin+isig*(hmax-hmin)/len(signatures)):hmin+int((isig+1)*(hmax-hmin)/len(signatures))-1,60:100
                ])
            # type_axes.append(ax)

            plotSignature(ax, combined_sigs[sig_type][sig],
                        np.array(combined_sigs[sig_type].index),
                        regexSearch("([A-Z]+)[0-9+]",sig_type, 1),
                        xticks=isig==len(signatures)-1,
                        fs=24, pad=0.3)
            # plt.ylim(0,max_val*1.05)
            # plt.text(-0.5, ax.get_ylim()[1]*0.97, sig, ha='left', va='top', fontsize=24,
            #          color=map_colors['novel'] if sig in list(relabel_map['novel'][sig_type].new) else map_colors['degasperi'] if sig in list(relabel_map['degasperi'][sig_type].new) else 'k')
            plt.text(-0.5, ax.get_ylim()[1]*0.97, sig, ha='left', va='top', fontsize=24,
                     color= map_colors['degasperi'] if sig in list(relabel_map['degasperi'][sig_type].old) else map_colors['novel'] if sig in list(relabel_map['novel'][sig_type].new) else 'k')
            if isig!=len(signatures)-1: ax.set_xticklabels([])


def groupSignatureActivities(sample_df, combined_acts, cohorts='tumour_group'):

    """
    Aggregate signatures by tumour group and return number of non-zero signatures and median non-zero activity
    """

    tables = {}
    tables['count'] = {}
    tables['median'] = {}
    signatures = {}
    count_threshold = 20
    for sig_type in combined_acts:
        acts = combined_acts[sig_type].copy()# if which_sigs=="combined" else decomposed_acts.copy()

        # if sig_type=="DBS78":
        if True:
            acts = acts[np.sum(acts, axis=1)>=count_threshold]

        signatures[sig_type] = acts.keys()[orderSignatures(acts.keys())]

        acts = pd.merge(sample_df[['sample_platekey', cohorts]],
                                acts,
                                left_on='sample_platekey', right_index=True)

        unique_tissues = np.unique(acts[cohorts])
        # Make plot grids:
        tables['count'][sig_type] = np.zeros((len(signatures[sig_type]), len(unique_tissues)), dtype=int)
        tables['median'][sig_type] = np.zeros((len(signatures[sig_type]), len(unique_tissues)))

        for isig, signature in enumerate(signatures[sig_type]):
            for it, tissue in enumerate(unique_tissues):
                tables['count'][sig_type][isig,it] = np.sum(acts[signature][acts[cohorts]==tissue]>0)
                tables['median'][sig_type][isig,it] = np.median(acts[signature][acts[cohorts]==tissue][
                    acts[signature][acts[cohorts]==tissue]>0
                ])

        tables['size'] = sample_df.groupby(cohorts).size().loc[unique_tissues]

    return tables, signatures



def plotCohortActivities(tables, signatures, orientation='horizontal', rename_sigs=None):

    """
    Plot signature activities by tumour group
    """

    width_ratios = [np.sum(np.sum(tables['count'][sig_type], axis=1)>0) for sig_type in tables['count']]
    width_ratios = np.insert((np.cumsum(width_ratios[1:])/np.sum(width_ratios[1:])*300).astype(int), 0,0)
    width_ratios[0]-=2
    width_ratios[-1]+=2

    if orientation=='horizontal':
        fig = plt.figure(constrained_layout=False, figsize=(50,15))
        gs = fig.add_gridspec(nrows=2, ncols=300, left=0.05, right=0.48,
                            wspace=0.1, hspace=0.4)
    else:
        fig = plt.figure(constrained_layout=False, figsize=(35,30))
        gs = fig.add_gridspec(nrows=300, ncols=2, left=0.05, right=0.48,
                            wspace=0.3, hspace=0.1)

    size_factor = lambda x: np.log10(x+1)**1.5 * 70

    unique_tissues = list(tables['size'].index)
    group_count = np.array(tables['size'])
    unique_tissue_labels = [tissue+f"({str(tables['size'].loc[tissue])})".rjust(7) for tissue in unique_tissues]

    for i, sig_type in enumerate(tables['count'].keys()):

        if i==0:
            ax = fig.add_subplot(gs[0,:]) if orientation=='horizontal' else fig.add_subplot(gs[:,0])
        else:
            ax = fig.add_subplot(gs[1, width_ratios[i-1]+3:width_ratios[i]-3]) if orientation=='horizontal' \
                                    else fig.add_subplot(gs[width_ratios[i-1]+2:width_ratios[i]-2,1])

        sig_subset = np.sum(tables['count'][sig_type], axis=1)>0

        xx_sig, yy_tis = np.meshgrid(np.arange(len(unique_tissues)), np.arange(np.sum(sig_subset)))

        vmax = 10**int(np.ceil(np.log10(np.nanmax(tables['median'][sig_type]))))

        count_table = tables['count'][sig_type][sig_subset].copy()
        if orientation=='horizontal':
            im = ax.scatter(yy_tis.flatten()+0.5, xx_sig.flatten()+0.5,
                        s=size_factor((100*count_table/group_count[None,:]).flatten()),
                        cmap='inferno_r', norm=LogNorm(vmin=1,vmax=vmax),
                        c=tables['median'][sig_type][sig_subset].flatten())

            ax.set_yticks(np.arange(len(unique_tissues)+1), minor=True)
            ax.set_xticks(np.arange(len(signatures[sig_type][sig_subset])+1), minor=True)
            ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.5)

            ax.set_xticks(np.arange(len(signatures[sig_type][sig_subset]))+0.5, minor=False)
            ax.set_yticks(np.arange(len(unique_tissues))+0.5, minor=False)
            ax.set_xticklabels(signatures[sig_type][sig_subset], fontsize=16, minor=False, rotation=90)

            if i<2:
                ax.set_yticklabels([label.replace("_"," ") for label in unique_tissue_labels],
                                fontsize=16, minor=False);
            else:
                ax.set_yticklabels([], fontsize=16, minor=False);

            ax.set_xlim(0,np.sum(sig_subset))
            ax.set_ylim(0,len(unique_tissues))
        else:
            im = ax.scatter(xx_sig.flatten()+0.5, yy_tis.flatten()+0.5,
                        s=size_factor((100*count_table[::-1]/group_count[None,:]).flatten()),
                        cmap='inferno_r', norm=LogNorm(vmin=1,vmax=vmax),
                        c=tables['median'][sig_type][sig_subset][::-1].flatten())

            ax.set_xticks(np.arange(len(unique_tissues)+1), minor=True)
            ax.set_yticks(np.arange(len(signatures[sig_type][sig_subset])+1), minor=True)
            ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.5)

            ax.set_yticks(np.arange(len(signatures[sig_type][sig_subset]))+0.5, minor=False)
            ax.set_xticks(np.arange(len(unique_tissues))+0.5, minor=False)

            ax.set_yticklabels(signatures[sig_type][sig_subset], fontsize=16, minor=False, rotation=0)

            if (i==0)|(i==4):
                ax.set_xticklabels([label.replace("_"," ") for label in unique_tissues],
                                fontsize=16, minor=False, rotation=90);
            else:
                ax.set_xticklabels([], fontsize=16, minor=False);

            ax.set_ylim(0,np.sum(sig_subset))
            ax.set_xlim(0,len(unique_tissues))

        divider = make_axes_locatable(ax)
        if orientation=='horizontal':
            cax = divider.append_axes('top', size='5%', pad=0.05)
            cbar = fig.colorbar(im, cax=cax, orientation='horizontal');
            cbar.ax.xaxis.set_ticks_position('top')
            cbar.ax.xaxis.set_label_position('top')
        else:
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cbar = fig.colorbar(im, cax=cax, orientation='vertical');
        cbar.set_label(r"Median non-zero signature", fontsize=20, labelpad=10)
        cbar.ax.tick_params(labelsize=20)

        for isig,sig in enumerate(signatures[sig_type][sig_subset]):
            if rename_sigs is not None:
                if sig in list(rename_sigs['degasperi'][sig_type].old):
                    ax.get_xticklabels()[isig].set_color(map_colors['degasperi']) if orientation=='horizontal' \
                    else ax.get_yticklabels()[isig].set_color(map_colors['degasperi'])
                elif sig in list(rename_sigs['novel'][sig_type].new):
                    ax.get_xticklabels()[isig].set_color(map_colors['novel']) if orientation=='horizontal' \
                    else ax.get_yticklabels()[isig].set_color(map_colors['novel'])
            else:
                raise ValueError("rename_sigs not provided")


        if i==0:
            l1 = plt.scatter([],[],s=size_factor(1), c='k', alpha=0.7)
            l2 = plt.scatter([],[],s=size_factor(10), c='k', alpha=0.7)
            l3 = plt.scatter([],[],s=size_factor(100), c='k', alpha=0.7)
            if orientation=='horizontal':
                legend = ax.legend([l1,l2,l3],[f'$1$',f'$10$',f'$100$'],loc='upper left', bbox_to_anchor=(1,1),
                                        fontsize=20, ncol=1, frameon=False, title="% of samples")
            else:
                legend = ax.legend([l1,l2,l3],[f'$1\%$',f'$10\%$',f'$100\%$'],loc='lower left', bbox_to_anchor=(0,1),
                                        fontsize=20, ncol=3, frameon=False, title="% of samples")
            legend.get_title().set_fontsize('25')


def fisherTests(combined_act_df):

    signature_present = ((np.array(combined_act_df)>0)&(np.array(combined_act_df)>np.median(combined_act_df, axis=0))).astype(int)

    TP = signature_present.T @ signature_present
    FP = signature_present.T @ (1-signature_present)
    FN = (1-signature_present.T) @ signature_present
    TN = (1-signature_present.T) @ (1-signature_present)
    contingency = np.moveaxis(np.array([[TP,FP],[FN,TN]]),(2,3),(0,1))

    fisher_test = pvalue_npy(contingency[:,:,0,0].flatten().astype(np.uint64),
            contingency[:,:,1,0].flatten().astype(np.uint64),
            contingency[:,:,0,1].flatten().astype(np.uint64),
            contingency[:,:,1,1].flatten().astype(np.uint64))
    fisher_pv = fisher_test[1].reshape(contingency.shape[:2])+1e-300

    # Don't need to run this
    fisher_pv = np.zeros(contingency.shape[:2])
    fisher_stat = np.zeros(contingency.shape[:2])
    for i in tqdm.tqdm(range(contingency.shape[0]), total=contingency.shape[0]):
        for j in range(i+1):
            fisher_test = scipy.stats.fisher_exact(contingency[i,j], alternative='two-sided')
            fisher_stat[i,j] = fisher_test[0]
            fisher_pv[i,j] = fisher_test[1]
            fisher_pv[i,j] += 1e-300
    fisher_stat = fisher_stat + fisher_stat.T - np.eye(fisher_stat.shape[0])
    fisher_pv = fisher_pv + fisher_pv.T - np.eye(fisher_pv.shape[0])

    return fisher_stat, fisher_pv


def clusterHierarchical(combined_act_df):

    X = np.log(np.array(combined_act_df) + 1)

    nan_subset = np.std(X, axis=0)==0
    X -= np.mean(X, axis=0)
    X /= np.std(X, axis=0)
    X = X[:,~nan_subset]

    # Selecting Annual Income and Spending Scores by index
    clusters = shc.linkage(pd.DataFrame(X.T, index=combined_act_df.keys()[~nan_subset]),
                method='ward',
                metric="euclidean")

    return clusters, nan_subset

def pearsonSpearmen(combined_act_df):

    # Log(act + 1) correlations
    log_act = np.log(np.array(combined_act_df)+1)[:,:]
    pearson_pv = np.zeros((log_act.shape[1],log_act.shape[1]))
    pearson_stat = np.zeros((log_act.shape[1],log_act.shape[1]))
    spearman_pv = np.zeros((log_act.shape[1],log_act.shape[1]))
    spearman_stat = np.zeros((log_act.shape[1],log_act.shape[1]))
    for i in tqdm.tqdm(range(log_act.shape[1]), total=log_act.shape[1]):
        for j in range(i+1):
            pearson_test = scipy.stats.pearsonr(log_act[:,i], log_act[:,j])
            pearson_stat[i,j] = pearson_test[0]
            pearson_pv[i,j] = pearson_test[1]
            pearson_pv[i,j] += 1e-300

            spearman_test = scipy.stats.spearmanr(log_act[:,i], log_act[:,j])
            spearman_stat[i,j] = spearman_test[0]
            spearman_pv[i,j] = spearman_test[1]
            spearman_pv[i,j] += 1e-300

    pearson_stat = pearson_stat + pearson_stat.T - np.eye(pearson_stat.shape[0])
    pearson_pv = pearson_pv + pearson_pv.T - np.eye(pearson_pv.shape[0])

    spearman_stat = spearman_stat + spearman_stat.T - np.eye(spearman_stat.shape[0])
    spearman_pv = spearman_pv + spearman_pv.T - np.eye(spearman_pv.shape[0])

    return pearson_stat, pearson_pv, spearman_stat, spearman_pv

def squish(x, s=1):

    z = 1/(1+np.exp(-(x-0.5)/s))
    z -= 0.5
    z *= 1/(1-2/(1+np.exp(0.5/s)))
    z += 0.5

    return z

def cluster_diagram(fig, ax, fisher_stat, fisher_pv, spearman_stat, spearman_pv,
                    PVMAX=30, PVMIN=0, fs=12, use_cbar=True,
                    cmap_sp=cm.bwr, cmap_f=cm.PRGn):

    plt.sca(ax)

    fisher_grid = (np.sign(np.log(fisher_stat))*-np.log10(fisher_pv)).T
    fisher_grid[fisher_grid<-PVMAX]=-PVMAX+1e-10
    fisher_grid[fisher_grid>PVMAX]=PVMAX
    full_mtx = np.tril(spearman_stat/2+0.5) - np.eye(np.sum(spearman_stat.shape[0])) + np.triu((fisher_grid+PVMAX)/(2*PVMAX)+1)
    full_mtx[np.arange(full_mtx.shape[0]),np.arange(full_mtx.shape[0])] = 1-1e-10

    # sample the colormaps that you want to use. Use 128 from each so we get 256
    # colors in total
    colors1 = cmap_sp(squish(np.linspace(0., 1, 127), s=0.2))
    colors2 = cmap_f(np.linspace(0.1, 0.9, 128))#plt.cm.bone_r(np.linspace(0., 1, 128))
    colors2[int(np.floor((1-PVMIN/PVMAX)/2 * 128)):            int(np.ceil((1+PVMIN/PVMAX)/2 * 128))] = [1.,1.,1.,1.]
    colors = np.vstack((colors1, np.array([[0,0,0,1]]), colors2))
    map1 = mcolors.LinearSegmentedColormap.from_list('my_colormap1', colors1)
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    im = plt.pcolor(full_mtx, cmap=mymap, vmin=0, vmax=2)

    plt.plot([0,full_mtx.shape[0]], [0,full_mtx.shape[1]], '-w', linewidth=6)

    if use_cbar:
        cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0+ax.get_position().height*0.02,
                             0.02,ax.get_position().height*0.46])
        cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=-PVMAX, vmax=PVMAX),
                            cmap=mcolors.LinearSegmentedColormap.from_list('my_colormap', colors2)),
                            cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=fs*1.5)
        cbar.set_label(fr"$\pm\log_{{10}}\left(\mathrm{{Fisher}}\,\,P-\mathrm{{value}}\right)$", fontsize=fs*2)

        cax2 = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0+ax.get_position().height*0.52,
                             0.02,ax.get_position().height*0.46])
        cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=-1, vmax=1), cmap=map1),
                            cax=cax2, orientation='vertical')
        cbar.ax.tick_params(labelsize=fs*1.5)
        cbar.set_label(r"Spearman correlation", fontsize=fs*2)

    else: return colors1, colors2


def fullClusterDiagram(combined_act_df, fisher_stat, fisher_pv, spearman_stat, spearman_pv):

    clusters, nan_subset = clusterHierarchical(combined_act_df)

    pv_sws = 0.05/(fisher_stat.shape[0]*(fisher_stat.shape[0]-1)/2)

    PrOr = mcolors.LinearSegmentedColormap.from_list('PrOr',
                                                    np.vstack((cm.Purples_r(np.linspace(0,0.7,100)),
                                                                cm.Oranges(np.linspace(0,0.7,100))))[::-1])

    group_labels = {"UV":["SBS7a","SBS7b","SBS7c","SBS7d","DBS1"],
                    "Smoking":["SBS4","DBS2","ID3"],
                    "MMR":["SBS44","SBS26","SBS15","SBS57"],#,"ID12","DBS12"],
                    "Clock":["SBS1","SBS5","ID1","ID2"],
                    r"$POLE$":["SBS10a","SBS10b","SBS28","DBS10","DBS3"],
                    "APOBEC":["SBS2","SBS13"],
                    "HRD":["SBS3","ID6","SV3"]}


    fig, axes = plt.subplots(2,1,figsize=(27,30), gridspec_kw={'height_ratios':[1,9]})
    plt.subplots_adjust(hspace=0.)

    plt.sca(axes[0])
    dendro = shc.dendrogram(Z=clusters, labels=combined_act_df.keys()[~nan_subset])
    # Reorder according to Dendrogram
    intersection = np.intersect1d(dendro['ivl'], combined_act_df.keys()[~nan_subset], return_indices=True)
    order = intersection[2][np.argsort(intersection[1])]
    plt.axis('off')

    plt.sca(axes[1])
    ax = axes[1]

    fs=12
    PVMAX=30
    PVMIN=-np.log10(pv_sws)

    cluster_diagram(fig, ax,
                    fisher_stat[~nan_subset][:,~nan_subset][order][:,order],
                    fisher_pv[~nan_subset][:,~nan_subset][order][:,order],
                    spearman_stat[~nan_subset][:,~nan_subset][order][:,order],
                    spearman_pv[~nan_subset][:,~nan_subset][order][:,order],
                    PVMIN=PVMIN, PVMAX=PVMAX, fs=fs, use_cbar=True,
                    cmap_sp=cm.bwr, cmap_f=cm.PRGn)

    ax.set_xticks(np.arange(len(order)+1), minor=True)
    ax.set_yticks(np.arange(len(order)+1), minor=True)
    ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.3)

    ax.set_xticks(np.arange(len(order))+0.5, minor=False)
    ax.set_xticklabels(combined_act_df.keys()[~nan_subset][order], rotation=90, fontsize=fs, minor=False);

    ax.set_yticks(np.arange(len(order))+0.5, minor=False)
    ax.set_yticklabels(combined_act_df.keys()[~nan_subset][order], fontsize=fs, minor=False);

    colour_id = np.unique(dendro['leaves_color_list'], return_inverse=True)[1]
    for i,col_i in enumerate(colour_id):
        ax.get_xticklabels()[i].set_color(default_colours[col_i+1])
        ax.get_yticklabels()[i].set_color(default_colours[col_i+1])

    plt.sca(axes[1])

    for label, gene_group in group_labels.items():
        intersection = np.intersect1d(dendro['ivl'], gene_group, return_indices=True)
        plt.text(np.mean(intersection[1])+0.5,-10, label, ha='center',va='center', fontsize=fs*1.5,
                 color=default_colours[colour_id[intersection[1]][0]+1], rotation=30)


def clusterGroups(combined_act_df, fisher_stat, fisher_pv, spearman_stat, spearman_pv):

    clusters, nan_subset = clusterHierarchical(combined_act_df)

    pv_sws = 0.05/(fisher_stat.shape[0]*(fisher_stat.shape[0]-1)/2)

    group_labels = {"UV":["SBS7a","SBS7b","SBS7c","SBS7d","DBS1"],
                    "Smoking":["SBS4","DBS2","ID3"],
                    "MMR":["SBS44","SBS26","SBS15"],
                    "Clock 1":["SBS1","SBS5","ID1","ID2"],
                    "Clock 2":["SBS40","ID5","ID8"],
                    r"$POLE$":["SBS10a","SBS10b","SBS28","DBS10","DBS3"],
                    "APOBEC":["SBS2","SBS13"],
                    "HRD":["SBS3","ID6","SV3", "CN17"],
                    "Sarcoma":["SV4","SV6","CN8"],
                    "Chromothripsis":["CN6","CN7","CN8"],
                    "CRC":["SBS17a","SBS17b","DBS7","SBS93","ID14","SBS18","DBS4","SV7"],
                    r"$POLG$":["SBS17a","SBS17b","SBS93","DBS4","DBS7","ID14"],
                    "Clustered Rearrangments":["SV4","SV6","SV13"]}#,"SV10","CN18","DBS6"]}

    fig, axes = plt.subplots(3,2,figsize=(16,24))

    dendro = shc.dendrogram(Z=clusters, labels=combined_act_df.keys()[~nan_subset], no_plot=True)

    PVMAX=30
    PVMIN=-np.log10(pv_sws)

    for i,group in enumerate(['UV', 'Smoking', 'MMR', r'$POLE$', 'HRD', r'$POLG$']):#'Clock 1']):

        subset_indices = np.intersect1d(group_labels[group], combined_act_df.keys(), return_indices=True)[2]
        subset_labels = np.array(combined_act_df.keys()[
                        np.sum((fisher_pv[subset_indices]<pv_sws)&\
                            (fisher_stat[subset_indices]>1), axis=0)==len(group_labels[group])])
        subset_indices = np.intersect1d(subset_labels, combined_act_df.keys(), return_indices=True)
        subset_indices = subset_indices[2][np.argsort(subset_indices[1])]

        print(subset_labels)

        # Reorder according to Dendrogram
        intersection = np.intersect1d(dendro['ivl'], subset_labels, return_indices=True)
        order = intersection[2][np.argsort(intersection[1])]

        ax = axes[i//2,i%2]
        fs=16
        plt.sca(ax)
        colors1, colors2 = cluster_diagram(fig, ax,
                                    fisher_stat[subset_indices][:,subset_indices][order][:,order],
                                    fisher_pv[subset_indices][:,subset_indices][order][:,order],
                                    spearman_stat[subset_indices][:,subset_indices][order][:,order],
                                    spearman_pv[subset_indices][:,subset_indices][order][:,order],
                                    PVMIN=-np.log10(pv_sws), fs=fs, use_cbar=False,
                                    cmap_sp=cm.bwr, cmap_f=cm.PRGn)

        ax.set_xticks(np.arange(len(order)+1), minor=True)
        ax.set_yticks(np.arange(len(order)+1), minor=True)
        ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.3)

        ax.set_xticks(np.arange(len(order))+0.5, minor=False)
        ax.set_xticklabels(np.array(subset_labels)[order], rotation=90, fontsize=fs, minor=False);

        ax.set_yticks(np.arange(len(order))+0.5, minor=False)
        ax.set_yticklabels(np.array(subset_labels)[order], fontsize=fs, minor=False);

        colour_id = np.unique(dendro['leaves_color_list'], return_inverse=True)[1][np.sort(intersection[1])]
        for i,col_i in enumerate(colour_id):
            ax.get_xticklabels()[i].set_color(default_colours[col_i+1])
            ax.get_yticklabels()[i].set_color(default_colours[col_i+1])

        col_i = np.unique(dendro['leaves_color_list'], return_inverse=True)[1][
                        np.argwhere(np.array(dendro['ivl'])==group_labels[group][0])[0][0]
                    ]
        plt.title(group, fontsize=24, color=default_colours[col_i+1])

    ax = axes[1,1]
    cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0-ax.get_position().height*1.1,
                        0.02,ax.get_position().height*1.5])
    cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=-PVMAX, vmax=PVMAX),
                        cmap=mcolors.LinearSegmentedColormap.from_list('my_colormap', colors2)),
                        cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(fr"$\pm\log_{{10}}\left(\mathrm{{Fisher}}\,\,P\mathrm{{-value}}\right)$", fontsize=fs*1.5)

    cax2 = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0+ax.get_position().height*0.6,
                        0.02,ax.get_position().height*1.5])
    cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=-1, vmax=1),
                                        cmap=mcolors.LinearSegmentedColormap.from_list('my_colormap', colors1)),
                        cax=cax2, orientation='vertical')
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(r"Spearman correlation", fontsize=fs*1.5)

    plt.subplots_adjust(hspace=0.25)

def clusterSets(combined_act_df, fisher_stat, fisher_pv, spearman_stat, spearman_pv,
                groups=['UV', r'dMMR and $POLE$', r'$POLE$', 'HRD', r'$POLG$']):

    clusters, nan_subset = clusterHierarchical(combined_act_df)

    pv_sws = 0.05/(fisher_stat.shape[0]*(fisher_stat.shape[0]-1)/2)

    group_labels = {"UV":["SBS7a","SBS7b","SBS7c","SBS7d","DBS1"],
                    "Smoking":["SBS4","DBS2","ID3"],
                    "MMR":["SBS44","SBS26","SBS15"],
                    "Clock 1":["SBS1","SBS5","ID1","ID2"],
                    "Clock 2":["SBS40","ID5","ID8"],
                    r"$POLE$":["SBS10a","SBS10b","SBS28","DBS10","DBS3"],
                    r"dMMR and $POLE$":["SBS44", "SBS26", "SBS10a","SBS10b","SBS28","DBS10","DBS3"],
                    "APOBEC":["SBS2","SBS13"],
                    "HRD":["SBS3","ID6","SV3", "CN17"],
                    "HRD and APOBEC":["SBS3","ID6","SV3", "CN17", "SBS2", "SBS13"],
                    "Sarcoma":["SV4","SV6","CN8"],
                    "Chromothripsis":["CN6","CN7","CN8"],
                    "CRC":["SBS17a","SBS17b","DBS7","SBS93","ID14","SBS18","DBS4","SV7"],
                    r"$POLG$":["SBS17a","SBS17b","SBS93","DBS4","DBS7","ID14"],
                    "Clustered Rearrangments":["SV4","SV6","SV13"]}#,"SV10","CN18","DBS6"]}

    nrow = int((len(groups)-0.1)//2+1)
    ncol = 2
    fig, axes = plt.subplots(nrow,ncol,figsize=(ncol*8,nrow*8))

    dendro = shc.dendrogram(Z=clusters, labels=combined_act_df.keys()[~nan_subset], no_plot=True)

    PVMAX=30
    PVMIN=-np.log10(pv_sws)

    for i,group in enumerate(groups):#'Clock 1']):

        # Reorder according to Dendrogram
        intersection = np.intersect1d(dendro['ivl'], group_labels[group], return_indices=True)

        # Get groups which are included and all of their signatures
        dendro_colors = np.unique(np.array(dendro['leaves_color_list'])[intersection[1]])
        subset_labels = np.array(dendro['ivl'])[np.in1d(dendro['leaves_color_list'], dendro_colors)]
        subset_indices = np.intersect1d(subset_labels, combined_act_df.keys(), return_indices=True)
        subset_indices = subset_indices[2][np.argsort(subset_indices[1])]

        # Reorder according to Dendrogram
        intersection = np.intersect1d(dendro['ivl'], subset_labels, return_indices=True)
        order = intersection[2][np.argsort(intersection[1])]

        ax = axes[i//ncol,i%ncol]
        fs=16
        plt.sca(ax)
        colors1, colors2 = cluster_diagram(fig, ax,
                                    fisher_stat[subset_indices][:,subset_indices][order][:,order],
                                    fisher_pv[subset_indices][:,subset_indices][order][:,order],
                                    spearman_stat[subset_indices][:,subset_indices][order][:,order],
                                    spearman_pv[subset_indices][:,subset_indices][order][:,order],
                                    PVMIN=-np.log10(pv_sws), fs=fs, use_cbar=False,
                                    cmap_sp=cm.bwr, cmap_f=cm.PRGn)

        ax.set_xticks(np.arange(len(order)+1), minor=True)
        ax.set_yticks(np.arange(len(order)+1), minor=True)
        ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.3)

        ax.set_xticks(np.arange(len(order))+0.5, minor=False)
        ax.set_xticklabels(np.array(subset_labels)[order], rotation=90, fontsize=fs, minor=False);

        ax.set_yticks(np.arange(len(order))+0.5, minor=False)
        ax.set_yticklabels(np.array(subset_labels)[order], fontsize=fs, minor=False);

        colour_id = np.unique(dendro['leaves_color_list'], return_inverse=True)[1][np.sort(intersection[1])]
        for i,col_i in enumerate(colour_id):
            ax.get_xticklabels()[i].set_color(default_colours[col_i+1])
            ax.get_yticklabels()[i].set_color(default_colours[col_i+1])

        col_i = np.unique(dendro['leaves_color_list'], return_inverse=True)[1][
                        np.argwhere(np.array(dendro['ivl'])==group_labels[group][0])[0][0]
                    ]
        plt.title(group, fontsize=24, color=default_colours[col_i+1])

    ax = axes[1,1]
    cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0-ax.get_position().height*1.1,
                        0.02, ax.get_position().height*1.5])
    cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=-PVMAX, vmax=PVMAX),
                        cmap=mcolors.LinearSegmentedColormap.from_list('my_colormap', colors2)),
                        cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(fr"$\pm\log_{{10}}\left(\mathrm{{Fisher}}\,\,P\mathrm{{-value}}\right)$", fontsize=fs*1.5)

    cax2 = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0+ax.get_position().height*0.6,
                        0.02,ax.get_position().height*1.5])
    cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=-1, vmax=1),
                                        cmap=mcolors.LinearSegmentedColormap.from_list('my_colormap', colors1)),
                        cax=cax2, orientation='vertical')
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(r"Spearman correlation", fontsize=fs*1.5)

    plt.subplots_adjust(hspace=0.25)

    for iax in range(len(groups), ncol*nrow):
        plt.sca(axes[iax//ncol,iax%ncol])
        plt.axis('off')


def plotBoxAndWhiskers(ax, samples, x, width, color='k', label=None):
    percentiles = np.percentile(samples, np.array([2.5,16,50,84,97.5]))
    patch = Rectangle((x-width*2/5, percentiles[1]), width*4/5, percentiles[3]-percentiles[1],
                      linewidth=1, edgecolor=color, facecolor=color, alpha=0.5)
    ax.add_patch(patch)

    plt.plot([x,x], [percentiles[0], percentiles[1]], linewidth=1, c=color)
    plt.plot([x,x], [percentiles[3], percentiles[4]], linewidth=1, c=color)
    plt.plot([x-width/3,x+width/3], [percentiles[2],]*2, linewidth=1, c=color)
    plt.plot([x-width/2,x+width/2], [percentiles[0],]*2, linewidth=1, c=color)
    plt.plot([x-width/2,x+width/2], [percentiles[4],]*2, linewidth=1, c=color, label=label)


def getMatrices(sample_df):
    # Load in full mutation matrices for each mutation type
    sig_totals = pd.DataFrame(index=sample_df.sample_platekey)

    cohorts = ["_".join(folder.split("_")[1:]) for folder in os.listdir("/re_gecip/cancer_pan/fprefect/botl/results/SIGmats/v3/matrices/") \
                    if "CANCERTYPE" in folder and not "PASS" in folder and not "_test" in folder]

    matrices = {"SBS":"SBS288", "DBS":"DBS78", "ID":"ID83"}
    drafts = {"SBS":"v2_draft", "DBS":"v2_draft", "ID":"v3"}
    max_sigs = {"SBS":30, "DBS":20, "ID":25}
    sig_mats = {}
    for mut_type in matrices:
        print(mut_type)
        for i,cohort in enumerate(cohorts):
            cohort_matrix = pd.read_csv(f"/re_gecip/cancer_pan/fprefect/botl/results/SIGmats/{drafts[mut_type]}/{matrices[mut_type]}/"\
                                        f"CANCERTYPE_{cohort}/{matrices[mut_type]}_random_max{max_sigs[mut_type]}sigs_500nmf_reps/Samples.txt",
                                            sep="\t", index_col=0)
            if i==0:
                sig_mats[matrices[mut_type]]=cohort_matrix.copy()
            else:
                sig_mats[matrices[mut_type]]=pd.merge(sig_mats[matrices[mut_type]], cohort_matrix,
                                            left_index=True, right_index=True, how='inner')
        sig_mats[matrices[mut_type]] = sig_mats[matrices[mut_type]].T
        if mut_type=='ID':
            sig_mats['ID83'] = pd.merge(sig_mats['ID83'], sample_df[['tumour_sample_platekey', 'sample_platekey']],
                                    left_index=True, right_on='tumour_sample_platekey', how='inner')\
                                    .drop('tumour_sample_platekey', axis=1).set_index('sample_platekey')
        sig_totals = pd.merge(sig_totals, pd.DataFrame(np.sum(sig_mats[matrices[mut_type]], axis=1), columns=[mut_type]),
                                                    left_index=True, right_index=True, how='left')

    sig_dir = f"{DATA_DIR}/copy_number_signatures/SigProfilerCNV48/"
    cohorts = [folder for folder in os.listdir(sig_dir) \
                    if not folder in ["eofiles", "Summary", "Pan"]]

    for i,cohort in enumerate(cohorts):
        cohort_matrix = pd.read_csv(f"{sig_dir}/{cohort}/sigs_1_to_15_bestrun_inclFailed/CNV48/Samples.txt",
                                    sep="\t", index_col=0)
        if i==0: sig_mat=cohort_matrix.copy()
        else: sig_mat=pd.merge(sig_mat, cohort_matrix, left_index=True, right_index=True, how='inner')
    sig_mats['CNV48'] = sig_mat.T
    sig_mats['CNV48'].set_index(sig_mats['CNV48'].index.map(lambda x: "_".join(x.split("_")[3:])), inplace=True)
    sig_mats['CNV48'] = pd.merge(sig_mats['CNV48'], sample_df[['tumour_sample_platekey', 'sample_platekey']],
                                left_index=True, right_on='tumour_sample_platekey', how='inner')\
                                .drop('tumour_sample_platekey', axis=1).set_index('sample_platekey')
    sig_totals = pd.merge(sig_totals, pd.DataFrame(np.sum(sig_mats['CNV48'], axis=1), columns=['CN']),
                          left_index=True, right_index=True, how='left')
    # sig_totals[f"TCNV"] = np.sum(sig_mats['CNV48'], axis=1)
    # sig_mat.set_index(sig_mat.index.map(lambda x: "_".join(x.split("_")[1:3])), inplace=True)

    sig_dir="/re_gecip/cancer_pan/atapinos/Results_SigProfiler_v2/"
    cohorts = [folder.split("_S1")[0] for folder in os.listdir(sig_dir) \
                    if "S1_S15_R500_PoisNoisRes_T" in folder and not "AllSamples" in folder]

    for i,cohort in enumerate(cohorts):
        cohort_matrix = pd.read_csv(f"{sig_dir}/{cohort}_S1_S15_R500_PoisNoisRes_T/SV32/Samples.txt",
                                    sep="\t", index_col=0)
        if i==0: sig_mat=cohort_matrix.copy()
        else: sig_mat=pd.merge(sig_mat, cohort_matrix, left_index=True, right_index=True, how='inner')
    sig_mats['SV32'] = sig_mat.T
    sig_totals = pd.merge(sig_totals, pd.DataFrame(np.sum(sig_mats['SV32'], axis=1), columns=['SV']),
                          left_index=True, right_index=True, how='left')

    sig_totals.fillna(0, inplace=True)

    return sig_mats, sig_totals

def plotMutationRates(sig_totals, groups=[], sig_types=['SBS', 'DBS', 'ID', 'CN', 'SV']):

    fig, ax = plt.subplots(1,1,figsize=(12,5))
    default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for j, sig_type in enumerate(sig_types):
    #     for j, sig_type in enumerate(['DBS']):
        for i, group in enumerate(groups):
    #         print(group)
            c = default_colours[j]
            plotBoxAndWhiskers(ax, sig_totals.loc[group][sig_type], i+0.5-(2-j)*0.15, 0.1, color=c,
                            label=sig_type if i==0 else None)

    plt.yscale('log')

    ax.set_xticks(np.arange(len(groups)+1), minor=True)
    ax.grid(True, which='minor', axis='x', linestyle='-', color='k')

    ax.set_xticks(np.arange(len(groups))+0.5, minor=False)
    ax.set_xticklabels([group.replace("_"," ") for group in groups], rotation=90, fontsize=16, minor=False);

    plt.ylabel("Total mutation rate")

    plt.xlim(0,len(groups))

    plt.legend(loc='lower right', bbox_to_anchor=(1,1), ncol=5)


if __name__=='__main__':

    publish_dir = FIGURE_DIR
    if not os.path.exists(publish_dir): os.mkdir(publish_dir)

    ### Load sample list
    samples_file = f"{DATA_DIR}/sample_lists_incl_SEGs/sample_list.tsv"
    sample_df = pd.read_csv(samples_file, usecols=['participant_id', 'tumour_sample_platekey', 'germline_sample_platekey',
                                                'tumour_group', 'signature_extraction_group'], sep="\t")
    sample_df['tumour_tissue'] = sample_df.tumour_group.map(lambda x: x.split("-")[0])
    sample_df['sample_platekey']=sample_df.participant_id.astype(str)+"_"\
                        +sample_df.tumour_sample_platekey.astype(str)+"_"\
                        +sample_df.germline_sample_platekey.astype(str)

    ### Load signatures and activities
    combined_acts, combined_sigs, relabel_map = loadSignatures()

    # Save relabel mappings
    relabel_dir = f"{COMBINED_SIGS_DIR}/combinedSignatures_relabel"
    if not os.path.exists(relabel_dir): os.mkdir(relabel_dir)
    for sig_set in combined_acts:
        relabel_set = pd.concat((relabel_map['reference'][sig_set],
                                relabel_map['degasperi'][sig_set],
                                relabel_map['novel'][sig_set]))
        relabel_set.to_csv(f"{relabel_dir}/{sig_set}.tsv", sep="\t")

    # Plot novel signature profiles
    plotNovelSignatures(combined_sigs, relabel_map)
    publish_fig("signature_distributions_alltypes", publish=publish_dir)

    # Load activity tables
    tables, signatures = groupSignatureActivities(sample_df, combined_acts, cohorts='signature_extraction_group')
    plotCohortActivities(tables, signatures, rename_sigs=relabel_map)
    publish_fig(f"signature_sizes_alltypes", publish=publish_dir)
    plotCohortActivities(tables, signatures, rename_sigs=relabel_map, orientation='vertical')
    publish_fig(f"signature_sizes_alltypes_vertical", publish=publish_dir)

    # Activities
    combined_act_df = pd.DataFrame(index=list(sample_df.sample_platekey))
    for sig_type in ["SBS288", "DBS78", "ID83", "CNV48", "SV32"]:
        combined_act_df = pd.merge(combined_act_df, combined_acts[sig_type],
                                    left_index=True, right_index=True, how='left').fillna(0)

    # Get fisher pvalues of signature activity associations
    fisher_stat, fisher_pv = fisherTests(combined_act_df)
    # Get pearson and spearman association coefficients
    pearson_stat, pearson_pv, spearman_stat, spearman_pv = pearsonSpearmen(combined_act_df)
    # Run hieriarchical clustering
    clusters = clusterHierarchical(combined_act_df)

    # Plot cluster diagram
    fullClusterDiagram(combined_act_df, fisher_stat, fisher_pv, spearman_stat, spearman_pv)
    publish_fig(f"signature_activity_clusters", publish=publish_dir)

    # Plot cluster diagram in signature type groups
    clusterGroups(combined_act_df, fisher_stat, fisher_pv, spearman_stat, spearman_pv)
    publish_fig(f"signature_activity_cluster_subgroups", publish=publish_dir)

    # Plot cluster diagram in signature type groups
    clusterSets(combined_act_df, fisher_stat, fisher_pv, spearman_stat, spearman_pv,
                groups=['UV', r'dMMR and $POLE$', r'$POLG$', 'HRD and APOBEC', 'Smoking'])
    publish_fig(f"signature_activity_cluster_subsets", publish=publish_dir)

    # Plot mutation rates across tumour groups
    sig_mats, sig_totals = getMatrices(sample_df)
    sig_totals = pd.merge(sig_totals, sample_df[['signature_extraction_group', 'sample_platekey']],
                     left_index=True, right_on='sample_platekey', how='inner').set_index('signature_extraction_group')
    groups = np.unique(sample_df.signature_extraction_group)
    plotMutationRates(sig_totals, groups=groups)
    publish_fig(f"total_mutation_rates", publish=publish_dir)
