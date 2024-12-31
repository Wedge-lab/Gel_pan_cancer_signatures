"combine signatures between cohorts - iteratively add cohort extracted signatures to pan-cancer COSMIC list"
import sys, os
import pandas as pd, numpy as np, re
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from signatures.utils import orderSignatures

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)

def plotSignature(ax, signature, mutation_labels, mutation_type="SBS", xticks=True, text="", fs=12,
                    ylim=None):

    plt.sca(ax)

    if mutation_type=='SBS':
        context_1 = signature.index.map(lambda x: re.search("([A-Z])\[([A-Z>]+)\]([A-Z])",x).group(1)) # type: ignore
        snv = signature.index.map(lambda x: re.search("([A-Z])\[([A-Z>]+)\]([A-Z])",x).group(2)) # type: ignore
        context_2 = signature.index.map(lambda x: re.search("([A-Z])\[([A-Z>]+)\]([A-Z])",x).group(3)) # type: ignore
        order = np.argsort(np.unique(snv, return_inverse=True)[1]*100 +\
                           np.unique(context_1, return_inverse=True)[1]*10 +\
                           np.unique(context_2, return_inverse=True)[1])
        signature = signature.iloc[order]
        mutation_labels = np.array(mutation_labels)[order]
        #colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:6], 16)
        #colors = np.repeat(np.array(['deepskyblue','k','red','silver','limegreen','pink']), 16)
        type_counts = np.repeat(16, 6)
        colors = np.repeat(np.array(['deepskyblue','k','red','silver','limegreen','pink']), type_counts)
        plt.xlim(-1,96)
        label_freq=1#4
    elif mutation_type=='DBS':
        type_counts = np.array([9,6,9,6,9,6,6,9,9,9])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        colors = np.repeat(np.array(['deepskyblue', 'steelblue','limegreen','forestgreen','salmon','firebrick','orange','darkorange', 'violet', 'indigo']), type_counts)
        plt.xlim(-1,78)
        label_freq=1#4
    elif mutation_type=='ID':
        type_counts = np.array([12,12,24,24,11])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        colors = np.repeat(np.array(['darkorange','forestgreen','firebrick','deepskyblue','darkviolet']), type_counts)
        label_freq=1#4
        plt.xlim(-1,83)
    elif mutation_type=='CNV':
        label_freq=1#2
        type_counts = np.array([3,5,5,5,5,5,5,5,5,5])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        plt.xlim(-1,48)
    elif mutation_type=='SV':
        label_freq=1
        type_counts = np.array([5,5,5,1,5,5,5,1])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        plt.xlim(-1,32)
    else: raise ValueError("mutation_type unknown")

    plt.bar(mutation_labels, signature, color=colors)

    if xticks:
        ax.xaxis.set_tick_params(which='both', labelbottom=True, labelsize=8 if len(mutation_labels)>60 else 12, rotation=90)
        for j,label in enumerate(ax.xaxis.get_ticklabels()):
            if (j+2)%label_freq!=0: label.set_visible(False)

    ylim = ax.get_ylim() if ylim is None else ylim
    for posx in np.cumsum(type_counts)[:-1]: plt.plot([posx-0.5,]*2,ylim,'-k', alpha=0.3)
    plt.ylim(ylim)

    plt.text(1, ylim[1]*0.9, text, ha='left', va='center', fontsize=fs)

def plotSignature2(ax, signature, mutation_labels, mutation_type="SBS", xticks=True, text="", fs=12,
                 pad=0.1):

    plt.sca(ax)

    if mutation_type=='SBS':
        context_1 = signature.index.map(lambda x: re.search("([A-Z])\[([A-Z>]+)\]([A-Z])",x).group(1)) # type: ignore
        snv = signature.index.map(lambda x: re.search("([A-Z])\[([A-Z>]+)\]([A-Z])",x).group(2)) # type: ignore
        context_2 = signature.index.map(lambda x: re.search("([A-Z])\[([A-Z>]+)\]([A-Z])",x).group(3)) # type: ignore
        order = np.argsort(np.unique(snv, return_inverse=True)[1]*100 +\
                           np.unique(context_1, return_inverse=True)[1]*10 +\
                           np.unique(context_2, return_inverse=True)[1])
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
        super_type = signature.index.map(lambda x: re.search("([A-Z]+)>([A-Z]+)",x).group(1)+">") # type: ignore
        sub_type = signature.index.map(lambda x: re.search("([A-Z]+)>([A-Z]+)",x).group(2)) # type: ignore
        mutation_labels = np.array(sub_type)

        type_counts = np.array([9,6,9,6,9,6,6,9,9,9])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        colors = np.repeat(np.array(['deepskyblue', 'steelblue','limegreen','forestgreen','salmon','firebrick','orange','darkorange', 'violet', 'indigo']), type_counts)
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

        plt.bar(np.arange(len(signature)), signature, color=colors)
        plt.xticks(np.arange(len(signature)), mutation_labels)

        pad=pad*6

    elif mutation_type=='SV':
        super_type = ["_".join(label.split("_")[:2])+"" for label in mutation_labels]
        sub_type = ["_".join(label.split("_")[2:]) for label in mutation_labels]
        mutation_labels = np.array(sub_type)

        label_freq=1
        type_counts = np.array([5,5,5,1,5,5,5,1])
        colors = np.repeat(plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(type_counts)], type_counts)
        plt.xlim(-1,32)

        pad=pad*4

    else: raise ValueError("mutation_type unknown")

    plt.bar(np.arange(len(signature)), signature, color=colors)
    plt.xticks(np.arange(len(signature)), mutation_labels)

    label_unique, label_counts = np.unique(super_type, return_counts=True)
    label_counts = label_counts[np.argsort(np.intersect1d(super_type, label_unique, return_indices=True)[1])]
    if xticks:
        ax.xaxis.set_tick_params(which='both', labelbottom=True, labelsize=fs*0.6, rotation=90)
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

def plot_cosine_similarity(matrix, xlabels, ylabels, fs=30):

    fig, ax = plt.subplots(1,figsize=(matrix.shape[1]*0.7,
                                      matrix.shape[0]*0.7),
                          sharey=True, constrained_layout=True)

    plt.sca(ax)
    breaks = np.arange(0,101,5)
    blue_cmap = mpl.colors.LinearSegmentedColormap.from_list("seg_blue",
                                                             plt.get_cmap("viridis")(np.linspace(0.1,1,len(breaks)-1)),
                                                             20)
    blue_cmap.set_under(plt.get_cmap("bwr")(0.5))
    im = plt.pcolor(matrix*100, cmap=blue_cmap, vmin=1e-10, vmax=100)

    def rect(pos, colour):
        r = plt.Rectangle(pos, 0.88,0.8, facecolor="none", edgecolor=colour, linewidth=2)
        plt.gca().add_patch(r)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if ylabels[i]==xlabels[j]:
                rect((j+0.06,i+0.1),'r')

    cax = fig.add_axes([ax.get_position().x0+0.05,ax.get_position().y1+0.01,
                        ax.get_position().width*9/10,0.02])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=breaks+1e-10)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    # cbar.ax.set_xticks(np.arange(0,101,20))
    cbar.ax.set_xticklabels(breaks, fontsize=fs*1.5)
    cbar.set_label(r"Cosine Similarity", fontsize=fs*1.8, labelpad=40)

    ax.set_xticks(np.arange(len(xlabels)+1), minor=True)
    ax.set_yticks(np.arange(len(ylabels)+1), minor=True)
    ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.3)

    ax.set_xticks(np.arange(len(xlabels))+0.5, minor=False)
    ax.set_xticklabels(xlabels, rotation=90, fontsize=fs, minor=False)

    ax.set_yticks(np.arange(len(ylabels))+0.5, minor=False)
    ax.set_yticklabels(ylabels, fontsize=fs, minor=False)

    plt.subplots_adjust(wspace=0.01)

    return fig, ax


if __name__=='__main__':

    dir_output, sig_type, samples_file = sys.argv[1:]
    if not os.path.exists(f'{dir_output}/figs'): os.mkdir(f'{dir_output}/figs')

    combined_sigs = pd.read_csv(f"{dir_output}/Combined_Solution_Signatures.tsv", sep="\t", index_col=0)
    combined_acts = pd.read_csv(f"{dir_output}/Combined_Solution_Activities.tsv", sep="\t", index_col=0)
    combined_maps = pd.read_csv(f"{dir_output}/Combined_Solution_Mapped.tsv", sep="\t", index_col=0)
    combined_maps.old_signatures = combined_maps.old_signatures.apply(lambda x: x[2:-2].split("', '"))

    ##### ------------ Combiner Assignment Figure ------------ ####
    input_sigs = pd.read_csv(f"{dir_output}/Decomposed_Solution_Signatures.tsv",
                             sep="\t", index_col=0).keys()
    try:
        denovo_map = pd.read_csv(glob.glob(f"{dir_output}/iter_0/Decompose_Solution/De_Novo_map_to_COSMIC_*.csv")[0],
                                 skipinitialspace=True)
        novel_sigs = input_sigs[denovo_map['De novo extracted']==denovo_map['Global NMF Signatures']]
    except IndexError:
        # If no input reference then it wouldn't have produced the file
        novel_sigs = input_sigs

    lost_map = {}
    for index, row in combined_maps.reset_index().iterrows():

        denovo_map = pd.read_csv(glob.glob(f"{dir_output}/iter_{index+1}/Decompose_Solution/De_Novo_map_to_COSMIC_*.csv")[0],
                                 skipinitialspace=True)

        remaining_sigs = input_sigs[denovo_map['De novo extracted']==denovo_map['Global NMF Signatures']]

        lost_map[index+1] = np.setxor1d(novel_sigs, np.concatenate((remaining_sigs, row.old_signatures)))
        novel_sigs = remaining_sigs

    plt.figure(figsize=(8,20))
    print(len(combined_maps))
    i=1
    plt.ylim(-1,1)

    for idx, row in combined_maps.reset_index().iterrows():
        old_sigs = row['old_signatures']
        lost_sigs = lost_map[idx+1]

        r = plt.Rectangle((0., 1.025-i*0.05-(len(old_sigs)+int(len(lost_sigs)>0))*0.05),
                          1., (len(old_sigs)+int(len(lost_sigs)>0))*0.05,
                          facecolor="none", edgecolor=plt.get_cmap("viridis")(0.5), linewidth=2)
        plt.gca().add_patch(r)

        plt.text(0.95, 1-i*0.05-(len(old_sigs[:-1])*0.05/2), row['new_signatures'],
                 fontsize=20, va='center', ha='right')
        for j, sig in enumerate(old_sigs):
            plt.text(0.05, 1-i*0.05, sig, fontsize=16, va='center')
            i+=1
        if len(lost_sigs)>0:
            plt.text(0.05, 1-i*0.05, f"({', '.join(lost_sigs[:3])}...)", fontsize=16,
                     va='center', color='r', alpha=0.5)
            i+=1
        i+=0.5
    plt.axis('off')

    plt.savefig(f'{dir_output}/figs/combined_maps_{sig_type}.png',
               bbox_inches='tight', dpi=200, facecolor='w', transparent=False)


    ##### ------------ New signature distribution figure ------------ ####
    decomposed_sigs = pd.read_csv(f"{dir_output}/Decomposed_Solution_Signatures.tsv", sep="\t", index_col=0)
    if sig_type=="SBS288": decomposed_sigs = decomposed_sigs.groupby(decomposed_sigs.index.str[2:9]).sum()

    mut_type = re.search("([A-Z]+)[0-9]+", sig_type).group(1)
    for row in combined_maps[combined_maps.old_signatures.map(lambda x: len(x)>0)].iterrows():
        new_signature = row[0]
        old_signatures = row[1].old_signatures

        fig, axes = plt.subplots(len(old_signatures), 2, figsize=(12, len(old_signatures)*4), sharex='col')

        if len(old_signatures)==1: ax=axes[0]
        else: ax=axes[0,0]
        plotSignature(ax, combined_sigs[new_signature], list(combined_sigs.index),
                      mutation_type=mut_type, text=new_signature, fs=16)
        for i in range(1,len(old_signatures)):
            if len(old_signatures)==1: ax=axes[0]
            else: ax=axes[i,0]
            plt.sca(ax)
            plt.axis('off')

        for i, sig in enumerate(old_signatures):
            if len(old_signatures)==1: ax=axes[1]
            else: ax=axes[i,1]
            cossim = np.sum(np.array(decomposed_sigs[sig])*np.array(combined_sigs[new_signature]))/\
                    np.sqrt(np.sum(decomposed_sigs[sig]**2)*np.sum(combined_sigs[new_signature]**2))
            plotSignature(ax, decomposed_sigs[sig], list(combined_sigs.index), mutation_type=mut_type,
                         xticks=(i==len(old_signatures)-1),
                          text=f"{sig}, cos-sim={cossim:.3f}", fs=16)

        plt.subplots_adjust(hspace=0.05)

        print(new_signature, old_signatures)

        plt.savefig(f'{dir_output}/figs/{new_signature}_combined_sigs.png',
                   bbox_inches='tight', dpi=200, facecolor='w', transparent=False)


    ##### ------------ New signature distribution figure ------------ ####
    decomposed_sigs = pd.read_csv(f"{dir_output}/Decomposed_Solution_Signatures.tsv", sep="\t", index_col=0)
    if sig_type=="SBS288": decomposed_sigs = decomposed_sigs.groupby(decomposed_sigs.index.str[2:9]).sum()

    mut_type = re.search("([A-Z]+)[0-9]+", sig_type).group(1)
    new_signatures = combined_maps.index#combined_sigs.keys()
    if len(new_signatures)>0:
        ncol = int(np.ceil(np.sqrt(len(new_signatures))))
        nrow = int(np.ceil(len(new_signatures)/ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*8, nrow*5), sharex='col')

        for i in range(ncol*nrow):
            if ncol*nrow==1: ax=axes
            elif nrow==1: ax=axes[i]
            else: ax=axes[i//ncol,i%ncol]

            if i<len(new_signatures):
                plotSignature(ax, combined_sigs[new_signatures[i]], list(combined_sigs.index),
                              mutation_type=mut_type, text=new_signatures[i], fs=16)
            else:
                plt.sca(ax)
                plt.axis('off')

        plt.subplots_adjust(hspace=0.3, wspace=0.1)

        plt.savefig(f'{dir_output}/figs/new_combined_sigs.png',
                   bbox_inches='tight', dpi=200, facecolor='w', transparent=False)


    ##### ------------ Activity Map ------------ ####
    index_col='tumour_sample_platekey'
    if sig_type in ['SBS288','DBS78','ID83','SV32','CNV48']: combined_acts.set_index(combined_acts.index.map(lambda x: "_".join(x.split("_")[1:3])), inplace=True)
    print(combined_acts.head())
    # if sig_type=="SV32": index_col='participant_id'
    sample_df = pd.read_csv(samples_file, usecols=[index_col, 'tumour_group', 'signature_extraction_group'], sep="\t")
    sample_df['tumour_tissue'] = sample_df.tumour_group.map(lambda x: x.split("-")[0])

    which_sigs = "combined"
    cohorts = 'signature_extraction_group'

    acts = combined_acts.copy() if which_sigs=="combined" else decomposed_acts.copy()

    sample_df[index_col] = sample_df[index_col].astype(str)
    acts.index = acts.index.map(str)

    signatures = acts.keys()[orderSignatures(acts.keys())]

    acts = pd.merge(sample_df[[index_col, cohorts]], acts,
                    left_on=index_col, right_index=True, how='inner')
    print(acts.head)
    print(acts.shape)

    unique_tissues = np.unique(acts[cohorts])
    # Make plot grids:
    tables = {}
    tables['count'] = np.zeros((len(signatures), len(unique_tissues)), dtype=int)
    tables['median'] = np.zeros((len(signatures), len(unique_tissues)))

    for isig, signature in enumerate(signatures):
        for it, tissue in enumerate(unique_tissues):
            tables['count'][isig,it] = np.sum(acts[signature][acts[cohorts]==tissue]>0)
            tables['median'][isig,it] = np.median(acts[signature][acts[cohorts]==tissue][
                acts[signature][acts[cohorts]==tissue]>0
            ])

    sig_subset = np.sum(tables['count'], axis=1)>0
    fig, ax = plt.subplots(1,1,figsize=(tables['count'].shape[1]*0.7,np.sum(sig_subset)*0.5))

    vmin = 10**np.floor(np.nanmin(np.log10(tables['median'][sig_subset])))
    vmax = 10**np.ceil(np.nanmax(np.log10(tables['median'][sig_subset])))

    xx_sig, yy_tis = np.meshgrid(np.arange(np.sum(sig_subset)), np.arange(len(unique_tissues)))
    im = plt.scatter(yy_tis.flatten()+0.5, xx_sig.flatten()+0.5,
                s=np.log10(tables['count'][sig_subset][::-1].T.flatten()+1)*100,
                # cmap='viridis', norm=LogNorm(vmin=vmin,vmax=1e5 if sig_type=="SBS288" else vmax),
                c=tables['median'][sig_subset][::-1].T.flatten())

    ax.set_xticks(np.arange(len(unique_tissues)+1), minor=True)
    ax.set_yticks(np.arange(len(signatures[sig_subset])+1), minor=True)
    ax.grid(True, which='minor', axis='both', linestyle='-', color='k', alpha=0.5)

    ax.set_yticks(np.arange(len(signatures[sig_subset]))+0.5, minor=False)
    ax.set_yticklabels([sig.replace("-","-") for sig in signatures[sig_subset][::-1]],
                                    fontsize=16, minor=False);
    ax.set_xticks(np.arange(len(unique_tissues))+0.5, minor=False)
    ax.set_xticklabels(unique_tissues, rotation=90, fontsize=16, minor=False);

    plt.ylim(0,np.sum(sig_subset))
    plt.xlim(0,len(unique_tissues))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('top', size=f"{50/tables['count'].shape[0]}%",
                                      pad=f"{50/tables['count'].shape[0]}%")
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.set_label(r"Median non-zero signature", fontsize=16)
    cbar.ax.tick_params(labelsize=16)

    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Scatter',
                              markerfacecolor='k', alpha=0.5, markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Scatter',
                              markerfacecolor='k', alpha=0.5, markersize=20),
                       Line2D([0], [0], marker='o', color='w', label='Scatter',
                              markerfacecolor='k', alpha=0.5, markersize=30)]
    l1 = plt.scatter([],[],s=100*np.log10(10), c='k', alpha=0.7)
    l2 = plt.scatter([],[],s=100*np.log10(100), c='k', alpha=0.7)
    l3 = plt.scatter([],[],s=100*np.log10(1000), c='k', alpha=0.7)
    ax.legend([l1,l2,l3],[f'$10$',f'$10^{2}$',f'$10^{3}$'],loc='upper left', bbox_to_anchor=(1,1),
             fontsize=16)

    plt.savefig(f"{dir_output}/figs/signature_sizes_{sig_type}_{which_sigs}_{cohorts}.png",
                bbox_inches='tight', dpi=200, facecolor='w', transparent=False)



    ##### ------------ Cosine Similarity Map ------------ ####
    sigmat = np.array(combined_sigs)/np.sqrt(np.sum(np.array(combined_sigs)**2, axis=0))[None,:]

    cosine_similarity_mtx = sigmat.T@sigmat

    XX,YY = np.meshgrid(np.arange(cosine_similarity_mtx.shape[0]),
                        np.arange(cosine_similarity_mtx.shape[0]))
    cosine_similarity_mtx[XX>=YY] = np.nan

    plot_cosine_similarity(cosine_similarity_mtx[::-1], combined_sigs.columns, combined_sigs.columns[::-1],
                          fs=16)

    plt.savefig(f"{dir_output}/figs/combined-combined_cosine_similarity_{sig_type}.png",
               bbox_inches='tight', dpi=200, facecolor='w', transparent=False)
