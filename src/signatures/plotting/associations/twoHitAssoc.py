import os
import scipy, scipy.stats
import pandas as pd, numpy as np, re
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import patches
from signatures.utils import BH_threshold, orderSignatures
from signatures.plotting.combinedSignatures import loadSignatures, signatureRenamer, map_colors
from dotenv import load_dotenv

load_dotenv()
RESULT_DIR = os.getenv("RESULT_DIR")
FIGURE_DIR = os.getenv("FIGURE_DIR")
DATA_DIR = os.getenv("DATA_DIR")

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)

default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']

group_labels = {'BileDuct-AdenoCA': 'BileDuct', 'Bladder-TCC': 'Bladder',
 'Breast-DuctalCA': 'DuctalCA', 'Breast-LobularCA': 'LobularCA',
 'CNS-Astro': 'Astro', 'CNS-GBM-IDHmut': 'GBM-IDHmut', 'CNS-GBM-IDHwt': 'GBM-IDHwt',
 'CNS-Menin': 'Menin', 'CNS-Oligo': 'Oligo',
 'ColoRect-AdenoCA': 'CRC',
 'Connective-Chondro': 'Chondro', 'Connective-Leiomyo': 'Leiomyo', 'Connective-Liposarc': 'Liposarcoma',
 'Connective-Myxofibro': 'Myxofibro', 'Connective-Osteosarc': 'Osteosarc', 'Connective-SCS': 'SCS',
 'Connective-SS': 'SS',
 'Eso-AdenoCA': 'Esophogeal',
 'Haem-ALL': 'ALL', 'Haem-AML': 'AML', 'Haem-CML': 'CML', 'Haem-MM': 'MM', 'Haem-MPN': 'MPN',
 'HeadNeck-SCC': 'Headneck',
 'Kidney-CCRCC': 'CCRCC', 'Kidney-ChRCC': 'ChRCC', 'Kidney-PRCC': 'PRCC',
 'Liver-HCC': 'HCC', 'Lung-AdenoCA': 'LungCA', 'Lung-LargeCell': 'LungLC', 'Lung-SCC': 'LungSCC',
 'Lung-SmallCell': 'LungSC',
 'Mes-Mesothelioma': 'Mesothelioma', 'Ovary-AdenoCA': 'Ovary', 'Panc-AdenoCA': 'Pancreas', 'Prost-AdenoCA': 'Prostate',
 'Skin-Melanoma': 'Melanoma', 'Stomach-AdenoCA': 'Stomach', 'Testis-GCT': 'Testis', 'Uterus-AdenoCA': 'Uterus'}


mock_match = np.vectorize(lambda x:bool(re.match('.*mock.*', x, re.IGNORECASE)))


def plotHitRates(sample_df, target_df, test_df, DNA_repair):

    # Get mutation rates
    # OncoKB
    somatic_df = pd.read_csv(f"{DATA_DIR}/cancGeneHits/somatic/output/OncoKB_somatic_cancer-gene_hits.tsv",
                            sep="\t", index_col=0)
    # LOH
    loh_df = pd.read_csv(f"{DATA_DIR}/cancGeneHits/somatic/output/lohfrac_gene_sample_matrix.tsv",
                            sep="\t").set_index("tumour_sample_platekey")
    # Add in platekeys which failed battenberg and assume no LoH
    failed_battenberg_platekeys = np.setxor1d(np.intersect1d(loh_df.index, somatic_df.index), somatic_df.index)
    loh_df = pd.concat((loh_df,
                        pd.DataFrame(np.zeros((len(failed_battenberg_platekeys), len(loh_df.keys()))),
                                    index=failed_battenberg_platekeys,
                                    columns=loh_df.keys())))
    # Only tumour suppressor genes
    group_df = pd.read_csv(f"{RESULT_DIR}/sample_lists_incl_SEGs/sample_list.tsv",
                        sep="\t", usecols=['participant_id', 'tumour_sample_platekey', 'germline_sample_platekey',
                                            'signature_extraction_group'])
    group_df['sample_platekey'] = group_df.participant_id.astype(str)+"_"+group_df.tumour_sample_platekey+"_"+\
                                group_df.germline_sample_platekey
    group_df.set_index('sample_platekey', inplace=True)

    sample_df = pd.merge(group_df[['signature_extraction_group']], sample_df,
            left_index=True, right_index=True, how='inner')

    group_col = 'signature_extraction_group'
    test_genes = np.intersect1d(np.unique(test_df.target), DNA_repair.Gene)
    test_genes = test_genes[[not bool(re.match("MOCK", gene)) for gene in test_genes]]
    column_subset = [re.search("^[A-Z0-9]+", key).group(0) in test_genes for key in target_df.keys()]
    target_df = target_df[target_df.keys()[column_subset]]
    target_df = pd.merge(sample_df[[group_col]], target_df,
                        left_index=True, right_index=True, how='inner')\
                    .replace({group_col:group_labels}).set_index(group_col)

    # Size
    group_size = target_df.groupby(group_col).size()

    hit_1 = (target_df[test_genes]>0).groupby(group_col).sum()
    hit_2 = (target_df[test_genes]>1).groupby(group_col).sum()

    hit_count = (target_df>0).groupby(group_col).sum()
    # Germline
    hit_germline = hit_count[test_genes+"_0"].rename(dict(zip(test_genes+"_0", test_genes)), axis=1)
    # Somatic
    hit_somatic = hit_count[test_genes+"_1"].rename(dict(zip(test_genes+"_1", test_genes)), axis=1)
    # LoH
    hit_loh = hit_count[test_genes+"_2"].rename(dict(zip(test_genes+"_2", test_genes)), axis=1)

    # Two hit
    hit_g2 = pd.DataFrame(np.array((target_df[test_genes+"_0"]>0))*np.array(target_df[test_genes]>1),
                index=target_df.index, columns=test_genes).groupby(group_col).sum()
    hit_s2 = pd.DataFrame(np.array((target_df[test_genes+"_1"]>0))*np.array(target_df[test_genes]>1),
                index=target_df.index, columns=test_genes).groupby(group_col).sum()
    hit_l2 = pd.DataFrame(np.array((target_df[test_genes+"_2"]>0))*np.array(target_df[test_genes]>1),
                index=target_df.index, columns=test_genes).groupby(group_col).sum()

    ncoh = len(group_size)
    ngene = len(test_genes)
    _ = plt.figure(constrained_layout=True, figsize=(50,20))

    color_maps = [cm.Greys, cm.Blues, cm.Reds, cm.Purples]

    for j,cohort in enumerate(group_size.index[:ncoh]):
        for i,gene in enumerate(test_genes[:ngene]):
            i_plot = j*(ngene+1) + i+1
            ax = plt.subplot(ncoh, ngene+1, i_plot)
            plt.sca(ax)

            single_hits = np.array([(hit_1/np.array(group_size)[:,None])[gene].loc[cohort],
                            (hit_germline/np.array(group_size)[:,None])[gene].loc[cohort],
                            (hit_somatic/np.array(group_size)[:,None])[gene].loc[cohort],
                            (hit_loh/np.array(group_size)[:,None])[gene].loc[cohort]])
            double_hits = np.array([(hit_2/np.array(group_size)[:,None])[gene].loc[cohort],
                            (hit_g2/np.array(group_size)[:,None])[gene].loc[cohort],
                            (hit_s2/np.array(group_size)[:,None])[gene].loc[cohort],
                            (hit_l2/np.array(group_size)[:,None])[gene].loc[cohort]])

            plt.bar(["T", "G", "S", "L"], double_hits, width=1.0,
                    color=[color_map(0.8) for color_map in color_maps]);
            plt.bar(["T", "G", "S", "L"], single_hits-double_hits, bottom=double_hits, width=1.0,
                    color=[color_map(0.5) for color_map in color_maps])

            plt.plot([-0.5,3.5], [0.02,0.02], '--k', alpha=0.7)
            plt.plot([-0.5,3.5], [0.1,0.1], '--k', alpha=0.7)
            plt.plot([-0.5,3.5], [0.5,0.5], '--k', alpha=0.7)

            if i==0: plt.ylabel(cohort.replace("_"," "),
                                fontsize=24, rotation=0, labelpad=20, horizontalalignment='right')
            else: ax.axes.yaxis.set_visible(False)
            if j==ncoh-1: ax.set_xlabel(fr"${gene}$", labelpad=20, fontsize=16, rotation=0)
            ax.set_xticks([])

            plt.ylim(0.005,0.999)
            plt.yscale('log')
            plt.subplots_adjust(hspace=0, wspace=0)
            plt.xlim(-0.5,3.5)

            if i+j==0:
                T_patch = patches.Patch(color=cm.Greys(0.8), label='All hits')
                G_patch = patches.Patch(color=cm.Blues(0.8), label='Germline')
                S_patch = patches.Patch(color=cm.Reds(0.8), label='Somatic')
                L_patch = patches.Patch(color=cm.Purples(0.8), label='LoH')
                plt.legend(loc="lower left", handles=[T_patch, G_patch, S_patch, L_patch],
                        bbox_to_anchor=(-0., 1.01), fontsize=24, ncol=4)

    ax2 = plt.subplot(1,ngene+1,ngene+1)
    plt.sca(ax2)
    plt.barh(group_size[:ncoh].index[::-1], np.array(group_size[:ncoh][::-1]), height=0.9, color='k', alpha=0.8);
    plt.xscale('log')
    plt.xlim(19,2500)
    plt.ylim(-0.5,ncoh-0.5)
    plt.xlabel("Cohort size", fontsize=24)
    ax2.set_yticks([])
    ax2.tick_params(axis='both', which='both', labelsize=16)



def getAssociationResults(results_dir, DNA_repair):
    # # Results

    results_zinb = pd.read_csv(f"{results_dir}/output/signature_target_assoc_nb.csv")
    results_log0 = pd.read_csv(f"{results_dir}/output/signature_target_assoc_logistic.csv")

    unique_groups = pd.DataFrame(results_zinb.groupby('group').size()).reset_index()
    unique_groups['tissue'] = unique_groups.group.map(lambda x: x.split("-")[0])

    unique_tissues = pd.DataFrame(unique_groups.groupby('tissue').size())

    unique_groups = pd.merge(unique_groups, unique_tissues.rename({0:'Count'}, axis=1).reset_index(),
            left_on='tissue', right_on='tissue', how='inner')

    group_labels = unique_groups[unique_groups['Count']==1].set_index('group')['tissue'].to_dict()

    results_combined = pd.merge(results_zinb, results_log0,
                                on=('target', 'signature', 'group'),
                                suffixes=("_zinb", '_log0')).replace({'group':group_labels})

    results_combined['zscore'] = results_combined.target_means_alt_zinb/np.sqrt(results_combined.target_covs_alt_zinb)
    results_combined['log10_pvalue_zinb'] = -np.log10(results_combined.resample_pvalue)
    results_combined['log10_pvalue_log0'] = -np.log10(results_combined.wilks_pvalue_log0)

    results_combined = pd.merge(results_zinb, results_log0,
                                on=('target', 'signature', 'group'),
                                suffixes=("_zinb", '_log0'))#.replace({'group':group_labels})

    results_combined['zscore'] = results_combined.target_means_alt_zinb/np.sqrt(results_combined.target_covs_alt_zinb)
    results_combined['log10_pvalue_zinb'] = -np.log10(results_combined.resample_pvalue)
    results_combined['log10_pvalue_log0'] = -np.log10(results_combined.wilks_pvalue_log0)
    print(len(results_combined))

    DNArepair = pd.read_csv("/re_gecip/shared_allGeCIPs/pancancer_signatures/data/human-dna-repair-genes.tsv", sep="\t")
    results_combined = pd.merge(results_combined, DNArepair, how='left', left_on='target', right_on='Gene')
    print(len(results_combined))

    results_combined = results_combined[results_combined.target.map(lambda x: x in np.array(DNA_repair.Gene))|
                    mock_match(results_combined.target)]
    print(len(results_combined))

    results_combined['exp_means_alt_zinb'] = np.exp(results_combined['target_means_alt_zinb'])
    results_combined["gain"] = np.exp(results_combined.target_means_alt_zinb)

    sample_df = pd.read_csv(f"{results_dir}/input/samples.tsv", sep="\t").set_index('sample_id')
    target_df = pd.read_csv(f"{results_dir}/input/targets.tsv", sep="\t").set_index('sample_id')
    test_df = pd.read_csv(f"{results_dir}/input/tests.tsv", sep="\t")

    group_col = 'group'
    test_genes = np.intersect1d(np.unique(test_df.target), DNA_repair.Gene)
    test_genes = test_genes[[not bool(re.match("MOCK", gene)) for gene in test_genes]]
    column_subset = [re.search("^[A-Z0-9]+", key).group(0) in test_genes for key in target_df.keys()]
    target_df = target_df[target_df.keys()[column_subset]]
    target_df = pd.merge(sample_df[[group_col]], target_df,
                        left_index=True, right_index=True, how='inner')\
                    .set_index(group_col)

    agg_df = pd.melt(target_df.groupby('group').sum().reset_index(), id_vars='group', var_name='target')
    agg_df = agg_df[(agg_df.value>0)&(agg_df.value<5)]
    agg_df['target'] = agg_df.target.map(lambda x: x.split("_")[0])
    agg_df.drop_duplicates(['group', 'target'], inplace=True)

    results_combined = pd.merge(results_combined, agg_df, how='outer', on=['group', 'target'], indicator=True)
    results_combined = results_combined[results_combined._merge.map(lambda x: x in ['both', 'left_only'])]
    print(len(results_combined))

    # Replace Connective with Sarcoma
    results_combined.group = results_combined.group.str.replace("Connective","Sarcoma")

    return results_combined



def plotTwoHitAssoc(results_combined,
                    highlight_set = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,10,2),
                                                ('POLE','Uterus-AdenoCA',1e110,8,2.2)],
                                'signature,target':[('SBS18','MUTYH',1e-10,2,1.5)]},
                    highlight_hit = [('POLG', 'DBS4', 'ColoRect-AdenoCA', -0.1, 2),
                                       ('POLG', 'SBS93', 'ColoRect-AdenoCA', -0.1, 5)],
                    xlim=[[-4,4],[-1,4]]):

    # Z-score - Q plot
    fig, axes = plt.subplots(1,2,figsize=(16,5))

    #passed_indices = BenjiminiHochberg(np.array(results_combined['resample_pvalue']), 0.001)
    p_threshold = BH_threshold(results_combined['resample_pvalue'], 0.01)

    for i in range(2):
        plt.sca(axes[i])

        # zscore = results_combined.target_means_alt/np.sqrt(results_combined.target_covs_alt)
        plt.scatter(results_combined.target_means_alt_zinb,
                    -np.log10(results_combined.resample_pvalue),
                    c=-np.log10(results_combined.wilks_pvalue_log0),
                    s=4, cmap='viridis_r', vmax=15)

        ylim = (0, axes[i].axes.get_ylim()[1])

        plt.plot(xlim[0],[-np.log10(p_threshold),-np.log10(p_threshold)], '--r')
        plt.plot([0,0],[-np.log10(p_threshold),ylim[1]], '--r')

        plt.xlim(xlim[0])#-np.max(np.abs(zscore)), np.max(np.abs(zscore)))
        plt.gca().set_ylim(bottom=0, top=ylim[1])

        cbar = plt.colorbar()
        cbar.set_label(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{logistic})$")

        plt.xlabel(r"$\beta$")
        plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{ZINB-dCRT})$")

    plt.sca(axes[1])
    plt.xlim(xlim[1])
    plt.ylim(p_threshold, 22)

    plt.sca(axes[0])

    for key,combos in highlight_set.items():
        for i,set_row in enumerate(combos):
            rset = results_combined[(results_combined[key.split(",")[0]]==set_row[0])&\
                                (results_combined[key.split(",")[1]]==set_row[1])&\
                                (results_combined.resample_pvalue<set_row[2])&\
                                (results_combined.signature!='CN1')]
            posy = np.mean(-np.log10(rset.resample_pvalue))+set_row[3]
            if np.max(-np.log10(rset.resample_pvalue))>20: ax = axes[0]
            else: ax = axes[1]
            if set_row[4]<0: ha='right'
            else: ha='left'
            plt.sca(ax)
            plt.text(set_row[4], posy, f"{set_row[0]}, \n{set_row[1]}", ha=ha, va='center')
            for index, row in rset.iterrows():
                plt.plot([set_row[4],row.target_means_alt_zinb],
                        [posy,-np.log10(row.resample_pvalue)], ':k', alpha=0.3)

    # Highlight hits
    for i, set_row in enumerate(highlight_hit):
        row = results_combined[(results_combined.target==set_row[0])&\
                                (results_combined.signature==set_row[1])&\
                                (results_combined.group==set_row[2])].iloc[0]
        pos = (row.target_means_alt_zinb, -np.log10(row.resample_pvalue))
        if pos[1]>20: ax = axes[0]
        else: ax = axes[1]
        plt.sca(ax)
        text_pos = (pos[0]+set_row[3], pos[1]+set_row[4])
        plt.text(text_pos[0], text_pos[1], f"{row.signature} (${row.target}$), \n{row.group}",
                ha='left' if set_row[3]>0 else 'right', va='bottom')
        plt.plot([pos[0],text_pos[0]],[pos[1], text_pos[1]], ':k', alpha=0.3)



def plotTwoHitAssocZoom(results_combined, p_threshold=0.05,
                        highlight_set = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,10,2),
                                                    ('POLE','Uterus-AdenoCA',1e110,8,2.2)],
                                    'signature,target':[('SBS18','MUTYH',1e-10,2,1.5)]},
                        highlight_hit = [('POLG', 'DBS4', 'ColoRect-AdenoCA', -0.1, 2),
                                        ('POLG', 'SBS93', 'ColoRect-AdenoCA', -0.1, 5)]):

    # Z-score - Q plot
    fig1, ax1 = plt.subplots(1,1,figsize=(8,5))
    fig2, ax2 = plt.subplots(1,1,figsize=(8,5))
    axes = [ax1,ax2]

    xlim=[-4,8]

    for i in range(2):
        plt.sca(axes[i])

        # zscore = results_combined.target_means_alt/np.sqrt(results_combined.target_covs_alt)
        plt.scatter(results_combined.target_means_alt_zinb,
                    -np.log10(results_combined.resample_pvalue),
                    c=-np.log10(results_combined.wilks_pvalue_log0),
                    s=4, cmap='viridis_r', vmax=15)

        ylim = (0, axes[i].axes.get_ylim()[1])

        plt.plot(xlim,[-np.log10(p_threshold),-np.log10(p_threshold)], '--r')
        plt.plot([0,0],[-np.log10(p_threshold),ylim[1]], '--r')

        plt.xlim(xlim)#-np.max(np.abs(zscore)), np.max(np.abs(zscore)))
        plt.gca().set_ylim(bottom=0, top=ylim[1])

        cbar = plt.colorbar()
        cbar.set_label(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{logistic})$")

        plt.xlabel(r"$\beta$")
        plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{ZINB-dCRT})$")

    plt.sca(axes[1])
    plt.xlim(-1,5)
    plt.ylim(p_threshold, 22)

    plt.sca(axes[0])

    for key,combos in highlight_set.items():
        for i,set_row in enumerate(combos):
            rset = results_combined[(results_combined[key.split(",")[0]]==set_row[0])&\
                                (results_combined[key.split(",")[1]]==set_row[1])&\
                                (results_combined.resample_pvalue<set_row[2])&\
                                (results_combined.signature!='CN1')]
            posy = np.mean(-np.log10(rset.resample_pvalue))+set_row[3]
            if np.max(-np.log10(rset.resample_pvalue))>20: ax = axes[0]
            else: ax = axes[1]
            if set_row[4]<0: ha='right'
            else: ha='left'
            plt.sca(ax)
            plt.text(set_row[4], posy, f"{set_row[0]},\n{set_row[1]}", ha=ha, va='center')
            for index, row in rset.iterrows():
                plt.plot([set_row[4],row.target_means_alt_zinb],
                        [posy,-np.log10(row.resample_pvalue)], ':k', alpha=0.3)


    rset = results_combined[(results_combined.resample_pvalue<1e-15)&\
                            (results_combined.signature=='CN1')]
    posy = np.mean(-np.log10(rset.resample_pvalue))+10
    if np.max(-np.log10(rset.resample_pvalue))>20: ax = axes[0]
    else: ax = axes[1]
    if -2<0: ha='right'
    else: ha='left'
    plt.sca(ax)
    plt.text(-2, posy, f"CN1", ha=ha, va='center')
    for index, row in rset.iterrows():
        plt.plot([-2,row.target_means_alt_zinb],
                [posy,-np.log10(row.resample_pvalue)], ':k', alpha=0.3)

    # Highlight hits
    for i, set_row in enumerate(highlight_hit):
        row = results_combined[(results_combined.target==set_row[0])&\
                                (results_combined.signature==set_row[1])&\
                                (results_combined.group==set_row[2])].iloc[0]
        pos = (row.target_means_alt_zinb, -np.log10(row.resample_pvalue))
        if pos[1]>20: ax = axes[0]
        else: ax = axes[1]
        plt.sca(ax)
        text_pos = (pos[0]+set_row[3], pos[1]+set_row[4])
        plt.text(text_pos[0], text_pos[1], f"{row.signature}, (${row.target}$) ,\n{row.group}",
                ha='left' if set_row[3]>0 else 'right', va='bottom')
        plt.plot([pos[0],text_pos[0]],[pos[1], text_pos[1]], ':k', alpha=0.3)

    return [fig1,fig2], axes


def plotZHist(results_combined):

    plt.figure(figsize=(8,5))

    Z = results_combined['target_means_alt_zinb']/np.sqrt(results_combined['target_covs_alt_zinb'])

    hist_kwargs = {"bins":np.linspace(-20,20,100), "density":True, "histtype":'step',
                "linewidth":2}
    plt.hist(Z[results_combined['signature']=="MUTmock"],
            label="Mock signature", color=default_colours[0], **hist_kwargs);
    plt.plot([np.nanmedian(Z[results_combined['signature']=="MUTmock"]),]*2, [0,0.5],
            c=default_colours[0], linestyle=":")

    plt.hist(Z[(results_combined['resample_pvalue']<1e-6)&(results_combined['resample_pvalue']>0)],
            label=f"$P-\mathrm{{value}} < 10^{{-6}}$", color=default_colours[1], **hist_kwargs)
    plt.plot([np.nanmedian(Z[(results_combined['resample_pvalue']<1e-6)&(results_combined['resample_pvalue']>0)]),]*2,
            [0,0.5], c=default_colours[1], linestyle=":")

    plt.hist(Z[(results_combined['resample_pvalue']>0.5)],
            color=default_colours[2], label=f"$P-\mathrm{{value}} > 0.1$", **hist_kwargs);
    plt.plot([np.nanmedian(Z[(results_combined['resample_pvalue']>0.5)]),]*2,
            [0,0.5], c=default_colours[2], linestyle=":")

    plt.plot(hist_kwargs['bins'], scipy.stats.norm.pdf(hist_kwargs['bins']), '--r')

    plt.legend(loc='upper left')
    plt.xlim(hist_kwargs['bins'][[0,-1]])
    plt.xlabel(r"Z")
    plt.ylabel("PDF")


def plotZBetaHist(results_combined, p_threshold=0.05):

    fig, axes = plt.subplots(1,2,figsize=(16,5))

    plt.sca(axes[0])

    hist_kwargs = {"bins":np.linspace(-2,2,100), "density":True, "histtype":'step',
                "linewidth":2}
    plt.hist(results_combined['target_means_alt_zinb'][results_combined['signature']=="MUTmock"],
            label="Mock signature", color=default_colours[0], **hist_kwargs);
    plt.plot([np.nanmedian(results_combined['target_means_alt_zinb'][
                    results_combined['signature']=="SBSmock"]),]*2, [0,3.5],
            c=default_colours[0], linestyle=":")

    plt.hist(results_combined['target_means_alt_zinb'][(results_combined['resample_pvalue']<p_threshold)&\
                                                    (results_combined['resample_pvalue']>0)],
            label=f"$P-\mathrm{{value}} < {p_threshold/10**np.floor(np.log10(p_threshold)):.2f}e^{{{int(np.floor(np.log10(p_threshold)))}}}$",
            color=default_colours[1], **hist_kwargs)
    sig_med = np.nanmedian(results_combined['target_means_alt_zinb'][
                (results_combined['resample_pvalue']<p_threshold)&(results_combined['resample_pvalue']>0)])
    plt.plot([sig_med,]*2,
            [0,3.5], c=default_colours[1], linestyle=":")
    print(sig_med)

    plt.hist(results_combined['target_means_alt_zinb'][(results_combined['resample_pvalue']>0.5)],
            color=default_colours[2], label=f"$P-\mathrm{{value}} > 0.1$", **hist_kwargs);
    plt.plot([np.nanmedian(results_combined['target_means_alt_zinb'][
                    (results_combined['resample_pvalue']>0.1)]),]*2,
            [0,3.5], c=default_colours[2], linestyle=":")

    plt.legend(loc='upper left')
    plt.xlim(hist_kwargs['bins'][[0,-1]])
    plt.xlabel(r"$\beta$")
    plt.ylabel("PDF")
    plt.ylim(0,3.5)

    plt.sca(axes[1])

    Z = results_combined['target_means_alt_zinb']/np.sqrt(results_combined['target_covs_alt_zinb'])

    hist_kwargs = {"bins":np.linspace(-20,20,100), "density":True, "histtype":'step',
                "linewidth":2}
    plt.hist(Z[results_combined['signature']=="MUTmock"],
            label="Mock Sig", color=default_colours[0], **hist_kwargs);
    plt.plot([np.nanmedian(Z[results_combined['signature']=="MUTmock"]),]*2, [0,0.5],
            c=default_colours[0], linestyle=":")

    sig_med = np.nanmedian(Z[(results_combined['resample_pvalue']<1e-6)&(results_combined['resample_pvalue']>0)])
    plt.hist(Z[(results_combined['resample_pvalue']<p_threshold)&(results_combined['resample_pvalue']>0)],
            label=f"$P-\mathrm{{value}} < 10^{{-6}}$", color=default_colours[1], **hist_kwargs)
    plt.plot([sig_med,]*2,
            [0,0.5], c=default_colours[1], linestyle=":")
    print(sig_med)

    plt.hist(Z[(results_combined['resample_pvalue']>0.5)],
            color=default_colours[2], label=f"$P-\mathrm{{value}} > 0.1$", **hist_kwargs);
    plt.plot([np.nanmedian(Z[(results_combined['resample_pvalue']>0.5)]),]*2,
            [0,0.5], c=default_colours[2], linestyle=":")

    plt.plot(hist_kwargs['bins'], scipy.stats.norm.pdf(hist_kwargs['bins']), '--r')

    #plt.legend(loc='upper left')
    plt.xlim(hist_kwargs['bins'][[0,-1]])
    plt.xlabel(r"Z")
    plt.ylabel("PDF")

    plt.ylim(0,0.5)


def plotCohortHits(results_combined, DNA_repair):
    # ## Cohort hits

    method = "zinb"

    fig_type='ivw'
    if fig_type=="most_hit":
        cohort_hit_table = pd.DataFrame(results_combined
                                        .sort_values('gain', ascending=False)\
                                        .drop_duplicates(['group', 'signature'], keep='first')\
                                        .groupby(['target', 'signature']).size())
        # Convert to grid
        cohort_hit_table = cohort_hit_table.pivot_table(index=['target'],
                                                        columns=['signature'], values=[0]).fillna(0).astype(int)
    elif fig_type=="ivw":
        results_test = results_combined.copy()[
            ['signature', 'target', f'target_means_alt_{method}', f'target_covs_alt_{method}']
        ]
        results_test['target_ivw'] = results_test[f'target_means_alt_{method}']/results_test[f'target_covs_alt_{method}']
        results_test['target_iv'] = 1/results_test[f'target_covs_alt_{method}']
        results_test = results_test.groupby(['signature','target']).sum().copy()
        results_test['target_z'] = results_test['target_ivw']/np.sqrt(results_test['target_iv'])
        results_test['target_mu'] = results_test['target_ivw']/results_test['target_iv']
        # Convert to grid
        cohort_hit_table = results_test.pivot_table(index=['target'],
                                                        columns=['signature'], values=['target_z'])

    # Reorder columns
    cohort_hit_table.columns = cohort_hit_table.columns.get_level_values(1)
    cohort_hit_table = cohort_hit_table[cohort_hit_table.columns[
        orderSignatures(cohort_hit_table.columns)]
                                    ]

    # Reorder rows
    DNA_repair_subset = pd.merge(pd.DataFrame(cohort_hit_table.index), DNA_repair,
                                left_on='target', right_on='Gene', how='inner')\
                        .replace({"Type":{"Chromatin Structure and Modification":"ChrMod",
                                    "Homologous recombination":"HRD",
                                    "DNA polymerases ":"POL",
                                    "Fanconi anemia":"FA",
                                    "Editing and processing nucleases":"Nuc",
                                    "Direct reversal of damage":"DRD"}})\
                        .sort_values("target").sort_values('Type', ascending=False)[::-1]
    cohort_hit_table = cohort_hit_table.loc[DNA_repair_subset.Gene]

    fig, ax = plt.subplots(1,1,figsize=(16,16*cohort_hit_table.shape[0]/cohort_hit_table.shape[1]))

    breaks = np.linspace(0,1,11)
    vmin=0; vmax=25#3
    im = plt.pcolor(cohort_hit_table, cmap='inferno_r', vmin=vmin, vmax=vmax)

    plt.scatter(np.arange(len(cohort_hit_table.columns))+0.5,
                np.nanargmax(np.array(cohort_hit_table), axis=0)+0.5,
                c='c', s=10)

    fs=16
    ax.set_yticks(np.arange(cohort_hit_table.shape[0]), minor=True)
    ax.set_xticks(np.arange(cohort_hit_table.shape[1]), minor=True)
    ax.grid(True, which='minor', axis='both', linestyle='-', color='k')

    ax.set_xticks(np.arange(cohort_hit_table.shape[1])+0.5, minor=False)
    ax.set_xticklabels(cohort_hit_table.columns, rotation=90, fontsize=fs, minor=False);
    ax.set_yticks(np.arange(cohort_hit_table.shape[0])+0.5, minor=False)
    ax.set_yticklabels([f"${gene}$" for gene in cohort_hit_table.index], rotation=0, fontsize=fs, minor=False);

    i=0
    for row in pd.DataFrame(DNA_repair_subset.groupby("Type").size()).iterrows():
        plt.text(-7,i+row[1]/2,row[0], ha='center', va='center', fontsize=fs*1.1)
        i+=row[1]
        plt.plot([-10,cohort_hit_table.shape[1]],[i,i],'-k', alpha=1, clip_on=False)
    plt.plot([-10,cohort_hit_table.shape[1]],[0,0],'-k', alpha=1, clip_on=False)
    plt.xlim(0,cohort_hit_table.shape[1])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='2%', pad='2%')
    # norm = mpl.colors.Normalize(vmin=0,vmax=10)
    cbar = fig.colorbar(im, cax=cax, orientation='vertical',
                        ticks=np.linspace(vmin,vmax,5))
    cbar.set_label("Pan-cancer Z-score", fontsize=fs)
    # cbar.ax.set_yticks()


def printMockStats(results_combined):
    # # Mock Results

    subset = (results_combined.signature.map(lambda x: x[0]!="T"))&\
            (results_combined.resample_pvalue>0)&\
            (results_combined._merge=='left_only')
    mock_gene = np.vectorize(lambda x:bool(re.match("MOCK[A-Z]+[0-9]+", x)))
    mock_size = np.vectorize(lambda x:bool(re.match("MOCK[0-9]+", x)))
    subset = subset&mock_size(results_combined.target)

    threshold = 0.05/np.sum(subset)
    threshold1 = BH_threshold(results_combined[subset].wilks_pvalue_zinb)
    threshold2 = BH_threshold(results_combined[subset].resample_pvalue)
    FDRw = 100*(np.sum(results_combined[subset]['wilks_pvalue_zinb']<threshold1))/np.sum(subset)
    FDRc = 100*(np.sum(results_combined[subset]['resample_pvalue']<threshold2))/np.sum(subset)

    subset = (~mock_match(results_combined.target))&\
            (~mock_match(results_combined.signature))&\
            (results_combined.signature.map(lambda x: x[0]!="T"))&\
            (results_combined.resample_pvalue>0)&\
            (results_combined.signature!='CN1')&\
            (results_combined._merge=='left_only')

    threshold = 0.05/np.sum(subset)
    threshold1 = BH_threshold(results_combined[subset].wilks_pvalue_zinb)
    threshold2 = BH_threshold(results_combined[subset].resample_pvalue)

    DRw = 100*(np.sum(results_combined[subset]['wilks_pvalue_zinb']<threshold1))/np.sum(subset)
    DRc = 100*(np.sum(results_combined[subset]['resample_pvalue']<threshold2))/np.sum(subset)

    print("False positives")
    print(f"Wilks: {FDRw:.2f}, "\
        f"Resample: {FDRc:.2f}")
    print("Positives")
    print(f"Wilks: {DRw:.2f}, "\
        f"Resample: {DRc:.2f}")
    print("True positives")
    print(f"Wilks: {DRw-FDRw:.2f}?, "\
        f"Resample: {DRc-FDRc:.2f}?")




def plotMockPP(results_combined):

    subset = (results_combined.resample_pvalue>0)&\
            (results_combined.model_alt_zinb.isin(['negbin', 'glm.binom', 'binomial']))&\
            (results_combined.group.map(lambda x: x.split("-")[0]!='Haem'))

    subset = (results_combined.signature.map(lambda x: x[0]!="T"))&\
            (results_combined.resample_pvalue>0)&\
            (results_combined._merge=='left_only')

    mock_gene = np.vectorize(lambda x:bool(re.match("MOCK[A-Z]+[0-9]+", x)))
    mock_size = np.vectorize(lambda x:bool(re.match("MOCK[0-9]+", x)))
    subsets = [subset&(results_combined.signature=="MUTmock"),
            subset&mock_size(results_combined.target),
            subset&mock_gene(results_combined.target)]
    titles = ["Signature mock", "Target resample", "Target perturbation"]

    fig, axes = plt.subplots(1,3,figsize=(15,5), sharey=True)

    ymax = 10

    for i, subset in enumerate(subsets):

        plt.sca(axes[i])

        n_test = np.sum(subset)
        k = np.arange(1,n_test+1)

        pv_expected = np.linspace(0,1-1/n_test,n_test)+1/(2*n_test)
        pv_expected = np.arange(1/n_test,1+1e-10,1/n_test)
        percentiles = scipy.stats.beta.ppf(np.array([0.05,0.5,0.95])[:,None], k[None,:],n_test+1-k[None,:])


        Y = np.sort(-np.log10(results_combined[subset]['wilks_pvalue_zinb']))
        Y[Y>ymax]=ymax+0.01
        plt.scatter(-np.log10(percentiles[1][::-1])[Y<ymax], Y[Y<ymax],
                    s=20, label=f"Wilks",
                    c=default_colours[3], marker=".")#, marker=markers[i])
        plt.scatter(-np.log10(percentiles[1][::-1])[Y>ymax], Y[Y>ymax],

                    s=20, c=default_colours[3], marker="^")

        plt.scatter(-np.log10(percentiles[1][::-1]),
                    np.sort(-np.log10(results_combined[subset]['resample_pvalue'])),
                    s=20, label=f"dCRT",
                    c=default_colours[0], marker=".")#, marker=markers[i])

        k = np.arange(1,n_test+1)
        percentiles = scipy.stats.beta.ppf(np.array([0.05,0.5,0.95])[:,None],
                                        k[None,:],
                                        n_test+1-k[None,:])
        plt.plot(-np.log10(percentiles[1][::-1]),-np.log10(percentiles[1][::-1]), c='k')
        plt.fill_between(-np.log10(percentiles[1][::-1]),
                        -np.log10(percentiles[0][::-1]),
                        -np.log10(percentiles[2][::-1]), alpha=0.3, color='k')

        plt.gca().set_xlim(left=0)
        plt.gca().set_ylim(bottom=0)
        plt.gca().set_ylim(top=ymax+0.1)

        plt.title(titles[i], fontsize=20, pad=10)

        if i==1:
            plt.xlabel(r"$-\log_{10}(P_\mathrm{expected})$")
        if i==0:
            plt.legend(loc="upper left")
            plt.ylabel(r"$-\log_{10}(P_\mathrm{observed})$")

    plt.subplots_adjust(wspace=0.01)


def correlationGrid(tgt_df, target_subset, p_threshold=0.05):

    # Generate target grid
    unique_targets, inverse_target = np.unique(tgt_df.target, return_inverse=True)
    unique_dept, inverse_dept = np.unique(tgt_df.dependent, return_inverse=True)

    target_grid = np.zeros((len(unique_targets), len(unique_dept)))

    target_grid[inverse_target, inverse_dept] = np.where(tgt_df.wilks_pvalue<p_threshold,
                                                         tgt_df.target_means_alt, 0)
    target_grid = pd.DataFrame(target_grid, index=unique_targets, columns=unique_dept)

    # Insert missing columns
    target_grid[np.setdiff1d(target_subset, target_grid.keys())] = np.nan
    # Insert missing rows
    target_grid = target_grid.T
    target_grid[np.setdiff1d(target_subset, target_grid.keys())] = np.nan
    target_grid = target_grid.T

    return target_grid

def plotTargetConfounding(tgt_df, target_subset, ordered_genes, p_threshold=0.05,
                          groups = ['Colorectal', 'Uterus', 'Breast']):

    fig, axes = plt.subplots(1,3,figsize=(50,17))
    for i,group in enumerate(groups):
        plt.sca(axes[i])
        target_grid = correlationGrid(tgt_df[tgt_df.group==group],
                                    target_subset,
                                    0.01)[ordered_genes].loc[ordered_genes]

        im = plt.pcolor(target_grid, vmin=-5, vmax=5, cmap='bwr')

        plt.title(f"Target-target associations in {group}", fontsize=36)

        axes[i].set_xticks(np.arange(len(ordered_genes))+0.5, minor=False)
        axes[i].set_yticks(np.arange(len(ordered_genes))+0.5, minor=False)

        axes[i].set_xticklabels([f"${gene}$" for gene in ordered_genes], rotation=90, fontsize=24, minor=False);
        axes[i].set_yticklabels([f"${gene}$" for gene in ordered_genes], fontsize=24, minor=False);


    fig.subplots_adjust(right=0.95, wspace=0.12, hspace=0.12)

    # add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
    cb_ax = fig.add_axes([0.96, 0.15, 0.01, 0.7])

    cbar = fig.colorbar(im, cax=cb_ax)
    cbar.set_label(r"$\beta$", fontsize=40)
    cbar.ax.tick_params(labelsize=30)


def getCovariateAssociations(results_combined):
    # # Covariate associations

    covariates = [re.search("([A-Z0-9_]+)_means_null_zinb",key,re.IGNORECASE).group(1) \
                for key in results_combined.keys() \
                if bool(re.search("[A-Z0-9_]+_means_null_zinb",key,re.IGNORECASE))]

    results_medians = results_combined[['group','signature','target']\
                                            +[f"{cov}_means_null_zinb" for cov in covariates[1:]]\
                                            +[f"{cov}_covs_null_zinb" for cov in covariates[1:]]\
                                            +[f"{cov}_means_null_log0" for cov in covariates[1:]]\
                                            +[f"{cov}_covs_null_log0" for cov in covariates[1:]]\
                                            +[f"Xsd_{cov}_zinb" for cov in covariates[1:]]]\
                    .groupby(['group','signature']).median()

    rename = dict(zip([f"{cov}_means_null_zinb" for cov in covariates]+\
                    [f"{cov}_zero_means_null" for cov in covariates]+\
                    [f"{cov}_covs_null_zinb" for cov in covariates]+\
                    [f"{cov}_means_null_log0" for cov in covariates]+\
                    [f"{cov}_covs_null_log0" for cov in covariates]+\
                    [f"Xsd_{cov}_zinb" for cov in covariates],
                    [f"mean-{cov}" for cov in covariates]+\
                    [f"zeromean-{cov}" for cov in covariates]+\
                    [f"cov-{cov}" for cov in covariates]+\
                    [f"meanlog0-{cov}" for cov in covariates]+\
                    [f"covlog0-{cov}" for cov in covariates]+\
                    [f"sd-{cov}" for cov in covariates]))
    results_medians.rename(rename, inplace=True, axis=1)

    results_medians = pd.wide_to_long(results_medians.reset_index(),
                    stubnames=["mean","cov","meanlog0","covlog0","zeromean","sd"],
                    i=['group','signature'],j='covariate',
                    sep="-", suffix=r'\w+').reset_index()

    # Inverse variance weighted mean
    results_medians['ivw'] = results_medians['mean']/results_medians['cov']
    results_medians['iv'] = 1/results_medians['cov']

    results_medians['ivwlog0'] = results_medians['meanlog0']/results_medians['covlog0']
    results_medians['ivlog0'] = 1/results_medians['covlog0']

    # Correct for standard-deviation normalistion
    results_medians['mean_sd'] = results_medians['mean']/results_medians['sd']
    results_medians['cov_sd'] = results_medians['cov']/results_medians['sd']**2
    results_medians['ivw_sd'] = results_medians['mean_sd']/results_medians['cov_sd']
    results_medians['iv_sd'] = 1/results_medians['cov_sd']

    results_medians['meanlog0_sd'] = results_medians['meanlog0']/results_medians['sd']
    results_medians['covlog0_sd'] = results_medians['covlog0']/results_medians['sd']**2
    results_medians['ivwlog0_sd'] = results_medians['meanlog0_sd']/results_medians['covlog0_sd']
    results_medians['ivlog0_sd'] = 1/results_medians['covlog0_sd']

    results_medians['tissue'] = results_medians.group.map(lambda x: x.split("-")[0])

    results_medians = results_medians.drop('group', axis=1)\
                                    .groupby(['tissue', 'signature', 'covariate']).sum().reset_index()

    sig_pos = dict(zip(
        np.unique(results_medians.signature)[orderSignatures(np.unique(results_medians.signature))],
        np.arange(len(np.unique(results_medians.signature))),))

    signatures = np.unique(results_medians.signature)[
        orderSignatures(np.unique(results_medians.signature))
    ]

    results_medians['z'] = results_medians['mean']/np.sqrt(results_medians['cov'])
    results_medians['z_log0'] = results_medians['meanlog0']/np.sqrt(results_medians['covlog0'])
    results_medians['log_odds'] = results_medians['mean']/results_medians['sd']
    results_medians['log_odds_zero'] = results_medians['zeromean']/results_medians['sd']
    results_medians['log_odds_log0'] = results_medians['meanlog0']/results_medians['sd']
    results_medians['logpv'] = scipy.stats.chi2.logsf(results_medians.z**2, df=1)
    results_medians['logpv_log0'] = scipy.stats.chi2.logsf(results_medians.z_log0**2, df=1)

    results_medians['z'] = results_medians['ivw']/np.sqrt(results_medians['iv'])
    results_medians['mu'] = results_medians['ivw']/results_medians['iv']
    results_medians['z_sd'] = results_medians['ivw_sd']/np.sqrt(results_medians['iv_sd'])
    results_medians['mu_sd'] = results_medians['ivw_sd']/results_medians['iv_sd']

    results_medians['zlog0'] = results_medians['ivwlog0']/np.sqrt(results_medians['ivlog0'])
    results_medians['mulog0'] = results_medians['ivwlog0']/results_medians['ivlog0']
    results_medians['zlog0_sd'] = results_medians['ivwlog0_sd']/np.sqrt(results_medians['ivlog0_sd'])
    results_medians['mulog0_sd'] = results_medians['ivwlog0_sd']/results_medians['ivlog0_sd']

    results_medians['zlog0'][np.abs(results_medians['zlog0'])>20]=20*np.sign(results_medians['zlog0'])

    return results_medians


def plotCovariateAssociations(results_medians, rename_sigs,
                              covs_ = np.array(['is_female', 'log_age', 'log_n_loh', 'pc1', 'pc2', 'pc3']),
                              cov_labels = {"log_age":f"$\log(\mathrm{{age}})$", "is_female":f"Sex (is female)",
                                            "pc1":f"PC1", "pc2":f"PC2", "pc3":f"PC3", "log_n_loh":f"$\log(N_\mathrm{{LoH}})$"}):

    z_var = 'zlog0'
    mu_var = 'mulog0_sd'

    z_min = np.sqrt(scipy.stats.chi2.isf(0.01, df=1))
    logpv_threshold = np.log(BH_threshold(np.exp(results_medians.logpv), 0.01))
    z_threshold = np.sqrt(scipy.stats.chi2.isf(np.exp(logpv_threshold), df=1))

    results_3d_grid = {}
    for variable in [z_var, mu_var]:

        results_pivoted = pd.pivot_table(results_medians,
                        values=variable, index=['tissue'], columns=['signature', 'covariate'])

        sigs_ = np.unique(np.array(list(results_pivoted.columns))[:,0])
        covs_ = np.unique(np.array(list(results_pivoted.columns))[:,1])
        tissues_ = np.array(list(results_pivoted.index))

        results_3d_grid[variable] = np.zeros((len(results_pivoted), len(sigs_), len(covs_)))
        for isig,sig in enumerate(sigs_):
            for icov,cov in enumerate(covs_):
                try: results_3d_grid[variable][:,isig,icov] = results_pivoted[sig][cov]
                except: pass

    z_min=0
    sig_subset = np.sum(np.abs(results_3d_grid[z_var])>z_min, axis=(0,2))
    sigs_ = sigs_[sig_subset>0]

    results_3d_grid[z_var] = results_3d_grid[z_var][:,sig_subset>0]
    results_3d_grid[mu_var] = results_3d_grid[mu_var][:,sig_subset>0]

    sigs_order_ = np.array(orderSignatures(sigs_))

    tissue_subset = np.nansum(np.abs(results_3d_grid[z_var][:,:,:])>z_min, axis=(1,2))>0
    tissue_subsets = [tissue_subset & (np.nansum(np.abs(results_3d_grid[z_var][:,:,icov])>0, axis=1)>0)\
                            for icov in range(len(covs_))]

    fig, axes = plt.subplots(len(covs_),1,figsize=(25,1.3*25*np.sum([np.sum(tss) for tss in tissue_subsets])\
                                                /results_3d_grid[z_var].shape[1]),
                            gridspec_kw={'height_ratios':[np.sum(tss) for tss in tissue_subsets]})
    plt.subplots_adjust(hspace=0.4)
    fs=16

    for icov, cov in enumerate(covs_):
        ax = axes[icov]
        plt.sca(ax)

        # Interquartile range
        dynamic_rng = np.diff(np.nanpercentile(np.array(results_3d_grid[mu_var][:,:,icov]).flatten(),
                                            np.array([25,75])))[0]

        tissue_subset = tissue_subsets[icov]
        z_table = results_3d_grid[z_var][:,:,icov][:,sigs_order_][tissue_subset]
        mu_table = results_3d_grid[mu_var][:,:,icov][:,sigs_order_][tissue_subset]

        tried_grid = 1-np.isnan(z_table).astype(int)
        tried_grid[np.abs(np.array(z_table))>z_threshold] = 5
        plt.pcolor(tried_grid, cmap='bone_r', vmin=0, vmax=10)

        size_factor=20
        XX,YY = np.meshgrid(np.arange(z_table.shape[0]), np.arange(z_table.shape[1]))
        im = plt.scatter(YY.flatten()+0.5, XX.flatten()+0.5,
                        c=mu_table.T.flatten(),
                        s=np.min(np.vstack((np.abs(z_table.T.flatten())*size_factor,
                                            np.ones(len(z_table.flatten()))*120)), axis=0),
                        cmap=cm.bwr, vmin=-dynamic_rng*2, vmax=dynamic_rng*2)

        ax.set_yticks(np.arange(np.sum(tissue_subset)), minor=True)
        ax.set_xticks(np.arange(len(sigs_)), minor=True)
        ax.grid(True, which='minor', axis='both', linestyle='-', color='k')

        ax.set_xticks(np.arange(len(sigs_))+0.5, minor=False)
        ax.set_xticklabels(sigs_[sigs_order_], rotation=90, fontsize=fs, minor=False);
        ax.set_yticks(np.arange(np.sum(tissue_subset))+0.5, minor=False)
        ax.set_yticklabels(tissues_[tissue_subset], rotation=0, fontsize=fs, minor=False);

        plt.text(0, np.sum(tissue_subset)+0.3, cov_labels[cov], ha='left', va='bottom',
                fontsize=fs*1.3)
        plt.xlim(0,len(sigs_))
        plt.ylim(0,np.sum(tissue_subset))

        # print(sigs_)
        if rename_sigs is not None:
            for isig,sig in enumerate(sigs_[sigs_order_]):
                for sig_type in rename_sigs['degasperi']:
                    if sig in list(rename_sigs['degasperi'][sig_type].old):
                        ax.get_xticklabels()[isig].set_color(map_colors['degasperi'])
                    elif sig in list(rename_sigs['novel'][sig_type].new):
                        ax.get_xticklabels()[isig].set_color(map_colors['novel'])

        cax = fig.add_axes([ax.get_position().x1*1.01,
                            ax.get_position().y0+ax.get_position().height*0,
                            0.015,
                            ax.get_position().height*1])
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=fs*1.3)
        #cbar.set_label(rf"$\partial\log\mathrm{{rate}}/\partial \mathrm{{X}}$", fontsize=fs*1.3)
        cbar.set_label(r"$\beta$", fontsize=fs*1.3)

    ax = axes[0]
    plt.sca(ax)
    l1 = plt.scatter([],[],s=size_factor*2, c='k', alpha=0.7)
    l2 = plt.scatter([],[],s=size_factor*z_threshold, c='k', alpha=0.7)
    l3 = plt.scatter([],[],s=size_factor*5, c='k', alpha=0.7)
    ax.legend([l1,l2,l3], [f'$Z=2$', fr'$Z={z_threshold:.2f}$', fr'$Z=5$'],
            loc='lower right', bbox_to_anchor=(1,0.95),
            fontsize=fs*1.3, ncol=3, frameon=False)




def plotCovariateAssociations2(results_medians,
                                cov_labels = {"log_age":f"$\log(\mathrm{{age}})$", "is_female":f"Sex (is female)",
                                            "pc1":f"PC1", "pc2":f"PC2", "pc3":f"PC3", "log_n_loh":f"$\log(N_\mathrm{{LoH}})$"}):

    z_var = 'z'
    mu_var = 'mu_sd'

    z_min = np.sqrt(scipy.stats.chi2.isf(0.01, df=1))
    logpv_threshold = BH_threshold(results_medians.logpv, 0.01)
    z_threshold = np.sqrt(scipy.stats.chi2.isf(np.exp(logpv_threshold), df=1))

    results_3d_grid = {}
    for variable in [z_var, mu_var]:

        results_pivoted = pd.pivot_table(results_medians,
                        values=variable, index=['tissue'], columns=['signature', 'covariate'])

        sigs_ = np.unique(np.array(list(results_pivoted.columns))[:,0])
        covs_ = np.unique(np.array(list(results_pivoted.columns))[:,1])
        tissues_ = np.array(list(results_pivoted.index))

        results_3d_grid[variable] = np.zeros((len(results_pivoted), len(sigs_), len(covs_)))
        for isig,sig in enumerate(sigs_):
            for icov,cov in enumerate(covs_):
                try: results_3d_grid[variable][:,isig,icov] = results_pivoted[sig][cov]
                except: pass

    sig_subset = np.sum(np.abs(results_3d_grid[z_var])>z_min, axis=(0,2))
    sigs_ = sigs_[sig_subset>0]

    results_3d_grid[z_var] = results_3d_grid[z_var][:,sig_subset>0]
    results_3d_grid[mu_var] = results_3d_grid[mu_var][:,sig_subset>0]

    sigs_order_ = np.array(orderSignatures(sigs_))

    tissue_subset = np.nansum(np.abs(results_3d_grid[z_var][:,:,:])>z_min, axis=(1,2))>0
    tissue_subsets = [tissue_subset & (np.nansum(np.abs(results_3d_grid[z_var][:,:,icov])>0, axis=1)>0)\
                            for icov in range(len(covs_))]
    fig, axes = plt.subplots(len(covs_),1,figsize=(20,1.3*20*np.sum([np.sum(tss) for tss in tissue_subsets])\
                                                /results_3d_grid[z_var].shape[1]),
                            gridspec_kw={'height_ratios':[np.sum(tss) for tss in tissue_subsets]})

    plt.subplots_adjust(hspace=0.4)

    fs=16

    for icov, cov in enumerate(covs_):
        ax = axes[icov]
        plt.sca(ax)

        # Interquartile range
        dynamic_rng = np.diff(np.nanpercentile(np.array(results_3d_grid[mu_var][:,:,icov]).flatten(),
                                            np.array([25,75])))[0]

        tissue_subset = tissue_subsets[icov]
        z_table = results_3d_grid[z_var][:,:,icov][:,sigs_order_][tissue_subset]
        mu_table = results_3d_grid[mu_var][:,:,icov][:,sigs_order_][tissue_subset]

        # im = plt.pcolor(z_table, cmap=cm.coolwarm, vmin=-dynamic_rng*3, vmax=dynamic_rng*3)
        tried_grid = 1-np.isnan(z_table).astype(int)
        tried_grid[np.abs(np.array(z_table))>z_threshold] = 5
        plt.pcolor(tried_grid, cmap='bone_r', vmin=0, vmax=10)

        size_factor=5
        XX,YY = np.meshgrid(np.arange(z_table.shape[0]), np.arange(z_table.shape[1]))
        size = np.abs(z_table.T.flatten())*size_factor
        im = plt.scatter(YY.flatten()+0.5, XX.flatten()+0.5,
                        c=mu_table.T.flatten(),
                        s=np.min(np.vstack((size, np.ones(len(size))*50)), axis=0),
                        cmap=cm.bwr, vmin=-dynamic_rng*2, vmax=dynamic_rng*2)

        #significant_loc = np.unravel_index(np.argwhere(np.abs(np.array(z_table))>z_threshold), z_table.shape)
        #plt.scatter(significant_loc[1][:,1]+0.5, significant_loc[1][:,0]+0.5, marker="*", s=10, color='k')

        ax.set_yticks(np.arange(np.sum(tissue_subset)), minor=True)
        ax.set_xticks(np.arange(len(sigs_)), minor=True)
        ax.grid(True, which='minor', axis='both', linestyle='-', color='k')

        ax.set_xticks(np.arange(len(sigs_))+0.5, minor=False)
        ax.set_xticklabels(sigs_[sigs_order_], rotation=90, fontsize=fs, minor=False);
        ax.set_yticks(np.arange(np.sum(tissue_subset))+0.5, minor=False)
        ax.set_yticklabels(tissues_[tissue_subset], rotation=0, fontsize=fs, minor=False);

        plt.text(0, np.sum(tissue_subset)+0.3, cov_labels[cov], ha='left', va='bottom',
                fontsize=fs*1.3)
        #plt.title(cov_labels[cov], fontsize=fs*2.5)
        plt.xlim(0,len(sigs_))
        plt.ylim(0,np.sum(tissue_subset))

        cax = fig.add_axes([ax.get_position().x1*1.01,
                            ax.get_position().y0+ax.get_position().height*0,
                            0.015,
                            ax.get_position().height*1])
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=fs*1.3)
        cbar.set_label(r"$\log(\mathrm{odds})$", fontsize=fs*1.3)


    ax = axes[0]
    plt.sca(ax)
    l1 = plt.scatter([],[],s=size_factor*2, c='k', alpha=0.7)
    l2 = plt.scatter([],[],s=size_factor*z_threshold, c='k', alpha=0.7)
    l3 = plt.scatter([],[],s=size_factor*5, c='k', alpha=0.7)
    ax.legend([l1,l2,l3], [f'$Z=2$', fr'$Z={z_threshold:.2f}$', fr'$Z=5$'],
            loc='lower right', bbox_to_anchor=(1,0.95),
            fontsize=fs*1.3, ncol=1, frameon=False)


def plotGeneHitsTreatment(DNA_repair, treatment_list={}, annotated=False):

    fig, axes = plt.subplots(1,1,figsize=(15,5))
    divider = make_axes_locatable(axes)
    fs=16

    for i,run_name in enumerate(["cancGeneSigs/twohit_dCRT_multibinom2_tumour-group_lognloh100",
                                "treatmentSigs/treatment_dCRT_tumour-group_InDelQC"]):

        sample_df = pd.read_csv(f"{RESULT_DIR}/{run_name}/input/samples.tsv", sep="\t").set_index('sample_id')
        target_df = pd.read_csv(f"{RESULT_DIR}/{run_name}/input/targets.tsv", sep="\t").set_index('sample_id')
        test_df = pd.read_csv(f"{RESULT_DIR}/{run_name}/input/tests.tsv", sep="\t")
        signature_df = pd.read_csv(f"{RESULT_DIR}/{run_name}/input/signatures.tsv", sep="\t").set_index('sample_id')

        results_zinb = pd.read_csv(
            f"{RESULT_DIR}/results/{run_name}/output/signature_target_assoc_nb.csv"
        )
        results_log0 = pd.read_csv(
            f"{RESULT_DIR}/results/{run_name}/output/signature_target_assoc_logistic.csv"
        )

        results_df = pd.merge(results_zinb, results_log0,
                                    on=('target', 'signature', 'group'),
                                    suffixes=("_zinb", '_log0'))
        results_df['zscore'] = results_df.target_means_alt_zinb/np.sqrt(results_df.target_covs_alt_zinb)
        results_df['log10_pvalue_zinb'] = -np.log10(results_df.resample_pvalue)
        results_df['log10_pvalue_log0'] = -np.log10(results_df.wilks_pvalue_log0)
        results_df['exp_means_alt_zinb'] = np.exp(results_df['target_means_alt_zinb'])
        if i==1:
            results_df = pd.merge(results_df,
                                    pd.DataFrame(treatment_list.keys(), columns=['target']),
                                    on='target', how='inner').replace({'target':treatment_list})
        else:
            results_df = results_df[results_df.target.map(lambda x: x in np.array(DNA_repair.Gene))]

        if i==0: ax = axes
        else: ax = divider.append_axes("right", size=f"100%", pad=0.8)
        plt.sca(ax)
        xlim=[-6,8]

        if i==0:
            plt.title("Gene inactivation", fontsize=fs)
            group_col = 'group'
            test_genes = np.intersect1d(np.unique(test_df.target), DNA_repair.Gene)
            test_genes = test_genes[[not bool(re.match("MOCK", gene)) for gene in test_genes]]
            column_subset = [re.search("^[A-Z0-9]+", key).group(0) in test_genes for key in target_df.keys()]
            target_df = target_df[target_df.keys()[column_subset]]
            target_df = pd.merge(sample_df[[group_col]], target_df,
                                left_index=True, right_index=True, how='inner')\
                            .set_index(group_col)

            agg_df = pd.melt(target_df.groupby('group').sum().reset_index(), id_vars='group', var_name='target')
            agg_df = agg_df[(agg_df.value>0)&(agg_df.value<5)]
            agg_df['target'] = agg_df.target.map(lambda x: x.split("_")[0])
            agg_df.drop_duplicates(['group', 'target'], inplace=True)

            results_df = pd.merge(results_df, agg_df, how='outer', on=['group', 'target'], indicator=True)
            results_df = results_df[results_df._merge.map(lambda x: x in ['both', 'left_only'])]
            subset = results_df._merge=='left_only'
        else:
            plt.title("Treatment", fontsize=fs)
            subset = np.ones(len(results_df))

        subset = (~mock_match(results_df.target))&\
                (~mock_match(results_df.signature))&\
                (results_df.signature.map(lambda x: x[0]!="T"))&\
                (results_df.resample_pvalue>0)&\
                (results_df.signature!='CN1')&\
                subset

        p_threshold = BH_threshold(np.array(results_df['resample_pvalue']), alpha=0.01)

        # zscore = results_df.target_means_alt/np.sqrt(results_df.target_covs_alt)
        im = plt.scatter(results_df.target_means_alt_zinb[subset],
                    -np.log10(results_df.resample_pvalue[subset]),
                    c=-np.log10(results_df.wilks_pvalue_log0[subset]),
                    s=5, cmap='viridis_r', vmax=12)

        ylim = (0, ax.axes.get_ylim()[1])

        plt.plot(xlim,[-np.log10(p_threshold),-np.log10(p_threshold)], '--r')
        plt.plot([0,0],[-np.log10(p_threshold),ylim[1]], '--r')

        plt.xlim(xlim)#-np.max(np.abs(zscore[subset])), np.max(np.abs(zscore[subset])))
        plt.gca().set_ylim(bottom=0, top=ylim[1])

        plt.xlabel(r"$\beta$")
        plt.ylabel(r"$-\log_{10}(P_\mathrm{NB/Poisson-dCRT})$")

        if i==1:
            cax = divider.append_axes("right", size=f"5%", pad=0.1)
            cbar = plt.colorbar(im, cax=cax, orientation='vertical')
            cbar.set_label(r"$-\log_{10}(P_\mathrm{logistic})$")
            # cbar.ax.tick_params(labelsize=fs*1.3)

        if annotated:

            if i==0:
                pos_sets = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,3,10),
                                            ('POLE','Uterus-AdenoCA',p_threshold,4,10)],
                            'signature,target':[('SBS18','MUTYH',p_threshold,4,0)],
                        }

                for key,combos in pos_sets.items():
                    for i,set_row in enumerate(combos):
                        rset = results_df[subset&(results_df[key.split(",")[0]]==set_row[0])&\
                                            (results_df[key.split(",")[1]]==set_row[1])&\
                                            (results_df.resample_pvalue<set_row[2])&\
                                            (results_df.signature!='CN1')]
                        posy = np.mean(-np.log10(rset.resample_pvalue))+set_row[4]
                        if set_row[3]<0: ha='right'
                        else: ha='left'
                        plt.sca(ax)
                        if key=='target,group':
                            plt.text(set_row[3], posy, f"${set_row[0]}$,\n{set_row[1]}", ha=ha, va='center')
                        else:
                            plt.text(set_row[3], posy, f"{set_row[0]} (${set_row[1]}$)", ha=ha, va='center')
                        for index, row in rset.iterrows():
                            plt.plot([set_row[3],row.target_means_alt_zinb],
                                    [posy,-np.log10(row.resample_pvalue)], ':k', alpha=0.3)


                rset = results_df[subset&(results_df.resample_pvalue<1e-15)&\
                                    (results_df.signature=='CN1')]
                posy = np.mean(-np.log10(rset.resample_pvalue))+10
                if -2<0: ha='right'
                else: ha='left'
                plt.sca(ax)
                plt.text(-2, posy, f"CN1", ha=ha, va='center')
                for index, row in rset.iterrows():
                    plt.plot([-2,row.target_means_alt_zinb],
                            [posy,-np.log10(row.resample_pvalue)], ':k', alpha=0.3)

                # Highlight hits
                select = [('POLG', 'SBS93', 'ColoRect-AdenoCA',-1.7, 21),
                        ('MSH6', 'CN25', 'ColoRect-AdenoCA', -2.2, 18),
                        ('POLG', 'DBS4', 'ColoRect-AdenoCA', -2, 10),
                        ('BRCA2', 'DBS2', 'Breast-DuctalCA', -2, 5)]
                for i, set_row in enumerate(select):
                    row = results_df[(results_df.target==set_row[0])&\
                                            (results_df.signature==set_row[1])&\
                                            (results_df.group==set_row[2])].iloc[0]
                    pos = (row.target_means_alt_zinb, -np.log10(row.resample_pvalue))
                    plt.sca(ax)
                    text_pos = (pos[0]+set_row[3], pos[1]+set_row[4])
                    plt.text(text_pos[0], text_pos[1], f"{row.signature} (${row.target}$),\n{row.group}",
                            ha='left' if set_row[3]>0 else 'right', va='bottom')
                    plt.plot([pos[0],text_pos[0]],[pos[1], text_pos[1]], ':k', alpha=0.3)


            elif i==1:
                plt.sca(ax)

                pos_sets = {'target,group':[],
                            'signature,target':[('DBS2','Radiotherapy',p_threshold,0,-1.2),
                                                ('ID5','Radiotherapy',p_threshold,0,4),
                                                ('ID8','Radiotherapy',p_threshold,-7,4.2)],
                            'signature,group':[('DBS5','ColoRect-AdenoCA',1e-10,3,-0.5)]}
                for key,combos in pos_sets.items():
                    for i,set_row in enumerate(combos):
                        print(set_row)
                        rset = results_df[subset&(results_df[key.split(",")[0]]==set_row[0])&\
                                            (results_df[key.split(",")[1]]==set_row[1])&\
                                            (results_df.resample_pvalue<set_row[2])]
                        posy = np.mean(-np.log10(rset.resample_pvalue))+set_row[3]
                        if set_row[4]<0: ha='right'
                        else: ha='left'
                        if key=='signature,group':
                            plt.text(set_row[4], posy, f"{set_row[0]}, \n{set_row[1]}", ha=ha, va='center')
                        else:
                            plt.text(set_row[4], posy, f"{set_row[0]} ({set_row[1]})", ha=ha, va='center')
                        for index, row in rset.iterrows():
                            plt.plot([set_row[4],row.target_means_alt_zinb],
                                    [posy,-np.log10(row.resample_pvalue)], ':k', alpha=0.3)

                # Highlight hits
                select = [('Oxaliplatin', 'DBS5', 'ColoRect-AdenoCA',0,-0.5),]
                for i, set_row in enumerate(select):
                    row = results_df[(results_df.target==set_row[0])&\
                                            (results_df.signature==set_row[1])&\
                                            (results_df.group==set_row[2])].iloc[0]
                    pos = (row.target_means_alt_zinb, -np.log10(row.resample_pvalue))
                    text_pos = (pos[0]+set_row[4],
                                pos[1]+set_row[3])
                    plt.text(text_pos[0], text_pos[1], f"{row.signature} ({row.target}),\n{row.group}",
                            ha='right' if set_row[4]<0 else 'left', va='center')
                    plt.plot([pos[0],text_pos[0]],[pos[1], text_pos[1]], ':k', alpha=0.3)

                ellipse = Ellipse((5.3,27), width=3, height=30, facecolor='none', edgecolor='grey', linestyle="--")
                ax.add_patch(ellipse)
                plt.text(5.3,44,"Breast-DuctalCA,\nSBS26/SBS44", ha='center', va='center')

    plt.text(0.087, 0.9, r"a)", fontsize=20, weight='bold', transform=plt.gcf().transFigure)
    plt.text(0.49, 0.9, r"b)", fontsize=20, weight='bold', transform=plt.gcf().transFigure)


def plotGermlineSomatic(results_1, results_2,
                        highlight_set = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,10,2),
                                                    ('POLE','Uterus-AdenoCA',1e110,8,2.2)],
                                    'signature,target':[('SBS18','MUTYH',1e-10,2,1.5)]},
                        highlight_hit = [('POLG', 'DBS4', 'ColoRect-AdenoCA', -0.1, 2),
                                         ('POLG', 'SBS93', 'ColoRect-AdenoCA', -0.1, 5)],
                        labels=['1','2']):

    # Merge the results tables
    results_merged = pd.merge(results_1, results_2,
                              how='outer', on=['group','signature', 'target'], suffixes=['_1', '_2'])

    # Z-score - Q plot
    fig, ax = plt.subplots(1,1,figsize=(8,8))

    pos_neg_sets = [(results_merged.target_means_alt_zinb_1<0)&(results_merged.target_means_alt_zinb_2<0),
              (results_merged.target_means_alt_zinb_1<0)&(results_merged.target_means_alt_zinb_2>0),
              (results_merged.target_means_alt_zinb_1>0)&(results_merged.target_means_alt_zinb_2<0),
              (results_merged.target_means_alt_zinb_1>0)&(results_merged.target_means_alt_zinb_2>0)]
    pos_neg_labels = [f"{labels[0]}<0, \n{labels[1]}<0",
              f"{labels[0]}<0, \n{labels[1]}>0",
              f"{labels[0]}>0, \n{labels[1]}<0",
              f"{labels[0]}>0, \n{labels[1]}>0"]
    pos_neg_colors = default_colours[:4]
    for i, pos_neg_set in enumerate(pos_neg_sets):
        plt.scatter(-np.log10(results_merged[pos_neg_set].resample_pvalue_1),
                    -np.log10(results_merged[pos_neg_set].resample_pvalue_2),
                    c=pos_neg_colors[i], s=4,
                    label=pos_neg_labels[i])

    ylim = (0, ax.axes.get_ylim()[1])
    xlim = (0, ax.axes.get_xlim()[1])
    plt.legend(loc='upper left', bbox_to_anchor=(1,1), frameon=False)

    p_threshold = BH_threshold(results_merged['resample_pvalue_2'], 0.01)
    plt.plot(xlim,[-np.log10(p_threshold),-np.log10(p_threshold)], '--r')

    p_threshold = BH_threshold(results_merged['resample_pvalue_1'], 0.01)
    plt.plot([-np.log10(p_threshold),-np.log10(p_threshold)], ylim, '--r')

    plt.xlim(xlim)#-np.max(np.abs(zscore)), np.max(np.abs(zscore)))
    plt.gca().set_ylim(bottom=0, top=ylim[1])

    plt.xlabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{ZINB-dCRT})$ "+labels[0])
    plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{ZINB-dCRT})$ "+labels[1])

    for key,combos in highlight_set.items():
        for i,set_row in enumerate(combos):
            rset = results_merged[(results_merged[key.split(",")[0]]==set_row[0])&\
                                  (results_merged[key.split(",")[1]]==set_row[1])&\
                                 ((results_merged.resample_pvalue_1<set_row[2])|\
                                  (results_merged.resample_pvalue_2<set_row[2]))&\
                                  (results_merged.signature!='CN1')]
            if set_row[4]<0: ha='right'
            else: ha='left'
            dx = -np.nanmean(-np.log10(rset.resample_pvalue_1) - set_row[3])/(xlim[1]-xlim[0])
            dy = -np.nanmean(-np.log10(rset.resample_pvalue_2) - set_row[4])/(ylim[1]-ylim[0])
            print(f"{set_row[0]}, \n{set_row[1]}", dx, dy)
            plt.text(set_row[3], set_row[4], f"{set_row[0]}, \n{set_row[1]}",
                     ha='center' if (np.abs(dy)>np.abs(dx)) else 'left' if (dx>0) else 'right',
                     va='center' if (np.abs(dy)<np.abs(dx)) else 'bottom' if (dy>0) else 'top')
            for index, row in rset.iterrows():
                plt.plot([set_row[3],-np.log10(row.resample_pvalue_1)],
                         [set_row[4],-np.log10(row.resample_pvalue_2)], ':k', alpha=0.3)

if __name__=='__main__':

    # Get renaming dictionary
    combined_acts, combined_sigs, relabel_map = loadSignatures()
    rename_dict = signatureRenamer(relabel_map)

    # ## Two hit rates
    run_name="twohit"
    run_dir = f"{RESULT_DIR}/results/associations/cancGeneHits/{run_name}"

    sample_df = pd.read_csv(f"{run_dir}/input/samples.tsv", sep="\t").set_index('sample_id')
    target_df = pd.read_csv(f"{run_dir}/input/targets.tsv", sep="\t").set_index('sample_id')
    test_df = pd.read_csv(f"{run_dir}/input/tests.tsv", sep="\t")
    signature_df = pd.read_csv(f"{run_dir}/input/signatures.tsv", sep="\t").set_index('sample_id')

    # Replace Connective with Sarcoma
    sample_df.group = sample_df.group.str.replace("Connective","Sarcoma")
    test_df.group = test_df.group.str.replace("Connective","Sarcoma")

    # Subset DNA repair gene list
    DNA_repair = pd.read_csv(f"{DATA_DIR}/human-dna-repair-genes.tsv", sep="\t")
    DNA_repair = DNA_repair[DNA_repair.Gene.map(lambda x: x in target_df.keys())]
    Types = ['BER','MMR','Homologous recombination','Fanconi anemia',
            'DNA polymerases ','NER','XPG', 'Editing and processing nucleases',
            'Direct reversal of damage', 'Chromatin Structure and Modification']
    DNA_repair = DNA_repair[DNA_repair.Type.map(lambda x: x in Types)]

    plotHitRates(sample_df, target_df, test_df, DNA_repair)
    publish_fig("two_hit_counts", publish=FIGURE_DIR)

    # ## Two hit rates
    for run_name in ["twohit", "twohit_lognloh"]:
        run_dir = f"{RESULT_DIR}/associations/cancGeneHits/{run_name}"

        sample_df = pd.read_csv(f"{run_dir}/input/samples.tsv", sep="\t").set_index('sample_id')
        target_df = pd.read_csv(f"{run_dir}/input/targets.tsv", sep="\t").set_index('sample_id')
        test_df = pd.read_csv(f"{run_dir}/input/tests.tsv", sep="\t")
        signature_df = pd.read_csv(f"{run_dir}/input/signatures.tsv", sep="\t").set_index('sample_id')


        results_combined = getAssociationResults(run_dir, DNA_repair)
        results_combined.replace({'signature':rename_dict}, inplace=True)

        subset = (~mock_match(results_combined.target))&\
                (~mock_match(results_combined.signature))&\
                (results_combined.signature.map(lambda x: x[0]!="T"))&\
                (results_combined.resample_pvalue>0)&\
                (results_combined.signature!='CN1')&\
                (results_combined._merge=='left_only')

        # Point association figures
        p_threshold = BH_threshold(results_combined['resample_pvalue'], 0.01)
        highlight_set = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,10,-1),
        #                                  ('BRCA1','Breast-DuctalCA',1e-15,10,2),
                                        ('POLE','Uterus-AdenoCA',1e-5,3,2.9)],
                        'signature,target':[('SBS18','MUTYH',p_threshold,2,2.2)]}
        highlight_hit = [('POLG', 'DBS4', 'ColoRect-AdenoCA', -0.01, 2),
                        ('POLG', 'SBS93', 'ColoRect-AdenoCA', -0.1, 7),
                        ('BRCA1', 'CN17', 'Breast-DuctalCA', 1, 7),
                        ('BRCA2', 'ID6', 'Breast-DuctalCA', 1, 7),
                        ('BRCA2', 'SBS3', 'Breast-DuctalCA', 0.5, 1.5)]
        plotTwoHitAssoc(results_combined[subset], highlight_set=highlight_set,
                                                            highlight_hit=highlight_hit)
        publish_fig(f"{run_name}_sig_bqplot", publish=FIGURE_DIR)

        figs, axes = plotTwoHitAssocZoom(results_combined[subset], p_threshold=p_threshold,
                                                    highlight_set=highlight_set, highlight_hit=highlight_hit)
        plt.sca(axes[0])
        publish_fig(f"{run_name}_sig_bqplot_full", publish=FIGURE_DIR)
        plt.sca(axes[1])
        publish_fig(f"{run_name}_sig_bqplot_zoom", publish=FIGURE_DIR)

        plotZHist(results_combined)
        publish_fig(f"{run_name}_sig_zscore_histograms", publish=FIGURE_DIR)

        plotZBetaHist(results_combined, p_threshold=p_threshold)
        publish_fig(f"{run_name}_sig_z-and-beta_histograms", publish=FIGURE_DIR)

        plotCohortHits(results_combined[subset&(results_combined.resample_pvalue<p_threshold)], DNA_repair)
        publish_fig(f"{run_name}_target-sig_meta-z", publish=FIGURE_DIR)

        print("Two hit association mock statistics")
        printMockStats(results_combined)

        plotMockPP(results_combined)
        publish_fig(f"{run_name}_mock_gene_qqplots_x3", publish=FIGURE_DIR)

        # # Target-target confounding
        tgt_df = pd.read_csv(f"{run_dir}/output/target_target_assoc.csv")
        tgt_df.wilks_pvalue[tgt_df.wilks_pvalue==0] = 1e-300

        # Replace Connective with Sarcoma
        tgt_df.group = tgt_df.group.str.replace("Connective","Sarcoma")

        # Get gene positions
        census_df = pd.read_csv(f"{DATA_DIR}/cancer_gene_census.csv")
        gene_loc = census_df['Genome Location'].str.split("[:-]",2,expand=True).rename({0:'chrom',1:'start',2:'end'}, axis=1)
        gene_loc['gene'] = census_df['Gene Symbol'].map(lambda x: x.replace('-','.'))

        # Subset of Genes to plot
        target_subset = np.unique(results_combined.target[subset&(results_combined.resample_pvalue<p_threshold)])
        target_subset = target_subset[~mock_match(target_subset)]
        target_subset = np.intersect1d(target_subset, gene_loc['gene'])

        gene_loc = gene_loc[(gene_loc.start!='')&(gene_loc.end!='')]

        gene_loc['chrom_num'] = gene_loc.chrom.replace('X','24').astype(int)
        gene_loc.start = gene_loc.start.astype(int)
        gene_loc.end = gene_loc.end.astype(int)

        gene_loc = gene_loc.set_index('gene').loc[target_subset]
        gene_loc = gene_loc.sort_values(['chrom_num','start'])
        ordered_genes = np.array(gene_loc.index)

        groups=['ColoRect-AdenoCA', 'Uterus-AdenoCA', 'Breast-DuctalCA']
        plotTargetConfounding(tgt_df, target_subset, ordered_genes,
                                        p_threshold=BH_threshold(tgt_df['wilks_pvalue']),
                                        groups =groups)
        publish_fig(f"{run_name}_target_grid_beta_{'_'.join(groups)}", publish=FIGURE_DIR)

        # Plot covariates
        results_medians = getCovariateAssociations(results_combined[subset])
        covs_ = np.array(['is_female', 'log_age', 'pc1', 'pc2', 'pc3'])
        plotCovariateAssociations(results_medians, relabel_map,
                                            covs_ = covs_)
        publish_fig(f"{run_name}_covariate_associations_betasdlog10_{'-'.join(covs_)}", publish=FIGURE_DIR)

    # # Gene and treatment
    annotated=False
    treatment_list = {'rtprescribeddose': 'Radiotherapy', 'RADIOTHERAPY':'RADIOTHERAPY',
                'OXALIPLATIN': 'OXALIPLATIN', 'CARBOPLATIN': 'CARBOPLATIN', 'CISPLATIN': 'CISPLATIN',
                'FLUOROURACIL': 'FLUOROURACIL', 'CYCLOPHOSPHAMIDE': 'CYCLOPHOSPHAMIDE', 'EPIRUBICIN': 'EPIRUBICIN',
                'BEVACIZUMAB': 'BEVACIZUMAB',
                'BOSUTINIB': 'BOSUTINIB',
                'CAPECITABINE': 'CAPECITABINE',
                'CETUXIMAB': 'CETUXIMAB',
                'DASATINIB': 'DASATINIB',
                'DOCETAXEL': 'DOCETAXEL',
                'DOXORUBICIN': 'DOXORUBICIN',
                'FOLINIC_ACID': 'FOLINIC_ACID',
                'GEMCITABINE': 'GEMCITABINE',
                'HYDROXYCARBAMIDE': 'HYDROXYCARBAMIDE',
                'IFOSFAMIDE': 'IFOSFAMIDE',
                'IMATINIB': 'IMATINIB',
                'IRINOTECAN': 'IRINOTECAN',
                'PACLITAXEL': 'PACLITAXEL',
                'PEMETREXED': 'PEMETREXED',
                'PERTUZUMAB': 'PERTUZUMAB',
                'TEMOZOLOMIDE': 'TEMOZOLOMIDE',
                'TRASTUZUMAB': 'TRASTUZUMAB'}
    treatment_list = {key:treatment_list[key][0].upper()+treatment_list[key][1:].lower() for key in treatment_list}
    plotGeneHitsTreatment(DNA_repair, treatment_list=treatment_list)
    publish_fig(f"sig_twohit-treatment_bqplot{'_annotated' if annotated else ''}", publish=FIGURE_DIR)


    # Germline
    run_dir = f"{RESULT_DIR}/associations/cancGeneHits/germline/"
    results_germline = getAssociationResults(run_dir, DNA_repair)
    results_germline.replace({'signature':rename_dict}, inplace=True)

    results_germline = pd.merge(results_germline.rename({'_merge':'_merge0'}, axis=1),
            results_combined[['group', 'signature', 'target', 'resample_pvalue']],
            on=['group', 'signature', 'target'], how='left', indicator=True,
            suffixes=("", "_twohit"))

    germline_subset = (~mock_match(results_germline.target))&\
                        (~mock_match(results_germline.signature))&\
                        (results_germline.signature.map(lambda x: x[0]!="T"))&\
                        (results_germline.resample_pvalue>0)&\
                        (results_germline.signature!='CN1')&\
                        (results_germline._merge0=='left_only')&\
                        (results_germline.resample_pvalue_twohit<p_threshold)

    germline_p_threshold = BH_threshold(results_germline[germline_subset]['resample_pvalue'], 0.01)
    highlight_set = {'target,group':[],
                    'signature,target':[('SBS18','MUTYH',germline_p_threshold,5,3),
                                        ('DBS2','BRCA2',germline_p_threshold,0.1,-0.2),
                                        ('SBS44','MSH6',germline_p_threshold,7,3),
                                        ('SBS26','MSH6',germline_p_threshold,2,3.5),
                                        ('SV5','BRCA2',germline_p_threshold,2.5,-0.2),
                                        ('ID6', 'BRCA2',germline_p_threshold,5,2),]}
    highlight_hit = [('POLD1', 'SBS23', 'Bladder-TCC', -2, -10),
                    ('BRCA2', 'SBS3', 'Breast-DuctalCA', -0.5, 3)]

    plotTwoHitAssoc(results_germline[germline_subset],
                    highlight_set=highlight_set, highlight_hit=highlight_hit)
    publish_fig("sig_germline_bqplot", publish=FIGURE_DIR)

    figs, axes = plotTwoHitAssocZoom(results_germline[germline_subset], p_threshold=germline_p_threshold,
                                     highlight_set=highlight_set, highlight_hit=highlight_hit)
    plt.sca(axes[0])
    publish_fig("sig_germline_bqplot_full", publish=FIGURE_DIR)
    plt.sca(axes[1])
    publish_fig("sig_germline_bqplot_zoom", publish=FIGURE_DIR)

    # Somatic
    run_dir = f"{RESULT_DIR}/associations/cancGeneHits/somatic/"
    results_somatic = getAssociationResults(run_dir, DNA_repair)
    results_somatic.replace({'signature':rename_dict}, inplace=True)

    results_somatic = pd.merge(results_somatic.rename({'_merge':'_merge0'}, axis=1),
            results_combined[['group', 'signature', 'target', 'resample_pvalue']],
            on=['group', 'signature', 'target'], how='left', indicator=True,
            suffixes=("", "_twohit"))

    somatic_subset = (~mock_match(results_somatic.target))&\
                        (~mock_match(results_somatic.signature))&\
                        (results_somatic.signature.map(lambda x: x[0]!="T"))&\
                        (results_somatic.resample_pvalue>0)&\
                        (results_somatic.signature!='CN1')&\
                        (results_somatic._merge0=='left_only')&\
                        (results_somatic.resample_pvalue_twohit<p_threshold)

    somatic_p_threshold = BH_threshold(results_somatic[somatic_subset]['resample_pvalue'], 0.01)
    highlight_set = {'target,group':[],
                    'signature,target':[]}
    highlight_hit = []

    plotTwoHitAssoc(results_somatic[somatic_subset],
                                highlight_set=highlight_set, highlight_hit=highlight_hit,
                                xlim=[[-4,10],[-4,10]])
    publish_fig("sig_somatic_bqplot", publish=FIGURE_DIR)

    figs, axes = plotTwoHitAssocZoom(results_somatic[somatic_subset], p_threshold=somatic_p_threshold,
                                                highlight_set=highlight_set, highlight_hit=highlight_hit)
    plt.sca(axes[0])
    publish_fig("sig_somatic_bqplot_full", publish=FIGURE_DIR)
    plt.sca(axes[1])
    publish_fig("sig_somatic_bqplot_zoom", publish=FIGURE_DIR)

    ##### ------ Germline vs Somatic associations ------ #####
    # ## Two hit rates
    result_dir = f"{RESULT_DIR}/associations/cancGeneHits/"
    all_results_combined = {}
    for run_name in ["germline", "somatic", "germline_30"]:
        run_dir = f"{result_dir}/{run_name}"

        sample_df = pd.read_csv(f"{run_dir}/input/samples.tsv", sep="\t").set_index('sample_id')
        target_df = pd.read_csv(f"{run_dir}/input/targets.tsv", sep="\t").set_index('sample_id')
        test_df = pd.read_csv(f"{run_dir}/input/tests.tsv", sep="\t")
        signature_df = pd.read_csv(f"{run_dir}/input/signatures.tsv", sep="\t").set_index('sample_id')

        results_combined = getAssociationResults(run_dir, DNA_repair)
        results_combined.replace({'signature':rename_dict}, inplace=True)

        subset = (~mock_match(results_combined.target))&\
                (~mock_match(results_combined.signature))&\
                (results_combined.signature.map(lambda x: x[0]!="T"))&\
                (results_combined.resample_pvalue>0)&\
                (results_combined.signature!='CN1')&\
                (results_combined._merge=='left_only')

        # Point association figures
        p_threshold = BH_threshold(results_combined['resample_pvalue'], 0.01)
        highlight_set = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,10,-1),
        #                                  ('BRCA1','Breast-DuctalCA',1e-15,10,2),
                                        ('POLE','Uterus-AdenoCA',1e-5,3,2.9)],
                        'signature,target':[('SBS18','MUTYH',p_threshold,2,2.2)]}
        highlight_hit = [('POLG', 'DBS4', 'ColoRect-AdenoCA', -0.01, 2),
                        ('POLG', 'SBS93', 'ColoRect-AdenoCA', -0.1, 7),
                        ('BRCA1', 'CN17', 'Breast-DuctalCA', 1, 7),
                        ('BRCA2', 'ID6', 'Breast-DuctalCA', 1, 7),
                        ('BRCA2', 'SBS3', 'Breast-DuctalCA', 0.5, 1.5)]
        all_results_combined[run_name] = results_combined[subset]


    # Germline CADD>20 vs CADD>30
    highlight_set = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,13,60),
                                    ('MSH6','Breast-DuctalCA',1e-15,8,100),
                                    ('POLG','Breast-DuctalCA',1e-15,2,70)],
                    'signature,target':[('ID6','BRCA2',p_threshold,9.5,25)]}
    plotGermlineSomatic(all_results_combined['germline'], all_results_combined['germline_30'],
                    highlight_set=highlight_set, highlight_hit=highlight_hit, labels=[r'$\mathrm{CADD}_{20}$',
                                                                                    r'$\mathrm{CADD}_{30}$'])
    publish_fig("sig_germline_20vs30_assoc", publish=FIGURE_DIR)

    # Germline vs Somatic
    highlight_set = {'target,group':[('MSH6','ColoRect-AdenoCA',1e-15,12,60),
                                    ('POLE','Uterus-AdenoCA',1e-20,10,100),
                                    ('POLE','ColoRect-AdenoCA',1e-20,10,120)
                                    ],
                    'signature,target':[('ID6','BRCA2',p_threshold,10,20),
                                        ('DBS4','BRCA2',p_threshold,15,20)]}
    plotGermlineSomatic(all_results_combined['germline'], all_results_combined['somatic'],
                        highlight_set=highlight_set, highlight_hit=highlight_hit,
                        labels=['Germline', 'Somatic'])
    publish_fig("sig_germlinevssomatic_assoc", publish=FIGURE_DIR)