# Genomics England pancancer signatures
Code and results for paper:

[Comprehensive repertoire of the chromosomal alteration and mutational signatures across 16 cancer types from 10,983 cancer patients](doi.org/10.1101/2023.06.07.23290970)

Disclaimer: This code was written inside the Genomics England Research Environment without github and it has not been tested outside the research environment. We provide it here to aid the transparency of the publication and make our analysis methods more accessible and reproducible. The pipelines provided will almost certainly not be of significant use outside Genomics England, however, we hope users find pieces of code helpful to understand our work and useful in their research.

# Getting started

Clone the repository
`git clone https://github.com/Wedge-lab/Gel_pan_cancer_signatures.git`

Create .env file
```
cd Gel_pan_cancer_signatures
cp .env-template .env
```
in `.env`, change all file names as required to your local files.

Run editable install on a fresh conda environment (python=3.9)
`pip install -e .`


# Contained in this repository


### Curate gene mutations
Bash scripts in `./scripts/data_prep` are for collecting and curating genetic mutations in WGS data available in Genomics England including inferring expected impact using CADD scores or OncoKB annotations.


### Combine signatures




### Association analysis

A series of scripts in `./scripts` run association ananlysis on mutation rates with various covariates including germline and somatic gene inactivations and treatment exposures:
```
germline_assoc.sh
treatment_assoc.sh
genotype_assoc.sh
somatic_assoc.sh
twohit_assoc.sh
```
Each script prepares data and runs the `signatureAssociations.nf` pipeline.

*This is the most interesting thing probably contained in this code*
Go to `pipelines/signaturessignatureAssociations.nf/README.md` for more details.
This pipeline contains the GLM association analysis with resampling which reduces the false positive rate of associations.


### Plots

This repository includes code for generating the figures in the publication. All plotting scripts
are contained in `src/signatures/plotting`.

#### Main plotting scripts

**`combinedSignatures.py`**
- `signature_distributions_alltypes` (`Figure 2`): Novel signature profiles across all mutation types (SBS, DBS, ID, CNV, SV)
- `signature_sizes_alltypes` (`Figure 1`): Signature activities by tumour group
- `signature_activity_clusters` (`Figure 3`): Full hierarchical clustering diagram of signature activities with Fisher's exact test and Spearman correlations
- `signature_activity_cluster_subgroups` (`Extended Data Figure 2`): Clustered signature groups (UV, Smoking, MMR, POLE, HRD, POLG)
- `signature_activity_cluster_subsets`: Selected signature subsets for specific pathways
- `total_mutation_rates` (`Extended Data Figure 1`): Total mutation rates across tumour groups for all signature types

**`signatureSelect.py`**
- `{cohort}_SV32_aic`: AIC plots for structural variant signature selection per cohort
- `{cohort_set}_SV32_aic` (`Supplementary Figure 13`): Combined AIC plots for multiple cohorts (e.g., Breast-Ovary-Uterus)

**`triangulation_plot.py`**
- `triangulation` (`Exteneded Data Figure 6`): Comprehensive overview figure showing signature associations across gene inactivation, genomic alterations, treatment, clinical status, survival, and timing

**`mutationTimeR.py`**
- `subclonal_fraction_rank_mutationtimer` (`Figure 7`): Subclonal fraction of signatures across tumour groups with statistical significance
- `lateclonal_fraction_rank_mutationtimer` (`Extended Data Figure 5`): Late clonal fraction of signatures (late vs early clonal) across tumour groups

#### Association plotting scripts (`associations/`)

**`clinicalAssoc.py`**
- `clinical_panel` (`Figure 6`): Multi-panel figure showing associations between signatures and clinical variables (grade, stage, tumour-specific characteristics)
- `histology_rate_comparison` (`Figure 4`): Histology subtype comparisons within tumour groups

**`genomicAlterations.py`**
- `genomic_alteration_assoc` (`Extended Data Figure 3`): Associations between signatures and genomic alterations (WGD, chromothripsis, chromoplexy, tandem duplications, kataegis) across all signature types

**`survivalFigs.py`**
- `survival_panel_survival_log_grade-4_unique` (`Figure 8`): Kaplan-Meier survival curves and CPH association scatter plots

**`treatmentFigs.py`**
- `sig_treatment_bqplot`: Beta-Q plot showing treatment-signature associations
- `treatment_mock_qqplots_x3`: Q-Q plots for mock associations (signature mock and target resample)
- `treatment_grid_beta_...` (`Supplementary Figure 6`): Treatment-treatment correlation grids for specific tumour groups
- `treatment_counts`: Heatmap of treatment sample counts across tumour groups

**`twoHitAssoc.py`**
- `two_hit_counts` (`Supplementary Figure 15`): Bar charts showing gene inactivation hit rates (total, germline, somatic, LoH) across cohorts
- `..._sig_bqplot`: Beta-Q plots for two-hit gene inactivation associations
- `..._sig_zscore_histograms`: Z-score histograms comparing mock vs significant associations
- `..._sig_z-and-beta_histograms`: Combined Z-score and beta histograms
- `..._target-sig_meta-z`: Pan-cancer meta-analysis Z-scores for gene-signature associations
- `..._mock_gene_qqplots_x3` (`Supplementary Figure 17`): Q-Q plots for mock gene associations
- `..._target_grid_beta_...` (`Supplementary Figure 5`): Gene-gene correlation grids showing confounding effects
- `twohit_covariate_associations_betasdlog10_is_female-log_age-pc1-pc2-pc3` (`Extended Data Figure 4`): Covariate association plots (sex, age, PCs)
- `sig_germline_bqplot`: Germline-only associations
- `sig_somatic_bqplot`: Somatic-only associations
- `sig_germline_20vs30_assoc`: Comparison of germline associations using CADD>20 vs CADD>30 thresholds
- `sig_germlinevssomatic_assoc`: Comparison of germline vs somatic associations
- `sig_twohit-treatment_bqplot` (`Figure 5`): Combined two-hit and treatment association plot
