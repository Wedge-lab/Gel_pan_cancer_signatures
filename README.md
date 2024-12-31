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


