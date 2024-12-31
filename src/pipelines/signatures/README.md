# SignatureAssociations.nf
Pipeline for running associations between signature activities and any target variables.

## Summary

Signature activities are counts of mutations.
Naively these could be modelled with a Poisson distribution however there are lots of un-measured covariates which determine the mean.
The Negative Binomial distribution is a Poisson distribution with mean drawn from a Gamma distribution. This enables excess variance over the pure Poisson.
The signature extraction also tends to set small signature activity counts to zero therefore a zero-inflated component to our model is also required.

This pipeline uses a distilled-CRT method on a zero-inflated Negative Binomial model to fit the signature activities in cohorts.
This should reduce the Type-I error induced by extreme outliers in the data.


## Pipeline

'''
nextflow run src/signatures/signatureAssociations.nf -with-trace -resume  \
        --filename_tests ${filename_tests} \
        --filename_signatures ${filename_signatures} \
        --filename_targets ${filename_targets} \
        --filename_samples ${filename_samples} \
        --max_signatures 1000 \
        --output_dir ${dir_output} \
        --run_name ${run_name}
'''

The pipeline requires a:
- filename_tests: TSV with signature, group, target, type - type is binary/continuous to determine whether the target is modelled with logistic/linear regression.
- filename_signatures: TSV of signature activities including 'sample_id' column and a column for each signature.
- filename_targets: TSV of target variables which we want to test associations including 'sample_id' column.
- filename_samples: TSV of sample covariates including 'sample_id' and 'group' columns. Group are the cohorts in which the fits are performed.
- output_dir: Directory to which results are saved.


## Output

glm_results.tsv

logistic_results.tsv
