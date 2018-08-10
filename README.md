## TF2Exp: regression models to predict the impact of altered TF binding on gene expression levels

Author: Wenqiang Shi.

Email:wqshi.nudt@gmail.com

## Introduction

TF2Exp is a gene-based framework to predict the impact of altered TF-binding events on gene expression levels based on regulatory variants. Using data from lymphoblastoid cell lines, TF2Exp models were applied successfully to predict the expression levels of ~3,100 genes. More details please see the citation.

In this repository, the provided pipeline would use the trained models to predict the impact of variants given by users.

The model training scripts of the TF2Exp paper are located in https://github.com/wqshi/TF2Exp_paper


## Example

Please see the s_test_tf2exp.R

Rscript s_tf2exp.R --chr_str chr22 \
                   --vcf_file ./data/test_data/test.vcf.gz \
                   --output_dir ./data/test_data \
                   --deepsea_dir ./data/tf2exp_data/deepsea/ \
                   --model_dir ./data/tf2exp_data/models/ 
                   

## Input parameters for TF2Exp pipeline:

TF2Exp models only require genotype data in vcf format as input. The vcf file should contain genotype data for at least one individual.

The output is the predicted expression levels for impacted genes.

In addition, TF2Exp further needs two types of data, pre-compiled deepsea data and trained models. The repository include the test data. They are located at ./data/tf2exp_data/models/ by default.

The complete set of the data is stored at XXXXX. We have pre-compiled all the variants for 358 European individuals in the 1000 Genomes Project. The pre-compiled set might miss some variants in the user vcf file. Alternatively, user can use deepsea to predict the impact of all the variants. A guideline will be provided soon.


## Software used in the testing:

* R 3.1.3

* OS: CentOS 5

* Other R packages include data.table.1.9.6, RUnit.0.4.29, plyr.1.8.3, stringi.0.5-5, dplyr0.5.0, stringr.1.0.0, futile.logger.1.4.1.

## To do:


## Citation:

Shi, Wenqiang, Oriol Fornes, and Wyeth W. Wasserman. "Altered transcription factor binding events predict personalized gene expression and confer insight into functional cis-regulatory variants." bioRxiv (2017): 228155.





















