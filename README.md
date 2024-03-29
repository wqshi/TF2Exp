## TF2Exp: regression models to predict the impact of altered TF binding on gene expression levels

Programed by: Wenqiang Shi.

Contact: wqshi.nudt@gmail.com

## Introduction

TF2Exp is a gene-based framework to predict the impact of altered TF-binding events on gene expression levels based on regulatory variants. Using data from lymphoblastoid cell lines, TF2Exp models were applied successfully to predict the expression levels of ~3,100 genes. More details please see the paper in citation.

In this repository, the provided pipeline would use the trained models to predict the impact of variants given by users.

The model training scripts of the TF2Exp paper are located in https://github.com/wqshi/TF2Exp_paper.


## Example

A test run with example data is provided in s_test_tf2exp.sh as following:

```
Rscript s_tf2exp.R --chr_str chr22 \
                   --vcf_file ./data/test_data/test.vcf.gz \
                   --output_dir ./data/test_data \
                   --deepsea_dir ./data/tf2exp_data/deepsea/ \
                   --model_dir ./data/tf2exp_data/models/
```


## Input parameters for TF2Exp pipeline

* --chr_str: chromosome to be processed, e.g. chr22

The script will only process variants in one chromosom each time to reduce memory usage.


* --vcf_file: input vcf file

TF2Exp models only require genotype data in vcf format as input. The vcf file should contain genotype data for different individuals. Either vcf file or gzipped vcf file (end with vcf.gz) is fine. The vcf file should use hg19 coordinates.


* --output_dir: path to output directory

Based on the input variants, TF2Exp will calculate the alteration of key TF-binding events in the TF2Exp model and further predict the alteration of gene expression in each individual. The predicted expression change of each individual is stored in a file named `predicted_exp.chr*.txt` in the output directory specified by `--output_dir`. The `predicted_exp.chr*.txt` presents a data matrix, in which a element is the predicted expression for a gene (indicated by column name) in a individual of input vcf file(indicated by row name). If the variants doesn't impact the expression of a gene across all the individuals, the predicted expression levels would be same for all individuals. This value is the mean expression value of that gene in the training data. GENCODE v19 transcriptome annotations are used for gene identifiers. In addition, altered features of these genes are stored in the `features` subfolder of the specified output directory.


* --deepsea_dir: path to deepsea prediction directory

* --model_dir: path to TF2Exp trained models

In addition, TF2Exp further needs two types of data, pre-compiled variant impact data (given by [DeepSEA](https://www.nature.com/articles/nmeth.3547), `--deepsea_dir` parameter) and trained models for each gene (`--model_dir` parameter). This repository include some example data. The example variant impact data are located in ./data/tf2exp_data/deepsea, and example trained models are in ./data/tf2exp_data/models. A complete set of impact and trained models is available at [tf2exp_data_full.tar.gz](https://zenodo.org/record/1343131). The directory organization is same as ./data/tf2exp_data/. In the complete data set, we have pre-compiled all the variants on chr1-chr22 for 358 European individuals in the 1000 Genomes Project. The pre-compiled set might miss some variants in the user vcf file. Alternatively, user can use scriptes in ./deepsea to predict the impact of all the input variants.


## Environment and packages used in the testing

* R 3.1.3

* OS: CentOS 5

* Other R packages include data.table.1.9.6, RUnit.0.4.29, plyr.1.8.3, stringi.0.5.5, dplyr0.5.0, stringr.1.0.0, futile.logger.1.4.1.


This package has been tested on Ubuntu 14.04.


## To do:
Compile all the variants data in 1000 Genomes project.

## Citation:

Shi, Wenqiang, Oriol Fornes, and Wyeth W. Wasserman. "Altered transcription factor binding events predict personalized gene expression and confer insight into functional cis-regulatory variants." bioRxiv (2017): 228155.






















