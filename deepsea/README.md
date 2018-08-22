## Use DeepSEA to calculate the impact of variants in TF-bound regions


Contact: wqshi.nudt@gmail.com

## Introduction

This directory contains the scripts to calculate the impact of variants on TF binding regions in TF2Exp. The scripts are modified from DeepSEA 0.93. 

## Install and data prepariation

First, DeepSEA should be [downloaded](http://deepsea.princeton.edu/media/code/deepsea.v0.94b.tar.gz) and installed.

Second, create a symbol link of 'Deepsea_install_directory/resources/' to the deepsea subfolder of git repository. Assuming DeepSEA is installed at /home/wqshi/packages/DeepSEA-v0.93/ and we are at the top level of the git repository, the link command would be:


```

ln -s /home/wqshi/packages/DeepSEA-v0.93/resources/ ./deepsea/

#Another necessary step to enable DeepSEA
cp /home/wqshi/packages/DeepSEA-v0.93/deepsea.cpu ./deepsea/

```

Thrid, download and unzip:

* TF-bould regions from the [encode_peaks.tar.gz](https://zenodo.org/record/1343131). 

* prepare the human reference fasta file [GRCh37/hg19](https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz)


## Example to run:

A test run with example data is like following:

```
python2.7 p_tf2exp_rundeepsea.py --vcf_file ./examples/deepsea/chr22.3k.vcf \                   
                   --out_dir out \
                   --tf_dir $PATH_TO_TF_DIR \
                   --hg19_file $PATH_TO_HG19
```

The PATH_TO_TF_DIR should be pointed to the [encode_peaks.tar.gz](https://zenodo.org/record/1343131). In the tf_dir, we expect to find narrowPeak files for different TFs in GM12878, e.g. haib-gm12878-atf2.narrowPeak. The output will be sorted in the ./out directory. For each tf, there will be three files, e.g. atf2.out.diff, atf2.out.ref, atf2.out.evalue. The output directory will be the input directory (--deepsea_dir) for TF2Exp pipeline. 


## Environment and packages used in the testing:

* Python 2.7

* OS: Ubuntu 14.04

* Pandas 0.12.0


## Citation:

Shi, Wenqiang, Oriol Fornes, and Wyeth W. Wasserman. "Altered transcription factor binding events predict personalized gene expression and confer insight into functional cis-regulatory variants." bioRxiv (2017): 228155.






