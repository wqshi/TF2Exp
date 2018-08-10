

library("optparse")

option_list = list(
    make_option(c("--chr_str"),      type="character", default=NULL, help="chromosome, e.g. chr22", metavar="character"),
    make_option(c("--vcf_file"), type="character", default='', help="Path to the input vcf file", metavar="character"),
    make_option(c("--output_dir"), type="character", default='./', help="Path for output directory", metavar="character"),
    make_option(c("--deepsea_dir"), type="character", default='./data/tf2exp_data/deepsea/', help="Path to deepsea prediction directory.", metavar="character"),
    make_option(c("--model_dir"), type="character", default='./data/tf2exp_data/models/', help="Path to TF2Exp trained models", metavar="character")
    
);



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source('s_prepare_predict_data.R')
source('s_predict_expression.R')
library(futile.logger)



deepsea_dir = sprintf('%s/%s', opt$deepsea_dir, opt$chr_str)
model_dir = sprintf('%s/%s', opt$model_dir, opt$chr_str)


f_parse_vcf_and_predict(opt$vcf_file, deepsea_dir, model_dir, opt$output_dir, opt$chr_str)


