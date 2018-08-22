
source('s_prepare_predict_data.R')
source('s_predict_expression.R')

t_test_one_gene <- function(){

    target_col= 'RNASEQ'

    gene = "ENSG00000100376.7"
    chr_str = 'chr22'
    test_batch_name = "800samples_regionkeepLow"

    batch_name = '358samples_regionkeepLow'

    output_dir = f_p('./data/%s/', batch_name)
    opt_name = "rm.histone_model.cv.glmnet_add.penalty_rm.YRI_population.None_batch.mode.TF_other.info.corR2keepZeroPopPeerRmdup"
    results_dir = f_p('%s/rnaseq/%s/%s/', output_dir, chr_str, opt_name )

    pred_result =f_load_test_data_and_predict(results_dir, gene, chr_str, test_batch_name, target_col, debug = T)


    checkEquals(pred_result$pred_cor[1,1], 0.2114608, 'Corr Wrong')
    
}


###Prepare the test vcf file
t_prepare_test_data <- function(){
    vcf_source_dir = '/homed/home/shi/expression_var/data/raw_data/wgs/1kg/'
    pred_dir = '/homed/home/shi/expression_var/data/test/predict_dir/'
    col_num = 25
    cmd0 = f_p("zcat %s/chr22.vcf.gz | head -n 1000 | grep -i '^#' | head -10 > %s/test.vcf", vcf_source_dir, pred_dir)
    cmd1 = f_p("zcat %s/chr22.vcf.gz | head -n 1000 | grep -i '^#' | tail -1 | cut -f1-%s >> %s/test.vcf", vcf_source_dir, col_num, pred_dir)
    cmd2 = f_p("zcat %s/chr22.vcf.gz | grep -v '#' | cut -f1-%s >> %s/test.vcf", vcf_source_dir, col_num, pred_dir)
    system(cmd0)
    system(cmd1)
    system(cmd2)
}


#t_prepare_test_data()

t_parse_vcf_and_predict <- function(){

    pred_dir = '/homed/home/shi/expression_var/data/test/predict_dir/'
    output_dir = pred_dir
    chr_str = 'chr22'
    test_vcf_file = f_p('%s/test.vcf', pred_dir)


    ##pred_dir = '/homed/home/shi/expression_var/data/raw_data/GTEx'
                                        #chr_str = 'chr1'
                                        #test_vcf_file = f_p('%s/processed.chr1.vcf.gz', pred_dir)

    deepsea_dir = sprintf('./data/358samples_regionkeepLow/deep_result/all/chrMergeTF/%s/', chr_str)
    opt_name = "rm.histone_model.cv.glmnet_rm.penalty_rm.YRI_population.None_new.batch.358samples.peer_batch.mode.TF_other.info.corR2Nzv5PopPeerRmdup"
    batch_name = '358samples_regionkeepLow'

    results_dir = f_p('./data/%s/rnaseq/%s/%s/', batch_name, chr_str, opt_name )

    f_parse_vcf_and_predict(test_vcf_file, deepsea_dir, results_dir, output_dir, chr_str, debug = T)

}


t_parse_vcf_and_predict_new_dir <- function(){

    pred_dir = './data/test_data'
    output_dir = pred_dir
    chr_str = 'chr22'
    test_vcf_file = f_p('%s/test.vcf.gz', pred_dir)

    deepsea_dir = sprintf('./data/tf2exp_data/deepsea/%s/', chr_str)
    
    results_dir = f_p('./data/tf2exp_data/models/%s/', chr_str)
    f_parse_vcf_and_predict(test_vcf_file, deepsea_dir, results_dir, output_dir, chr_str, debug = F)

}


t_test_GT_GQ_vcf_type <- function(){

    pred_dir = './data/test_data'
    output_dir = pred_dir
    chr_str = 'chr22'
    test_vcf_file = f_p('%s/test_GT_GQ.vcf', pred_dir)

    deepsea_dir = sprintf('./data/tf2exp_data/deepsea/%s/', chr_str)
    
    results_dir = f_p('./data/tf2exp_data/models/%s/', chr_str)
    f_parse_vcf_and_predict(test_vcf_file, deepsea_dir, results_dir, output_dir, chr_str, debug = F)

}


#t_test_GT_GQ_vcf_type()
t_parse_vcf_and_predict_new_dir()
#t_test_one_gene()






  
