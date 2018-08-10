


#setwd('../R/')
source('s_project_funcs.R')
source('r_bedtools.R')
source('s_predict_expression.R')
library(stringr)
options(warn=2)

f_genotype_to_numbers <- function(vcf_data){

    het_df = vcf_data
    genotype_df = het_df
    het_df[,] = 0
    
    for(loc_indiv in colnames(vcf_data)){
        het_selection = grepl('(0.1|1.0)', genotype_df[, loc_indiv])
        het_df[het_selection, loc_indiv] = 0.5
        homo_selection = grepl('1.1', genotype_df[, loc_indiv])
        het_df[homo_selection, loc_indiv] = 1
    }

    return (list(het_df=het_df, genotype_df = genotype_df))

}



f_transform_tf_to_deepsea_file_name <- function(input_tfs){
    if (any(grepl('ENSG', input_tfs))){
        input_tfs = str_replace(input_tfs, 'ENSG.*[|]', '')
    }
    tmp_name = str_replace(str_replace(input_tfs, 'promoter.|enhancer.', ''), pattern = '[.][0-9]+$','')
    tmp_name=str_replace(tmp_name, 'USF$', 'USF.1')
    tmp_name=str_replace(tmp_name, 'EGR$', 'EGR.1')
    tmp_name=str_replace(tmp_name, 'PU$', 'PU.1')
    tmp_name = str_replace(tolower(tmp_name), '[.]', '')
    return (tmp_name)
}


f_get_vcf_header <- function(input_vcf_file){
    library(data.table)
    if (grepl('.*vcf.gz', input_vcf_file)){
        header_cmd = f_p("zcat %s | head -n 10000 | grep -i '^#' | tail -1", input_vcf_file)
    }else{
        header_cmd = f_p("cat %s | head -n 10000 | grep -i '^#' | tail -1", input_vcf_file)

    }
    
    vcf_colnames = colnames(fread(header_cmd))
    return (vcf_colnames)
}




f_process_one_tf_key_features <- function(diff_file, feature_vcf_combine, vcf_samples, debug){

    if(f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
        
    }
    
    target_tf = str_replace( basename(diff_file), '.out.diff', '')

    impact_data_raw = read.table(diff_file, header = T, sep = ',', row.names = 1)
    impact_data = impact_data_raw[, c('chr', 'pos', 'name', 'ref', 'alt', grep('GM12878',colnames(impact_data_raw), value = T))]
    
    if(ncol(impact_data) > 6){        
        flog.warn('Impact file has >6 columns')
        impact_data = impact_data[,1:6]
    }
    colnames(impact_data) = c('chr', 'pos', 'name', 'ref', 'alt', 'tf_impact')

    

    
    #input_data = feature_vcf_combine %>% select(-one_of(c('qual', 'filter','info', 'format', 'id', 'overlap'))) %>% 
    #    filter(grepl(target_tf, feature_vcf_combine$tf_name, ignore.case = T))

    input_data = feature_vcf_combine %>% select(one_of(c('chr', 'tf_start', 'tf_end', 'tf_name', 'pos', 'ref', 'alt', vcf_samples))) %>% 
        filter(grepl(target_tf, feature_vcf_combine$tf_name, ignore.case = T))
    
    head(input_data)

    merge_data <- input_data %>% 
        left_join(impact_data, by =c('chr', 'pos', 'ref', 'alt'))


    na_impact_genotypes <- merge_data %>% filter(is.na(tf_impact)) %>% dplyr::select(one_of(vcf_samples))
    if (nrow(na_impact_genotypes) > 0) 
        f_assert(all(na_impact_genotypes == 0), sprintf('No match tf impact of the genotype %s out of %s', nrow(na_impact_genotypes), nrow(input_data) ))

    
    merge_data[,vcf_samples] = merge_data[, vcf_samples] * merge_data$tf_impact
    #tf_alteration <- merge_data %>% filter(!is.na(tf_impact)) %>% group_by(chr, tf_start, tf_end, tf_name)  %>%
    #    summarise_at(.cols = vcf_samples, .funs = c(sum = "sum")) %>% as.data.frame()

    vcf_samples_new = str_replace(vcf_samples,'-', '__')

    new_col_name = str_replace(colnames(merge_data),'-', '__')

    colnames(merge_data) = new_col_name

    tf_alteration <- merge_data %>% filter(!is.na(tf_impact)) %>% group_by(chr, tf_start, tf_end, tf_name)  %>%
        summarise_at(.cols = vcf_samples_new, .funs = c(sum = "sum")) %>% as.data.frame()

    colnames(tf_alteration) = c('chr', 'tf_start', 'tf_end', 'tf_name', vcf_samples) 

    return (tf_alteration)


}





f_key_features_alteration_given_vcf <- function(test_vcf_file, deepsea_dir, key_features, debug = T){

    if(f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
        browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
        
    }
    ##Read feature data#####
    if (class(key_features) == 'character'){
        feature_data = read.table(file = key_features, header = T)
    }else{
        feature_data = key_features
    }
    
    feature_bed = feature_data[,c('chr','feature_start', 'feature_end', 'tf')] %>%
    filter(tf != '(Intercept)') %>% unique %>% arrange(chr, feature_start)
    head(feature_bed)


    ##Read the vcf data
    vcf_colnames = f_get_vcf_header(test_vcf_file)
    match_index = grep("(chrom|start|end|chr|pos|id|ref|alt|qual|filter|info|format)", vcf_colnames, ignore.case = T)
    vcf_samples = grep("(chrom|chr|start|end|pos|id|ref|alt|qual|filter|info|format)", vcf_colnames, ignore.case = T, value = T, invert = T)
    vcf_colnames[match_index] = tolower(vcf_colnames[match_index])
    vcf_colnames[1] = 'chr_vcf'



    ##Overlap the vcf data with TF key features
    out = f_bedtools(input_A = feature_bed, input_B = test_vcf_file, fun = 'intersect', paras = '-wo -sorted', debug = F)
    out = as.data.frame(out)
    colnames(out) = c('chr', 'tf_start', 'tf_end', 'tf_name', vcf_colnames, 'overlap')

    out[vcf_samples]= f_genotype_to_numbers(out[vcf_samples])$het_df


    ##Create TF impact matrix??
    pattern_tf = paste(unique(out$tf_name), collapse = '|')
    file_list = grep(f_p('.*(%s).*diff$', pattern_tf), list.files(deepsea_dir, full.names = T), value = T, ignore.case = T)

    f_assert(length(unique(out$tf_name)) == length(file_list),'Missing TFs')
    diff_file = file_list[1]
    tf_alteration = data.frame()
    for (diff_file in file_list){
        flog.info('Process %s', diff_file)

        if (grepl('.*dnase.*xxx', diff_file)){
##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
            

        }
        one_tf_impact  <- f_process_one_tf_key_features(diff_file, out, vcf_samples, debug = debug)
        tf_alteration = rbind(tf_alteration, one_tf_impact)
    }

    return (list(tf_alteration = tf_alteration, feature_bed = feature_bed, vcf_out = out))

}



f_parse_vcf_and_predict <- function(test_vcf_file, deepsea_dir, results_dir, output_dir, chr_str, debug = F){


    if(f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
        
    }

    ##deepsea_dir = sprintf('%s/%s/', pred_dir, chr_str)
    feature_file = f_p('%s/features', results_dir)
    feature_data = read.table(file = feature_file, header = T)


    feature_data$chr = chr_str
    feature_data$loc_tf = str_replace(feature_data$name, 'ENSG.*[|]', '')
    feature_data$tf = f_transform_tf_to_deepsea_file_name(feature_data$loc_tf)

    feature_data = subset(feature_data, !is.na(feature_data))


    key_features = feature_data[,c('chr','feature_start', 'feature_end', 'tf')] %>% filter(!is.na(feature_start))
    process_list <- f_key_features_alteration_given_vcf(test_vcf_file, deepsea_dir, key_features, debug = F)

    
    ##Write the overlapped vcf variants for DeepSEA to process
    if (FALSE){
        vcf_data = process_list$vcf_out[,c('chr_vcf', 'end', 'id', 'ref', 'alt')]
        write.table(vcf_data, sprintf('%s/chr1.vcf', output_dir), quote = F, sep = '\t', row.names = F, col.names = F)
    }


    tf_alteration = process_list$tf_alteration
    tmp_cols = str_replace(colnames(tf_alteration), '_sum', '')
    tmp_cols[2:4] = c('feature_start', 'feature_end', 'feature')
    colnames(tf_alteration) = tmp_cols



    table(tf_alteration$feature)
    sample_cols = colnames(tf_alteration)[5:ncol(tf_alteration)]

    
    gene_list = unique(feature_data$gene)
    #Predict expression for each gene
    #
    loc_gene = gene_list[1]
    pred_df = data.frame( name =  sample_cols)

    closeAllConnections()

    for (loc_gene in gene_list){


        #flog.debug('Process %s', loc_gene)
        
        
        
        output_file = f_p('%s/%s.enet',  results_dir, loc_gene)

        result <-try(
        fit<- readRDS(f_p('%s.model', output_file)),
        silent = T
        )

        if (class(result) == "try-error"){
            
            next
        }else{

            flog.info('Process %s', loc_gene)
            
        }
        feature_df = read.table(f_p('%s.features.gz', output_file), header = T)

       
        feature_df$feature = f_transform_tf_to_deepsea_file_name(feature_df$name)
        feature_df$loc_tf = str_replace(feature_df$name, '.*[|]', '')

        table(feature_df$feature)

        table(tf_alteration$feature)

        target_data<-feature_df[,c('feature','feature_start', 'feature_end','name', 'loc_tf')] %>%
            left_join(tf_alteration, by = c("feature", "feature_start", "feature_end"))

        

        if (!all(rowSums(is.na(target_data[,sample_cols])) == length(sample_cols))){

            flog.debug('None empty prediction for %s', loc_gene)

        }


        rownames(target_data) = target_data$loc_tf
        
        predict_data = as.data.frame(t(target_data[,sample_cols]))
        rownames(predict_data)=sample_cols

        non_zero_cols=grep('Intercept',names(fit$key_features)[fit$key_features != 0],invert = T, value = T)

        non_zero_cols
        
                                        #f_assert(all(!is.na(predict_data[,non_zero_cols])), 'Key features are not zero')
        target_col = ''
        pred_result = f_predict_exp(fit, predict_data, target_col, debug = F)

        

        pred_df[loc_gene] = pred_result$pred_value[sample_cols,1]

        
        #flog.info('Perf: %s', pred_result$pred_value[1,1]    )
    }


    sum(apply(pred_df[,2:ncol(pred_df)],2,sd, na.rm=T) != 0)

    rownames(pred_df) = pred_df$name 
    pred_df$name = NULL

    write.table(pred_df, sprintf('%s/pred_full.txt', output_dir), quote = F, sep = '\t')

}

















