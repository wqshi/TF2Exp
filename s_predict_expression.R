

#test_data_raw = final_data[test_samples, input_features ]
#fit = fit

source('s_gene_data_class.R')

f_predict_exp <- function(fit, test_data_raw, target_col = '', debug = F){

    input_features = grep('Intercept', names(fit$key_features), value =T, invert = T)
    
    if (length(setdiff(input_features, colnames(test_data_raw))) == 0){
        if (f_judge_debug(debug)){
            ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        }

        test_data_raw[is.na(test_data_raw)] = 0

        if (is.null(fit$scale_mean)){
            test_data = scale(as.matrix(test_data_raw[,input_features]))
        }else{
            test_data = scale(as.matrix(test_data_raw[,input_features ]), center = fit$scale_mean[input_features], scale = fit$scale_sd[input_features])
        }

        
        
                                        #test_data = (as.matrix(final_data[test_samples, grep('Intercept', names(fit$key_features), value =T, invert = T) ]))
        
        pred_value = glmnet::predict.cv.glmnet(fit$finalModel, newx = test_data )
        
        if (target_col %in% colnames(test_data_raw) & sd(pred_value) != 0){
            pred_cor=cor(pred_value, test_data_raw[, target_col], method = 'spearman')
        }else{
            pred_cor = NA
        }
        return (list(pred_value = pred_value, pred_cor = pred_cor))
        
    }else{
        return (list(pred_value = NA, pred_cor = NA))
    }
}


f_load_test_data_and_predict <- function(results_dir, gene, chr_str, test_batch_name, target_col= '', debug = F){

    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        

    }

    output_file = f_p('%s/%s.enet',  results_dir, gene)
    fit = readRDS(f_p('%s.model', output_file))
    
    feature_df = read.table(f_p('%s.features.gz', output_file), header = T)
    feature_df$loc_tf = str_replace(feature_df$name, '.*[|]', '')
    feature_df$row_id = str_replace(feature_df$name, '[.][.0-9]*$', '')


    head(feature_df)

    pred_gene <- GENE(data = data.frame(), gene_name = str_replace(gene, '[.].*', ''), chr_str = chr_str, batch_name = test_batch_name)
    pred_gene$read_data()
    pred_gene$data$gene = gene


    non_sample_cols =setdiff(colnames(pred_gene$data), pred_gene$get_samples())



    target_data<- pred_gene$data %>%
        dplyr::arrange(gene, -cor ) %>%
        distinct(feature, feature_start, feature_end, gene, .keep_all = T) %>%
        mutate(row_id = str_replace(paste0(gene, '|', type,'.',feature), '[.][.0-9]*$', '')) %>% 
        left_join(x=feature_df[,c('loc_tf','row_id','feature_start', 'feature_end')])


    rownames(target_data) = target_data$loc_tf

    f_assert(length(unique(target_data$loc_tf)) == nrow(feature_df), 'Unique map')

    test_samples = pred_gene$get_samples()

    predict_data = as.data.frame(t(target_data[,test_samples]))
    rownames(predict_data)=test_samples

    predict_data[[target_col]] = t(pred_gene$data[pred_gene$data$feature == 'RNASEQ',test_samples])

    pred_result = f_predict_exp(fit, predict_data, target_col, debug = F)


    
    return (pred_result)
}





























