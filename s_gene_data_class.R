#setwd('~/projects/expression_var/R/')
library(methods)
#source('~/R/s_function.R', chdir = T)
library(futile.logger)
source('s_project_funcs.R')
library(stringr)
library(plyr)
#library(tidyr)
GENE <- setRefClass("GENE",
                     fields = list( data = "data.frame",
                                   gene_name = "character",
                                   batch_name = 'character',
                                   chr_str = 'character'),
                    methods = list(

                        get_samples = function(input_data = NULL, invert = FALSE){
                            if (is.null(input_data)){
                                sample_cols = sort(grep('(NA|HG)[0-9]+', colnames(data),value = T, invert = invert))
                            }else{
                                sample_cols = sort(grep('(NA|HG)[0-9]+', colnames(input_data),value = T, invert = invert))
                            }
                            return (sample_cols)
                        },

                        subset_features_to_snp_contain_region = function(batch_mode, debug = F){
                                        #Subset the features only to the SNP regions
                            if (f_judge_debug(debug)){
                                ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
                                
                            }
                            set_feature_id()
                            if ( !grepl('TFsnpMatch', batch_mode)) return (NULL)

                            if (batch_name %in% c('445samples_sailfish')){
                                other_batch_name = '445samples_snpOnly'
                            }else{
                                        #if (grepl('Conserve', batch_name)){
                                        #    other_batch_name = str_replace(batch_name, 'Conserve[0-9]+', 'Conserve0')
                                        #}else{
                                other_batch_name = f_p('%s_snpOnly', str_replace(batch_name, '(_rareVar|Conserve.*|RareConserve.*)', ''))
                                        #}
                                
                            }
                            snp_gene <- GENE(data = data.frame(), gene_name = gene_name, chr_str = chr_str, batch_name = other_batch_name)
                            snp_gene$read_data()

                            snp_gene$rmdup_based_on_features()
                            snp_gene$set_feature_id()
                            
                            
                            snp_tf_features = grep('enhancer|promoter', rownames(snp_gene$data), value =T)
                            loc_tf_features = grep('enhancer|promoter', rownames(data), value = T)

                            rare_tf_features=setdiff(loc_tf_features, snp_tf_features)
                            if (!grepl('rareVar', batch_name)){
                                f_ASSERT( length(setdiff(f_simple_locs(snp_tf_features), f_simple_locs(loc_tf_features))) == 0, 'SNP TF features are bigger than All TF features' )
                            }
                            flog.info('Remove rare features: %s out of %s', length(rare_tf_features), length(loc_tf_features))
                            
                            selected_data = data[! (rownames(data) %in% rare_tf_features),]
                            data <<- selected_data 
                            
                        },
                        shuffle_rare_features_out_of_snp_contain_region = function(batch_mode, debug = F){
                                        #Subset the features only to the SNP regions
                            if (f_judge_debug(debug)){
                                ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
                                
                            }
                            set_feature_id()
                            if ( !grepl('TFrareShuffle|TFrarePromoter', batch_mode)) return (NULL)

                            if (batch_name %in% c('445samples_sailfish')){
                                other_batch_name = '445samples_snpOnly'
                            }else{
                                        #if (grepl('Conserve', batch_name)){
                                        #    other_batch_name = str_replace(batch_name, 'Conserve[0-9]+', 'Conserve0')
                                        #}else{
                                other_batch_name = f_p('%s_snpOnly', str_replace(batch_name, '(_rareVar|Conserve.*|RareConserve.*)', ''))
                                        #}
                                
                            }
                            snp_gene <- GENE(data = data.frame(), gene_name = gene_name, chr_str = chr_str, batch_name = other_batch_name)
                            snp_gene$read_data()

                                        #snp_gene$rmdup_based_on_features()
                            snp_gene$set_feature_id()

                            snp_gene$print_types()

                            
                            snp_tf_features = grep('enhancer|promoter', rownames(snp_gene$data), value =T)
                            loc_tf_features = grep('enhancer|promoter', rownames(data), value = T)

                            rare_tf_features=setdiff(loc_tf_features, snp_tf_features)
                            if (!grepl('rareVar', batch_name)){
                                f_ASSERT( length(setdiff(f_simple_locs(snp_tf_features), f_simple_locs(loc_tf_features))) == 0, 'SNP TF features are bigger than All TF features' )
                            }
                            set.seed(2)
                            sample_cols = grep('t_',get_samples(), value = T, invert = TRUE)
                            if (length(rare_tf_features) > 0){
                                new_data = data
                                if (batch_mode == 'TFrareShuffle'){
                                    flog.info('Shuffle rare features: %s out of %s', length(rare_tf_features), length(loc_tf_features))

                                    

                                    permutate_names = sample(x =sample_cols, size = length(sample_cols))

                                    
                                        #sample_data = new_data[, sample_cols]
                                        #sample_data[is.na(sample_data)] = 0
                                        #new_data[,sample_cols] = sample_data
                                    new_data[rare_tf_features, sample_cols] = data[rare_tf_features, permutate_names]
                                    
                                }else if(batch_mode == 'TFrarePromoter'){
                                    
                                    rare_enh = grep('enhancer', rare_tf_features, value = T)
                                    flog.info('Enhancer rare features: %s out of %s', length(rare_enh), length(loc_tf_features))
                                    print(table(data$type))
                                    new_data = data[!(rownames(new_data) %in% rare_enh), ]
                                    print(table(new_data$type))

                                }
                                data <<- new_data
                                return  (rare_tf_features)
                             }
                             
                         },

set_na_to_zero = function(){
    new_data = data
    sample_cols = get_samples()
    sample_data = new_data[, sample_cols]
    sample_data[is.na(sample_data)] = 0
    new_data[,sample_cols] = sample_data
    
    data <<- new_data
    
},
rmdup_based_on_features = function(){
    new_data <- data %>% dplyr::arrange(gene, -cor ) %>% distinct(feature, feature_start, feature_end, gene, .keep_all = T) %>% as.data.frame
    flog.info('New: %s out of old %s', nrow(new_data), nrow(data))
    data <<- new_data
},
set_feature_id = function(){
    index = paste(data$type, data$feature, data$hic_fragment_id, data$pair, data$feature_start, data$feature_end ,sep=':')
    duplicated_index = index[duplicated(index)]
    merge_cols = c('feature','feature_start', 'feature_end', 'hic_fragment_id', 'type', 'pair')
    
    f_ASSERT( all(grepl('SNP', duplicated_index)), 'Duplicated snps error')
    flog.info('Number of duplicated SNPs %s', length(duplicated_index))
    data <<- data[!duplicated(data[,merge_cols]),] #Remove the SNPs with same rs names. Potential problem
    rownames(data) <<- paste(data$type, data$feature, data$hic_fragment_id, data$pair, data$feature_start, data$feature_end ,sep=':')
},

get_validation_samples = function(other_info){
    validation_samples = read.table('./data/raw_data/CEU47/validating_samples.txt', header = T)
    if (grepl('ctcf', other_info, ignore.case = T)){
        return (validation_samples$ctcf)
    }else if(grepl('pu1', other_info, ignore.case = T)){
        return (validation_samples$pu1)
    }else{
        return (NULL)
    }
    

},

hic_within_1Mb = function(return_pos = FALSE){
    
    all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes.gene_start',sep='\t', header = TRUE, quote = "")
    row.names(all.entrezgene) = all.entrezgene$ensembl_transcript_id
    loc_gene_name = str_replace(gene_name, '[.][0-9]+', '')
    gene_start = all.entrezgene[loc_gene_name, 'transcript_start'] - 1000000
    gene_end = all.entrezgene[loc_gene_name, 'transcript_end'] + 1000000

    
    hic_ids = as.numeric(data$hic_fragment_id)

    hic_ids[is.na(hic_ids)] = all.entrezgene[loc_gene_name, 'transcript_start'] #For the SNPs only overlap with hic-id.
    
    subset_data = data[ hic_ids > gene_start & hic_ids < gene_end,]
    flog.info('Filtered %s features out of 1Mb regions', nrow(data) - nrow(subset_data))
    data <<- subset_data

    if (return_pos == TRUE){
        return (c(gene_start, gene_end))
    }
    
},
read_test_data = function(test_batch, test_individuals = NULL, debug = FALSE) {
    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
                                 browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
    }

                                        #Seperate part of samples for final validation
    if (!is.null(test_individuals)){
        
        part_data = data[, setdiff(colnames(data), test_individuals)]
        shared_data = data[, intersect(colnames(data), test_individuals)]
        colnames(shared_data) = paste0('t_', colnames(shared_data))
        data <<- cbind(part_data, shared_data)
        return (TRUE)
    }
    

                                        #Read additional data.
    output_dir = f_p('./data/%s/', test_batch)
    expression_file = f_p('%s/rnaseq/%s/%s.txt', output_dir, chr_str, test_gene_name())                             

    if ( !file.exists(expression_file) ){
        flog.info('Missing testing data! %s', expression_file)
        return (FALSE)
    }
    test_data <- read.table(expression_file, header = T, na.strings = 'NA')
    
    sample_cols = get_samples(test_data)
    
    flog.info('Test samples: %s', length(sample_cols))

                                        #colnames(test_data) = paste0('t_', sample_cols)
    if ('cor'  %in% colnames(data)){
        test_data <- test_data[,c('chr', 'start', 'end', 'gene','feature',
                                  'feature_start', 'feature_end', 'hic_fragment_id' ,
                                  'type', 'cor', 'pair', sample_cols )]
    }else{

        test_data <- test_data[,c('chr', 'start', 'end', 'gene','feature',
                                  'feature_start', 'feature_end', 'hic_fragment_id' ,
                                  'type', sample_cols )]
    }

    merge_cols = c('feature','feature_start', 'feature_end', 'hic_fragment_id', 'type', 'pair')
    test_merge_data = test_data[!duplicated(test_data[, merge_cols]),c(merge_cols,sample_cols)]
    test_samples = paste0('t_', sample_cols)
    colnames(test_merge_data) = c(merge_cols, test_samples )
    
    
    merge_data = merge(data, test_merge_data, by = merge_cols, all.x = TRUE)

    f_assert(nrow(merge_data) == nrow(data), 'Merge error')
    
    inserted_data = merge_data[, test_samples]

    duplicated_rows=f_duplicated_cols(data[, merge_cols])
    missing_data = merge_data[rowSums(is.na(inserted_data)) == length(test_samples),merge_cols]

    flog.info('Missing features in the test data: %s', nrow(missing_data))

    data <<- merge_data
    return (TRUE)
},


test_gene_name = function(){
    return (str_replace(gene_name, '.[0-9]+$',''))
},

read_data = function(debug=FALSE) {

    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
        
    }
    
    output_dir = f_p('./data/%s/', batch_name)
    expression_file = f_p('%s/rnaseq/%s/%s.txt', output_dir, chr_str, gene_name)

    
    result <- try(
        data <<- read.table(f_p('%s.gz', expression_file), header = T, na.strings = 'NA'),    
        silent = TRUE
        
    )
    if (class(result) == "try-error"){

        cat(f_p('Read %s.gz failed ', expression_file), '\n')
        data <<- read.table(expression_file, header = T, na.strings = 'NA')    
        
    }
    
    
    sample_cols = get_samples()
    flog.info('Number of samples: %s, features %s', length(sample_cols), nrow(data))


    if ('cor'  %in% colnames(data)){
        data <<- data[,c('chr', 'start', 'end', 'gene','feature',
                         'feature_start', 'feature_end', 'hic_fragment_id' ,
                         'type', 'cor', 'pair', sample_cols )]
    }else{

        data <<- data[,c('chr', 'start', 'end', 'gene','feature',
                         'feature_start', 'feature_end', 'hic_fragment_id' ,
                         'type', sample_cols )]
    }               
    
},                         

change_expression =function(new_batch_name, batch_mode = 'TF', debug = FALSE) {
    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }
                                        #Change the mRNA values to another normalization method,
                                        #e.g. quantile to peer normalization.
    if (!is.null(new_batch_name)){
        if (batch_mode != 'random'){
            new_data = read.table(f_p('./data/%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', new_batch_name), header = TRUE)
        }else{
            new_data = shuffle_gene_expression(new_batch_name)
        }
        rownames(new_data) = new_data$gene
        samples = get_samples()

                                 mRNA_index = which(data$type == 'gene')
                                 name_shared = intersect(samples , colnames(new_data))
                                 name_diff = setdiff(samples , colnames(new_data))
                                 flog.info('Set diff: %s', paste(name_diff, sep = '', collapse = '.'))
                                 if (grepl('randomN[0-9]+', batch_mode)){
                                     seed_num=f_extract_pattern_nums(batch_mode,'randomN')
                                     set.seed(seed_num)
                                     flog.info('Set random seed: %s', seed_num)
                                     permutate_names = sample(x =name_shared, size = length(name_shared))
                                     data[mRNA_index, name_shared] <<- new_data[gene_name, permutate_names]
                                 }else{
                                     data[mRNA_index, name_shared] <<- new_data[gene_name, name_shared]
                                 }
                                 

                                 
                                 
                             }else{
                                 flog.info('New batch_name is null: %s', new_batch_name)
                             }
                         },

                         add_snp_location = function(loc_batch_name)
                         {
                             if(!grepl('SNP', loc_batch_name, ignore.case = T)){
                                 return (NULL)

                             }


                             promoter_hic = unique(subset(data, type %in% c( 'promoter'))$hic_fragment_id)

                             snp_data = data %>% filter(type == 'SNP')

                             promoter_snps = data$type == 'SNP' & data$hic_fragment_id %in% promoter_hic
                             enhancer_snps = (data$type == 'SNP' &!( data$hic_fragment_id %in% promoter_hic))

                             new_data = data
                             new_data[promoter_snps,'feature'] = paste0(data[promoter_snps,'feature'],'.pro')
                             new_data[enhancer_snps,'feature'] = paste0(data[enhancer_snps,'feature'],'.enh')

                             data <<- new_data
                         },


                         shuffle_gene_expression = function(loc_batch_name){
                             shuffle_file = f_p('./data/%s/rnaseq/GEUVADIS.GeneRandom.DATA_MATRIX', loc_batch_name)
                             if(!file.exists(shuffle_file)){
                                 loc_data = read.table(f_p('./data/%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', loc_batch_name), header = TRUE)
                                 set.seed(897)
                                 permutate_names = sample(x =1:nrow(loc_data), size = nrow(loc_data))
                                 random_data = loc_data
                                 random_data$gene = loc_data[permutate_names, 'gene']
                                 write.table(random_data, file = shuffle_file, quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)
                                 flog.info('Write random expression file: %s', shuffle_file)
                             }else{
                                 random_data = read.table(shuffle_file, header = TRUE)
                             }
                             return (random_data)
                         },
                         print_types = function(){
                             print(table(data$type))

                         },


                       subset_snps_in_tf_regions = function(batch_mode, debug = F){
                           if (f_judge_debug(debug)){
                               ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
                               
                           }
                           
                           if ( !grepl( 'SNPinTF', batch_mode)) return (NULL)
                           
                           regulatory_regions = subset(data, type == 'promoter' | type == 'enhancer')[, c('chr', 'feature_start', 'feature_end', 'feature')]
                           tf_binding_regions = regulatory_regions[grep('H3K', regulatory_regions$feature, invert = T),]

                           SNP_regions = subset(data, type == 'SNP')[, c('chr', 'feature_start', 'feature_end', 'feature')]
                           
                           tf_bed = makeGRangesFromDataFrame(tf_binding_regions)
                           SNP_bed = makeGRangesFromDataFrame(SNP_regions)

                           matches = as.data.frame( findOverlaps(tf_bed, SNP_bed) )
                           removed_SNPs = setdiff(rownames(SNP_regions), rownames(SNP_regions[unique(matches[,2]),]))

                           if (batch_mode == 'randomSNPinTF'){
                               removed_SNPs = sample(rownames(SNP_regions), size = length(removed_SNPs))
                           }

                           flog.info('%s of %s SNPs are filtered out of the TF binding regions', length(removed_SNPs) , nrow(SNP_regions))
                           data <<- data[setdiff(rownames(data), removed_SNPs),]
                           
                       }
                         

                         
                         
                     ))

.DollarNames.GENE <- function(x, pattern){
    my_funcs = grep('callSuper|show|copy|export|field|getClass|getRefClass|import|trace|init|usingMethods',getRefClass(class(x))$methods(), value = T, invert = T) 
    grep(pattern, c( my_funcs, names(getRefClass(class(x))$fields())), value=TRUE)
}

f_simple_locs <- function(input_features){

    #simple_locs<-ldply( str_split(input_features, ':') ) %>%
    #unite(col=new, c(V1,V2,V5,V6),sep=':') %>% select(new) %>%
    #unique %>% unlist

    simple_locs<-ldply( str_split(input_features, ':') ) %>%
    mutate(new=paste(V1,V2,V5,V6,sep=':')) %>% select(new) %>%
    unique %>% unlist

    return (simple_locs)

}


library(RUnit)
options(run.main=FALSE)
if (getOption('run.main', default=TRUE)) {
    runTestFile('./test/t_gene_data_class.R',testFuncRegexp = '^t_.*')
}



#runTestFile('./test/t_gene_data_class.R',testFuncRegexp = '^t_.*')
