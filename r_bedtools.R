#setwd('~/projects/expression_var/R/')
library(methods)
#source('~/R/s_function.R', chdir = T)
library(futile.logger)
source('s_project_funcs.R')

options(scipen=999)


t_test_bedtools <- function(){
    input_A = data.frame(chr = c('chr1', 'chr1', 'chr2'),
                         start = c(1000, 100000, 100020),
                         end =c(2000, 200000, 200020),
                         name = c('g1','g2', 'g3'))


    input_B = data.frame(chr = c('chr1', 'chr1', 'chr2'),
                         start = c(3000, 200000, 200020),
                         end =c(4000, 200001, 200021),
                         name = c('g4','g5', 'g6'))

    fun = 'intersect'
    paras = '-wao'

    f_bedtools( input_A, input_B, fun = 'intersect', paras = '')
    
}



f_df_to_bed <- function(loc_file, input_loc){
    write.table(input_loc, file = loc_file, quote = F, sep = '\t', row.names = F, col.names = F)
}

require(stringi)

f_bedtools <- function(input_A, input_B, fun = 'intersect', paras = '', debug = FALSE){
    if( f_judge_debug(debug) ){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }
    
    my_tmpdir = f_p('tmp.%s', stri_rand_strings(n=1, length=8, pattern="[A-Za-z0-9]"))
    dir.create(my_tmpdir)

    A_file = input_A
    if (class(input_A) == 'data.frame'){
        A_file = tempfile(pattern = 'tmp.bed.', tmpdir = my_tmpdir)
        f_df_to_bed(A_file, input_A)
    }


    B_file = input_B
    if (class(input_B) == 'data.frame'){
        B_file = tempfile(pattern = 'tmp.bed.', tmpdir = my_tmpdir)
        f_df_to_bed(B_file, input_B)
    }

    cmd = f_p('bedtools %s -a %s -b %s %s | tr -s \'\t\' \'\t\'', fun, A_file, B_file, paras)
    print(cmd)
    library(data.table)
    df <- fread(cmd)
    unlink(my_tmpdir, recursive = T)
    return (df)
}



f_extract_bed_from_features <- function(feature_df, chr_str =NULL, pattern = 'promoter|enhancer'){
    if (!is.null(chr_str)){
        feature_df$chr = chr_str
    }
    PE_df = feature_df[grep(pattern, feature_df$name), c('chr', 'feature_start', 'feature_end', 'name') ]
    colnames(PE_df ) = c('chr', 'start', 'end', 'name')
    return (PE_df)
}







