
#loire='shi@loire.cmmt.ubc.ca:/homed/home/shi/'

loire='..'

rsync -avu $loire/expression_var/R/s_predict_expression.R .
rsync -avu $loire/expression_var/R/t_predict_expression.R .
rsync -avu $loire/expression_var/R/s_prepare_predict_data.R .
rsync -avu $loire/expression_var/R/*tf2exp* .
rsync -avu $loire/expression_var/R/s_project_funcs.R .
rsync -avu $loire/expression_var/R/r_bedtools.R .
rsync -avu $loire/expression_var/R/s_gene_data_class.R .
rsync -avu $loire/expression_var/R/README.rd .
#rsync -avu $loire/expression_var/R/.git/* ./.git/


###data####
rsync -avu $loire/expression_var/data/tf2exp_data/deepsea/chr22/*.diff ./data/tf2exp_data/deepsea/chr22/
rsync -avu $loire/expression_var/data/tf2exp_data/models/chr22 ./data/tf2exp_data/models/
rsync -avu $loire/expression_var/data/test_data/test.vcf.gz ./data/test_data/
