
import os
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
project_dir = '%s/expression_var/' % home_dir
deepsea_dir = '%s/deepsea/' % project_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, deepsea_dir)
import pandas as pd
import p_mymodule as my
from p_pd import *
import pybedtools
import argparse
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
import re

import p_region_table as p_region_table
reload(p_region_table)
from p_region_table import *


hic_file = '%s/data/raw_data/hic/high_interaction_hic.addChr.txt' % project_dir
#Change the default id file to new file.
hic_id_file = '%s/data/raw_data/hic/fragmentStartEnd.addChr.addMissing.txt' % project_dir
#hic_id_file = '%s/data/raw_data/hic/fragmentStartEnd.addChr.txt' % project_dir
#gene_expression_file = '%s/data/rnaseq/transcript_data.bed' % project_dir

def log(title,content):
    logging.debug(title + ":")
    logging.debug(content)


def f_get_batch_output_dir(batch_name):
    return '%s/data/%s/' % (project_dir, batch_name)


def f_get_tf_peak_list(project_dir, version = 'processed', tf_dir=None):

    if tf_dir is None:
        tf_dir = '%s/data/raw_data/tf/encode_peaks/%s/' % (project_dir, version)
    
    peak_list_raw = my.f_shell_cmd( "find %s -name '*gm12878-*.narrowPeak'"%(tf_dir), quiet = True).split('\n')
    black_list = my.grep_list(".*(--|Rep[1-9]|-myc|xyy1|test|pax5n19|embl|encode-)", peak_list_raw)
    duplicate_list =['uta-gm12878-ctcf.narrowPeak', 'uw-gm12878-ctcf.narrowPeak', 'sydh-gm12878-yy1.narrowPeak',
                    'sydh-gm12878-rad21.narrowPeak', 'haib-gm12878-p300.narrowPeak', 'ut-gm12878-cmyc.narrowPeak',
                    'haib-gm12878-pol24h8.narrowPeak', 'sydh-gm12878-pol2.narrowPeak', 'uta-gm12878-pol2.narrowPeak']
    peak_list = list(set(peak_list_raw) - set(['']) - set(black_list) - set( my.grep_list( '.*(%s)'%'|'.join(duplicate_list), peak_list_raw)))
    logging.info('Length of peaks: %s' % len(peak_list))
    return peak_list

chr_list = []
for i in range(1,23):
    chr_list.append('%s' % i)
#for sample_id in sample_list:
#    print '===========%s============' % sample_id
#    #my.f_ensure_make_dir('%s/chr_vcf_files/%s/' % (kg_dir, sample_id))
#    for chr_num in chr_list:
#       bcf_subset_cmd = "~/packages/bcftools/bcftools-1.2/bcftools view -c1 -Ov -s %s %s | grep -v -i 'MULTI_ALLELIC\|ALU' | gzip > %s/chr_vcf_files/chr%s/%s.vcf.gz " % (sample_id, node_chr_vcf, output_dir, chr_num, sample_id)
#        my.f_shell_cmd(bcf_subset_cmd)
chr_list.append('X')
chr_list.append('Y')
#chr_list = ['22', '21', '20'] # It has to be str not int.


def f_check_loc_cols(input_data):
    
    matched_cols=my.grep_list('start|end|pos', input_data.columns.tolist())
    for loc_col in matched_cols:
        if input_data[loc_col].dtype != 'int64':
            logging.error('Error data type in %s', loc_col)
            input_data[loc_col] = input_data[loc_col].astype(int)
    return input_data


def f_stop_script():
    import sys
    sys.exit("Stop the stript on purpose")


def f_get_peak_file_df_rmdup(project_dir, version='', file_path = None, tf_dir = None):
    if file_path != None or version == '1k_extension':
        file_path = './data/raw_data/tf/encode_peaks/AWG/selected_tf_list.txt'
        peak_file_df_rmdup = pd.read_csv(file_path, header=0,sep='\t',)
        peak_file_df_rmdup['deepsea_tf'] = peak_file_df_rmdup.deepsea.str.replace('[|]','.')
        peak_file_df_rmdup['file_path'] = peak_file_df_rmdup.file_path.str.replace('narrowPeaks','processed')
        peak_file_df_rmdup.index = peak_file_df_rmdup.deepsea_tf
    
        return peak_file_df_rmdup
    
    peak_list = f_get_tf_peak_list(project_dir, version, tf_dir)
    tf_list = [[re.split('[.]|-',os.path.basename(peak_file))[2], os.path.getsize(peak_file)] for peak_file in peak_list ]
    peak_file_df = pd.DataFrame(data = tf_list, columns = ['tf', 'file_size'])
    peak_file_df['file_path'] = peak_list

    ##Keep the large peak file.
    peak_file_df_rmdup = peak_file_df.sort(columns = 'file_size', ascending = False ).drop_duplicates(['tf'])

    peak_file_df_rmdup['deepsea_tf'] = peak_file_df_rmdup.tf
    peak_file_df_rmdup.index = peak_file_df_rmdup.tf
    
    peak_file_df_rmdup.ix['pu1','deepsea_tf'] = 'pu.1'
    peak_file_df_rmdup.ix['pax5c20','deepsea_tf'] = 'pax5-c20'
    peak_file_df_rmdup.ix['pax5','deepsea_tf'] = 'pax5-n19'
    peak_file_df_rmdup.ix['pax5c20','deepsea_tf'] = 'pax5-c20'
    peak_file_df_rmdup.ix['egr1','deepsea_tf'] = 'egr-1'
    peak_file_df_rmdup.ix['nfyb','deepsea_tf'] = 'nf-yb'
    peak_file_df_rmdup.ix['usf1','deepsea_tf'] = 'usf-1'
    peak_file_df_rmdup.ix['cmyc','deepsea_tf'] = 'c-myc'
    peak_file_df_rmdup.ix['cfos','deepsea_tf'] = 'c-fos'
    peak_file_df_rmdup.ix['nfya','deepsea_tf'] = 'nf-ya'
    peak_file_df_rmdup.ix['nfe2','deepsea_tf'] = 'nf-e2'

    missing_tfs = set(['znf384','cdp','irf3', 'spt20','mafk','gcn5','erra','srebp2','srebp1']).intersection(set(peak_file_df_rmdup.tf))
    peak_file_df_rmdup.ix[list(missing_tfs), 'deepsea_tf'] = None

    return peak_file_df_rmdup.ix[~peak_file_df_rmdup.deepsea_tf.isnull(),:]
    #return peak_file_df_rmdup


def f_judge_debug(debug):
    import socket
    server_name = socket.gethostname()
    return debug == True and  server_name in ['loire', 'wqshi']

def f_add_break_point():
    import socket
    server_name = socket.gethostname()
    if server_name == 'wqshi':
        import pdb; pdb.set_trace()
    else:
        import ipdb; ipdb.set_trace()

def f_add_debug_point_in_loire(debug = False):
    if f_judge_debug(debug):
        import ipdb; ipdb.set_trace()

def f_add_suffix_on_duplicates(dup_list):
    from collections import Counter
    counter = Counter()
    deduped = []
    for name in dup_list:
        new = name  + ('.' + str(counter[name])) if counter[name] else name
        counter.update({name: 1})
        deduped.append(new)
    return deduped


def f_sync_scripts_to_run_server(target_server = 'clust'):
    import subprocess
    print ''
    print 'Sync scripts ...........'
    print ''

    if target_server == 'clust':
        target_str = 'shi@clustdell.cmmt.ubc.ca:/home/shi/'
    else:
        target_str = 'wenqiang@orcinus.westgrid.ca:/home/wenqiang/'
    
    
    rsync_cmd1  = "rsync -rav --include '*.py' --exclude '*' /homed/home/shi/projects/expression_var/python/ %s/projects/expression_var/python/" % target_str
    subprocess.call(rsync_cmd1,shell=True)









