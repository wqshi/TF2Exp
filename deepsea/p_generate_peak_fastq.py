import os
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, '%s/expression_var/python/' % home_dir)
import pandas as pd
import p_mymodule as my
from p_project_metadata import *


import argparse
parser = argparse.ArgumentParser(description='Extract the deepsea predictions of one tf.Deepsea prediction is sample based. The output of this script is TF based.')
if __doc__ is None and __name__ == "__main__":
 
    parser.add_argument("--peak_file", help="TF binding regions", default=None)
    parser.add_argument('--vcf_file', help = "vcf file of the chromose", default =None)
    parser.add_argument('--hg19_file', help = "referecne genome fasta file for hg19", default =None)
    parser.add_argument('--tmp_dir', help = "tmp_dir for the intermdiate results", default =None)
    args = parser.parse_args()
 
    peak_file = args.peak_file
    vcf_file = args.vcf_file
    hg19_file = args.hg19_file

    if hg19_file == 'None':
        hg19_file = None
    
    tmp_dir = args.tmp_dir

class chipseq_region():
    peak_file = None
 
    binding_df = None
    
    def __init__(self, file_path, compression=None):
        self.peak_file = file_path
        self.binding_df = None
        self.compression = compression
        self.check_peakMax()
        
    def check_peakMax(self, debug = False):
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()
    
        self.binding_df = pd.io.parsers.read_csv(self.peak_file, sep="\t", header=None, compression = self.compression)
        self.binding_df.columns = ['chr', 'start', 'end', 'name', 'col4', 'col5', 'col6', 'col7', 'col8', 'peakMax']
        self.binding_df.sort(['chr', 'start', 'end', 'peakMax'], inplace = True)
        self.binding_df.index = range(self.binding_df.shape[0])
        if any(self.binding_df.peakMax < 0) or any( self.binding_df.end - self.binding_df.start < self.binding_df.peakMax ):
            logging.info('Wrong state of the Peak max position ')
            self.binding_df.peakMax = (self.binding_df.end.astype(float) - self.binding_df.start)/2
            self.binding_df.peakMax = self.binding_df.peakMax.astype(int)

    def merge_single_file(self, debug = False):
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()

        #self.binding_df = pd.io.parsers.read_csv(self.peak_file, sep="\t", header=None, compression = self.compression)
        #self.binding_df.columns = ['chr', 'start', 'end', 'name', 'col4', 'col5', 'col6', 'col7', 'col8', 'peakMax']
        #self.binding_df.sort(['chr', 'start', 'end', 'peakMax'], inplace = True)
        #self.binding_df.index = range(self.binding_df.shape[0])
        
        
        bed_obj = my.f_pd_to_bed_based_on_file(self.binding_df)
        bed_merge = bed_obj.merge(c=1, o='count')

        merged_df = my.f_bed_to_pd2(bed_merge)
        merged_df.columns = ['chr', 'start', 'end', 'overlap_count']
        self.binding_df = merged_df
        
    def merge_overlapped_peaks(self, debug = False):
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()

        #self.binding_df = pd.io.parsers.read_csv(self.peak_file, sep="\t", header=None, compression = self.compression)
        #self.binding_df.columns = ['chr', 'start', 'end', 'name', 'col4', 'col5', 'col6', 'col7', 'col8', 'peakMax']
        #self.binding_df.sort(['chr', 'start', 'end', 'peakMax'], inplace = True)
        #self.binding_df.index = range(self.binding_df.shape[0])
        #if any(self.binding_df.peakMax < 0) or any( self.binding_df.end - self.binding_df.start < self.binding_df.peakMax ):
        #    logging.info('Wrong state of the Peak max position ')
        #    self.binding_df.peakMax = (self.binding_df.end.astype(float) - self.binding_df.start)/2
        #    self.binding_df.peakMax = self.binding_df.peakMax.astype(int)
        
        bed_obj = my.f_pd_to_bed_based_on_file(self.binding_df)
        bed_merge = bed_obj.merge(c=1, o='count')

        merged_bed = bed_merge.intersect(bed_obj, wao = True)

        merged_df = my.f_bed_to_pd2(merged_bed)
        merged_df.columns = ['chr', 'start', 'end', 'overlap_count' ,'chr2', 'start2', 'end2', 'name', 'col4', 'col5', 'col6', 'col7', 'col8', 'peakMax', 'overlap']
        merged_df['peak_max_new'] = merged_df.start2 + merged_df.peakMax - merged_df.start
        print merged_df.ix[merged_df.overlap_count == 2].head()
        
        assert merged_df.shape[0] == self.binding_df.shape[0], 'Merge error '
        assert all(merged_df.peak_max_new > 0), 'Merge error'
        self.binding_df = merged_df.ix[:, ['chr', 'start', 'end', 'name', 'col4', 'col5', 'col6', 'col7', 'col8', 'peak_max_new'] ]

        
    def split_peaks_with_multiple_peakMax(self, debug = False):
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()

        if self.binding_df is None:
            self.binding_df = pd.io.parsers.read_csv(self.peak_file, sep="\t", header=None)

        self.binding_df.columns = ['chr', 'start', 'end', 'name', 'col4', 'col5', 'col6', 'col7', 'col8', 'peakMax']
        self.binding_df.sort(['chr', 'start', 'end', 'peakMax'], inplace = True)

        self.binding_df.index = range(self.binding_df.shape[0])

        
        #handle the cases where peak max is -1

        #Multiple peakMax in the same region
        alt_data = self.binding_df.copy()
        alt_data['name'] = alt_data.chr + '_' + alt_data.start.map(str) + '_' + alt_data.end.map(str)
        multiple_rows = self.binding_df.index[ (alt_data.name.duplicated()) |  (alt_data.name.duplicated(take_last = True))]
        logging.info('Duplicated rows %s out of %s' %(len(multiple_rows), alt_data.shape[0]))
        for row_id in multiple_rows[range(len(multiple_rows) -1)]:
            #if alt_data.ix[row_id,'start'] == 43044464:
            #    import ipdb; ipdb.set_trace()
                
            if alt_data.ix[row_id, 'name'] == alt_data.ix[row_id + 1, 'name']:
                split_pos = self.binding_df.ix[row_id, 'start'] + int((float(self.binding_df.ix[row_id, 'peakMax']) + self.binding_df.ix[row_id +1, 'peakMax'])/2)
                alt_data.ix[row_id, 'end'] = split_pos
                alt_data.ix[row_id + 1, 'peakMax'] = self.binding_df.ix[row_id + 1, 'peakMax'] + self.binding_df.ix[row_id + 1, 'start'] - split_pos 
                alt_data.ix[row_id + 1, 'start'] = split_pos
        assert all(alt_data.start < alt_data.end), 'Start <= End'
        assert all(alt_data.peakMax >= 0), 'PeakMax > 0'
        self.binding_df = alt_data
        
        
    def bed_trim_binding_regions(self, peak_max_col = 9, debug = False):
        
        #shutil.copyfile(self.peak_file, self.peak_file + '.trim_bak')
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()
            
        if self.binding_df is None:
            self.binding_df = pd.io.parsers.read_csv(self.peak_file, sep="\t", header=None)
            self.check_peakMax()
            
        alt_data = self.binding_df.copy().ix[:,range(3) + [9]]
        alt_data.columns = ['chr', 'start', 'end', 'name']
        #PeakMax might be bigger than 550
        alt_data['alt_start'] = alt_data.start + alt_data.name - 550
        alt_data['alt_end'] = alt_data.start + alt_data.name + 550
        change_start = (alt_data.alt_start > alt_data.start) | (alt_data.name > 550)
        self.binding_df.ix[change_start,1] = alt_data.ix[change_start, 'alt_start']
        self.binding_df.ix[change_start,9] = 550
        self.binding_df.ix[alt_data.alt_end < alt_data.end, 2 ] = alt_data.ix[ alt_data.alt_end < alt_data.end, 'alt_end' ]

        assert all(self.binding_df.ix[:,9] <= 550), 'Peak Max > 550'
        #self.binding_df.to_csv(self.peak_file.replace('.trim_bak',''), header=None, index = None, sep='\t')



    def bed_trim_binding_1k(self, peak_max_col = 9, debug = False):
        
        #shutil.copyfile(self.peak_file, self.peak_file + '.trim_bak')
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()
            
        if self.binding_df is None:
            self.binding_df = pd.io.parsers.read_csv(self.peak_file, sep="\t", header=None)
            self.check_peakMax()
            
        alt_data = self.binding_df.copy().ix[:,range(3) + [9]]
        alt_data.columns = ['chr', 'start', 'end', 'name']
        alt_data['alt_start'] = alt_data.start + alt_data.name - 550
        alt_data['alt_end'] = alt_data.start + alt_data.name + 550

        alt_data.ix[alt_data.alt_start < 0,'alt_end'] = 1100
        alt_data.ix[alt_data.alt_start < 0,'alt_start'] = 1
        
        change_start = self.binding_df.index
        self.binding_df.ix[change_start,1] = alt_data.ix[change_start, 'alt_start']
        self.binding_df.ix[change_start,9] = 550
        self.binding_df.ix[change_start, 2 ] = alt_data.ix[ change_start, 'alt_end' ]

        self.binding_df.sort(['chr', 'start', 'end', 'peakMax'], inplace = True)
        
        assert all(self.binding_df.ix[:,9] <= 550), 'Peak Max > 550'
        assert all(self.binding_df.ix[:,1] >= 0), 'Peak Start > 0'
        
    def bed4_format_peak(self, other_col):
        #import ipdb; ipdb.set_trace()
        self.binding_df = pd.io.parsers.read_csv(self.peak_file, sep="\t", header=None, compression = self.compression).ix[:,range(3) + [other_col]]

        assert all(self.binding_df.ix[:,9] <= 550), 'Peak Max > 550'
        

            
    def vcf_to_tmp_bed(self, vcf_file, tmp_dir, debug = False):
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()
        if vcf_file.endswith('gz'):
            compress_format = 'gzip'
        else:
            compress_format = None
        vcf_df_raw = pd.io.parsers.read_csv(vcf_file, sep="\t", header=None, compression = compress_format).ix[:,0:5]
        vcf_df_raw.columns = ['chr', 'pos', 'name', 'ref', 'alt']

        vcf_df = vcf_df_raw.ix[vcf_df_raw.alt.str.contains('^[ATGCatgc]+$'),:]
        if vcf_df.shape[0] != vcf_df_raw.shape[0]:
            logging.info('Filter %s illegal alt variants' % (vcf_df_raw.shape[0] - vcf_df.shape[0]))
            
        vcf_df['start'] = vcf_df.ix[:,'pos'] -1
        vcf_df['name'] = vcf_df['chr'] + '_' + vcf_df['pos'].map(str) + '_' + vcf_df['ref'] + '_'  +vcf_df['alt']
        tmp_bed_file = tmp_dir + '/' + my.f_generate_tmp_file_name('bed')
        vcf_df.ix[:,['chr', 'start', 'pos', 'name']].to_csv(tmp_bed_file , header = False, index = False ,sep = '\t')
        return tmp_bed_file
        
    def overlap_with_other_bed(self, other_bed, debug = True):
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()
        bed_data = my.f_pd_to_bed_based_on_file(self.binding_df)
        overlap_regions = bed_data.intersect(other_bed, wo=True)
        if overlap_regions.count() == 0:
            raise Exception('Empty overlap between features and vcf file')
        
        overlap_db = my.f_bed_to_pd2(overlap_regions).ix[:, [0,1,2,3,7]]
        print overlap_db.head()
        overlap_db.columns = ['chr','start','end','info','overlap_name']
        assert overlap_db.duplicated('overlap_name').sum() == 0, 'Duplicated vcf positions'
        return overlap_db

    def save_bed(self, output_file):
        self.binding_df.to_csv(output_file, header=False, index=False, sep="\t")
        
    def bed_to_fastq(self, overlap_db, tmp_dir, hg19_file = None, debug = False):
        if f_judge_debug(debug):
            import ipdb; ipdb.set_trace()

        #Get the 1100bp centered at the peak max position
        overlap_db['fastq_start'] = overlap_db['start'].astype(int) + overlap_db['info'].astype(int) - 550
        overlap_db['fastq_end'] = overlap_db['start'].astype(int) + overlap_db['info'].astype(int) + 550
        overlap_db['name'] = overlap_db['overlap_name'] + '_' + overlap_db['fastq_start'].map(str)
        fastq_bed=my.f_pd_to_bed_based_on_file(overlap_db.ix[:,['chr','fastq_start', 'fastq_end', 'name']])
        #fastq_bed = peak_bed.slop(b = 550, genome='hg19')
        if hg19_file is None:
            hg19_file = my.f_get_reference_genome()
            
        fasta = pybedtools.example_filename(hg19_file)
        a = fastq_bed.sequence(fi=fasta, name = True)
        fasta_file= os.path.join(tmp_dir, 'infile.vcf.wt1100.fasta')
        import shutil
        shutil.copyfile(a.seqfn, fasta_file)
        return fasta_file

    def interval_sum(self):
        return (self.binding_df.end - self.binding_df.start).sum()

    def preprocess_1k_peaks(self):
        #import ipdb; ipdb.set_trace()
        self.bed_trim_binding_1k()
        self.merge_overlapped_peaks()
        self.split_peaks_with_multiple_peakMax(debug = False)
        

    
def f_prepare_deepsea_fastq_based_on_vcf(peak_file, vcf_file, tmp_dir, hg19_file=None, debug = False):
    if f_judge_debug(debug):
        import ipdb; ipdb.set_trace()

    tf_region = chipseq_region(file_path = peak_file)
    tf_region.bed4_format_peak(other_col = 9)
    vcf_bed = tf_region.vcf_to_tmp_bed(vcf_file, tmp_dir, debug =False)
    overlap_pd = tf_region.overlap_with_other_bed(vcf_bed, debug = False)
    fastq_file = tf_region.bed_to_fastq(overlap_pd, tmp_dir, hg19_file, debug = False)
    return fastq_file
















import unittest
class TestDatabaseTable(unittest.TestCase):
    def setUp(self):
        a =0

    def test_basic(self):
        import ipdb; ipdb.set_trace()
            
        other_col = 9
        peak_file = '%s/deepsea/tests/data/yy1.sorted.bed' % project_dir        

        chr_str = 'chr22'
        vcf_file = '%s/deepsea/tests/data/%s.merge.head.vcf.gz'%(project_dir, chr_str)
        tmp_dir = '%s/deepsea/tmp/%s/'%(project_dir, my.f_generate_tmp_file_name('t'))
        tmp_dir = '%s/deepsea/tmp/'%(project_dir)
        my.f_ensure_make_dir(tmp_dir)
        
        fastq_file=f_prepare_deepsea_fastq_based_on_vcf(peak_file, vcf_file, tmp_dir)
    
    def test_dnase(self):
        other_col = 9
        

        peak_file = '/homed/home/shi/expression_var/data/raw_data/tf/encode_peaks/processed/uw-gm12878-dnase.narrowPeak'
        #peak_file = '/homed/home/shi/expression_var/data/raw_data/tf/encode_peaks/processed/haib-gm12878-runx3.narrowPeak'
        #peak_file = '/homed/home/shi/expression_var/data/raw_data/tf/encode_peaks/processed/sydh-gm12878-ctcf.narrowPeak'
        chr_str = 'chr22'
        vcf_file = '%s/deepsea/tests/data/%s.merge.head.vcf.gz'%(project_dir, chr_str)
        
        tmp_dir = '%s/deepsea/tmp/'%(project_dir)
        my.f_ensure_make_dir(tmp_dir)
        #import ipdb; ipdb.set_trace()
        fastq_file=f_prepare_deepsea_fastq_based_on_vcf(peak_file, vcf_file, tmp_dir)
        a = 0
    def test_duplicate_pos(self):
        chr_str = 'chr22'
        #import ipdb; ipdb.set_trace()
        peak_file = '/homed/home/shi/expression_var/data/raw_data/tf/encode_peaks/haib-gm12878-stat5a.narrowPeak'
        tmp_dir = '%s/deepsea/tmp/'%(project_dir)
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.merge_overlapped_peaks(debug = False)
        self.assertTrue( (tf_region.binding_df.start == 3590297 ).sum() == 2, 'One case wrong')
        vcf_file = '%s/deepsea/examples/deepsea/example.vcf' % project_dir
        #f_prepare_deepsea_fastq_based_on_vcf(peak_file, vcf_file, tmp_dir, debug =True)
        
    def test_empty_vcf_overlap_with_bed(self):
        #import ipdb; ipdb.set_trace()
        other_col = 9
        peak_file = '/homed/home/shi/expression_var//data/raw_data/tf/encode_peaks/uw-gm12878-dnase.narrowPeak'
        chr_str = 'chr22'
        vcf_file = '%s/deepsea/examples/deepsea/example.vcf' % project_dir
        tmp_dir = '/tmp/tmpegec75'
        my.f_ensure_make_dir(tmp_dir)
        with self.assertRaises(Exception) as context:
            fastq_file=f_prepare_deepsea_fastq_based_on_vcf(peak_file, vcf_file, tmp_dir, debug = False)
        
        self.assertTrue('Empty overlap between features and vcf file' in context.exception)

    def test_ref_allele_not_match(self):
        other_col = 9
        peak_file = '/homed/home/shi/expression_var//data/raw_data/tf/encode_peaks/haib-gm12878-pol2.narrowPeak'

        #import ipdb; ipdb.set_trace()
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.bed_trim_binding_regions()

        peak_file2 = '/homed/home/shi/expression_var//data/raw_data/tf/encode_peaks/haib-gm12878-pol24h8.narrowPeak'
        tf_region = chipseq_region(file_path = peak_file2)
        tf_region.bed_trim_binding_regions()
        #fastq_file=f_prepare_deepsea_fastq_based_on_vcf(peak_file, vcf_file, tmp_dir, debug = False)

    def test_trim_region(self):
        other_col = 9
        peak_file='/homed/home/shi/expression_var/data/raw_data/tf/encode_peaks/AWG/narrowPeaks/wgEncodeAwgDnaseDukeFibroblUniPk.narrowPeak'

        
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.bed_trim_binding_1k(debug = False)


    def test_trim_region_start_error(self):
        other_col = 9
        peak_file='/homed/home/shi/expression_var/data/test_data/python/c.txty'
        
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.bed_trim_binding_1k(debug = True)

        
        print tf_region.interval_sum()
        tf_region.merge_overlapped_peaks()
        print tf_region.interval_sum()

    def test_sum_merge(self):
        other_col = 9
        peak_file='/homed/home/shi/expression_var/data/test_data/python/b.bed'
        
        tf_region = chipseq_region(file_path = peak_file)

        import ipdb; ipdb.set_trace()
        print tf_region.interval_sum()
        tf_region.merge_single_file()
        print tf_region.interval_sum()

        
    def test_split_multiple_peakmax(self):
        #import ipdb; ipdb.set_trace()
        peak_file = '/homed/home/shi/expression_var//data/raw_data/tf/encode_peaks/haib-gm12878-nfic.narrowPeak' #
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.split_peaks_with_multiple_peakMax()

    def test_split_multiple_peakmax_error(self):
        
        peak_file = '/homed/home/shi/expression_var//data/raw_data/tf/encode_peaks/haib-gm12878-pol2.narrowPeak' #
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.split_peaks_with_multiple_peakMax()
        #import ipdb; ipdb.set_trace()
        print tf_region.binding_df.ix[ tf_region.binding_df.start== 43044464,:]
        print tf_region.binding_df.ix[18342,:]
        print my.f_get_reference_genome()

    def test_trim_1k_region(self):
        
        peak_file='/homed/home/shi/expression_var/data/test_data/python/d.bed'
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.preprocess_1k_peaks()
        import ipdb; ipdb.set_trace()
        print tf_region.binding_df.head()
    def test_preprocess_error(self):
        peak_file='./data/raw_data/tf/encode_peaks/AWG//narrowPeaks/wgEncodeAwgTfbsHaibGm12878Atf2sc81188V0422111UniPk.narrowPeak'
        tf_region = chipseq_region(file_path = peak_file)
        tf_region.preprocess_1k_peaks()
        
        print tf_region.binding_df.head()
        
if __name__ == "__main__":

    if __doc__ is not None:
        testname = "test_sum_merge"

        testname = 'test_trim_1k_region'
        testname = 'test_preprocess_error'
        testcase = TestDatabaseTable(testname)
        #testcase.test_sum_merge()
        
        #suite = unittest.TestLoader().loadTestsFromTestCase( TestDatabaseTable )
        unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( testcase )
    else:
        fastq_file=f_prepare_deepsea_fastq_based_on_vcf(peak_file, vcf_file, tmp_dir, hg19_file)








