
import numpy as np
import h5py
import pandas as pd
import gzip
import sys
import os
from subprocess import *
from statsmodels.distributions import ECDF
import joblib


header = np.loadtxt('./resources/predictor.names',dtype=np.str)

if __doc__ is None:
    vcf_prefix = sys.argv[1]
    ref_pred_file = sys.argv[2]
    alt_pred_file = sys.argv[3]
    target_tf = sys.argv[4]
elif False:
    vcf_prefix = 'infile.vcf'
    ref_pred_file = 'tests/t_h5ToOutput/infile.vcf.wt1100.fasta.ref.h5.pred.h5'
    alt_pred_file = 'tests/t_h5ToOutput/infile.vcf.mut1100.fasta.ref.h5.pred.h5'
    target_tf = 'dnase'
else:
    tmp_dir = './tmp/aaa'
    vcf_prefix = '%s/infile.vcf' % tmp_dir
    ref_pred_file = '%s/infile.vcf.wt1100.fasta.ref.h5.pred.h5' % tmp_dir
    
    alt_pred_file = '%s/infile.vcf.mut1100.fasta.ref.h5.pred.h5' % tmp_dir
    target_tf = 'GM12878.CTCF.None'
    

print sys.argv
print vcf_prefix
print ref_pred_file
print 'test'

#import pdb; pdb.set_trace()
#load  column headers
if vcf_prefix.endswith('.vcf'):
    vcf_file = vcf_prefix+'.wt1100.fasta.ref.vcf'
    print vcf_file
    coordata = pd.read_csv( vcf_file, header=None,delimiter='\t').iloc[:,np.asarray([0,1,2,3,4])]
    coordata.columns = ['chr','pos','name','ref','alt']
    coordata.pos=coordata.pos.astype(int)
    wfile1 = open(vcf_prefix+'.out.ref','w')
    #wfile2 = open(vcf_prefix+'.out.alt','w')
    #wfile3 = open(vcf_prefix+'.out.logfoldchange','w')
    wfile4 = open(vcf_prefix+'.out.diff','w')
    wfile6 = open(vcf_prefix+'.out.evalue','w')    

    data1=np.asarray(h5py.File(ref_pred_file,'r')['pred'])
    data2=np.asarray(h5py.File(alt_pred_file,'r')['pred'])
    #compute relative diffrence and absolute difference of chromatin feature prediciton between reference and alternative alleles
    data=np.hstack([np.log2(data2/(1-data2+1e-12))-np.log2(data1/(1-data1+1e-12)),data2-data1])
    data=data[:(data.shape[0]/2),:]/2.0+data[(data.shape[0]/2):,:]/2.0
    
    #compute E-values for chromatin effects
    ecdfs=joblib.load('./resources/ecdfs/ecdf_difflogdiff.rescaled.pkl')
    datae=np.ones((data.shape[0],919))
    for i in range(919):
        #datae[:,i]=1-ecdfs[i](np.abs(data[:,i+919]*data[:,i]))
        datae[:,i]=np.abs(data[:,i+919])*data[:,i]
    #datae[datae==0]=1e-6

    SORTIND=np.arange(data.shape[0])

    #write E-values for chromatin effects
    temp=pd.DataFrame(datae[:,:919])
    temp.columns = header
    import re
    #selected_columns=header #filter(lambda x:re.search(r'GM12878', x), header)
    #selected_columns=filter(lambda x:re.search(r'GM12878[|]%s[|]'% target_tf, x, re.IGNORECASE), header)
    selected_columns=filter(lambda x:re.search(r'%s'% target_tf, x, re.IGNORECASE), header)
    if len(selected_columns) > 1:
        selected_columns = list(set(selected_columns))
    print 'Grep TF: %s' % target_tf, selected_columns
    datae=pd.concat([coordata,temp.ix[:,selected_columns]],axis=1)
    datae=datae.iloc[SORTIND,:]
    datae.to_csv(wfile6,float_format='%.4e')
    
    #write reference allele prediction, alternative allele prediction, relative difference and absolution difference files
    data1=data1[:(data1.shape[0]/2),:]/2.0+data1[(data1.shape[0]/2):,:]/2.0
    #data2=data2[:(data2.shape[0]/2),:]/2.0+data2[(data2.shape[0]/2):,:]/2.0
    temp = pd.DataFrame(data1)
    temp.columns = header
    data1=pd.concat([coordata,temp.ix[:,selected_columns]],axis=1)
    data1=data1.iloc[SORTIND,:]
    data1.to_csv(wfile1, float_format='%.4e')
    #temp = pd.DataFrame(data2)
    #temp.columns = header
    #data2=pd.concat([coordata,temp],axis=1)
    #data2=data2.iloc[SORTIND,:]
    #data2.to_csv(wfile2, float_format='%.4e')

    #temp=pd.DataFrame(data[:,:919])
    #temp.columns = header
    #data3=pd.concat([coordata,temp.ix[:,selected_columns]],axis=1)
    #data3=data3.iloc[SORTIND,:]
    #data3.to_csv(wfile3, float_format='%.4e') #Name it as *.diff. unify the file name.
    
    temp=pd.DataFrame(data[:,919:])
    temp.columns = header
    data4=pd.concat([coordata,temp.ix[:,selected_columns]],axis=1)
    data4=data4.iloc[SORTIND,:]
    data4.to_csv(wfile4, float_format='%.4e')

    wfile1.close()
    wfile6.close()
    wfile4.close()
    





















