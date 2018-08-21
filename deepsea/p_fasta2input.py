#Example usage: python 1_fasta2input.py infile.fasta window_size

#Convert sequences to hdf5 format. Used for processing variant format.
#It only supports fasta file produced by 0_coor2fasta.R, which includes ref and alt
#allele information in sequence names. 

from Bio import SeqIO
import numpy as np
import sys
import math
from os.path import basename,dirname
import os
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, '%s/expression_var/python/' % home_dir)
import pandas as pd
import p_mymodule as my
from p_project_metadata import *



writeFasta=True
writeh5=True
writebed=True
writevcf=True

#if vcf_original_allele_check is False, always use the ref and alt alleles user specified. 
#Otherwise ref allele must match the reference genome
vcf_original_allele_check=True

import argparse
parser = argparse.ArgumentParser(description='Extract the deepsea predictions of one tf.Deepsea prediction is sample based. The output of this script is TF based.')
if __doc__ is None:
    parser.add_argument("--fasta_file", help="fasta file", default=None)
    args = parser.parse_args()
    input_fasta = args.fasta_file
else:
    input_fasta = '%s/deepsea/tests/data/input.wt.fasta' %(project_dir)
    #input_fasta = '/tmp/tmpQ7pYXi/infile.vcf.wt1100.fasta'
#mutpos=inputwindow/2-1


fasta_sequences = SeqIO.parse(open(input_fasta),'fasta')


np.random.seed(1)




#write h5 file

def writeh5(seqs,filename,offset_values=None):
    seqsnp=np.zeros((len(seqs),4,1000),np.bool_)

    mydict={'A':np.asarray([1,0,0,0]),'G':np.asarray([0,1,0,0]),'C':np.asarray([0,0,1,0]),'T':np.asarray([0,0,0,1]),'N':np.asarray([0,0,0,0]),'H':np.asarray([0,0,0,0]),'a':np.asarray([1,0,0,0]),'g':np.asarray([0,1,0,0]),'c':np.asarray([0,0,1,0]),'t':np.asarray([0,0,0,1]),'n':np.asarray([0,0,0,0])}
    n=0
    for line in seqs:
        cline = line[ int(math.floor(( (len(line)-1000)/2.0))):int(math.floor(len(line)-(len(line)-1000)/2.0))]
        for c,i in zip(cline,range(len(cline))):
            try:
                seqsnp[n,:,i]=mydict[c]
            except:
                print line
                print c,i
                raise Exception('Error in dict')
        n=n+1
    
    #get the complementary sequences
    dataflip=seqsnp[:,::-1,::-1];
    seqsnp=np.concatenate([seqsnp, dataflip],axis=0)


    seqsnp = seqsnp.astype(np.uint8)
    f=h5py.File(filename,'w')
    f.create_dataset('testxdata', data= seqsnp,compression="gzip")
    f.close()

    
    
seqs=[str(fasta.seq) for fasta in fasta_sequences]

oris=[]
muts=[]
chrs=[]
poss=[]
annos=[]
names=[]
mutpos=[]
bed_start=[]
fasta_sequences = SeqIO.parse(open(input_fasta),'fasta')

#import pdb; pdb.set_trace()

for fasta in fasta_sequences:
    anno = fasta.name.split('_')
    annos.append(fasta.name)
    oris.append(anno[2])
    muts.append(anno[3])
    chrs.append(anno[0])
    bed_start.append(anno[4])
    #print int(anno[4])
    #chr, variant position (1-based), ref, alt, fastq_start.
    #print anno
    
    mutpos.append( int(anno[1]) - 1 - int(anno[4]) )

    poss.append(int(anno[1]))
    if len(anno)>5:
        names.append('_'.join(anno[5:]))

oris=np.asarray(oris)
muts=np.asarray(muts)
chrs=np.asarray(chrs)
seqs=np.asarray(seqs)
poss=np.asarray(poss)
annos=np.asarray(annos)
names=np.asarray(names)
mutpos = np.asarray(mutpos)
bed_start = np.asarray(bed_start)
print len(mutpos)
print len(oris)

seqsmut=[]
inds=[]

for i in range(len(seqs)):
    
    
    if vcf_original_allele_check:
        if type(oris[i]) is not np.string_ or seqs[i][mutpos[i]:(mutpos[i]+len(oris[i]))].upper()!=oris[i]:
            print annos[i], 'Seq ref', seqs[i][mutpos[i]:(mutpos[i]+len(oris[i]))].upper()
            continue
    else:
        if type(oris[i]) is not np.string_:
            continue
    
    #assert seqs[i][mutpos[i]:(mutpos[i]+len(oris[i]))].upper()==oris[i], "Reference allele doesn't match"
    
    inds.append(i)
    seqsmut.append(seqs[i][:mutpos[i]].lower()+ muts[i]+seqs[i][(mutpos[i]+len(oris[i])):].lower())
    seqs[i]=seqs[i][:mutpos[i]].lower()+ oris[i]+seqs[i][(mutpos[i]+len(oris[i])):].lower()

print "Number of input variants:"
print(len(inds))
assert float(len(inds))/len(seqs) > 0.98, 'High ref mismatch: %s out of %s' %(len(inds), len(seqs))

print "Number of valid variants:"
print(len(seqs))
inds = np.asarray(inds)
seqsmut = np.asarray(seqsmut)
seqs = seqs[inds]
oris = oris[inds]
muts = muts[inds]
chrs = chrs[inds]
poss = poss[inds]
annos = annos[inds]
bed_start = bed_start[inds]
if len(names)>0:
    names = names[inds]

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


if writeFasta:
    allrecwt=[]
    allrecmut=[]
    for i in range(len(seqsmut)):
        allrecmut.append(SeqRecord(Seq(seqsmut[i]),id=annos[i], description = ''))
        allrecwt.append(SeqRecord(Seq(seqs[i]),id=annos[i], description = ''))
    SeqIO.write(allrecwt,open(input_fasta+'.ref.fasta','w'),'fasta')
    SeqIO.write(allrecmut,open(input_fasta.replace('wt','mut')+'.ref.fasta','w'),'fasta')

    
if my.f_get_server_name() == 'loire':
    f_stop_script()
    
import h5py
if writeh5:
    writeh5(seqs,input_fasta+'.ref.h5')
    writeh5(seqsmut,dirname(input_fasta)+'/'+basename(input_fasta).replace('wt','mut')+'.ref.h5')

if writebed:
    myfile=open(input_fasta+'.ref.bed','w')
    for i in range(len(seqs)):
        myfile.write(chrs[i]+'\t'+str(poss[i]-mutpos[i])+'\t'+str(poss[i]+mutpos[i]+1)+'\t.\t0\t*\n')
    myfile.close()

if writevcf:
    myfile=open(input_fasta+'.ref.vcf','w')
    if len(names)>0:
        for i in range(len(seqs)):
            myfile.write(chrs[i]+'\t'+str(poss[i])+'\t'+names[i]+'\t'+oris[i]+'\t'+muts[i]+'\t'+bed_start[i]+'\n')
    else:
        for i in range(len(seqs)):
            myfile.write(chrs[i]+'\t'+str(poss[i])+'\t1\t'+oris[i]+'\t'+muts[i]+'\t'+bed_start[i]+'\n')
    myfile.close()
    
        
    
