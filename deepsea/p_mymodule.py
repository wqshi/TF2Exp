#!
import re
import subprocess
import os
import sys
import socket
import pandas as pd
import logging


server_name = socket.gethostname()
if (server_name == "loire"):
    sys.path.append("/homed/home/shi/packages/TFFM")
elif "bugaboo" in server_name:
    sys.path.append("/home/wenqiang/packages/TFFM")
else:
    sys.path.append("/home/shi/packages/TFFM")

def grep_list(pattern, word_list, ignore_case = True ):
    if (ignore_case == True):
        expr = re.compile(pattern,re.I)
    else:
        expr = re.compile(pattern)
    return [elem for elem in word_list if expr.match(elem)]


def grep_file(pat,file):
    """ finds patt in file - patt is a compiled regex
        for negative match, we can use (?!pattern)
        returns all lines that match patt (case insensitive) """

    patt=re.compile(pat,re.I)
    matchlines = []
    filetxt = open(file)
    lines = filetxt.readlines()
    for line in lines:
        #print line
        match = patt.search(line)
        if match:
           matchlines.append(line.replace('\n',''))
    if matchlines:
        return matchlines
    else:
        return None


def f_server_name_test():
    #
    server_name = subprocess.check_output(["uname","-a"])
    return server_name

def f_server_name():
    #This won't give you loire string on loire. See f_get_server_name
    server_name = subprocess.check_output(["uname","-a"])
    return server_name


def f_shell_pipe(cmd):
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    output = ps.communicate()[0]
    return output

def f_shell_pipe_bash(cmd):
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT, executable='/bin/bash')
    output = ps.communicate()[0]
    return output


def f_shell_fun_pipe(shell_fun):
    cmd = "source ~/library2/s_vcf_function2.sh;source ~/library2/s_function.sh;" + shell_fun
    ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    output = ps.communicate()[0]
    return output

def f_call_shell_fun(shell_fun):
    cmd = "source ~/library2/s_vcf_function2.sh;source ~/library2/s_function.sh;" + shell_fun
    print subprocess.check_output(cmd,shell=True)


def f_shell_cmd(cmd,quiet=False):
    
    if quiet == False:
        print cmd
        print subprocess.check_output(cmd,shell=True)
    else:
        return subprocess.check_output(cmd,shell=True)
    
    
def f_download(urlbase, filename):
    if not os.path.isfile(filename):
        url = urlbase + filename
        try:
            urllib.urlretrieve(url, filename=filename)
        except IOError, e:
            raise IOError('No Internet connection')
        # Check if downloaded file is an HTML file, which
        # is what github.com returns if the URL is not existing
        f = open(filename, 'r')
        if 'DOCTYPE html' in f.readline():
            raise IOError('URL %s does not exist' % url)


def f_import_file(full_path_to_module):
    try:
        import os
        module_dir, module_file = os.path.split(full_path_to_module)
        module_name, module_ext = os.path.splitext(module_file)
        save_cwd = os.getcwd()
        os.chdir(module_dir)
        module_obj = __import__(module_name)
        module_obj.__file__ = full_path_to_module
        globals()[module_name] = module_obj
        os.chdir(save_cwd)
    except:
        raise ImportError


def f_grep_wget_from_given_file(index_file, pattern,  data_url, output_dir, prefix, test_flag=False, quiet=False, debug = False):
    if debug == True:
        import ipdb; ipdb.set_trace()
        
    import urllib
    matched_lines=grep_file(pattern, index_file)
    if matched_lines == None:
        if quiet == False:
            print "-----------------------Warning--------------------------\nNo matching for the pattern %s in %s\n"%(pattern, index_file)
        return "failed"
    
    file_names=[re.split('\t',line)[0] for line in matched_lines]
    #print file_names

    for file_name in file_names:

        file_suffix=re.match(r"[a-zA-Z0-9_]*\.(.*)",file_name).group(1)
        #print file_suffix
    
        match_object=re.match(r".*(Rep[1-9]).*",file_name,flags=re.IGNORECASE)
        if match_object:
            output_name=  prefix + "-"+match_object.group(1) + "." + file_suffix
        else:
            output_name = prefix + "." + file_suffix

        output_file = output_dir + "/" + output_name

        if test_flag == False:
            urllib.urlretrieve(url=data_url + "/" + file_name, filename= output_file )
            
        print "Downlaod " + file_name + " " + data_url
        
        match_object=re.match(r".*\Peak.gz$",file_name)
        if match_object:
            if test_flag == False:
                f_unzip_targz(output_file)
            print "Unzip " + output_name
        
    return "success"

        
def f_create_file_prefix(cell, tf, lab=None,rep=None):
    if lab == None:
        file_prefix = cell + "-" + tf
    else:
        file_prefix = lab + "-" + cell + "-" + tf

    if rep != None:
        file_prefix = file_prefix + "-" + rep
    
    return file_prefix

def f_create_file_name(data_dir,cell, tf, suffix ,lab=None,rep=None):
    return "%s/%s.%s"%(data_dir,f_create_file_prefix(cell, tf, lab,rep),suffix)



def f_unzip_targz(gz_file):

    cmd=""
    if re.match(".*.tar.gz",gz_file):
        cmd="tar -xvzf " + gz_file
    elif re.match(".*\.gz",gz_file):
        cmd="gzip -d -f " + gz_file
    else:
        print  "Unknown file" + gz_file
    print cmd
    f_shell_pipe(cmd)



def f_cat_files(data_dir, file_list, output_file, test_flag = False):
    old_dir=os.getcwd()
    
    os.chdir(data_dir)
    input_files=" ".join(file_list)
    if test_flag == True:
        cmd = "cat %s | head -n 10000 > %s" % (input_files, output_file)
    else:
        cmd = "cat %s > %s" % (input_files, output_file)
    f_shell_cmd(cmd)
    os.chdir(old_dir)


def f_ensure_make_dir(dnase_dir):
        if not os.path.exists(dnase_dir):
                os.makedirs(dnase_dir)
        else:
            print 'DIR: %s exists' % dnase_dir

def f_send_list_para(input_list):
    if input_list == []:
        return 'None'
    else:
        return ":".join(input_list)

def f_recieve_list_para(input_string):
    if input_string == 'None':
        return []
    else:
        return input_string.split(":")

def f_copy_and_remove_dir(input_dir, output_dir):
    cmd = "cp %s/* %s/; rm -r %s" % (input_dir, output_dir, input_dir)
    f_shell_cmd(cmd) 

def f_copy_to_dir(input_dir, pattern ,output_dir, option=""):
    try:
        if input_dir != output_dir:
            cmd = "cp %s %s/%s %s/" % (option,input_dir, pattern ,output_dir)
            f_shell_cmd(cmd)
        else:
            print "f_copy_to_dir: input and output dir are same"
    except:
        logging.error('Copy cmd failed: %s ' % cmd)

def f_copy_to_bugaboo_dir(input_dir, pattern ,output_dir, option=""):

    if "^/homed/home/shi" in output_dir:
        bgb_output_dir = output_dir.replace("/homed/home/shi", "/home/wenqiang")
    elif "^/home/shi" in output_dir:
        bgb_output_dir = output_dir.replace("/home/shi", "/home/wenqiang")
    else:
        bgb_output_dir = output_dir
        
    cmd = "scp %s %s/%s wenqiang@bugaboo.westgrid.ca:/%s/" % (option,input_dir, pattern ,output_dir)
    f_shell_cmd(cmd)


def f_get_server_name():
    import socket
    return socket.gethostname()

def f_remove_dir(data_dir):
    if os.path.exists(data_dir):
        cmd = "rm -r %s" % data_dir
        f_shell_cmd(cmd)


def f_remove_empty_dir(data_dir):
    if not os.listdir(data_dir):
        os.rmdir(data_dir)
    else:
        print "Not empty dir"
    
def f_get_prefix(input_string, path=True):
    if path == True:
        match_object=re.match("([a-zA-X\-_0-9/]*)\..*", input_string)
    else:
        match_object=re.match("([a-zA-X\-_0-9/]*)\..*", os.path.basename(input_string))
    prefix=match_object.group(1)

    return prefix

def f_get_dir_name_from_file(file_name):
    return os.path.dirname(file_name) + "/"
    

def f_scp_to_loire(from_dir, file_pattern, target_dir):
    cmd="scp %s/%s shi@loire.cmmt.ubc.ca:%s"%(from_dir, file_pattern, target_dir)
    print cmd
    f_shell_cmd(cmd)


#Test cases
#input_dir="/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/histone"
#output_dir='/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/'
#f_copy_and_remove_dir(input_dir, output_dir)    



def f_vcf_to_bed(vcf_file, output_file):
    cmd = """sed '/^chr[ ,\t]/d' %s | awk -F"\\t" -v OFS="\\t" '{print $1,$2-1,$2,$1"-"$2}' > %s""" % (vcf_file, output_file)
    print cmd
    f_shell_pipe(cmd)


def f_vcf_to_allele(vcf_file, output_file):
    cmd = """echo "chr\tstart\tend\tref\talt" > %s; sed '/^chr[ ,\t]/d' %s | awk -F"\\t" -v OFS="\\t" '{print $1,$2,$2+1,$3,$4}' >> %s""" % (output_file,vcf_file, output_file)
    print cmd
    f_shell_pipe(cmd)



def f_read_tffm(tffm_dir, tf_name):
    import tffm_module
    from constants import TFFM_KIND
    
    file_list=os.listdir(tffm_dir)

    if tf_name.lower() != "ctcf" and tf_name.lower() != "ebf1":
        print "TFFM not found for %s"%tf_name
        print "[Warning: change tf to CTCF motif]"
        tf_name = "ebf1"
    
    pattern="%s.*detail.*xml"%(tf_name)

    matched_file=grep_list(pattern, file_list)
    if len(matched_file) == 1:
        tffm_file=tffm_dir + grep_list(pattern, file_list)[0]

        #logging.debug(tffm_file)

        #tffm_first_order = tffm_module.tffm_from_xml(tffm_file, TFFM_KIND.FIRST_ORDER)
        tffm_detailed = tffm_module.tffm_from_xml(tffm_file, TFFM_KIND.DETAILED)
        return tffm_detailed
    else:
        return None


def f_remove_files(input_dir, pattern ):
    cmd = "rm %s/%s" % (input_dir, pattern)
    f_shell_cmd(cmd)


def f_list_server_dir(target_server,data_dir, file_suffix):
    if "clustdell" in target_server :
        cmd="ssh shi@clustdell.cmmt.ubc.ca /home/shi/s_list_file.sh %s \"\*.%s\""%(data_dir, file_suffix)
    elif "loire" in target_server:
        cmd="cd %s; ls *.%s"%(data_dir, file_suffix)
    else:
        cmd="ssh wenqiang@bugaboo.westgrid.ca /home/wenqiang/s_list_file.sh %s \"\*.%s\""%(data_dir, file_suffix)
    file_list_raw=f_shell_pipe(cmd)
    file_list=re.split("\n",file_list_raw)
    return file_list

def f_list_dir(data_dir, pattern="*",path=True, quiet = False):
    if path == False:
        cmd="cd %s; ls %s"%(data_dir, pattern)
    else:
        cmd="ls %s/%s"%(data_dir, pattern)
    file_list_raw=f_shell_pipe(cmd)

    if "No such file" in file_list_raw:
        return None
    
    file_list=re.split("\n",file_list_raw)
    del file_list[-1]
    return file_list
    
    
def f_copy_from_bugagoo_to_clustdell(input_dir, pattern ,output_dir, option=""):

    if "/home/wenqiang" in output_dir:
        clustdell_output_dir = output_dir.replace("/home/wenqiang", "/home/shi")
    else:
        clustdell_output_dir = output_dir

    import time
    cur_time=time.strftime("%Y%m%d_%H%M%S")
    tmp_dir="/homed/home/shi/data/clustdell/%s" % cur_time
    cmd="ssh shi@loire.cmmt.ubc.ca 'mkdir %s'"% tmp_dir
    f_shell_cmd(cmd)

    import time
    time.sleep(5)
    file_list=f_grep_files_from_dir(input_dir, pattern)
    print file_list
    cmd = "scp %s %s shi@loire.cmmt.ubc.ca:%s/" % (option, " ".join(file_list) ,tmp_dir)
    f_shell_cmd(cmd)
    time.sleep(5)
    cmd="ssh shi@loire.cmmt.ubc.ca /homed/home/shi/s_scp_to_clustdell.sh %s %s"%(clustdell_output_dir, tmp_dir)
    f_shell_cmd(cmd)



def f_grep_files_from_dir(data_dir, pattern, path=True, debug=False):
    file_list=os.listdir(data_dir)
    match_files=grep_list(pattern, file_list)

    if debug == True:
        print "\n ====debug in  f_grep_files_from_dir()===="
        print data_dir, pattern, file_list
    
    if path == True:
        return [ data_dir + "/" + f for f in match_files]
    else:
        return match_files


def f_grep_cell_tf_file_from_dir(base_dir, cell_name, feature, suffix, labs=None):
    #import ipdb; ipdb.set_trace()
    cell_dir = os.path.join(base_dir, cell_name)
    pattern=f_create_pattern([cell_name],[feature],suffix, labs, strict = True)
    matched_files=f_grep_files_from_dir(cell_dir, pattern)
    matched_file=f_unique_element_in_list(matched_files)
    return matched_file




def f_parse_file_name(file_name):
    
    style_num=file_name.count("-")

    if style_num == 2:
        matched_obj=re.match("([a-zA-Z0-9]+)-([a-zA-Z0-9]+)-([a-zA-Z0-9]*)[.]", file_name)
        lab=matched_obj.group(1)
        cell=matched_obj.group(2)
        tf   =matched_obj.group(3)
    elif style_num == 1:
        matched_obj=re.match("([a-zA-Z0-9]+)-([a-zA-Z0-9]+)[.]", file_name)
        lab=None
        cell=matched_obj.group(1)
        tf   =matched_obj.group(2)
    else:
        print "Invalid Input file name"
        lab=None
        cell=None
        tf=None

    return [lab, cell, tf]


def f_parse_rep_file_name(file_name):
    
    style_num=file_name.count("-")

    if style_num == 3:
        matched_obj=re.match("([a-zA-Z0-9]+)-([a-zA-Z0-9]+)-([a-zA-Z0-9]*)-(Rep[0-9])[.]", file_name)
        lab=matched_obj.group(1)
        cell=matched_obj.group(2)
        tf   =matched_obj.group(3)
        rep = matched_obj.group(4)
    elif style_num == 2:
        matched_obj=re.match("([a-zA-Z0-9]+)-([a-zA-Z0-9]*)-(Rep[0-9])[.]", file_name)
        lab=None
        cell=matched_obj.group(1)
        tf   =matched_obj.group(2)
        rep = matched_obj.group(3)
    else:
        print "Invalid Input file name"

    return [lab, cell, tf, rep]

def f_gzip_file(file_name):
    cmd="gzip -f %s" % file_name
    print f_shell_cmd(cmd)



def f_grep_and_copy(input_dir, pattern, output_dir, option=""):
    file_list=f_grep_files_from_dir(input_dir, pattern)
    print file_list
    if file_list != []:
        cp_cmd= "cp %s %s %s" % (option, " ".join(file_list) ,output_dir)
        #print cp_cmd
        f_shell_cmd(cp_cmd)
    else:
        print "File not find: %s %s"%(input_dir, pattern)
        
def f_grep_and_rm(input_dir, pattern, option="", quiet = False):
    file_list=f_grep_files_from_dir(input_dir, pattern)
    if quiet == False:
        print file_list
    if file_list != []:
        rm_cmd = "rm %s %s " % (option, " ".join(file_list) )
        #print rm_cmd
        f_shell_cmd(rm_cmd, quiet)
    else:
        print "File not find: %s %s"%(input_dir, pattern)


def f_write_locer_file(cur_time, server_name):
    import time
    if server_name != "loire":
        cmd="ssh shi@loire.cmmt.ubc.ca 'echo '=' > $HOME/locker_dir/%s'"%(cur_time)
    else:
        cmd="echo '=' > $HOME/locker_dir/%s"%(cur_time)
    f_shell_cmd(cmd)


def f_check_locer(cur_time):
    
    import os.path
    print "Check the locker:" + cur_time
    while True:
        if os.path.isfile("/homed/home/shi/locker_dir/%s"%(cur_time)):
            print "Get the locker"
            break
            

        import time
        time.sleep(10)

def f_grep_and_scp_to_loire(input_dir, pattern, output_dir):
    file_list=f_grep_files_from_dir(input_dir, pattern)
    print file_list
    
    f_scp_to_loire("", " ".join(file_list), output_dir)

def f_create_pattern(cell_list, tf_list, suffix, labs='', strict = False):
    if strict == True:
        pattern="(%s).*(%s).*(%s)[.]?%s"%('|'.join(labs), "|".join(cell_list), "|".join(tf_list), suffix)
    else:
        pattern="(%s).*(%s).*(%s).*%s"%('|'.join(labs), "|".join(cell_list), "|".join(tf_list), suffix)
    return pattern

def f_unique_element_in_list(input_list, notice_message=""):
    if len(input_list) == 1:
        return input_list[0]
    else:
        print "==Unique element faild=="
        print notice_message;
        print "None unique list: ", input_list
        return None
        #sys.exit()


def f_scp_python_script_to_clustdell(script_name, script_type = 'python'):
    from os.path import expanduser
    home = expanduser("~")
    cmd="scp shi@loire.cmmt.ubc.ca:/homed/home/shi/%s/%s %s/%s/" %(script_type, script_name, home, script_type)
    f_shell_cmd(cmd)

def f_get_cell_dir_name(cell_name):
        
    if "gm1" in cell_name and "gm12878" not in cell_name:
        #dir_name="gm12xxx"
        dir_name = cell_name
    else:
        dir_name=cell_name

    return dir_name

def f_parse_tf_list_file(tf_list_file):
    tf_list=list(pd.io.parsers.read_csv(tf_list_file).ix[:,0])
    tf_list_unique=list(set([tf.split("_")[0].lower().replace("egfp-","e").replace("(","").replace(")","").replace("-","") for tf in tf_list]))
    tf_list_unique.sort()
    return tf_list_unique


def f_id_generator(size=15):
    import string
    import random
    chars=string.ascii_uppercase + string.digits
    return ''.join(random.choice(chars) for _ in range(size))



def f_generate_tmp_file_name(discription):
    return "tmp.%s.%s"%(discription, f_id_generator(15))



def f_duplicated_index(data, keys):
    return data.duplicated(keys) | data.duplicated(keys, take_last = True)

def f_bed_to_pd(bed_data):
    #import ipdb; ipdb.set_trace()
    pd_data=[]
    import pandas as pd
    for record in bed_data:
        #print list(record)
        pd_data.append(list(record))

    pd_frame=pd.DataFrame(pd_data)
    return pd_frame

def f_bed_to_pd2(bed_data):
    #import ipdb; ipdb.set_trace()
    import pandas as pd
    #bed_data.to_csv(bed_data.fn, header=None, index = None, sep='\t')
    pd_frame= pd.read_table(bed_data.fn, header = None, sep='\t')
    return pd_frame

def f_pd_to_bed(pd_data):
    import pybedtools
    if pd_data.shape[1] >10:
        logging.error("Input data exceed 10 cols for bed")
    bed_str=pd_data.to_string(header=False, index=False)
    bed_file=pybedtools.BedTool(bed_str, from_string=True)
    return bed_file

def f_pd_to_bed_based_on_file(pd_data):
    import pybedtools
    if pd_data.shape[1] >10:
        logging.error("Input data exceed 10 cols for bed")

    bed_data=pybedtools.BedTool('chrX 1 100', from_string=True)
    pd_data.to_csv(bed_data.fn, header=None, index = None, sep='\t')
    return bed_data



def f_get_reference_genome():
    server_name = f_get_server_name()
    import os
    if server_name == "loire":
        hg_file="/shared2/RAW/shi/data/wgs/hg19.fa"
    elif server_name == 'wqshi':
        
        home_dir = os.path.expanduser('~')
        hg_file="%s/projects/wgs/hg19.fa" % home_dir
    else:
        #hg_file="/raid2/local/exome/data/hg19/hg19-orderMito.fasta"
        hg_file="./data/raw_data/wgs/hg19.fa"

    if os.path.exists(hg_file) == False:
        logging.error('Reference genome file missing, %s' % hg_file)
        
    return hg_file



def f_pd_replace_missing_value(data_frame, missing_value = ".", replace_value = "0"):
    #import ipdb; ipdb.set_trace()
    for col_name in data_frame.columns:
        print col_name
        data_frame.ix[data_frame[col_name] == missing_value,col_name] = replace_value

    return data_frame



def f_contain_and_replace(input_vector, pattern_list, replace_list):
    #import ipdb; ipdb.set_trace()
    output_vector  = input_vector
    for i in range(0, len(pattern_list)):
        output_vector[input_vector.str.contains(pattern_list[i])] = replace_list[i]
    return output_vector


def f_print_list(input_list):
    print '\n'.join(input_list)



def f_cat_multiple_files(tf_peak_list, output_basename):
    output_file = os.path.dirname(tf_peak_list[0]) + '/' + output_basename
    cat_cmd = 'cat %s | cut -f1-9 > %s' %(' '.join(tf_peak_list), output_file)
    f_shell_cmd(cat_cmd)
    return output_file



def f_create_cell_pair_database_name(target_dir, cell, guest_cells, loc_tf):
    return f_create_file_name(target_dir, '%s'%('-'.join( [cell] + guest_cells)) , loc_tf,"database")



def f_get_tf_known_cofactors(input_tf):

    targeted_tf = input_tf.upper()
    cofactor_file = os.path.expanduser('~/projects/wgs/cofactor_list.txt')
    cofactors_pd=pd.read_csv(cofactor_file, sep = '\t')

    cofactors_pd.index = cofactors_pd.tf_name

    if targeted_tf in cofactors_pd.tf_name:
        cofactors = cofactors_pd.ix[targeted_tf, 'cofactors'].split(',')
    else:
        cofactors = []

    return cofactors


def f_sync_scripts_to_run_server(target_server = 'clust'):
    import subprocess
    print ''
    print 'Sync scripts ...........'
    print ''

    if target_server == 'clust':
        target_str = 'shi@clustdell.cmmt.ubc.ca:/home/shi/'
    else:
        target_str = 'wenqiang@orcinus.westgrid.ca:/home/wenqiang/'
    
    
    rsync_cmd1  = "rsync -rav --include '*.py' --exclude '*' /homed/home/shi/python/ %s/python/" % target_str
    rsync_cmd2 = "rsync -rav --include '*.sh' --exclude '*' /homed/home/shi/library2/ %s/library2/" % target_str
    rsync_cmd3  = "rsync -rav --include '*.py' --exclude '*' /homed/home/shi/python/deepASB/ %s/projects/chipseq_snp/python/deepASB" % target_str
    rsync_cmd4  = "rsync -rav --include '*.*' --exclude '*' /homed/home/shi/python/deepASB/cnn/ %s/python/deepASB/cnn/" % target_str
    subprocess.call(rsync_cmd1,shell=True)
    subprocess.call(rsync_cmd2,shell=True)
    subprocess.call(rsync_cmd3,shell=True)
    subprocess.call(rsync_cmd4,shell=True)


def f_get_file_size(input_file):
    import os
    if os.path.isfile(input_file):
        file_size=os.path.getsize(input_file)/float(1000)
    else:
        file_size = 0

    return file_size

def f_debug(debug_flag = True):
    if f_get_server_name() == 'loire' and debug_flag == True:
        return True
    else:
        return False


