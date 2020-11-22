#!/usr/bin/python3
import os
import pandas as pd
import hashlib

from glob import glob


def getFastqPrefixes(fastqfiles, delim="_",test_num_delim=None, explicit_delim=None, getmaxlen=False, dohex=False, hexfun=hashlib.md5):
    fastqfiles_base = [os.path.basename(filename) for filename in fastqfiles]
    test_num_delim = fastqfiles[0].count(delim) if test_num_delim is None else test_num_delim


    fastqfiles_series = pd.Series(fastqfiles_base,index=fastqfiles)
    current_delim = 1
    top_counts = 0
    for i in range(1, test_num_delim):
        counts = sum( pd.Series([ delim.join(f.split(delim)[:i]) for f in fastqfiles_base]).value_counts() == 2)
        counts = counts - sum(pd.Series([ delim.join(f.split(delim)[:i]) for f in fastqfiles_base]).value_counts() > 2)
        if getmaxlen:
            if counts >= top_counts:
                top_counts = counts
                current_delim = i
        else:
            if counts > top_counts:
                top_counts = counts
                current_delim = i



    fastq_prefixmap = pd.DataFrame(dict(fastq=fastqfiles, prefix = [delim.join(fastqbase.split(delim)[:current_delim]) for fastqbase in fastqfiles_base])).set_index("prefix")
    
    fastq_prefixmap["counts"] = pd.Series(fastq_prefixmap.index).value_counts()
    fastq_prefixmap["sample"] = fastq_prefixmap.index

    if dohex:
        fastq_prefixmap["hash"] = fastq_prefixmap.fastq.apply(lambda x:hexfun(open(x,"rb").read(1000000)).hexdigest())

    return fastq_prefixmap





def makeinputtable(folder):
	
    from glob import glob
    
    #Only accept gz files for now
    files = glob("%s/*.gz"%folder)

    #sample_table = getFastqPrefixes(glob("/home/werner/Sources/Kamini/fastqs/*")).sort_values("fastq")
    sample_table = getFastqPrefixes(files).sort_values("fastq")
    return sample_table
