# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

import subprocess
import pandas as pd
import argparse
import os
import sys
import warnings

## frequency

def freq_check(in_head):
    """
    reads frequency checker by extracting the infromation from fasta head lines
    Args:
        in_head: input header for extraction
    """
    head = in_head.rstrip()
    head = head.split('=')[1]
    return head

## length
def length_check(in_seq):
    """
    reads length checker
    Args:
        in_seq: input sequence
    """
    length = (len(in_seq))
    return length

## stop codon
def stopcount(seq_record, table, frame = (1,2,3)):
    """
    number of stop codons checker, developed by Dr. Thomas Creedy
    Args: 
        seq_record: input sequence for checking
        table: translation table, which follows the NCBI numbering convention
        frame: ORF 1,2,3
    """
    # Check input types
    run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
    # Run counting
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        counts = [seq_record.seq[(i-1):].translate(table = table).count("*") for i in run_frame]
    # Return string or list depending on length
    if(len(counts) > 1):
        return counts
    else:
        return counts[0]

## main ##

# now open the file and do something :P
def check_property(input,stop_table):
    global uniqs
    global out_df
    with open(input) as infasta : 
        uniqs = [] # to save row names
        freq_list = [] # to save frequency info
        len_list = []# to save length info
        stop_list = [] # to save stop codon
        #
        # LOOP (is the best!)
        print('#Start checking freq and len#')
        # size = len([head for head, seq in SimpleFastaParser(infasta)]) # progress bar
        # step = 0 # progress bar
        # bar = 20 # progress barlength
        for head, seq in SimpleFastaParser(infasta):
            uniqs.append(head.split(';')[0].replace('>',''))## split by ; symbol, keep the first part then remove the > symbol
            freq_list.append(freq_check(head))# add frequency info
            len_list.append(length_check(seq))# add length info
            stop_codon = min(stopcount(SeqRecord(Seq(seq)),stop_table, frame = (1,2,3))) # check stop codon and retain the minimum
            stop_list.append(stop_codon)# append stop codon info
            # step += 1 # progress bar
            # percen = step/size # progress bar
            # heases = '#'* int(percen*bar)# progress bar
            # spaces = '-'* (bar - len(hashes))# progress bar
            # percen = round(percen*100,2)# progress bar
            # sys.stdout.write("\r %d%% |%s| %d/%d lines"%(percen,heases+spaces,step, size))# progress bar
            # sys.stdout.flush()# progress bar
        #
        #
        out_df = pd.DataFrame(data = uniqs, columns =['name'])# reate dataframe for data
        out_df['freq'] = freq_list # add list to dataframe 
        out_df['length'] = len_list
        out_df['stop_codon'] = stop_list
    infasta.close()


# to check chimera
def check_chimera(input):
    global out_df
    print('#Using Vsearch Filtering Chimera#')
    #
    # now let's do the chimera filtering!
    subprocess.Popen("vsearch --uchime3_denovo %s --chimeras %s_chimeras.fa"%(input,filename),shell = True).wait()# chimeras will be exported
    chime_list = ['False']*len(uniqs)# make all the cell be 'false'
    print('#Start checking Chimeras output#')
    with open(filename+'_chimeras.fa') as chimes:
        # size = len([head for head, seq in SimpleFastaParser(chimes)]) # progress bar
        # step = 0 # progress bar
        # bar = 30 # progress bar
        for head ,seq in SimpleFastaParser(chimes):
                index = head.split(';')[0].replace('uniq','')
                chime_list[int(index)-1] = 'True'# check the name id and locate on the list, turn False to True
                # step += 1 # progress bar
                # percen = step/size # progress bar
                # heases = '#'* int(percen*bar)# progress bar
                # spaces = '-'* (bar - len(hashes))# progress bar
                # percen = round(percen*100,2)# progress bar
                # sys.stdout.write("\r %d%% |%s| %d/%dlines"%(percen,heases+spaces,step,size))# progress bar
                # sys.stdout.flush()# progress bar
        out_df['chime'] = chime_list
    chimes.close()


def out(output):
    # output export
    out_df.to_csv(filename+'_info.csv',index= False, header = True)  
    print('done with info table')


if __name__ == "__main__":
    #
    parser = argparse.ArgumentParser(description = "This is a tool that extract frequency, length, stop codon and chimera information from QUALITY FILTERED and DEREPLICATED input fasta file")
    # input
    parser.add_argument("input", help = "input file path", metavar = "INPUT.fasta/INPUT.fa")
    # extract filename
    # translation table numer
    parser.add_argument("-t","--table", help = "the number referring to the translation table to use which follows the NCBI numbering convention (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)", metavar = "TABLE_NUM")
    # output file name with default
    parser.add_argument("-o","--output", help = "output directory (default is current directory)", default = "./", metavar = "OUTPUTFILENAME")# output
    # qulity filter parameter
    parser.add_argument("-q","--quality_maxee", help = "customise --fastq_maxee option in quality filtering with vsearch)", default = "1", metavar = "NUM")# 
    parser.add_argument("-l","--derep_label", help = "customise --relabel option in dereplication with vsearch)", default = "uniq", metavar = "string")
    args = parser.parse_args()
    filename = os.path.splitext(os.path.basename(args.input))[0]
    ## potentially can add switch to each function?
    #check the inputfile and options:
    if os.path.getsize(args.input) == 0:# in input have contents
        sys.exit("Error: input file is empty")
    # maybe add checker to table?
    #
    # quality filter, and dereplication:
    print("###filtering the quality of input file###")
    subprocess.Popen("vsearch --fastx_filter %s --fastq_maxee %s --fastaout %s_qf.fasta"%(args.input,args.quality_maxee,filename),shell=True).wait()
    #
    # dereplication
    print("###Doing dereplication ###")
    subprocess.Popen("vsearch --derep_fulllength %s --output %s_dereped.fasta --sizeout --relabel %s"%(filename+"_qf.fasta",filename,args.derep_label),shell=True).wait()
    #
    # maping derep with raw
    print('##mapping raw with derep##')
    subprocess.Popen("vsearch --search_exact %s -db %s -otutabout %s" %(filename+"_qf.fasta", filename+"_dereped.fasta", args.output+filename+"_origin_map.tsv"), shell = True).wait()
    print("done with info table")
    #
    check_property(filename+"_dereped.fasta",args.table)
    # if (chimera switch)
    check_chimera(filename+"_dereped.fasta")
    out(args.output+"[OUT]"+filename+"_info.csv")
    #
    # OTU delimitation_before:
    print("###Delimitation with data before  fiiltering###")
    subprocess.Popen("vsearch --cluster_size %s --id %s --centroids %s --sizein --relabel %s"%(filename+"_dereped.fasta","0.97",filename+"_OTU_before.fasta","OTU"),shell=True).wait()
    #
    # OTU_before mapping
    subprocess.Popen("vsearch --usearch_global %s -db %s -id %s -otutabout %s"%(filename+"_dereped.fasta",filename+"_OTU_before.fasta","0.97",filename+"_OTU_bef_map.tsv"),shell=True).wait()
    # 
    #OTU_after mapping
    