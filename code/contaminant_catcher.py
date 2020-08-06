# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
import argparse
import os
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
# import subprocess

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "test code")
    # parser.add_argument("-db" ,"--blast_datdbase",help= "blast searching database directory",,metavar = '<directory>/nt' default= '~/db/blast-db-2019-10-21/nt')
    parser.add_argument("input" , help = "input barcodes files dircetory")
    parser.add_argument("-bs" ,"--bc_selected", help= "input NAPselected barcode files")
    parser.add_argument("-fa" ,"--filtered_ASVs", help= "input filtered ASVs fasta file")
    parser.add_argument("-mb" ,"--mb_file", help= "input metabarcoding concatenated files")
    parser.add_argument("-mn" ,"--mb_file_names", help= "input metabarcoding concatenated files names")
    parser.add_argument("-cs","--conta_score_out", help = "output filename for contaminant score in mb (with default) ", default = "[INFO]contaminant_score.csv", metavar = "SCOREFILENAME")# score output
    parser.add_argument("-cr","--conta_reads_out", help = "output filename for bc contaminat reads that occur in mb (with default)", default = "[INFO]contaminant_reads.fasta", metavar="READFILENAME")    
    args = parser.parse_args()

    ## extract ref dict
    ref_dic = {}
    with open (args.bc_selected) as ref:
        for head , seq in SimpleFastaParser(ref):
            ref_dic[head.strip()] = seq.strip()
    
    ##extract contamenant dict
    conta_dic = {}# store contaminat sequences
    conta_count_dic = {}# store scores
    for path, dirs ,f in os.walk(args.input):
        for ref_head in ref_dic.keys():
            if ref_head.strip()+".fasta" in f:#if the file be found in input folder 
                with open(path+ref_head.strip()+".fasta") as fq:# open it, then record from the start 
                    
                    tmp_seq_list = []#empty list to store conta seq
                    for head , seq in SimpleFastaParser(fq):
                        
                        if seq == ref_dic[ref_head]:
                            # bar = ref_head.strip() +":"+ head
                            if tmp_seq_list == []:
                                conta_count_dic[ref_head] = 0
                                pass
                            else:
                                conta_count_dic[ref_head]=len(tmp_seq_list)
                                conta_dic[ref_head] = tmp_seq_list
                                break
                        else:
                            tmp_seq_list.append(seq)
            else:
                pass
    
    ## filter with cleaned ASVs:
    # first extract dict
    ASV_list = []
    with open (args.filtered_ASVs) as ASV_file:
        for head , seq in SimpleFastaParser(ASV_file):
            ASV_list.append(seq.rstrip())
    
    # then filter concat dict:
    for key, values in conta_dic.items():
        # rm = 0 # record how many been removed
        for seq in values:
            if seq not in ASV_list:
                conta_dic[key].remove(seq)
        conta_count_dic[key] = len(values)
    
    # extract mb file names
    with open(args.mb_file_names,"r") as mb_names:
        mb_names_list = mb_names.readlines()
        mb_names_list = [s.rstrip() for s in mb_names_list]
    
    ## loop through all mb files
    score = {}
    with open(args.mb_file,"r") as input:
        for mb_f in mb_names_list:
            score[mb_f] = {}
            mb_tmp = []
            for head, seq in SimpleFastaParser(input) :
                if re.search(mb_f,head):
                    mb_tmp.append(seq)
            for ref_key , ref_seq in ref_dic.items():
                score[mb_f][ref_key] = None
                if ref_seq in mb_tmp:
                    tmp_score = 0
                    if ref_seq in conta_dic.keys() and len(conta_dic[ref_key]) != 0 :#check if barcoding exist in contaminant dict and contaminat is not empty
                        for con_seq in conta_dic[ref_key]:
                            if con_seq in mb_tmp:
                                tmp_score += 1 
                        tmp_score = tmp_score/ len(conta_dic[ref_key])
                    score[mb_f][ref_key] = tmp_score
    # score = {}
    # for  path,dirs , files in os.walk(args.mb_dir):
    #     for f in files:
    #         score[f] = {}
            
    #         tmp = {}
    #         with open (path+f,'r') as mb:#open each mb file
    #             for head, seq in SimpleFastaParser(mb):#extract a tmp dict
    #                 tmp[head] = seq
                    
    #         for ref_key , ref_seq in ref_dic.items(): #loop through ref_dict
    #             score[f][ref_key] = None #inf retaining score
    #             if ref_seq in tmp.values():
    #                 tmp_score = 0
    #                 if ref_key in conta_dic.keys() and len(conta_dic[ref_key]) != 0 :#check if barcoding exist in current file
    #                     for con_seq in conta_dic[ref_key]:
    #                         if con_seq in tmp.values():
    #                             tmp_score += 1 
    #                     tmp_score = tmp_score/ len(conta_dic[ref_key])
    #                 score[f][ref_key] = tmp_score
    
    out_df = pd.DataFrame(score).T
    out_df.to_csv(args.conta_score_out,index = True)
    with open(args.conta_reads_out, "w") as conta_reads:
        for key, value in conta_dic.items():
            for reads in value:
                conta_reads.write(key+"_"+str(value.index(reads)+1)+"\n")
                conta_reads.write(reads+"\n")
