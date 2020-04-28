# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
import pandas as pd
import argparse

def extract_multiple_dictionary(input_map):
    """create one key multiple values dictionary"""
    with open(input_map,"r") as input_blast_map:
        uniq_dic = {}
        for l in input_blast_map: ## from each lib
            key = l.split()[0].split(';')[1].split('=')[1]## the filename
            value = l.split()[1].split(';')[0]##the uniq
            uniq_dic.setdefault(key,[]).append(value)
    return uniq_dic

def extract_dictionary(input_map):
    """create one key one value dictionary"""
    with open(input_map,"r") as input_blast_map:
        dic = {}
        for l in input_blast_map:
            key = l.split()[0].split(';')[0]
            value = l.split()[1]
            dic[key]=value
    return dic


def get_len(uniq_dic):
    """convert dictionary to length of it"""
    for key in uniq_dic.keys():
        uniq_dic[key] = len(set(uniq_dic[key]))
    return uniq_dic

def get_unique(uniq_dic):
    """convert dictionary to unique contents"""
    for key in uniq_dic.keys():
        uniq_dic[key] = list(set(uniq_dic[key]))
    return uniq_dic

def replace(in_list,dictionary):
    """replace values in the list with given dictionary"""
    rep_list = [dictionary[x] if x in dictionary else x for x in in_list]
    return rep_list

if __name__ == "__main__":
    # parser settings
    #
    parser = argparse.ArgumentParser(description = "This tool match the OTUs for filtered and unfiltered unique reads")
    # dereped input
    parser.add_argument("-s","--sample_dereped_input", help = "sample vs dereplicated, quality filtered map input file", metavar = "INPUT.txt//tsv")
    # otu before input
    parser.add_argument("-b","--otu_bef_map", help = "uniqs vs OTUs for unfiltered(but quality filtered and dereplcated) reads; geerated with vsearch in blast6out method", metavar = "OTU_bef.txt")
    # output file path with default
    parser.add_argument("-o","--output", help = "output directory (default is current directory)", default = "./", metavar = "OUTPUTFILENAME")# output
    # otu after input
    parser.add_argument("-a","--otu_aft_map", help = "uniqs vs. OTUs for filtered reads", metavar = "OTU_aft.txt")# 
    # parser.add_argument("-id","--id_otu", help = "identity for OTU clustering, e.g. for 3%% clustering should input 0.97",metavar = "float", default = "0.97" )
    args = parser.parse_args()
    #
    # FOR EACH SAMPLE:
    # 0. HOW MANY READS match to a uniq? i.e. len of value of sample_uniqs_dic
    # 1. How many total uniqs?
    # 2. How many uniqs match to bef_otus?
    # 3. How many bef_otus?
    # 4. How many uniqs match to aft_otus?
    # 5. How many aft_otus?
    # 6. How many uniqs match to bef_otus and/or aft_otus?
    # 7. How many uniqs match to bef_otus AND aft_otus?
    # 8. How many uniqs do not match to either?
    # 9. How many unique OTUs are there in total (i.e. need to use a mapping of aft_otus to bef_otus)   
    
    
    #extract dic from origin map
    # oringin_map = open(args.sample_dereped_input,"r")
    sample_uniqs_dic = extract_multiple_dictionary(args.sample_dereped_input)
    #
    #extract dic from bef_OTU map
    # bef_map = open(args.otu_bef_map,"r")
    bef_dic = extract_dictionary(args.otu_bef_map)
    #
    #extract dic from aft_OTU map
    # aft_map = open(args.otu_aft_map,"r")
    aft_dic = extract_dictionary(args.otu_aft_map)
    #
    #replace uniqs in sample_uniq_dic with OTU_bef
    sample_befOTU_dic = {}
    for k in sample_uniqs_dic:
        sample_befOTU_dic[k] = replace(sample_uniqs_dic[k],bef_dic) # What happens if a uniq in the list isn't a key in the bef_dic
    #
    #replace uniqs in sample_uniq_dic with OTU_aft
    sample_aftOTU_dic = {}
    for k in sample_uniqs_dic:
        sample_aftOTU_dic[k] = replace(sample_uniqs_dic[k],aft_dic)
    #
    #convert info dictionary to len
    sample_uniqs_dic = get_len(sample_uniqs_dic)
    sample_befOTU_dic = get_len(sample_befOTU_dic)
    sample_aftOTU_dic = get_len(sample_aftOTU_dic)
    #
    out_df = pd.DataFrame([sample_uniqs_dic,sample_befOTU_dic,sample_aftOTU_dic],index = ["uniqs","bef_OTU","aft_OTU"])
    #
    out_df.T.to_csv(args.output+'[INFO]'+'uniq_and_OTU_info.csv',index=True)
