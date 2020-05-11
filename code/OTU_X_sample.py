# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
import pandas as pd
import argparse

def extract_OTU_uniqs(input_map):
    with open(input_map,"r") as in_map:
        out_dic = {}
        for l in in_map:
            key = l.split()[1]
            value = l.split()[0].split(';')[0]
            out_dic.setdefault(key,[]).append(value)
    return out_dic

def extract_sample_uniqs(input_map):
    with open(input_map,"r") as in_map:
        out_dic = {}
        for l in in_map:
            key = l.split()[0].split('=')[1].split(';')[0]
            value = l.split()[1].split(';')[0]
            out_dic.setdefault(key,[]).append(value)
    return out_dic



if __name__ == "__main__":
    # parser settings
    #
    parser = argparse.ArgumentParser(description = "This tool count uniqs that clustered in same OTUs in each sample")
    # dereped input
    parser.add_argument("-s","--sample_dereped_input", help = "sample vs dereplicated, quality filtered map input file", metavar = "INPUT.txt//tsv")
    parser.add_argument("-o","--output", help = "output filename ", default = "[INFO]sample_X_OTU.csv", metavar = "OUTPUTFILENAME")# output
    # otu after input
    parser.add_argument("-a","--otu_aft_map", help = "uniqs vs. OTUs for filtered reads", metavar = "OTU_aft.txt")#
    args = parser.parse_args()
    sample_dic = extract_sample_uniqs(args.sample_dereped_input)
    otu_dic = extract_OTU_uniqs(args.otu_aft_map)
    sample_X_dic = {}
    ##generate nested dic
    for samp in sample_dic:
        sample_X_dic[samp]= {}
        for OTU in otu_dic:
            sample_X_dic[samp][OTU] =  len(set(sample_dic[samp]).intersection(set(otu_dic[OTU])))
    df  = pd.DataFrame(sample_X_dic).T
    df.fillna(0,inplace = True)
    df.to_csv(args.output,index = True)
    exit()

