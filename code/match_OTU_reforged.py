# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

import subprocess
import pandas as pd
import argparse
# import sys
# import warnings

def map_to_dic(input_map,transpose = False):
    """convert frequency map to dictionary without preserving frequency info"""
    print("READing the map...")
    the_map =  open(input_map,"r")
    dic = {}
    for l in the_map:
        key = l.split()[0].split(';')[0]
        value = l.split()[1]
        dic[key] = value
    return dic

def match_otu(bef_dic,aft_dic,otus_dic):
    """match OTU before and after in a dictionary and identify if they are same OTU with otus_dic"""
    id_otu = {}
    for key in bef_dic.keys():
        before = bef_dic[key]
        if key in aft_dic.keys():
            after = aft_dic[key]
            if before in otus_dic.keys() and after == otus_dic[before]:
                id_otu[key] = "True"
            else:
                id_otu[key] = "False"
        else:
            id_otu[key] = "False"
    return id_otu


if __name__ == "__main__":
    # parser settings
    #
    parser = argparse.ArgumentParser(description = "This tool match the OTUs for filtered and unfiltered unique reads")
    # dereped input
    parser.add_argument("-d","--dereped_input", help = "dereplicated, quality filtered input file", metavar = "INPUT.fasta/INPUT.fa")
    # extract filename
    # otu before input
    parser.add_argument("-b","--otu_bef", help = "OTUs for unfiltered(but quality filtered and dereplcated) reads", metavar = "OTU_bef.fasta")
    # output file name with default
    parser.add_argument("-o","--output", help = "output directory (default is current directory)", default = "./", metavar = "OUTPUTFILENAME")# output
    # otu after input
    parser.add_argument("-a","--otu_aft", help = "OTUs for filtered reads", metavar = "OTU_aft.fasta")# 
    parser.add_argument("-id","--id_otu", help = "identity for OTU clustering, e.g. for 3%% clustering should input 0.97",metavar = "float", default = "0.97" )
    args = parser.parse_args()

    # filename = os.path.splitext(os.path.basename(args.input))[0]
    ## potentially can add switch to each function?
    #check the inputfile and options:
    # if os.path.getsize(args.dereped_input) == 0:# in input have contents
    #     sys.exit("Error: dereped_input file is empty")
    #
    print("START generating maps for analysis")
    # OTU_before mapping
    subprocess.Popen("vsearch --usearch_global %s -db %s -id %s -blast6out %s"%(args.dereped_input,args.otu_bef,args.id_otu,"[tmp]OTU_bef_map.txt"),shell=True).wait()
    #
    bef_map = "[tmp]OTU_bef_map.txt"
    #
    # 
    # OTU_after mapping
    subprocess.Popen("vsearch --usearch_global %s -db %s -id %s -blast6out %s"%(args.dereped_input,args.otu_aft,args.id_otu,"[tmp]OTU_aft_map.txt"),shell=True).wait()
    #
    aft_map = "[tmp]OTU_aft_map.txt"
    #
    # 
    # OTUs_match_map
    subprocess.Popen("vsearch --search_exact %s -db %s -blast6out %s"%(args.otu_bef,args.otu_aft,"[tmp]OTUs_map.txt"),shell=True).wait()
    #
    otu_map = "[tmp]OTUs_map.txt"
    # 
    # parse the map and match otus
    print("Parsing and matchin OTUs")
    bef_dic = map_to_dic(bef_map)
    aft_dic = map_to_dic(aft_map)
    otu_dic = map_to_dic(otu_map)
    id_otu = match_otu(bef_dic,aft_dic,otu_dic)
    #
    #
    print('writing otu table')
    #merge the dictionaries to daraframe:
    otu_df = pd.DataFrame([bef_dic,aft_dic,id_otu],index = ["bef","aft","id_OTU"])
    #
    #
    #save the output
    otu_df.T.to_csv(args.output+'[INFO]'+'OTU_info.csv',index=True)

