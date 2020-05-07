# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
import pandas as pd
import argparse

def extract_multiple_dictionary(input_map):
    """create one key multiple values dictionary"""
    count = 0
    with open(input_map,"r") as input_blast_map:
        uniq_dic = {}
        for l in input_blast_map: ## from each lib
            count += 1
            key = l.split()[0].split(';')[1].split('=')[1]## the filename
            value = l.split()[1].split(';')[0]##the uniq
            uniq_dic.setdefault(key,[]).append(value)
    return uniq_dic , count

def extract_dictionary(input_map):
    """create one key one value dictionary"""
    count = 0
    with open(input_map,"r") as input_blast_map:
        dic = {}
        for l in input_blast_map:
            count += 1
            key = l.split()[0].split(';')[0]
            value = l.split()[1]
            dic[key]=value
    return dic, count


def get_u_len(uniq_dic):
    """convert dictionary to length of it"""
    dic = {}
    for key in uniq_dic.keys():
        dic[key] = len(set(uniq_dic[key]))
    return dic

def get_len(uniq_dic):
    """convert dictionary to length of it"""
    dic = {}
    for key in uniq_dic.keys():
        dic[key] = len(uniq_dic[key])
    return dic

def get_unique(uniq_dic):
    """convert dictionary to unique contents"""
    for key in uniq_dic.keys():
        uniq_dic[key] = list(set(uniq_dic[key]))
    return uniq_dic

def replace_ignore(in_list,dictionary):
    """replace values in the list with given dictionary
    and ignore any thing not in dictionary"""
    rep_list = [dictionary[x] for x in in_list if x in dictionary]
    return rep_list

def replace_missing(in_list,dictionary):
    """retain values that are not in the given dictionary"""
    rep_list = [x for x in in_list if x not in dictionary]
    return rep_list

def replace_retain(in_list,dictionary):
    """replace values in the list with given dictionary
    and keep any thing not in dictionary in the list"""
    rep_list = [dictionary[x] if x in dictionary else x for x in in_list]
    return rep_list

def check_exist(in_list,dictionary):
    """check if values are in the given dictionary
    and ignore any thing not in dictionary"""
    out_list = [x for x in in_list if x in dictionary]
    return out_list

# def check_or_exist(index_dic,check_dic_1,check_dic_2):
#     out_dic = {}
#     for k in index_dic:
#         if k in check_dic_1.keys or k in check_dic_2.keys():
#             out_dic[k] = True
#         else :
#             out_dic[k] = False
#     return out_dic

# def check_and_exist(index_dic,check_dic_1,check_dic_2):
#     out_dic = {}
#     for k in index_dic:
#         if k in check_dic_1.keys and k in check_dic_2.keys():
#             out_dic[k] = True
#         else :
#             out_dic[k] = False
#     return out_dic

if __name__ == "__main__":
    # parser settings
    #
    parser = argparse.ArgumentParser(description = "This tool match the OTUs for filtered and unfiltered unique reads")
    # dereped input
    parser.add_argument("-s","--sample_dereped_input", help = "sample vs dereplicated, quality filtered map input file", metavar = "INPUT.txt//tsv")
    # otu before input
    parser.add_argument("-b","--otu_bef_map", help = "uniqs vs OTUs for unfiltered(but quality filtered and dereplcated) reads; geerated with vsearch in blast6out method", metavar = "OTU_bef.txt")
    # output file path with default
    parser.add_argument("-o","--output", help = "output filename ", default = "[INFO]sample_uniqs_OTU_count.csv", metavar = "OUTPUTFILENAME")# output
    # otu after input
    parser.add_argument("-a","--otu_aft_map", help = "uniqs vs. OTUs for filtered reads", metavar = "OTU_aft.txt")#
    parser.add_argument("-m","--otus_b_a_map", help = "OTU_bef vs OTU_aft map with blast6out method", metavar = "OTUs_map.txt")#
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
    #
    #
    #extract dic from origin map
    # oringin_map = open(args.sample_dereped_input,"r")
    sample_uniqs_dic, reads_count = extract_multiple_dictionary(args.sample_dereped_input)#
    print("Stats:")
    print("%d reads in %d samples were match to a uniq"%(reads_count, len(sample_uniqs_dic)))#0
    uniqs_total = len(set(sum(sample_uniqs_dic.values(),[])))
    print("%d uniqs detected"%(uniqs_total))#1
    #
    #extract dic from bef_OTU map
    # bef_map = open(args.otu_bef_map,"r")
    bef_dic, u_bef_count = extract_dictionary(args.otu_bef_map)
    print("%d uniqs were match to a bef_otu"%(len(bef_dic)))#2
    bef_otu = len(set(bef_dic.values()))
    print("%d bef_otu detected"%(bef_otu))#3
    #
    #extract dic from aft_OTU map
    # aft_map = open(args.otu_aft_map,"r")
    aft_dic, u_aft_count = extract_dictionary(args.otu_aft_map)
    print("%d uniqs were match to a aft_otu"%(len(aft_dic)))#4
    aft_otu = len(set(aft_dic.values()))
    print("%d aft_otu detected"%(aft_otu))#5
    print("%d uniqs match to at least one of the OTU"%(len(bef_dic)+ len(aft_dic) -len(set(bef_dic.keys())&set(aft_dic.keys()))))#6
    print("%d uniqs got both OTU_bef&aft"%(len(set(bef_dic.keys())&set(aft_dic.keys()))))#7
    print("%d uniqs do not match any OTU"%(uniqs_total- len(bef_dic)- len(aft_dic) +len(set(bef_dic.keys())&set(aft_dic.keys()))))#8
    # print("%d unique OTUs detected"%(bef_otu+aft_otu-)#9
    otu_match_dic, c = extract_dictionary(args.otus_b_a_map)
    #
    reads_count_dic = get_len(sample_uniqs_dic)#0
    sample_uniqs_dic = get_unique(sample_uniqs_dic)
    #replace uniqs in sample_uniq_dic with OTU_bef
    sample_befOTU_dic = {}
    for k in sample_uniqs_dic:
        sample_befOTU_dic[k] = replace_ignore(sample_uniqs_dic[k],bef_dic) # What happens if a uniq in the list isn't a key in the bef_dic
    #
    #replace uniqs in sample_uniq_dic with OTU_aft
    sample_aftOTU_dic = {}
    for k in sample_uniqs_dic:
        sample_aftOTU_dic[k] = replace_ignore(sample_uniqs_dic[k],aft_dic)
    ### for OTU nor
    # before move on:
    sample_m_befOTU_dic = {}
    for k in sample_uniqs_dic:
        sample_m_befOTU_dic[k] = replace_missing(sample_uniqs_dic[k],bef_dic)
    #
    sample_m_aftOTU_dic = {}
    for k in sample_uniqs_dic:
        sample_m_aftOTU_dic[k] = replace_missing(sample_uniqs_dic[k],aft_dic)
    #output data:
    sample_uniqs_count_dic = get_u_len(sample_uniqs_dic)#1
    u_bef_count_dic = get_len(sample_befOTU_dic)#2
    sample_bef_count_dic = get_u_len(sample_befOTU_dic)#3
    u_aft_count_dic = get_len(sample_aftOTU_dic)#4
    sample_aft_count_dic = get_u_len(sample_aftOTU_dic)#5
    #6 7 8
    u_OTU_and = {}
    u_OTU_or = {}#6 empty dic
    # u_OTU_and = {}
    u_OTU_nor = {}#8
    #
    for k in sample_uniqs_dic:## fill in dictionaries
        tmp = check_exist(sample_uniqs_dic[k],bef_dic) + check_exist(sample_uniqs_dic[k],aft_dic)
        u_OTU_and[k] = len(set([x for x in tmp if tmp.count(x) > 1 ]))
        u_OTU_or[k] = len(set(tmp))
        u_OTU_nor[k] = sample_uniqs_count_dic[k]- u_OTU_or[k]
        # if sample_befOTU_dic[k] != [] and sample_aftOTU_dic[k] != []:
        #     u_OTU_and[k] = len(sample_befOTU_dic[k]+sample_aftOTU_dic[k])[]
    #9
    master_OTUs = {}
    for k in sample_befOTU_dic:
        sample_befOTU_dic[k] = replace_retain(sample_befOTU_dic[k],otu_match_dic)
        master_OTUs[k] = len(set(sample_befOTU_dic[k]+sample_aftOTU_dic[k]))
    out_df = pd.DataFrame([reads_count_dic, sample_uniqs_count_dic, u_bef_count_dic, sample_bef_count_dic, u_aft_count_dic, sample_aft_count_dic, u_OTU_or, u_OTU_and, u_OTU_nor, master_OTUs ],index = ["reads","uniqs","uniqs_in_bef_OTU","bef_OTU","uniqs_in_aft_OTU","aft_OTU","Got_a_OTU","Got_both","Got_nor","master_OTU"])
    #
    out_df.T.to_csv(args.output,index=True)
    # exit()