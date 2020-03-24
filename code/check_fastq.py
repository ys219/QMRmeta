# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-


import subprocess
import pandas as pd
from Bio import SeqIO
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

import argparse
import os
import sys
import warnings

f = open("04_RNMB_dereped_test.fasta")
f.close()


## frequency

# for l in f : ## to match the index
#     if l.startswith('>') :
#         l = l.rstrip()
#         index = l.split(';')[0].replace('>','')
#         freq = l.split('=')[1]
#         out_df.loc[index, 'freq'] = freqquit
def freq_check(in_head):
        head = in_head.rstrip()
        head = head.split('=')[1])
    return head



## length
def length_check(in_seq):
    length = (len(in_seq))
    return length


# len(l)
# for l in f:
#     if l.startswith('>uniq') and block :# if it's the name line and block got contents
#         seq = l# merge all the contents in the block
#         block = [] # empty the block
#     elif seq != None : 
#         seq = seq.join()
#             # len_list.append(len(seq))
#             # seq = None
#     else:
#         block.append(l.strip())
# out_df.index

# len_list








# for l in f : ## to make list of data then append to df
#     if not l.startswith('>'):
#         l = l.rstrip()
#         freq_list.append(len(l))


# out_df['len'] = len_list





##  stop codons
subprocess.Popen(['filtertranslate.py','04_RNMB_dereped_test.fasta','5','-f_failed']) # do the stop codon detection with `filtertranslate.py`

# do the results extraction
f = open("04_RNMB_dereped_test_transfiltered.fa")
translate = []
for l in f :
    if l.startswith('>'):
        if l.endswith('_failed\n'):
            translate.append('FAIL')
        else:
            translate.append('P')

f.close()

out_df['filter_trans'] = translate


# chimera
subprocess.Popen(['vsearch', '--uchime3_denovo ', '04_RNMB_dereped_test.fasta', '--chimeras', '06_chime.fasta', '--nonchimeras', '06_non_chim.fasta'
])
chim_inf = []
# with openC as chime_f :
#     for l in chime_f :
#         if l.startswith('>uniq'):
#             chim_inf.append(l.rstrip()+'_chimera')

# with open('06_non_chim.fasta') as non_chime_f :
#     for l in non_chime_f :
#         if l.startswith('>uniq'):
#             chim_inf.append(l.rstrip())

out_df['chime'] = 'False'
f = open('06_chime.fasta')
for l in f : 
    if l.startswith('>') :
        index = l.split(';')[0].replace('>','')
        out_df.loc[out_df.name == index ,'chime'] = 'True'


# out_df.loc[out_df.name == index,'chime'] = 'True'

out_df.head



out_df.to_csv(r'08_check_fasta.csv',index= False, header = True)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "")

    parser.add_argument("input", help = "input file path", metavar = "FASTA")#input

    parser.add_argument("-o","--output_directory", help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")# output
    
    args = parser.parse_args()
    
    filename = os.path.splitext(os.path.basename(args.input))[0]

    with open(args.input) as infasta : 
        uniqs = [] # to save row names
        for l in infasta : #get row names
            if l.startswith('>'):
                uniqs.append(l.split(';')[0].replace('>','')) # split by ; symbol, keep the first part then remove the > symbol
        out_df = pd.DataFrame(data = uniqs, columns = ['name'])

        freq_list = []
        len_list = []

        for head, seq in SimpleFastaParser(infasta):

            freq_list.append(freq_check(head))

            len_list.append(length_check(seq))
        
        out_df['freq'] = freq_list
        out_df['length'] = len_list

        subprocess.Popen(['vsearch','--uchime3_denovo',args.input,'--nonchimeras',filename+'_chimeras.fa'])
        out_df['chime'] = 'False'
        with open(filename+'_chimeras.fa') as chimes:
            for l in chimes : 
                if l.startswith('>') :
                    index = l.split(';')[0].replace('>','')
                    out_df.loc[out_df.name == index ,'chime'] = 'True'
    


        


            






## create dataframe:




# append output

