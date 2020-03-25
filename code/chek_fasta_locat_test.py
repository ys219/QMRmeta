with open('05_dereped_filtered.fasta') as infasta : 
    filename = '05_dereped_filtered'
    uniqs = [] # to save row names
    freq_list = []
    len_list = []
    stop_list = []
    for head, seq in SimpleFastaParser(infasta):
        uniqs.append(head.split(';')[0].replace('>',''))# split by ; symbol, keep the first part then remove the > symbol
        freq_list.append(freq_check(head))
        len_list.append(length_check(seq))
        stop_codon = min(stopcount(SeqRecord(Seq(seq)),5, frame = (1,2,3)))
        stop_list.append(stop_codon)
    out_df = pd.DataFrame(data = uniqs, columns =['name'])
    out_df['freq'] = freq_list
    out_df['length'] = len_list
    out_df['stop_codon'] = stop_list

subprocess.Popen(['vsearch','--uchime3_denovo','05_dereped_filtered.fasta','--chimeras',filename+'_chimeras.fa'])
out_df['chime'] = 'False'
with open(filename+'_chimeras.fa') as chimes:
    for head, seq in SimpleFastaParser(chimes):
        index = head.split(';')[0].replace('>','')
        out_df.loc[out_df.name == index ,'chime'] = 'True'

print('done')



