#!/bin/bash

#NOTE:bash run in the directory where data folders stored(mbc/ or sandbox/)

echo "######################################################" >> $output
echo "This is the START of Data Filtering" `date` >> $output

output=../output_log.txt #where to store the ouput logs

START=$(date +%s) #timming tool to record the running time

echo "######################################################" >> $output
echo "check No. reads in 03_mbc_concat.fastq @ " `date` >> $output

grep -c "^@" 03_mbc*.fastq | tee -a $output

echo "######################################################" >> $output
echo "STEP4_Quality_Filtering @ " `date` >> $output
echo "######################################################" >> $output

vsearch --fastx_filter 03_mbc_concat.fastq --fastq_maxee 1 --fastaout 04_error_filtered.fasta | tee -a $output

echo "######################################################" >> $output
echo "check No. reads in 04_error_filtered.fasta @ " `date` >> $output

grep -c "^>" 04_error_filtered.fasta | tee -a $output
NOW1=$(date +%s)
echo "######################################################" >> $output
echo "STEP4_Quality_filtering took $(($NOW1 - $START))s" >> $output

echo "######################################################" >> $output
echo "STEP5_Dereplication @ " `date` >> $output
echo "######################################################" >> $output

vsearch --derep_fulllength 04_error_filtered.fasta --output 05_dereped.fasta --sizeout --relabel uniq | tee -a $output


echo "######################################################" >> $output
echo "check No. reads in 05_dereped.fasta @ " `date` >> $output

grep -c "^>" 05_dereped.fasta | tee -a $output
NOW2=$(date +%s)
echo "######################################################" >> $output
echo "STEP5_Dereplication took $(($NOW2 - $NOW1))s" >> $output


echo "######################################################" >> $output
echo "STEP5.1_unwrapped the data @ " `date` >> $output
echo "######################################################" >> $output

perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' 05_dereped.fasta > 05_derep_unwrapped.fasta ## perl exacutively do: if[?] current line number[$.] is greater than 1 and match ^> [/^>/] then print new line symbol["\n"] else [:] chop off the new line symbol.

echo "######################################################" >> $output
echo "STEP6_indel_filtering @ " `date` >> $output
echo "######################################################" >> $output

vsearch --fastx_filter 05_derep_unwrapped.fasta --fastq_minlen 418 --fastq_maxlen 418 --fastaout 06_indel_418.fasta

##vsearch --fastx_filter 05_derep_unwrapped.fasta --fastq_minlen 415 --fastq_maxlen 421 --fastaout 06_indel_418.fasta


echo "######################################################" >> $output
echo "check No. reads in 06_indel_418.fasta @ " `date` >> $output

grep -c "^>" 06_indel_418.fasta | tee -a $output
NOW3=$(date +%s)
echo "######################################################" >> $output
echo "STEP6_indel_filtering took $(($NOW3 - $NOW2))s" >> $output

echo "######################################################" >> $output
echo "STEP7_denoising @ " `date` >> $output
echo "######################################################" >> $output

vsearch --cluster_unoise 06_indel_418.fasta --minsize 4 --unoise_alpha 2 --centroids 07_denoise.fasta

echo "######################################################" >> $output
echo "check No. reads in 07_denoise.fasta @ " `date` >> $output

grep -c "^>" 07_denoise.fasta | tee -a $output
NOW4=$(date +%s)
echo "######################################################" >> $output
echo "STEP7_denoising took $(($NOW4 - $NOW3))s" >> $output

echo "######################################################" >> $output
echo "Number of sequences and size dirtribution now is" >> $output
grep "^>" 07_denoise.fasta | sed -e "s/size=\([^;]\)/\1/" | tee -a $output

echo "######################################################" >> $output
echo "STEP8_point_error_filtering @ " `date` >> $output
echo "######################################################" >> $output

filtertranslate.py 07_denoise.fasta 5 | tee -a $output

mv 07_denoise_transpass.fa 08_point_error_filtered.fa 

echo "######################################################" >> $output
echo "check No. reads in 08_point_error_filtered @ " `date` >> $output

grep -c "^>" 08_point_error*.fa | tee -a $output

echo "######################################################" >> $output
echo "check No. reads in 07_denoise_transfail.fa @ " `date` >> $output

grep -c "^>" 07_denoise_transfail.fa | tee -a $output

rm 07_denoise_transfail.fa


NOW5=$(date +%s)
echo "######################################################" >> $output
echo "STEP8_point_error_filtering took $(($NOW5 - $NOW4))s" >> $output
echo "######################################################" >> $output
echo "STEP9_chimera_filtering @ " `date` >> $output
echo "######################################################" >> $output


vsearch --uchime3_denovo 08_point_error_filtered.fa --nonchimeras 09_ALL_filtered.fasta

echo "######################################################" >> $output
echo "check No. reads in 09_ALL_filtered.fasta @ " `date` >> $output
echo "######################################################" >> $output

grep -c "^>" 09_ALL_filtered.fasta | tee -a $output

NOW6=$(date +%s)
echo "######################################################" >> $output
echo "STEP9_chimera_filtered took $(($NOW6 - $NOW5))s" >> $output
echo "######################################################" >> $output

echo "######################################################" >> $output
echo "This is the END of Data Filtering" `date` >> $output
echo "######################################################" >> $output

END=$(date +%s)
echo "Data Filtering took $(($END - $START))s" >>$output
echo "######################################################" >> $output
## NAPselect -output . -blasttarget Coleoptera -seqlength 418 -base_var 3 -var_by_codon 09_ALL_filtering_done.fasta


