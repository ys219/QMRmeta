#!/bin/bash

#NOTE:bash run in the directory where data folders stored(mbc/ or sandbox/)
hashline="################################################################" 
echo $hashline 
echo "This is the START of Data Filtering" `date` 



START=$(date +%s) #timming tool to record the running time

cd 02_trimmed

for f in *; do sed -e "s/\(^@.*\) .*$/\1;sample=${f%.*};/" $f >> ../03_mbc_concat.fastq; done

cd ../
echo $hashline 
echo "check No. reads in 03_concat.fastq @ " `date`
echo $hashline 

grep -c "^@" 03_mbc_concat*


NOW0=$(date +%s)
echo $hashline 
echo "STEP3_Concatenation took $(($NOW4 - $NOW3))s" 


echo $hashline 
echo "check No. reads in 03_mbc_concat.fastq @ " `date` 

grep -c "^@" 03_mbc*.fastq

echo $hashline 
echo "STEP4_Quality_Filtering @ " `date` 
echo $hashline 

vsearch --fastx_filter 03_mbc_concat.fastq --fastq_maxee 1 --fastaout 04_error_filtered.fasta

echo $hashline 
echo "check No. reads in 04_error_filtered.fasta @ " `date` 

grep -c "^>" 04_error_filtered.fasta
NOW1=$(date +%s)
echo $hashline 
echo "STEP4_Quality_filtering took $(($NOW1 - $NOW0))s" 

echo $hashline 
echo "STEP5_Dereplication @ " `date` 
echo $hashline 

vsearch --derep_fulllength 04_error_filtered.fasta --output 05_dereped.fasta --sizeout --relabel uniq


echo $hashline 
echo "check No. reads in 05_dereped.fasta @ " `date` 

grep -c "^>" 05_dereped.fasta
NOW2=$(date +%s)
echo $hashline 
echo "STEP5_Dereplication took $(($NOW2 - $NOW1))s" 


echo $hashline 
echo "STEP5.1_unwrapped the data @ " `date` 
echo $hashline 

perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' 05_dereped.fasta > 05_derep_unwrapped.fasta ## perl exacutively do: if[?] current line number[$.] is greater than 1 and match ^> [/^>/] then print new line symbol["\n"] else [:] chop off the new line symbol.

echo $hashline 
echo "STEP6_indel_filtering @ " `date` 
echo $hashline 

vsearch --fastx_filter 05_derep_unwrapped.fasta --fastq_minlen 418 --fastq_maxlen 418 --fastaout 06_indel_418.fasta

##vsearch --fastx_filter 05_derep_unwrapped.fasta --fastq_minlen 415 --fastq_maxlen 421 --fastaout 06_indel_418.fasta


echo $hashline 
echo "check No. reads in 06_indel_418.fasta @ " `date` 

grep -c "^>" 06_indel_418.fasta
NOW3=$(date +%s)
echo $hashline 
echo "STEP6_indel_filtering took $(($NOW3 - $NOW2))s" 

echo $hashline 
echo "STEP7_denoising @ " `date` 
echo $hashline 

vsearch --cluster_unoise 06_indel_418.fasta --minsize 4 --unoise_alpha 2 --centroids 07_denoise.fasta

echo $hashline 
echo "check No. reads in 07_denoise.fasta @ " `date` 

grep -c "^>" 07_denoise.fasta
NOW4=$(date +%s)
echo $hashline 
echo "STEP7_denoising took $(($NOW4 - $NOW3))s" 

echo $hashline 
echo "Number of sequences and size dirtribution now is" 
grep "^>" 07_denoise.fasta | sed -e "s/size=\([^;]\)/\1/"

echo $hashline 
echo "STEP8_point_error_filtering @ " `date` 
echo $hashline 

filtertranslate.py 07_denoise.fasta 5

mv 07_denoise_transpass.fa 08_point_error_filtered.fa 

echo $hashline 
echo "check No. reads in 08_point_error_filtered @ " `date` 

grep -c "^>" 08_point_error*.fa

echo $hashline 
echo "check No. reads in 07_denoise_transfail.fa @ " `date` 

grep -c "^>" 07_denoise_transfail.fa

rm 07_denoise_transfail.fa


NOW5=$(date +%s)
echo $hashline 
echo "STEP8_point_error_filtering took $(($NOW5 - $NOW4))s" 
echo $hashline 
echo "STEP9_chimera_filtering @ " `date` 
echo $hashline 


vsearch --uchime3_denovo 08_point_error_filtered.fa --nonchimeras 09_ALL_filtered.fasta

echo $hashline 
echo "check No. reads in 09_ALL_filtered.fasta @ " `date` 
echo $hashline 

grep -c "^>" 09_ALL_filtered.fasta

NOW6=$(date +%s)
echo $hashline 
echo "STEP9_chimera_filtered took $(($NOW6 - $NOW5))s" 
echo $hashline 

echo $hashline 
echo "This is the END of Data Filtering" `date` 
echo $hashline 

END=$(date +%s)
echo "Data Filtering took $(($END - $START))s" >>$output
echo $hashline 
## NAPselect -output . -blasttarget Coleoptera -seqlength 418 -base_var 3 -var_by_codon 09_ALL_filtering_done.fasta


