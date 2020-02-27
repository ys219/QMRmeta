#!/bin/bash

#Author: Yige Sun
#
#NOTE:bash run in the directory where data folders stored(mbc/ or sandbox/)


hashline="################################################################" 
echo $hashline 
echo "This is the START of Data Filtering" `date` 

START=$(date +%s) #timming tool to record the running time


echo $hashline 
echo "STEP4_Quality_Filtering_by_library @ " `date` 
echo $hashline 

mkdir 04_error_filtered 

while read f ; do vsearch --fastx_filter 02_trimmed/$f --fastq_maxee 1 --fastaout 04_error_filtered/$f & done< <(ls 02_trimmed)

echo $hashline 
echo "check No. reads in 04_error_filtered.fasta @ " `date` 

wait

grep -c "^>" 04_error_filtered/*
NOW1=$(date +%s)
echo $hashline 
echo "STEP4_Quality_filtering took $(($NOW1 - $START))s" 

echo $hashline 
echo "STEP5_Dereplication @ " `date` 
echo $hashline 

mkdir 05_dereped

while read f ; do vsearch --derep_fulllength 04_error_filtered/$f --output 05_dereped/$f --sizeout --relabel uniq & done< <(ls 04_error_filtered)

wait

echo $hashline 
echo "check No. reads in 05_dereped@ " `date` 

grep -c "^>" 05_dereped/*

NOW2=$(date +%s)
echo $hashline 
echo "STEP5_Dereplication took $(($NOW2 - $NOW1))s" 


echo $hashline 
echo "STEP5.1_unwrapped the data @ " `date` 
echo $hashline 

mkdir 05_dereped_unwrapped ## have to be in a individual folder/ or everything got removed

while read f ; do perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' 05_dereped/$f > 05_dereped_unwrapped/$f & done< <(ls 05_dereped) ## perl exacutively do: if[?] current line number[$.] is greater than 1 and match ^> [/^>/ i.e.start with > symbol]  then print new line symbol["\n"] else [:] chop off the new line symbol if exist.

wait

echo $hashline 
echo "STEP6_indel_filtering @ " `date` 
echo $hashline 

mkdir 06_indel
while read f ; do vsearch --fastx_filter 05_dereped_unwrapped/$f --fastq_minlen 418 --fastq_maxlen 418 --fastaout 06_indel/$f & done< <(ls 05_dereped_unwrapped)

##vsearch --fastx_filter 05_derep_unwrapped.fasta --fastq_minlen 415 --fastq_maxlen 421 --fastaout 06_indel_418.fasta ## to keep an insertion/deletion of codon

wait

echo $hashline 
echo "check No. reads in 06_indel_418.fasta @ " `date` 

grep -c "^>" 06_indel/*
NOW3=$(date +%s)
echo $hashline 
echo "STEP6_indel_filtering took $(($NOW3 - $NOW2))s" 

echo $hashline 
echo "STEP7_denoising @ " `date` 
echo $hashline 


mkdir 07_denoise
while read f ; do vsearch --cluster_unoise 06_indel/$f --minsize 4 --unoise_alpha 2 --centroids 07_denoise/$f & done< <(ls 06_indel)

echo $hashline 
echo "check No. reads in 07_denoise.fasta @ " `date` 

wait

grep -c "^>" 07_denoise/*
NOW4=$(date +%s)
echo $hashline 
echo "STEP7_denoising took $(($NOW4 - $NOW3))s" 

echo $hashline 
echo "Number of sequences and size dirtribution now is" 
grep "^>" 07_denoise/* | sed -e "s/size=\([^;]\)/\1/"


echo $hashline 
echo "STEP8_point_error_filtering @ " `date` 
echo $hashline 
mkdir  08_point_error_filtered
while read f ; do filtertranslate.py 07_denoise/$f 5 & done< <(ls 07_denoise)

wait

mv *transpass* 08_point_error_filtered/ 

echo $hashline 
echo "check No. reads in 08_point_error_filtered @ " `date` 

grep -c "^>" 08_point_error_filtered/*

echo $hashline 
echo "check No. reads in 07_denoise_transfail.fa @ " `date` 

grep -c "^>" *transfail.fa

rm 07_denoise/*transfail.fa


NOW5=$(date +%s)
echo $hashline 
echo "STEP8_point_error_filtering took $(($NOW5 - $NOW4))s" 
echo $hashline 
echo "STEP9_chimera_filtering @ " `date` 
echo $hashline 
mkdir 09_ALL_filtered

while read f ; do vsearch --uchime3_denovo 08_point_error_filtered/$f --nonchimeras 09_ALL_filtered/$f & done< <(ls 08_point_error_filtered)

echo $hashline 
echo "check No. reads in 09_ALL_filtered.fasta @ " `date` 
echo $hashline 

wait

grep -c "^>" 09_ALL_filtered/*

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


