#!/bin/bash

##NOTE:bash run in the directory where data folders stored(mbc/ or sandbox/)
echo "######################################################" >> $output
echo "This is the START of Read Processing" `date` >> $output

mkdir 01_trimmed 02_merged #create directories to store data

output=../output_log.txt #where to store the ouput logs

START=$(date +%s) #timming tool to record the running time
echo "######################################################" >> $output
echo "check No. reads in 0_demux @ " `date` >> $output
echo "######################################################" >> $output

grep -c "^@" 0_demux/*.fq | tee -a $output

echo "######################################################" >> $output
echo "STEP1_Primer_Removal @ " `date` >> $output
echo "######################################################" >> $output 

while read s; do cutadapt -g ^CCNGAYATRGCNTTYCCNCG -G ^TANACYTCNGGRTGNCCRAARAAYCA -o 01_trimmed/${s}_R1.fq -p 01_trimmed/${s}_R2.fq --discard-untrimmed 0_demux/${s}_R1.fq 0_demux/${s}_R2.fq; done < <(ls 0_demux | rev | cut -d_ -f 2- | rev | sort | uniq) | tee -a $output



echo "######################################################" >> $output
echo "check No. reads in 01_trimmed @ " `date` >> $output
echo "######################################################" >> $output

grep -c "^@" 01_trimmed/*.fq | tee -a $output
NOW1=$(date +%s)
echo "######################################################" >> $output
echo "STEP1_STEP1_Primer_Removal took $(($NOW1 - $START))s" >> $output

echo "######################################################" >> $output
echo "STEP2_Pair_Merging @ " `date` >> $output
echo "######################################################" >> $output

while read l; do pear -f 01_trimmed/${l}_R1.fq -r 01_trimmed/${l}_R2.fq -o 02_merged/$l -q 26 -v 100; done < <(ls 01_trimmed | rev | cut -d_ -f 2- | rev | sort | uniq) | tee -a $output


echo "######################################################" >> $output
echo "check No. reads in 02_merged @ " `date` >> $output
echo "######################################################" >> $output

grep -c "^@" 02_merged/* | tee -a $output

rm 02_merged/*discarded* 02_merged/*unassembled* && rename -e "s/assembled\.//" 02_merged/*

echo "######################################################" >> $output
echo "check No. reads in 02_merged-after cleaning up @ " `date` >> $output
echo "######################################################" >> $output

grep -c "^@" 02_merged/* | tee -a $output

NOW2=$(date +%s)
echo "######################################################" >> $output
echo "STEP2_Pair_Merging took $(($NOW2 - $NOW1))s" >>$output

echo "######################################################" >> $output
echo "STEP3_Concatenation @ " `date` >> $output
echo "######################################################" >> $output


cd 02_merged

for f in *; do sed -e "s/\(^@.*\) .*$/\1;sample=${f%.*};/" $f >> ../03_mbc_concat.fastq; done 

cd ../
echo "######################################################" >> $output
echo "check No. reads in 03_concat.fastq @ " `date` >> $output
echo "######################################################" >> $output

grep -c "^@" 03_mbc_concat* | tee -a $output


NOW3=$(date +%s)
echo "######################################################" >> $output
echo "STEP3_Concatenation took $(($NOW3 - $NOW2))s" >>$output

echo "######################################################" >> $output
echo "This is the END of Read Processing" `date` >> $output
echo "######################################################" >> $output

END=$(date +%s)
echo "Read Processing took $(($END - $START))s" >>$output
echo "######################################################" >> $output

