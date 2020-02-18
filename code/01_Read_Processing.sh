#!/bin/bash

hashline="#################################################################"
##NOTE:bash run in the directory where data folders stored(mbc/ or sandbox/)
## this code has been editted from primer removal -> pair merging script, original codes were commented
 
echo $hashline 
echo "This is the START of Read Processing" `date` 

# mkdir 01_trimmed 02_merged #create directories to store data

mkdir 01_merged 02_trimmed 

output=../output_log.txt #where to store the ouput logs

START=$(date +%s) #timming tool to record the running time
echo $hashline 
echo "check No. reads in 0_demux @ " `date` 
echo $hashline 

grep -c "^@" 0_demux/*.fq 

echo $hashline 
# echo "STEP1_Primer_Removal @ " `date` 
echo "STEP1_Pair_Merging @ " `date` 
echo $hashline  

# while read s; do cutadapt -g ^CCNGAYATRGCNTTYCCNCG -G ^TANACYTCNGGRTGNCCRAARAAYCA -o 01_trimmed/${s}_R1.fq -p 01_trimmed/${s}_R2.fq --discard-untrimmed 0_demux/${s}_R1.fq 0_demux/${s}_R2.fq; done < <(ls 0_demux | rev | cut -d_ -f 2- | rev | sort | uniq) 

while read f; do pear -f 0_demux/${f}_R1.fq -r 0_demux/${f}_R2.fq -o 01_merged/${f} -q 26 -v 100 & done < <(ls 0_demux | rev | cut -d_ -f2- | rev | sort | uniq) # & before done make it parallel



wait ## ensure all the parallel commands finished then do the next step
echo $hashline 
# echo "check No. reads in 01_trimmed @ " `date`
echo "check No. reads in 01_merged @ " `date`
echo $hashline 

 
# grep -c "^@" 01_trimmed/*.fq 
#grep -c "^@" 01_merged/*
#NOW1=$(date +%s)
#echo $hashline 

rm 01_merged/*discarded* 01_merged/*unassembled* && rename -e "s/assembled\.//" 01_merged/*

grep -c "^@" 01_merged/*
NOW1=$(date +%s)
echo $hashline
echo "STEP1_Pair_merging took $(($NOW1 - $START))s" 
echo $hashline 
echo "STEP2_Primer_removal @ " `date` 
echo $hashline 

# while read l; do pear -f 01_trimmed/${l}_R1.fq -r 01_trimmed/${l}_R2.fq -o 02_merged/$l -q 26 -v 100; done < <(ls 01_trimmed | rev | cut -d_ -f 2- | rev | sort | uniq) 


while read f; do cutadapt -j 10 -g CCNGAYATRGCNTTYCCNCG...TGRTTYTTYGGNCAYCCNGARGTNTA -o 02_trimmed/${f} --discard-untrimmed 01_merged/${f}; done < <(ls 01_merged) 


wait 
echo $hashline 
echo "check No. reads in 02_primer_removal @ " `date` 
echo $hashline 

grep -c "^@" 02_trimmed/* 

# rm 02_merged/*discarded* 02_merged/*unassembled* && rename -e "s/assembled\.//" 02_merged/*

# echo $hashline 
# echo "check No. reads in 02_merged-after cleaning up @ " `date` 
# echo $hashline 

# grep -c "^@" 02_merged/* 


NOW2=$(date +%s)
echo $hashline 
echo "STEP2_Primer_Removal took $(($NOW2 - $NOW1))s" 
echo $hashline
# echo "checking failed files, following files have 0 reads in it"
# echo $hashline
# echo < <(grep -c "^@" 02_merged/* | rev | while read n; do if [[ "${n%:*}" = "0" ]] ; then echo $n; fi ; done | rev) 
##print the failed file names< <grep the file name with no. of reads in 02_merged/ directory | turn the file name around | read through the informations; if the no. of reads("number" before:) is 0 which means merge failed; then print the file name; finish if condition; done in while loop | turn the name back
# mkdir 1_merge_first

# echo $hashline
# echo "Doing pair merge followed by primer removal"
# echo $hashline

# while read f; do pear -f 0_demux/${f}_R1.fq -r 0_demux/${f}_R2.fq -o 1_merge_first/${f} -q 26 -v 100; done < <(grep -c "^@" 02_merged/* | rev | while read n; do if [[ "${n%:*}" = "0" ]] ; then echo $n; fi ; done | rev | cut -d. -f 1 | cut -d/ -f 2)
##do the pair merging with pear for all piped in files(failed files)< <grep the file name with no. of reads in 02_merged/ directory | turn the file name around | read through the informations; if the no. of reads("number" before:) is 0 which means merge failed; then print the file name; finish if condition; done in while loop | turn the name back | cut and get contents before  .fastq

# echo $hashline
# echo "checking files in 1_merge_first"

# grep -c "^@" 1_merge_first/*

# echo $hashline
# rm 1_merge_first/*discarded* 1_merge_first/*unassembled* && rename -e "s/assembled\.//" 1_merge_first/*

# echo $hashline
# echo "Check no. of reads in 1_merge_first after cleaning"
# grep -c "^@" 1_merge_first/*
# echo $hashline

##now cutadapt check if the primer are correct!! save to 02_merged
#while read f; do cutadapt -g CCNGAYATRGCNTTYCCNCG...TGRTTYTTYGGNCAYCCNGARGTNTA -o 02_merged/${f} --discard-untrimmed 1_merge_first/${f}; done < <(ls 1_merge_first) 
#echo $hashline
#echo "check No.reads in 02_merged @ " `date`
#echo $hashline

#grep -c "^@" 02_merged/*

#NOW3=$(date +%s)
#echo $hashline
#echo "STEP1&2 took $(($NOW3 - $NOW2))s"

#echo $hashline 
#echo "STEP3_Concatenation @ " `date` 
#echo $hashline 

#### maybe I should move concatenation to bulk filtering script
#cd 02_merged

#for f in *; do sed -e "s/\(^@.*\) .*$/\1;sample=${f%.*};/" $f >> ../03_mbc_concat.fastq; done 

#cd ../
#echo $hashline 
#echo "check No. reads in 03_concat.fastq @ " `date` 
#echo $hashline 

#grep -c "^@" 03_mbc_concat* 


#NOW4=$(date +%s)
#echo $hashline 
#echo "STEP3_Concatenation took $(($NOW4 - $NOW3))s" 

echo $hashline 
echo "This is the END of Read Processing" `date` 
echo $hashline 

END=$(date +%s)
echo "Read Processing took $(($END - $START))s" 
echo $hashline 

