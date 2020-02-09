#!/bin/bash

## otu delimitation
vsearch --cluster_size 09_ALL_filtered.fasta -id 0.97 --centroids 10_otus_3pc.fasta --sizein --relabel otu


## map the results
vsearch --usearch_global 09_ALL_filtered.fasta -db 10_otus_3pc.fasta -id 0.97 -otutabout 11_mapresults.tsv


