#!/bin/bash

vsearch --cluster_size 09_ALL_filtered.fasta -id 0.97 --centroids 10_otus_3pc.fasta --sizein --relabel otu
