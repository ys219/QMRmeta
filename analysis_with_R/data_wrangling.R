# Remove all objects from the environment
rm(list = ls())
# Set global options
options(stringsAsFactors = F)
# packages ----
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

#import data ----
taxonomy = read.table("../../data/10_midori_sintax_1.tsv",sep = "\t") # 
colnames(taxonomy) <- c("uniq","taxonomy", "strand", "selected") #
#sintax output

info = read.table("../../data/03_bcmb_concat_info.csv", sep = ",", header = TRUE)
names(info)[1] = "uniq";rownames(info) = info$uniq
# info is the output of check fasta

otu = read.table("../../data/10_des_otus_table.csv", sep = ",", header = TRUE)
names(otu)[1] = "uniq";rownames(otu) = otu$uniq
# match OTU output
# taxonomy ----
taxonomy$uniq = sub(";.*","",taxonomy$uniq) # cleaning: rm size= in row names
taxdetailed = separate(taxonomy,taxonomy,into = c('k','p','c','o','f','g','s'),sep = ",")%>%# split by taxon
  gather(.,n,inf,2:8)%>%#put them by row
  select(uniq,inf)%>%
  separate(.,inf,into = c('level','taxon','score',NA),sep = "[:()]")%>%# split woth special symmbols
  unique(.)# rm repeated rows
# which were selected


taxdetailed$selected = taxdetailed$score == "1.0000"
taxselect = separate(taxonomy,selected,into = c('k','p','c','o','f','g','s'),sep = ",",fill = 'left')%>%# split by taxon
  select(uniq,s)%>%# keep first and last col
  separate(s,into = c('sel_level','sel_taxon'),sep = ":")
rownames(taxselect) = taxselect$uniq

# Rename the taxonomic levels
taxdetailed$level = mapvalues(taxdetailed$level,c('k','p','c','o','f','g','s'),c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
taxselect$sel_level = mapvalues(taxselect$sel_level,c('k','p','c','o','f','g','s'),c("kingdom", "phylum", "class", "order", "family", "genus", "species"))

#  overwrite taxonomy with taxon inf
taxonomy = taxdetailed[1:3] %>% spread(level,taxon)
# Set the row names
# row.names(taxonomy) = taxonomy$uniq
# rm last col
taxonomy=taxonomy[,-ncol(taxonomy)]
# sort col and row:
taxonomy = taxonomy[,c("uniq","kingdom", "phylum", "class", "order", "family", "genus", "species")]


# master table ----
df_list <- list(info,otu,taxonomy,taxselect)
master_table <- Reduce(function(x,y) merge(x,y,by = "uniq"),df_list)

# write.csv(master_table,"../data/11_analysis.csv",row.names = FALSE)

rm(list = setdiff(ls(),"master_table"))#clear everything other than coord

####### load other data #######
metadata = read.table("../../data/metadata_info_v2.csv",header = T, sep = ',')

samp_unfil_uniq = read.table("../../data/03_bcmb_concat_origin_map.tsv", header = T, sep = "\t", row.names = 1, comment.char = '')
