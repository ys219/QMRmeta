options(stringsAsFactors = F)
rm(list = ls())# clear workspace
d1=read.csv('../data/ruie_demultiplexing_secondrun_metabarcode_-well_informations06_22_2020.csv',header = T,stringsAsFactors = FALSE)# import data
d2=read.csv('../data/meta_full_for_RN.csv', header = T, stringsAsFactors = FALSE)# import data
d1 = subset(d1,select = -c(plate,well,Tag.number,tagF,tagR,unrareads,OUT.No.))
d2 = subset(d2,select = -c(seqset,plate,well,unrarereads))

colnames(d1)= colnames(d2)# rename 
d2 = na.omit(d2)
d3= plyr::rbind.fill(d1,d2)
sp_fun = function(input){
  out = strsplit(input,"_")[[1]][2]
  return(out)
}


contents = sapply(d3$seqname,sp_fun)
plate_contents = ifelse(contents == 'Mix', "Samples from multiple sources", contents)
plate_contents = ifelse(startsWith(plate_contents,"FX"), "Feng samples only", paste(plate_contents , "samples only"))

d3$plate_contents = plate_contents

write.csv(d3,"../data/metadata_info_v2.csv",row.names = FALSE)


# d_d = d3[which(d3$seqname == d3$seqname[duplicated(d3$seqname)][1]),]
# 
# for(i in d3$seqname[duplicated(d3$seqname)][2:12]){
#   d_d = rbind(d_d,d3[which(d3$seqname == i),])
#   }
# write.csv(d_d,"../data/metadata_info_duplicated.csv",row.names = FALSE)
# dup_label=d1[duplicated(d1$label)[1],]
# for(i in d1$label[duplicated(d1$label)]){
#   dup_label = rbind(dup_label,d1[which(d1$label == i),])
# }
