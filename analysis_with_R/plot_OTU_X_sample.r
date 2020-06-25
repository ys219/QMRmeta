
rm(list = ls())
library("ggplot2")
options(stringsAsFactors = F)
# sample_X_OTU = read.csv("../../data/[INFO]sample_X_OTU.csv")

sample_X_OTU = read.csv("../data/[INFO]sample_X_OTU_cleaned_ASVs.csv")

info = read.csv("../data/metadata_info_v2.csv")

rnbc = sample_X_OTU[startsWith(sample_X_OTU$X,"RNBC"),]
rnmb = sample_X_OTU[startsWith(sample_X_OTU$X,"RNMB"),]

row.names(rnbc) = rnbc$X 
row.names(rnmb) = rnmb$X

# plot1 = data.frame(bio_div = c(rowSums(rnbc[,-1]),rowSums(rnmb[,-1])),
#            samp = c(rnbc$X,rnmb$X))
# row.names(plot1)= plot1$samp

# ggplot(plot1,aes(x = samp ,fill = samp, y = bio_div))+
#          scale_y_log10()+
#          geom_violin()+
#          geom_boxplot(position = position_dodge(width = 0.9),width = 0.3)

plot2 = data.frame(bio_div = rowSums(rnmb[,-1]),
                   samp = rnmb$X)
row.names(plot2) = plot2$samp

mb_names = unique(info[,1:2])[startsWith(unique(info[,1:2])[,1],"RNMB"),]
mb_names$location[which(mb_names$location == "Feng city")] = "Fengcity"## unify the names
mb_names = mb_names[-which(mb_names$location == "Foursites"),] # drop 4sites rows
mb_names = unique(mb_names); colnames(mb_names)[1] = "samp"


plot2 = merge(plot2, mb_names, by = "samp"
)

# ggplot(plot2,aes(x = location ,fill = location, y = bio_div))+ 
#   scale_y_log10()+
#   geom_violin()+
#   geom_boxplot(position = position_dodge(width = 0.9),width = 0.2) +
#   labs(y = "unfiltered_uniqs")

ggplot(plot2,aes(x = location ,fill = location, y = bio_div))+ 
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.2) +
  labs(y = "ASVs")
