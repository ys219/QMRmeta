# Remove all objects from the environment
rm(list = ls())
# Set global options
options(stringsAsFactors = F)
#
library("ggplot2")

#load data
info = read.csv("../data/[INFO]count_uniqs_otus.txt")# numbers of reads, uniqs and pre/post OTUs in each samples 

sample_X_OTU = read.csv("../data/[INFO]sample_X_OTU_cleaned_ASVs.csv")# unlike the map generated with vsearch, this is the one with umber of ASVs in each cell,

metadata = read.csv("../data/metadata_info_v2.csv")# sampling information

info_cole = read.csv("../data/[INFO]count_uniqs_otus_coleop.csv")# same as info but coleoptra OTUs only

sample_X_OTU_cole = read.csv("../data/[INFO]sample_X_OTU_coleop_ASVs.csv")# same as sample_X_OTU but coleoptra OTUs only
metadata = metadata[,-5]#remove labels column

# standarise data 
metadata[which(metadata$location == "Feng city"),]$location = "Fengcity"

metadata = metadata[-which(metadata$location == "Foursites"),]

info = info[startsWith(info$X , "RNMB"),]

sample_X_OTU = sample_X_OTU[startsWith(sample_X_OTU$X, "RNMB"),]

sample_X_OTU = data.frame(X = sample_X_OTU$X, ASVs_in_OTU = rowSums(sample_X_OTU[,-1]))

info = info[,-4]; info = info[,-5] ; info = info[,-6:-9]
# uniqs in pre_filtered OTU; uniqs in pre_filtered OTU; and information among pry/post filtered OTUs are removed. 

# unify colnames
colnames(metadata)[1] = "X"

plot_data = merge(info, sample_X_OTU, by = "X")
plot_data = merge(plot_data, metadata, by = "X")

info_cole = info_cole[startsWith(info_cole$X , "RNMB"),]

sample_X_OTU_cole = sample_X_OTU_cole[startsWith(sample_X_OTU_cole$X, "RNMB"),]

sample_X_OTU_cole = data.frame(X = sample_X_OTU_cole$X, ASVs_in_OTU = rowSums(sample_X_OTU_cole[,-1]))

info_cole = info_cole[,-4]; info_cole = info_cole[,-5] ; info_cole = info_cole[,-6:-9]

plot_data_cole = merge(info_cole, sample_X_OTU_cole, by = "X")
plot_data_cole = merge(plot_data_cole, metadata, by = "X")

# write.csv(plot_data,"../../data/all_mb_data.csv",row.names = F)
# write.csv(plot_data_cole,"../../data/coleop_data.csv",row.names = F)

# unique vs reads
uniq_vs_reads = ggplot(data = plot_data, aes(x = reads, y = uniqs))+
  geom_point(size = 1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Number of reads(quality filtered)", y= "uniqs")+
  geom_smooth(method = 'lm')

uniq_vs_reads + facet_grid(cols = vars(location))

# uniq_vs_reads = ggplot(data = plot_data, aes(x = reads, y = uniqs))+
ASV_vs_reads =ggplot(data = plot_data, aes(x = reads, y = ASVs))+
  geom_point(size = 1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Number of reads(quality filtered)", y= "ASV")+
  geom_smooth(method = 'lm')

ASV_vs_reads + facet_grid(cols = vars(location))

# OTU vs reads 

otu_vs_reads =ggplot(data = plot_data, aes(x = reads, y = aft_OTU))+
  geom_point(size = 1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Number of reads(quality filtered)", y= "OTU")+
  geom_smooth(method = 'lm')

otu_vs_reads + facet_grid(cols = vars(location))
# correlations 
cor(plot_data$reads,plot_data$uniqs)
cor.test(plot_data$reads,plot_data$bef_OTU)
cor.test(plot_data$reads,plot_data$ASVs)
cor.test(plot_data$reads,plot_data$aft_OTU)

cor.test(plot_data$nparts,plot_data$uniqs)
cor.test(plot_data$nparts,plot_data$bef_OTU)
cor.test(plot_data$nparts,plot_data$ASVs)
cor.test(plot_data$nparts,plot_data$aft_OTU)

# gather data
plot_data = gather(plot_data,"data_type","frequency",c(3:6))
# all vs reads
all_vs_reads =ggplot(data = plot_data, aes(x = reads, y = frequency))+
  geom_point(size = 1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Number of reads(quality filtered)", y= "count")+
  geom_smooth(method = 'lm')

all_vs_reads + facet_grid(cols = vars(data_type))


# all vs nparts
all_vs_nparts =ggplot(data = plot_data, aes(x = nparts, y = frequency))+
  geom_point(size = 1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Number of specimens", y= "count")+
  geom_smooth(method = 'lm')

all_vs_nparts + facet_grid(cols = vars(data_type))
#
