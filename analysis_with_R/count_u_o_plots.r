# Remove all objects from the environment
rm(list = ls())
# Set global options
options(stringsAsFactors = F)
#
library("ggplot2")
library("tidyr")
library("cowplot")

info = read.csv("../data/[INFO]count_uniqs_otus.txt")
info$type = substring(info$X,3,4)
info = info[order(info$X),]
info$index = seq(1:nrow(info))

info = info%>%
  gather(.,"OTU_type","no_OTU",c(5,7,11)) %>%
  gather(.,"uniq_status","no_uniqs",6:8)

## need to check normalisation
nor_info = info
nor_info$reads = mean(info$reads)/nor_info$reads
nor_info$uniqs = nor_info$reads*nor_info$uniqs
nor_info$no_uniqs = nor_info$reads*nor_info$no_uniqs
nor_info$no_OTU = nor_info$reads*nor_info$no_OTU
# info = gather(info,"d_type","No",2:3)
# nor_info = gather(nor_info, "d_type", "No",2:3)



###### original data
z1= ggplot(info,aes(x = type, fill = type, y = reads))+
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.3)+
  labs(x = "", y = "Frequency", title = "reads")+
  scale_fill_grey(name = "", labels = c("1 individual/sample",'~50 individual/sample') , start = 1, end = 0.7)+
  theme_bw()+
  scale_x_discrete(labels = c('BC'="",'MB'=''))+
  theme(legend.position="bottom")

a1= ggplot(info,aes(x = type, fill = type, y = uniqs))+
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.3)+
  labs(x = "", y = "Frequency", title = "uniqs")+
  scale_fill_grey(name = "", labels = c("1 individual/sample",'~50 individual/sample') , start = 1, end = 0.7)+
  theme_bw()+
  scale_x_discrete(labels = c('BC'="",'MB'=''))+
  theme(legend.position="bottom")


b1= ggplot(info,aes(x = OTU_type, fill = type, y = no_OTU))+
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.3)+
  labs(x = "", y = "OTU", title = "")+
  scale_fill_grey(name = "", labels = c("1 individual/sample",'~50 individual/sample') ,start = 1, end = 0.7)+
  theme_bw()+
  scale_x_discrete(labels = c("aft_OTU" = "after filtering", "bef_OTU" = "before filtering","master_OTU" ="overall"))+
  theme(legend.position="bottom")

c1= ggplot(info,aes(x = uniq_status, fill = type, y = no_uniqs))+
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.15)+
  labs(x = "", y = "uinqs", title = "")+
  scale_fill_grey(name = "", labels = c("1 individual/sample",'~50 individual/sample'),start = 1, end = 0.7)+
  theme_bw()+
  scale_x_discrete(labels = c("mapped in a OTU", "mapped in both b&a OTU","mapped in none" ))+
  theme(legend.position="bottom")
  # theme(axis.text.x = element_text(angle=345))


plot_grid(z1,a1,b1,c1,labels = c("a","b","c"))

#### norm data


a2= ggplot(nor_info,aes(x = d_type, fill = type, y = No))+
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.1)+
  labs(x = "", y = "log Frequency", title = "")+
  scale_fill_grey(start = 1, end = 0.7)+
  theme_bw()


b2= ggplot(nor_info,aes(x = OTU_type, fill = type, y = no_OTU))+
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.3)+
  labs(x = "", y = "log OTU", title = "")+
  scale_fill_grey(start = 1, end = 0.7)+
  theme_bw()

c2= ggplot(nor_info,aes(x = uniq_status, fill = type, y = no_uniqs))+
  scale_y_log10()+
  geom_violin()+
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.15)+
  labs(x = "", y = "log uinqs", title = "")+
  scale_fill_grey(start = 1, end = 0.7)+
  theme_bw()

plot_grid(a2,b2,c2,labels = c("a","b","c"))


#Hiseq vs Miseq plot #checking as requested
HM_info = read.csv("../data/Hi_Mi_info.csv",header = FALSE)
colnames(HM_info) = c("X","platform") 
info_with_HM = merge(info,HM_info,by = "X")
info_with_HM = unique(info_with_HM[,c("X","reads","uniqs","platform")]) 
info_with_HM$type = NA
info_with_HM[startsWith(info_with_HM$X,"RNBC"),]$type = rep("BC",length(info_with_HM[startsWith(info_with_HM$X,"RNBC"),]))
info_with_HM[startsWith(info_with_HM$X,"RNMB"),]$type = rep("MB",length(info_with_HM[startsWith(info_with_HM$X,"RNMB"),]))

HM = ggplot(data = info_with_HM, aes(x = reads, y = uniqs, col = platform))+
  geom_point(aes(shape = platform), size = 1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  labs(x = "Number of reads(quality filtered)", y= "Number of uniqs")+
  geom_smooth(method = 'lm')

#facet by BC MB
HM + facet_grid(cols = vars(type))
  # p0 =ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = reads))+
#   scale_y_log10()
# 
# p1=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs))+
#   scale_y_log10()
# 
# p2=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs_in_bef_OTU))+
#   scale_y_log10()
# 
# p3=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = bef_OTU))+
#   scale_y_log10()
# 
# p4=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs_in_aft_OTU))+
#   scale_y_log10()
# 
# p5=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = aft_OTU))+
#   scale_y_log10()
# 
# p8=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs_got_nor))+
#   scale_y_log10()
# 
# p9=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = master_OTU))+
#   scale_y_log10()
# 
# grid.arrange(p0, p1, p2, p3, p4, p5, p8, p9,nrow = 3)
# 
# ## not logged
# p0 =ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = reads))
# 
# p1=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs))
# 
# p2=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs_in_bef_OTU))
# 
# p3=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = bef_OTU))
# 
# p4=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs_in_aft_OTU))
# 
# p5=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = aft_OTU))
# 
# p8=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = uniqs_got_nor))
# 
# p9=ggplot(data = info,aes(x=index,color = type))+
#   geom_point(aes(y = master_OTU))
# 
# grid.arrange(p0, p1, p2, p3, p4, p5, p8, p9,nrow = 3)
# 
##violin
# p0 =ggplot(data = info,aes(x= type,y = reads,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# p1=ggplot(data = info,aes(x= type,y = uniqs,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# p2=ggplot(data = info,aes(x= type,y = uniqs_in_bef_OTU,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# p3=ggplot(data = info,aes(x= type,y = bef_OTU,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# p4=ggplot(data = info,aes(x= type,y = uniqs_in_aft_OTU,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# p5=ggplot(data = info,aes(x= type,y = aft_OTU,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# p8=ggplot(data = info,aes(x= type,y = uniqs_got_nor,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# p9=ggplot(data = info,aes(x= type,y = master_OTU,fill = type))+
#   geom_violin()+
#   geom_boxplot(width = 0.2,na.rm = TRUE)+
#   scale_y_log10()+
#   theme(legend.position="none")
# 
# plot_grid(p0, p1, p2, p3, p4, p5, p8, p9,nrow = 3)
