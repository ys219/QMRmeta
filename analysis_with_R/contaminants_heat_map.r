rm(list = ls())
library(tidyr)
library(tidyverse)
library(ggplot2)



map = read.csv("../data/[INFO]conta_score.csv",row.names = 1)
collapsed_map = map %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
# collapsed_map[is.na(collapsed_map)] = 0

ggplot(collapsed_map, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()+
  scale_fill_distiller(palette = "GnBu")+
  scale_x_discrete(breaks="")+
  scale_y_discrete(breaks="")+
  labs(x = "MB samples", y = "BC samples",title = "original")

omit_na_map = map[apply(map,1, function(x)any(!is.na(x))),]
omit_na_map = omit_na_map[,apply(omit_na_map,2, function(x)any(!is.na(x)))]
omit_coll_map = omit_na_map %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

omit_coll_map[is.na(omit_coll_map)] = 0

ggplot(omit_coll_map, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()+
  scale_fill_distiller(palette = "GnBu")+
  scale_x_discrete(breaks="")+
  scale_y_discrete(breaks="")+
  labs(x = "MB samples", y = "BC samples", title = "NA omitted")
