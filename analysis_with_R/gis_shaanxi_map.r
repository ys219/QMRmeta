# rm(list= ls())
# library(maps)
# library(mapdata)
# library(maptools)
# library(sp)
# library(sf)
library(rgdal)
# library(ggplot2)

## 
# china_map = st_read("../data/china_map/CHN_adm0.shp")
plot(china_map[1],col = 'white')
plot(shaanxi_map[1],col = 'white')
test_map = st_read("../data/china_map/CHN_adm3.shp")
shaanxi_map = test_map[which(test_map$NAME_1 == "Shaanxi"),9]
local_map = shaanxi_map[which(shaanxi_map$NAME_3 == "Feng" |shaanxi_map$NAME_3 == "Mei" | shaanxi_map$NAME_3 == "Foping"|shaanxi_map$NAME_3 == "Zhouzhi"| shaanxi_map$NAME_3 == "Foping"|shaanxi_map$NAME_3 == "Taibai"),]


test_coord = coord[41:44,2:3] ; test_coord$lat = as.numeric(test_coord$lat); test_coord$long = as.numeric(test_coord$long)
# st_crs(test_coord)
test_coord = st_as_sf(test_coord, coords = c('long','lat'),crs = 4326)#

plot(local_map, asp = 1, col = 'white',add = T)
par(new = TRUE)
plot(st_geometry(test_coord), pch = 20, col = 'red',cex = 1)



## rgdal 

china_map = readOGR(dsn = "../data/china_map",layer = "CHN_adm3")
head(china_map@data, n = 2)
