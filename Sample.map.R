library(ggplot2)
library(dplyr)
library(stringr)
library(devtools)
library(maps)
library(mapdata)
library(ggthemes)
library(mapproj)
require(dplyr); require(RColorBrewer); require(ggplot2)
require(mapdata); require(maptools)

theme_set(theme_bw())

setwd("/Users/lorishapiro/Dropbox/Matrix_forsharing_16S/Bee-beetle-only-beta")

samples <- read.table("samples.txt", sep="\t", header=TRUE, fileEncoding="UTF-16")

usa <- map_data("usa")
canada <- map_data("worldHires", "Canada")
mexico <- map_data("worldHires", "Mexico")

NAmap <- ggplot() + geom_polygon(data = usa, 
                                 aes(x=long, y = lat, group = group), 
                                 fill = "grey88", 
                                 color="black") +
  ggtitle("Insect Collection Locations") +
  geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
               fill = "grey88", color="black") + 
  geom_polygon(data = mexico, aes(x=long, y = lat, group = group), 
               fill = "grey88", color="black") +
  coord_fixed(xlim = c(-125, -65),  ylim = c(18, 50), ratio = 1.3) +
  geom_point(data=samples, aes(x=Longitude, y=Latitude),
             fill=samples$Color, color = "black", shape=21, size=5.0) #+
  theme(line = element_blank(),
        text = element_blank(), 
        rect = element_blank()) 
NAmap 
NAmap + labs(colour = "Type") ## Can't figure out how to get the legend to print

dev.off()

## To find RGB colors for Illustrator
#http://research.stowers.org/mcm/efg/R/Color/Chart/
  
colors()[grep("sky",colors())]

#The function col2rgb can be used to extract the RGB (red-green-blue) components of a color, e.g.,
col2rgb("orange")
col2rgb("darkolivegreen")
col2rgb("darkolivegreen2")
