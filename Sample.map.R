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

setwd("~/Dropbox/Matrix_forsharing_16S/Bee-beetle-only-beta")

samples <- read.table("samples.txt", sep="\t", header=TRUE, fileEncoding="UTF-16")
str(samples)
samples <- as_tibble(samples) 

#samples <- filter(samples, Type != "Peponapis pruinosa")## Subset data

usa <- map_data("usa")
canada <- map_data("worldHires", "Canada")
mexico <- map_data("worldHires", "Mexico")

NAmap <- ggplot() + geom_polygon(data = usa, 
                                 aes(x=long, y = lat, group = group), 
                                 fill = "grey88", 
                                 color="black") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Insect Collection Locations") +
  geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
               fill = "grey88", color="black") + 
  geom_polygon(data = mexico, aes(x=long, y = lat, group = group), 
               fill = "grey88", color="black") +
  coord_fixed(xlim = c(-125, -65),  ylim = c(18, 50), ratio = 1.3) +
  geom_point(data=samples, aes(x=Longitude, y=Latitude),
             fill=samples$Color, color = "black", shape=21, size=samples$sampleSize) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
            axis.text.x  = element_text(face="bold", color = "black", hjust=0.5, size=10),
            axis.title.y = element_text(face="bold", color = "black", size=12),
            axis.title.x = element_text(face="bold", color = "black", size=12),
            axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
            legend.title = element_text(size=12, face="bold", hjust = 0.5),
            legend.text = element_text(size = 10, face = "bold") +
theme(legend.position="top"))

NAmap  ## Can't figure out how to get the legend to print

ggsave("sample.map.pdf", height=6, width=6)

dev.off()

## To find RGB colors for Illustrator
#http://research.stowers.org/mcm/efg/R/Color/Chart/
  
colors()[grep("sky",colors())]

#The function col2rgb can be used to extract the RGB (red-green-blue) components of a color, e.g.,
col2rgb("orange")
col2rgb("darkolivegreen")
col2rgb("darkolivegreen2")

-------------------------------------------------------------------
  ## Map for grant
  NAmap <- ggplot() + geom_polygon(data = usa, 
                                   aes(x=long, y = lat, group = group), 
                                   fill = "grey88", 
                                   color="black") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Erwinia tracheiphila epidemic region") +
  geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
               fill = "grey88", color="black") + 
  geom_polygon(data = mexico, aes(x=long, y = lat, group = group), 
               fill = "grey88", color="black") +
  coord_fixed(xlim = c(-125, -65),  ylim = c(18, 50), ratio = 1.3) #+
#  geom_point(data=samples, aes(x=Longitude, y=Latitude),
#             fill=samples$Color, color = "black", shape=21, size=samples$sampleSize) +
NAmap + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=12),
        axis.text.x  = element_text(face="bold", color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=10),
        axis.title.x = element_text(face="bold", color = "black", size=10),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title = element_text(size=12, face="bold", hjust = 0.5),
        legend.text = element_text(size = 10, face = "bold") +
          theme(legend.position="top"))

NAmap  ## Can't figure out how to get the legend to print

ggsave("sample.map.pdf", height=6, width=6)

dev.off()
