################# 
### Peponapis ###
#################

## Subset top 20 otus in Peponapis
# Start with Peponapis
# Proportional barchart
insect.peponapis <- subset_samples(insect2, Genus=="Peponapis")
TopNOTUs = names(sort(taxa_sums(insect.peponapis), TRUE)[1:20])

insect_b <- insect2 %>%
  subset_samples(
      Genus=="Peponapis" & State2!="Guanajuato"
                 ) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(TopNOTUs, .) %>%
  tax_glom(., taxrank="Family") %>% 
  transform_sample_counts(., function(x) 100 * x/sum(x))
  

p <- plot_bar(insect_b, x="State_detail", fill="Family") + 
  ggtitle("Family breakdown of top 50 OTUS in in Peponapis (unrarefied)") +
  ylab("Proportion of Reads per Taxonomic Group") +
  xlab("Individual Bee") +
  coord_flip()
p + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
      axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
      axis.title.y = element_text(face="bold", color = "black", size=12),
      axis.title.x = element_text(face="bold", color = "black", size=12),
      axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
      legend.title = element_text(size=12, face="bold", hjust = 0.5),
      legend.text = element_text(size = 10, face = "bold"))

ggsave("Peponapis.proportional.barchart.pdf", height = 6, width = 8)

dev.off()

################################ 
## Acalymma - non Time Course ##
################################

## Subset top 20 otus in Acalymma

insect.acalymma <- subset_samples(insect2, Genus=="Acalymma")
TopNOTUs = names(sort(taxa_sums(insect.acalymma), TRUE)[1:20])

insect_b <- insect2 %>% 
  subset_samples(
    Genus=="Acalymma" &
    Time_Course=="No"
    ) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(TopNOTUs, .) %>%
  tax_glom(., taxrank="Family") %>% 
  transform_sample_counts(., function(x) 100 * x/sum(x))

acalymma.all <- plot_bar(insect_b, x="State_detail", fill="Family") + 
  ggtitle("Family breakdown of top 10 OTUS\nin Acalymma spp. (unrarefied)") +
  ylab("Proportion of Reads per Taxonomic Group") +
  xlab("Individual Striped Beetle") +
  coord_flip()
acalymma.all + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
          axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
          axis.title.y = element_text(face="bold", color = "black", size=12),
          axis.title.x = element_text(face="bold", color = "black", size=12),
          axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
          legend.title = element_text(size=12, face="bold", hjust = 0.5),
          legend.text = element_text(size = 10, face = "bold"))

ggsave("Acalymma.proportional.barchart.pdf", height = 8, width = 12)
dev.off()


############################ 
## Acalymma - Time Course ##
############################


## Subset top 20 otus in Acalymma

insect.acalymma <- subset_samples(insect2, Genus=="Acalymma")
TopNOTUs = names(sort(taxa_sums(insect.acalymma), TRUE)[1:20])

insect_b <- insect2 %>% 
  subset_samples(
    Genus=="Acalymma" &
      Time_Course=="Yes"
  ) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(TopNOTUs, .) %>%
  tax_glom(., taxrank="Family") %>% 
  transform_sample_counts(., function(x) 100 * x/sum(x))

acalymma.timeCourse <- plot_bar(insect_b, x="State_detail", fill="Family") + 
  ggtitle("Family breakdown of top 10 OTUS\nin Acalymma spp. in Massachussetts (unrarefied)") +
  ylab("Proportion of Reads per Taxonomic Group") +
  xlab("Individual Striped Beetle") +
  coord_flip()
acalymma.timeCourse + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
                     axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
                     axis.title.y = element_text(face="bold", color = "black", size=12),
                     axis.title.x = element_text(face="bold", color = "black", size=12),
                     axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
                     legend.title = element_text(size=12, face="bold", hjust = 0.5),
                     legend.text = element_text(size = 10, face = "bold"))

ggsave("Acalymma.timeCourse.proportional.barchart.pdf", height = 8, width = 12)

dev.off()







###########################
### Merge samples by category of insect species
#insectRm.Species = merge_samples(insect_b, "sample_Species")

### Merge samples by category of state 
insectRm.State2 = merge_samples(insect.50, "State2")

# Repair the merged values associated with each surface after merge.
sample_data(insectRm.State2)$sample_Species <- levels(sample_data(insect_b)$sample_Species)
insectRm.Species.Class <- tax_glom(insectRm.State2, taxrank="Class")

# Transform to percentages of total available.
insectRm = transform_sample_counts(insectRm.Species.Class, function(x) 100 * x/sum(x))

# 6x3 inches plot saved to pdf
title = "Class proportion in Acalymma by state"
Striped.State <-plot_bar(insectRm, fill = "Class", title = title) + coord_flip() 
Striped.State + 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5))

############################ 
## Acalymma - Time Course ##
############################

insect_b <- insect2 %>%
  subset_samples(
      Genus=="Acalymma" &
      Time_Course=="Yes"
    )  %>%
  prune_samples(sample_sums(.) > 0, .)

## Subset top 100 otus
TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:50])
insect.50 = prune_taxa(TopNOTUs, insect_b)

insect_prop = transform_sample_counts(insect.50, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 100 OTUs in Acalymma spp."
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2

### Merge samples by category of insect species
#insectRm.Species = merge_samples(insect_b, "sample_Species")

### Merge samples by category of state 
insectRm.State2 = merge_samples(insect.50, "State2")

# Repair the merged values associated with each surface after merge.
sample_data(insectRm.State2)$sample_Species <- levels(sample_data(insect_b)$sample_Species)
insectRm.Species.Class <- tax_glom(insectRm.State2, taxrank="Class")

# Transform to percentages of total available.
insectRm = transform_sample_counts(insectRm.Species.Class, function(x) 100 * x/sum(x))

# 6x3 inches plot saved to pdf
title = "Class proportion in Massachusetts Acalymma spp. over the season"
Striped.State <-plot_bar(insectRm, fill = "Class", title = title) + coord_flip() 
Striped.State + 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5))


################ 
## Diabrotica ##
################

# This should be faceted by species 
insect_b <- insect2 %>%
  subset_samples(
      Genus=="Diabrotica" &
      Time_Course=="No"
    ) %>%
  prune_samples(sample_sums(.) > 0, .)

## Subset top 100 otus
TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:50])
insect.50 = prune_taxa(TopNOTUs, insect_b)

insect_prop = transform_sample_counts(insect.50, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 100 OTUs in Acalymma spp."
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2

### Merge samples by category of insect species
#insectRm.Species = merge_samples(insect_b, "sample_Species")

### Merge samples by category of state 
insectRm.State2 = merge_samples(insect.50, "State2")

# Repair the merged values associated with each surface after merge.
sample_data(insectRm.State2)$sample_Species <- levels(sample_data(insect_b)$sample_Species)
insectRm.Species.Class <- tax_glom(insectRm.State2, taxrank="Class")

# Transform to percentages of total available.
insectRm = transform_sample_counts(insectRm.Species.Class, function(x) 100 * x/sum(x))

# 6x3 inches plot saved to pdf
title = "Class proportion in Diabrotica by state"
Striped.State <-plot_bar(insectRm, fill = "Class", title = title) + coord_flip() 
Striped.State + 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5))



