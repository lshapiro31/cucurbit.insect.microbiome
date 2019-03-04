################ 
### Barplots ###
################

# Select data for analysis
# Start with Peponapis
insect_b <- insect2
insect_b <- subset_samples(insect_b, State2!="Guanajuato")
insect_b <- subset_samples(insect_b, Genus=="Peponapis")
insect_b <- prune_samples(sample_sums(insect_b) > 1, insect_b) ## Restrict to samples that have at least 100 reads

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:50])
insect_b = prune_taxa(TopNOTUs, insect_b)
insect.Peponapis.Family <- tax_glom(insect_b, taxrank="Family")

insect_prop = transform_sample_counts(insect.Peponapis.Family, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Family breakdown of top 50 OTUS in in Peponapis (unrarefied)"
p <- plot_bar(insect_prop, x="State_detail", fill="Family", title=title) + coord_flip()
p + theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5))

############################# 
### Proportional Barplots ###
#############################

################################ 
## Acalymma - non Time Course ##
################################

insect_b <- insect2
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_samples(insect_b, Time_Course=="No")
insect_b <- prune_samples(sample_sums(insect_b) > 1, insect_b) ## Restrict to samples that have at least 100 reads

## Subset top 50 otus
TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:50])
insect.50 = prune_taxa(TopNOTUs, insect_b)
insect.Acalymma.Family <- tax_glom(insect.50, taxrank="Family")

insect_prop = transform_sample_counts(insect.Acalymma.Family, function(x) 100 * x/sum(x))

title="Normalized counts of top 50 OTUs in Acalymma spp. (unrarefied)"
p2 <- plot_bar(insect_prop, x="State_detail", fill="Family", title=title) + coord_flip()
p2 + theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5))

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

################################ 
## Acalymma - Time Course ##
################################

insect_b <- insect2
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")
insect_b <- prune_samples(sample_sums(insect_b) > 1, insect_b) ## Restrict to samples that have at least 100 reads

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
insect_b <- insect2
insect_b <- subset_samples(insect_b, Genus=="Diabrotica")
insect_b <- subset_samples(insect_b, Time_Course=="No")

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



