################ 
### Barplots ###
################

# Select data for analysis
insect_b <- insect2
insect_b <- subset_samples(insect2, State2!="Guanajuato")
insect_b <- subset_samples(insect_b, Genus=="Peponapis")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")
#insect_b <- subset_samples(insect2, State=="California")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)


p <- plot_bar(insect_b, x="State_detail", fill="OTU")
p

############################# 
### Proportional Barplots ###
#############################

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 20 OTUs in Peponapis"
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2

TopNOTUs = names(sort(taxa_sums(insectR), TRUE)[1:100])
insectR100 = prune_taxa(TopNOTUs, insectR)

title = "Top 100 OTUs in the dataset"
plot_bar(insectR100, "sample_Species", fill = "Phylum", title = title) + coord_flip()

## Merge samples by category of insect species
insectRm = merge_samples(insectR100, "sample_Species")

# Repair the merged values associated with each surface after merge.
sample_data(insectRm)$sample_Species <- levels(sample_data(insectR)$sample_Species)

# Transform to percentages of total available.
insectRm = transform_sample_counts(insectRm, function(x) 100 * x/sum(x))

title = "Class proportion in each insect species"
plot_bar(insectRm, fill = "Class", title = title) + coord_flip() 
