##################
## Erwinia OTUS ##
##################

## The otus that are probably Et are
denovo577182
denovo590454
denovo634643

insect2 # phyloseq object with no filtering
insectR #phyloseq object where wolbachia + chloroplasts have been removed, normalized to 1000 reads per sample

insect_b <- insectR
insect_b <- insect2 ## insectR is rarefied to 1000

#insect_b <- subset_samples(insect_b, Genus=="Acalymma") ## Subset to Acalymma beetles
insect_b <- subset_samples(insect_b, Time_Course=="Yes") ## Subset to beetles from the time course in Cambridge
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia") ## Subset to Erwinia genus
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b) ## Restrict to taxa that are present at least once

## Make a list of Et OTUs
## Subset the phyloseq object to only Et OTUs
Et <- c("denovo577182", "denovo590454", "denovo634643")
insect.Et = prune_taxa(Et, insect_b)
insect.Et.merged = merge_taxa(insect.Et, taxa_names(insect.Et)[1:3])

sample_sums(insect.Et.merged)  # How many Erwinia reads are in each sample?
sample_sums(insect.Et)

Top.10.erwinia.OTUs = names(sort(taxa_sums(insect_b), TRUE)[1:10])
insect_b = prune_taxa(Top.10.erwinia.OTUs, insect_b)

## This is less than with unrarified data
title="Erwinia spp. in Acalymma"
p2 <- plot_bar(insect.Et, x="State_detail", fill="OTU", title=title)
p2

#####################
## Unrarefied data ##
#####################

insect_b <-insect2 # Unrarefied data
insect_b <- subset_samples(insect_b, Time_Course=="Yes") ## Subset to beetles from the time course in Cambridge
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia") ## Subset to Erwinia genus
insect_b <- prune_taxa(taxa_sums(insect_b) > 0, insect_b) ## Restrict to taxa that are present at least once

ntaxa(insect_b)
sample_sums(insect_b)  # How many Erwinia reads are in each sample?


## This is less Erwinia than with unrarified data

title="Erwinia spp. in Acalymma"
p2 <- plot_bar(insect_b, x="State_detail", fill="OTU", title=title)
p2
## Very closely related to a plant associated Erwinia! Present in Western species as well

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 20 OTUs in Peponapis"
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2
