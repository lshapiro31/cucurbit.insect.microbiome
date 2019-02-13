insect_b <- insectR ## insectR is rarefied to 1000
insect_b <- subset_samples(insect_b, Genus=="Acalymma") ## Subset to Acalymma beetles
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia") ## Subset to Erwinia genus
insect_b <- prune_taxa(taxa_sums(insect_b) > 0, insect_b) ## Restrict to taxa that are present at least once

ntaxa(insect_b)
sample_sums(insect_b)  # How many Erwinia reads are in each sample?
