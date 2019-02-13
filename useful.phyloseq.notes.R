## Getting to know your phyloseq object 
## Using filtered and rarefied phyloseq object (insectR)

# Are there any OTUs included in this dataset that have no counted reads at all prior to preprocessing? If so, how many?
any(taxa_sums(insectR) == 0)
sum(taxa_sums(insectR) == 0)

otu_table(insectR)[1:5, 1:10]
tax_table(insectR)[1:5, 1:4]
myTaxa = names(sort(taxa_sums(insectR), decreasing = TRUE)[1:10])
myTaxa
str(insectR)
nsamples(insectR)
ntaxa(insectR)
sample_names(insectR)[1:5]
sample_names(insectR)
sample_variables(insectR)
rank_names(insectR)

# Number of distinct phyla
tax_table(insectR)[,"Phylum"] %>% unique %>% na.exclude %>% length
tax_table(insectR)[,"Family"] %>% unique %>% na.exclude %>% length

## Number of samples per variable
table(sample_data(insectR)$State)
table(sample_data(insectR)$Species)
table(sample_data(insectR)$Genus)

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
insect_b = prune_taxa(TopNOTUs, insect_b)