library(biomformat)
library(vegan)
library(plotly)
library(phyloseq)
library(ggplot2)
library(genefilter)
#library(impute)
library(PMA)
library(ade4)
library(ggrepel)
library(metagenomeSeq)
library(ape)
library(dplyr)
library(scales)
library(grid)
library(reshape2)

## To remove attached packages 
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

theme_set(theme_bw())

source("miseqR.R", local = TRUE) # From here http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

## https://github.com/edamame-course/2015-tutorials/wiki
## Picrust tutorial https://github.com/LangilleLab/microbiome_helper/wiki/STAMPS-2016-PICRUSt-Tutorial
## Phyloseq Tutorial http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
## AWESOME EXAMPLES http://joey711.github.io/phyloseq-demo/Restroom-Biogeography
## Metagenomics explanations: http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity
## https://rpubs.com/collnell/manova
setwd("/Users/lorishapiro/Dropbox/Matrix_forsharing_16S/Bee-beetle-only-beta")

biom_file = read_biom("otu_filter_meta__Type2_insect__.biom")
map_file = ("New_mapping_file2__Type2_insect__.txt")

trefile = ("rep_set.tre")
insect.tree <- read.tree(trefile)

## This is unrarified data!!

biomot = import_biom(biom_file, parseFunction = parse_taxonomy_default)
bmsd = import_qiime_sample_data(map_file)
class(bmsd)
dim(bmsd)

#tax_table <- ("seqs_rep_set_tax_assignments.txt")
insect = merge_phyloseq(biomot, bmsd, insect.tree)

colnames(tax_table(insect)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(insect)

# Are there any OTUs included in this dataset that have no counted reads at all prior to preprocessing? If so, how many?
any(taxa_sums(insect) == 0)
sum(taxa_sums(insect) == 0)

insect2 <- insect

### Many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate
#insect3 <- prune_taxa(taxa_sums(insect2) > 0, insect2)
#insect3 <- prune_taxa(taxa_sums(insect2) > 1, insect2)
insect3 <- prune_taxa(taxa_sums(insect2) > 2, insect2)
insect2 <- insect3
ntaxa(insect2)

levels(sample_data(insect2)$sample_Species)
levels(sample_data(insect2)$State_detail)

## Remove samples with less than 500 reads (there shouldn't be any)
# insect2 = prune_samples(sample_sums(insect) > 500, insect)

## Number of sequencing reads per sample, pre-filtering
sample_sums(insect)

##################
## Otus to omit ##
##################

# https://github.com/joey711/phyloseq/issues/652
badTaxa = read.table("Wolbachia_Chloroplast.otus", header=F) ## Wolbachia; should go back and look at distribution of htis later
badTaxa<-as.character(badTaxa[[1]])
allTaxa = taxa_names(insect2)
length(allTaxa)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
length(allTaxa)
insect_filt = prune_taxa(allTaxa, insect2)
# new phyloseq object with just the taxa you kept.
str(insect_filt)
insect3 <- insect_filt  ## Reassign filtered taxa to insect3 phyloseq object
insect2 <- insect3

ntaxa(insect2)
ntaxa(insect3)

## Rarefy data
set.seed(28132)
sample_sums(insect2)
insectR = rarefy_even_depth(insect2, sample.size = 1000)
sample_sums(insectR)

otu_table(insectR)[1:5, 1:10]
tax_table(insectR)[1:5, 1:4]
myTaxa = names(sort(taxa_sums(insectR), decreasing = TRUE)[1:10])
myTaxa
insectR
nsamples(insectR)
ntaxa(insectR)
sample_names(insectR)[1:5]
sample_names(insectR)
sample_variables(insectR)
rank_names(insectR)

## Select data for analyses
insect_b <- insectR
insect_b <- insect2
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Genus=="Diabrotica")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")
insect_b <- subset_samples(insect2, State2=="Iowa")
ntaxa(insect_b)
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")
insect_b <- subset_samples(insect2, State2!="Guanajuato")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
insect_R <- subset_samples(insectR, Specialist=="Yes")

## 
# melt to long format (for ggploting) 
# prune out phyla below 2% in each sample

insect_phylum <- insect_b %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

# Plot 
ggplot(insect_phylum, aes(x = sample_Species2, y = Abundance, fill = Phylum)) + 
#  facet_grid(State2~.) + ## To ask: How to facet so that there are no empty bars?
  geom_bar(stat = "identity", position="fill") +  # position="fill" brings bars to 100%
  coord_flip() +
  scale_fill_manual(values = phylum_colors) +
#  scale_x_discrete(
#    breaks = c("7/8", "8/4", "9/2", "10/6"),
#    labels = c("Jul", "Aug", "Sep", "Oct"), 
#    drop = FALSE
#  ) +
  # Remove x axis title
#  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Bug Guts \n by Sampling Location") 

# One of the best exploratory analyses for amplicon data is 
# unconstrained ordinations. Here we will look at ordinations of our 
# full community samples. We will use the scale_reads() function in miseqR.R 
# to scale to the smallest library size, which is the default. 
# If you want to scale to another depth, you can do so by setting the “n” argument.

# Scale reads to even depth 
#insectb_scale <- insect_b %>%
#  scale_reads(round = "round") 

## Use rarefied data 
insect_b <- insectR
insect_b <- subset_samples(insect_b, State2!="Guanajuato")
insect_b <- subset_samples(insect_b, Specialist=="Yes")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
ntaxa(insect_b)

# Ordinate
insectR_ord <- ordinate(
  physeq = insect_b, 
  method = "PCoA", 
  distance = "bray"
)

# Plot 
plot_ordination(
  physeq = insect_b,
  ordination = insectR_ord,
  color = "State2",
  shape = "insect_Type",
  title = "PCoA of Bacterial Communities within a Season"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19","black", "blue",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta",
                                "pink", "orange")
  ) +
  geom_point(aes(color = State2), alpha = 0.7, size = 5) +
  geom_point(colour = "grey90", size = 0.15) 

# Let’s try an NMDS instead. 
# For NMDS plots it’s important to set a seed since the starting positions 
# of samples in the alogrithm is random.
# Important: if you calculate your bray-curtis distance metric 
# “in-line” it will perform a square root transformation and 
# Wisconsin double standardization. If you don’t want this, 
# you can calculate your bray-curtis distance separately

## Use rarefied data 
insect_b <- insectR
insect_b <- subset_samples(insect_b, State2!="Guanajuato")
insect_b <- subset_samples(insect_b, Specialist=="Yes")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
ntaxa(insect_b)

set.seed(1)

# Ordinate
insect_nmds <- ordinate(
  physeq = insect_b, 
  method = "NMDS", 
  distance = "bray"
)

# Plot 
plot_ordination(
  physeq = insect_b,
  ordination = insect_nmds,
  color = "State2",
  shape = "insect_Type",
  title = "PCoA of Bacterial Communities within a Season"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19","black", "blue",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta",
                                "pink", "orange")
  ) +
  geom_point(aes(color = State2), alpha = 0.7, size = 5) +
  geom_point(colour = "grey90", size = 0.15) 

# Here is an example of how to run a permanova test using the 
# adonis function in vegan. 
# In this example we are testing the hypothesis that the 
# three stations we collected samples from have different centroids

set.seed(1)

# Calculate bray curtis distance matrix
insect_bray <- phyloseq::distance(insect_b, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(insect_b))

# Adonis test
adonis(insect_bray ~ insect_Type, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(insect_bray, sampledf$insect_Type)
permutest(beta)

# This output tells us that our adonis test is 
# significant so we can reject the null hypothesis that 
# our three sites have the same centroid.

# Additionally, our betadisper results are not significant, 
# meaning we cannot reject the null hypothesis that our groups 
# have the same dispersions. This means we can be more confident 
# that our adonis result is a real result, and not due to differences 
# in group dispersions
# There is a lot more analysis that can be done here. 
# We could use a distance metric other than Bray-curtis, 
# we could test different grouping variables, or we could 
# create a more complex permanova by testing a model that 
# combines multiple variables. 
# Unfortunately, there are currently no post-hoc tests 
# developed for adonis.

insect_b <- insectR
insect_b <- subset_samples(insect_b, Time_Course=="Yes")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)

set.seed(1)

# Calculate bray curtis distance matrix
insect_bray <- phyloseq::distance(insect_b, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(insect_b))

# Adonis test
time.course.results <- adonis(insect_bray ~ State2, data = sampledf)
time.course.results
summary(time.course.results)

# Homogeneity of dispersion test
beta <- betadisper(insect_bray, sampledf$State2)
permutest(beta)

## http://microbiome.github.io/microbiome/Core.html
# https://bioconductor.org/packages/release/bioc/html/microbiome.html

library(microbiome)

insect_b <- insectR
insect_b <- subset_samples(insect_b, Time_Course=="Yes")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)

# Relative population frequencies; at 1% compositional abundance threshold:
insect_b.rel <- microbiome::transform(insect_b, "compositional")
head(prevalence(insect_b.rel, detection = 1/100, sort = TRUE))

# Absolute population frequencies (sample count):
head(prevalence(insect_b.rel, detection = 1/100, sort = TRUE, count = TRUE))

# Core microbiota analysis
core.taxa.standard <- core_members(insect_b.rel, detection = 0, prevalence = 50/100)
head(core.taxa.standard)

# A full phyloseq object of the core microbiota is obtained as follows:
insect_b.core <- core(insect_b.rel, detection = 0, prevalence = .5)

# Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(insect_b.core)

# Total core abundance in each sample (sum of abundances of the core members):
core.abundance <- sample_sums(core(pseq.rel, detection = .01, prevalence = .95))






