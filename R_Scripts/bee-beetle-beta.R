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

theme_set(theme_bw())

## Picrust tutorial https://github.com/LangilleLab/microbiome_helper/wiki/STAMPS-2016-PICRUSt-Tutorial
## Phyloseq Tutorial http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
## AWESOME EXAMPLES http://joey711.github.io/phyloseq-demo/Restroom-Biogeography
## Metagenomics explanations: http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity

setwd("/Users/lorishapiro/Dropbox/Matrix_forsharing_16S/Bee-beetle-only-beta")

#biom_file = read_biom("Beetle_Bee.biom")
#map_file = ("Beetle_Bee.txt")

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
insect3 <- prune_taxa(taxa_sums(insect2) > 0, insect2)
insect3 <- prune_taxa(taxa_sums(insect2) > 1, insect2)
insect3 <- prune_taxa(taxa_sums(insect2) > 2, insect2)
insect2 <- insect3
ntaxa(insect2)

levels(sample_data(insect2)$sample_Species)
levels(sample_data(insect2)$State_detail)

## Display total reads per sample, 
## before chloroplast removal

readsumsdf = data.frame(nreads = sort(taxa_sums(insect2), TRUE), sorted = 1:ntaxa(insect2), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(insect2), TRUE), sorted = 1:nsamples(insect2), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
  
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

## Display total reads per sample, 
## after chloroplast removal

readsumsdf = data.frame(nreads = sort(taxa_sums(insect2), TRUE), sorted = 1:ntaxa(insect2), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(insect2), TRUE), sorted = 1:nsamples(insect2), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

## Number of sequencing reads per sample, post-filtering
## Prune samples that have less than 50 bacterial reads
sample_sums(insect2)
insect.prune <- prune_samples(sample_sums(insect2) > 50, insect2)
sample_sums(insect.prune)
insect2 <-insect.prune

## Transform to equal sample sum
## an alternative to random subsampling, a simple proportional transformation
insectP = transform_sample_counts(insect.prune, function(x) 100 * x/sum(x))
#insect2<-insectP
sample_sums(insectP)


# For sanity-check, replot the sample sums of each of these new data objects, 
## to convince ourselves that all of the samples now sum to 1500.

par(mfrow = c(1, 2))
title = "Sum of reads for each sample, insectR"
plot(sort(sample_sums(insectP), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 2000))
title = "Sum of reads for each sample, insectP"
plot(sort(sample_sums(insectR), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 2000))

## Visualize OTU variance
#hist(log10(apply(otu_table(insect2), 1, var)), xlab = "log10(variance)", main = "A large fraction of OTUs have very low variance")

## Filter out otus with high variance
#varianceThreshold = 5
#keepOTUs = apply(otu_table(insect2), 1, var) > varianceThreshold
#insect3 = prune_taxa(keepOTUs, insect3)
#insect3
#ntaxa(insect3)
#insect2 <- insect3

## Remove features with ambiguous phylum annotations
# insect2 = subset_taxa(insect2, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
# ntaxa(insect2)

### Initial exploration of data

otu_table(insect2)[1:5, 1:10]
tax_table(insect2)[1:5, 1:4]
myTaxa = names(sort(taxa_sums(insect2), decreasing = TRUE)[1:10])
myTaxa
insect2
nsamples(insect2)
ntaxa(insect2)
sample_names(insect2)[1:5]
sample_names(insect2)
sample_variables(insect2)
rank_names(insect2)

# Number of distinct phyla
tax_table(insect2)[,"Phylum"] %>% unique %>% na.exclude %>% length
tax_table(insect2)[,"Family"] %>% unique %>% na.exclude %>% length

## Number of samples per variable
table(sample_data(insect2)$State)
table(sample_data(insect2)$Species)
table(sample_data(insect2)$Genus)

### Remove taxa not seen more than 3 times in at least 1% of the samples

#insect = filter_taxa(insect2, function(x) sum(x > 3) > (0.01*length(x)), TRUE)
#ntaxa(insect)

###  Change order of bars
# https://github.com/joey711/phyloseq/issues/291
# sample_data(insect)$State <- factor(sample_data(insect)$State, levels = c("Cambridge_MA_Late", "Cambridge_MA_Mid", "Cambridge_MA_early"))
# levels(sample_data(insect)$State)
# insect_m <- merge_samples(insect, "State")

########################## 
### Most Abundant taxa ###
##########################
## Use this to decide if any should be filtered out

insect_b <- insect2
insect_b <- subset_taxa(insect_b, Class=="c__Chloroplast")
insect_b <- subset_samples(insect2, Specialist=="Yes")
insect_b <- subset_samples(insect2, Genus=="Peponapis")
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Time_Course=="Yes")

insect_b <- subset_samples(insect2, Genus=="Diabrotica")

insect_b <- prune_taxa(taxa_sums(insect_b) > 0, insect_b)

source("../taxa_summary.R", local = TRUE)
mdt = fast_melt(insect_b)

prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Sort by total counts in ascending order 
prevdt <- prevdt[order(-prevdt$TotalCounts),]
prevdt

write.table(prevdt, "remaining_Chloroplast.txt", sep="\t")
write.table(prevdt, "16S-otus_Acalymma_time_course_Wolbachia_or_Chloroplast.txt", sep="\t")
write.table(prevdt, "16S-otus_in_Acalymma_no_Wolbachia_or_Chloroplast.txt", sep="\t")
write.table(prevdt, "16S_deNovo_otus_in_Diabrotica_no_Wolbachia_or_Chloroplast.txt", sep="\t")
write.table(prevdt, "16S_deNovo_otus_in_bees_no_Wolbachia_or_Chloroplast.txt", sep="\t")

# Sort by total counts in descending order 
prevdt[order(-prevdt$TotalCounts),]

# Sort by prevalence in ascending order 
prevdt[order(prevdt$Prevalence),]

# Sort by prevalence in descending order
prevdt[order(-prevdt$Prevalence),]

### Rarefy data
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
#########################################
### Everything - Acalymma + Peponapis ###
#########################################

## Everything - Acalymma + Peponapis
insect_b <- insectP
insect_b <- insectR
insect_b <- insect2
insect_b <- subset_samples(insect2, Genus=="Peponapis")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Order=="o__Lactobacillales")

insect_R <- insectR
insect_R <- subset_samples(insectR, Specialist=="Yes")

## Alpha richness
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(insect_b, "sample_Species", "State2", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=State2, y=value, color=sample_Species), alpha=0.1)

##### Simple ordination of everything

## Configuration 1
title = "PCoA of Bray-Curtis distance, everything"
insectR_ord = ordinate(insectR, "PCoA", "bray")
p = plot_ordination(insect_R, insectR_ord, shape = "sample_Species", color = "State2")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
p + facet_wrap(~sample_Species)

#############################
## Proportional bar charts ##
#############################

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

##### Simple ordination of everything

## Configuration 1
title = "PCoA of Bray-Curtis distance, everything"
insectR_ord = ordinate(insectR, "PCoA", "bray")
p = plot_ordination(insectR, insectR_ord, color = "sample_Species", shape = "Region")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
p + facet_wrap(~State2)

## Configuration 2
title = "PCoA of Bray-Curtis distance, everything"
insectR_ord = ordinate(insectR, "PCoA", "bray")
p = plot_ordination(insectR, insectR_ord, color = "sample_Species", shape = "State2")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
p + facet_wrap(~sample_Species)

### Initial exploration of data

otu_table(insectR100)[1:5, 1:10]
tax_table(insectRm)[1:5, 1:4]
myTaxa = names(sort(taxa_sums(insectRm), decreasing = TRUE)[1:10])
myTaxa
insectRm
nsamples(insectRm)
ntaxa(insectRm)
sample_names(insectRm)[1:5]
sample_names(insectRm)
sample_variables(insectRm)
rank_names(insectRm)

### 

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
insect_b = prune_taxa(TopNOTUs, insect_b)

title="Distribution of 30 most abundant OTUs in Acalymma and Peponapis"
plot_tree(insect_b, nodelabf=nodeplotblank, ladderize="left", label.tips="taxa_names", shape="Phylum", color="sample_Species",title=title)


dev.off()

###########
## bees! ##
###########

insect_b <- insectP
insect_b <- insectR
insect_b <- insect2
insect_b <- subset_samples(insect2, State2!="Guanajuato")
insect_b <- subset_samples(insect_b, Genus=="Peponapis")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
ntaxa(insect_b)
nsamples(insect_b)

title="PCoA (bray) from all bacteria in Peponapis"
insect.ord <- ordinate(insect_b, "PCoA", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State", shape="State")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)

## Alpha richness
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(insect_b, "State", "State", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=State, y=value, color=State), alpha=0.1)

### Barplots
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

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 20 OTUs in Peponapis"
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2

## Trees 

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
insect_b = prune_taxa(TopNOTUs, insect_b)

title="Distribution of 30 most abundant OTUs in Acalymma and Peponapis"
plot_tree(insect_b, nodelabf=nodeplotblank, ladderize="left", label.tips="taxa_names", shape="Phylum", color="State2",title=title)


### Subset to proteobacteria in peponapis

insect_b <- subset_samples(insect2, Genus=="Peponapis")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")

title="PCoA (bray) from Proteobacteria in Peponapis"
insect.ord <- ordinate(insect_b, "PCoA", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State", shape="State")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)

insect_b <- subset_samples(insect2, Genus=="Peponapis")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top Proteobacteria OTUs in Peponapis (grouped by Phylum)"
p2 <- plot_bar(insect_prop, x="State_detail", fill="Order", title=title)
p2

########################
## Beetles - Acalymma ##
########################

######### PCoA all bacteria in Acalymma
insect_b <- insectR ## Rarefied to 1000
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
ntaxa(insect_b)
sample_sums(insect_b)

title="PCoA (bray) from all bacteria in Acalymma"
insect.ord <- ordinate(insect_b, "PCoA", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State", shape="Species")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)

## Alpha richness
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(insect_b, "State", "State", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=State, y=value, color=State), alpha=0.1)

#############
insect_b <- insectR
insect_b <- insect2
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Genus=="Diabrotica")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")

insect_b <- subset_samples(insect2, Time_Course=="Yes")
insect_b <- subset_samples(insect2, State2=="Iowa")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
ntaxa(insect_b)


title="PCoA (bray) from all bacteria in Acalymma"
insect.ord <- ordinate(insect_b, "PCoA", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State", shape="Species")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)

insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")

title="PCoA (bray) from Proteobacteria in Acalymma"
insect.ord <- ordinate(insect_b, "PCoA", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State", shape="Species")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)

## Alpha richness
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(insect_b, "State", "State", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=State, y=value, color=State), alpha=0.1)

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

## Barplots
insect_b <- insect2
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect_b, State_detail!="Cambridge_MA_Late-2")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")


insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")
insect_b <- prune_taxa(taxa_sums(insect_b) > 1, insect_b)
ntaxa(insect_b)

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

#p <- plot_bar(insect_b, x="Description", fill="OTU")
#p

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 20 OTUs in Acalymma"
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2

## Trees
TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

title="Distribution of 30 most abundant OTUs in Acalymma and Peponapis"
plot_tree(insect_b, nodelabf=nodeplotblank, ladderize="left", label.tips="taxa_names", shape="sample_Species", color="Order",title=title)

#######################################
## Beetles - Acalymma and Diabrotica ##
#######################################

## Create a dummy variable representing both Acalymma and Diabrotica

insect_R <- insectR
insect_R <- subset_samples(insect_R, Spotted_Or_Striped!="Neither")

insect_R <- subset_samples(insect_R, Genus=="Acalymma")
insect_R <- subset_samples(insect_R, Type=="beetle")

insect_b <- insect2
insect_b <- subset_samples(insect2, Type=="beetle")
insect_b <- subset_samples(insectR, State2=="Iowa")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
ntaxa(insect_b)
nsamples(insect_b)

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)


### Acalymma and diabrotica in Iowa PCoA
insect_R <- insectR
insect_R <- subset_samples(insect_R, Spotted_Or_Striped!="Neither")
insect_b <- subset_samples(insectR, State2=="Iowa")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)
ntaxa(insect_b)
nsamples(insect_b)

title="PCoA (bray)  in Acalymma + Diabrotica from Iowa"
insect.ord <- ordinate(insect_b, "PCoA", "bray")
p = plot_ordination(insect_R, insect.ord, shape = "Species", color="Species")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)

## Alpha richness
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(insect_b, "Species", "Species", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=Species, y=value, color=Species), alpha=0.1)

p <- plot_bar(insect_b, x="Description", fill="OTU")
p

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 20 OTUs in Acalymma + Diabrotica"
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2

## Trees
TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
insect_b = prune_taxa(TopNOTUs, insect_b)

## Ordination 
## Configuration 1
title = "PCoA of Bray-Curtis distance, everything"
insectR_ord = ordinate(insect_R, "PCoA", "bray")
p = plot_ordination(insect_R, insectR_ord, color = "State3", shape = "sample_Species")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
p + facet_wrap(~State3)

#######################################################################


biom_file = read_biom("otu_filter_time-course.biom")
map_file = ("time_course-map.txt")
trefile = ("rep_set.tre")
insect.tree <- read.tree(trefile)

biomot = import_biom(biom_file, parseFunction = parse_taxonomy_default)
bmsd = import_qiime_sample_data(map_file)
class(bmsd)
dim(bmsd)

#tax_table <- ("seqs_rep_set_tax_assignments.txt")
time.course = merge_phyloseq(biomot, bmsd, insect.tree)

colnames(tax_table(time.course)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")

### Many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate
time.course2 <- prune_taxa(taxa_sums(time.course) > 0, insect2)
time.course <- time.course2
ntaxa(time.course)

### Initial exploration of data

otu_table(insect2)[1:5, 1:10]
tax_table(insect2)[1:5, 1:4]
myTaxa = names(sort(taxa_sums(insect2), decreasing = TRUE)[1:10])
myTaxa
insect2
nsamples(insect2)
ntaxa(insect2)
sample_names(insect2)[1:5]
sample_names(insect2)
sample_variables(insect2)
rank_names(insect2)

# Number of distinct phyla
tax_table(insect2)[,"Family"] %>% unique %>% na.exclude %>% length

#############
### Frass ###
#############

insect
insect_b <- subset_samples(insect2, Type=="beetle")
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Genus=="Diabrotica")
insect_b <- subset_samples(insect2, Genus=="Peponapis")

insect_b <-insect2
insect_b <- subset_samples(insect2, State2=="New_Mexico")
insect_b <- subset_samples(insect2, State2=="Guanajuato")
insect_b <- subset_samples(insect2, State2=="Iowa")
insect_b <- subset_samples(insect2, State2=="Vermont")
insect_b <- subset_samples(insect2, State2=="California")
insect_b <- subset_samples(insect2, State2=="Massachussetts")

insect_b <-insect2
insect_b <- subset_samples(insect2, Species=="Acalymma_vittatum")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")

insect_b <- subset_samples(insect2, Species=="Acalymma_trivittatum")
insect_b <- subset_samples(insect2, Species=="Diabrotica_undecimpunctata")
insect_b <- subset_samples(insect2, Spotted_Or_Striped=="Striped")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

title="PCoA (bray) from all bacteria in Peponapis"
insect.ord <- ordinate(insect_b, "PCoA", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State", shape="State")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)

title="Proteobacteria in striped beetles collected early, mid, late season"
insect.ord <- ordinate(insect_b, "MDS", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State", shape="Species")
p + geom_point(size = 5, alpha = 0.8) + ggtitle(title)


p + geom_polygon(aes(fill=State)) + geom_point(size=5) + ggtitle("samples")


insect.ord <- ordinate(insect_b, "NMDS", "bray")
p1 = plot_ordination(insect_b, insect.ord, type="taxa", color="Family", title="taxa")
print(p1) 
p1 + facet_wrap(~Phylum, 3)

## Plot samples

title="PCoA of Peponapis"
insect.ord <- ordinate(insect_b, "MDS", "bray")
p = plot_ordination(insect_b, insect.ord, color = "State2", shape="Region")
p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)
#p + geom_polygon(aes(fill=Region)) + geom_point(size=5) + ggtitle("samples")

################
## Plot phlya ##
################

## To use the unsubsetted dataset
insect_b <-insect2

insect
insect_b <- subset_samples(insect2, Type=="beetle")
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Genus=="Diabrotica")
insect_b <- subset_samples(insect2, Genus=="Peponapis")

insect_b <- subset_samples(insect2, Species=="Acalymma_vittatum")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")

insect_b <- subset_samples(insect2, Species=="Acalymma_trivittatum")
insect_b <- subset_samples(insect2, Species=="Diabrotica_undecimpunctata")
insect_b <- subset_samples(insect2, Spotted_Or_Striped=="Striped")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Class=="c__Bacilli")
insect_b <- subset_taxa(insect_b, Order=="o__Lactobacillales")
insect_b <- subset_taxa(insect_b, Family=="f__Lactobacillaceae")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

insect_b <- subset_samples(insect2, State2=="Vermont")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

kingdom <- plot_bar(insect_b, "sample_Species",fill="Phylum", title="Lactobacilliaceae genera in Acalymma Guts") + coord_flip() + 
  xlab("Species")  
kingdom + theme(axis.title.y = element_text(face="bold", size=14),
                axis.text.y  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(lineheight=1.0, face="bold")) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, siz=18))


## Plot genera in proteobacteria
insect_proteo = subset_taxa(insect2, Phylum=="p__Proteobacteria")

TopNOTUs = names(sort(taxa_sums(insect_proteo), TRUE)[1:15])
insect_proteo = prune_taxa(TopNOTUs, insect_proteo)


Proteo <- plot_bar(insect_proteo, "sample_Genus", fill="Order", title="Proteobacteria in Beetle and Bee Guts") + coord_flip() + 
  xlab("Species") 
Proteo + theme(axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(lineheight=1.0, face="bold")) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, siz=18))

## Plot genera in Enterobacteriaceae
insect_entero = subset_taxa(insect2, Family=="f__Enterobacteriaceae")
insect_entero <- subset_taxa(insect2, Genus=="g__Erwinia")

TopNOTUs = names(sort(taxa_sums(insect_entero), TRUE)[1:15])
insect_proteo = prune_taxa(TopNOTUs, insect_entero)

Entero <- plot_bar(insect_entero, "sample_Genus", fill="Genus", title="Enterobacteriaceae in Beetle and Bee Guts") + coord_flip() + 
  xlab("Species") 
Entero + theme(axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(lineheight=1.0, face="bold")) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, siz=18))

##########################
## Barplots with facets ##
##########################

## To use the unsubsetted dataset
insect_b <-insect2

insect
insect_b <- subset_samples(insect2, Type=="beetle")
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Genus=="Diabrotica")
insect_b <- subset_samples(insect2, Genus=="Peponapis")

insect_b <- subset_samples(insect3, State2=="Vermont")

insect_b <- subset_samples(insect3, Species=="Acalymma_vittatum")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")

insect_b <- subset_samples(insect3, Species=="Acalymma_trivittatum")
insect_b <- subset_samples(insect3, Species=="Diabrotica_undecimpunctata")
insect_b <- subset_samples(insect3, Spotted_Or_Striped=="Striped")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Class=="c__Bacilli")
insect_b <- subset_taxa(insect_b, Order=="o__Lactobacillales")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:10])
ent10 = prune_taxa(TopNOTUs, insect_b)

# title="Top Proteobacterial genera in Massachusetts Acalymma\ncollected early, middle, and late season"
title ="bugs"
#p <- plot_bar(ent10, "Sample", fill = "Phylum",facet_grid=~Region, title=title) 
p <- plot_bar(insect_b, "State", fill = "Phylum", title=title) + facet_wrap(~Phylum, 3)
p + theme(plot.title = element_text(hjust = 0.5, siz=18, lineheight=1.0, face="bold"))

p <- plot_bar(insect2, "Sample", fill = "Genus", title=title)
p + facet_wrap(~State2, 1)

p = plot_bar(ent10, "sample_Species", fill="Genus", facet_grid=sample_Species~Region)
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

#############
## Heatmap ##
#############
## https://github.com/joey711/phyloseq/issues/230

## To use the unsubsetted dataset
insect_b <-insect2

insect
insect_b <- subset_samples(insect2, Type=="beetle")
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Genus=="Diabrotica")
insect_b <- subset_samples(insect2, Genus=="Peponapis")

insect_b <- subset_samples(insect2, Species=="Acalymma_vittatum")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")

insect_b <- subset_samples(insect2, Species=="Acalymma_trivittatum")
insect_b <- subset_samples(insect2, Species=="Diabrotica_undecimpunctata")
insect_b <- subset_samples(insect2, Spotted_Or_Striped=="Striped")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

rank_names(insect_b)
title="Erwinia in Acalymma beetle guts"
heat <- plot_heatmap(insect_b, sample.label="State", sample.order="State", taxa.label="Family", method="Bray", low="#000033", high="#FF3300", title=title)
heat + theme(axis.title.x = element_text(face="bold", size=12),
             axis.text.x  = element_text(vjust=0.5, size=12, face="bold"))


heatmap(otu_table(insect_b), Colv = Type)

###################
## Plot Network. ##
###################
# I'm not really sure what this means...

insect_b <-insect2

insect
insect_b <- subset_samples(insect2, Genus=="Acalymma")
insect_b <- subset_samples(insect2, Genus=="Diabrotica")
insect_b <- subset_samples(insect2, Genus=="Peponapis")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

ig = make_network(insect_b, max.dist = .9)
plot_network(ig, insect_b, color = "State2", shape = "Species", line_weight = 0.4, 
             label = NULL)

ig = make_network(insect2, max.dist = 0.9)
plot_network(ig, insect2, color = "Genus", shape = "Genus", line_weight = 0.4, 
             label = NULL)


#####################
## Alpha diversity ##
#####################

## To use the unsubsetted dataset
insect_b <-insect2

insect
insect_b <- subset_samples(insect3, Type=="beetle")
insect_b <- subset_samples(insect3, Genus=="Acalymma")
insect_b <- subset_samples(insect3, Genus=="Diabrotica")
insect_b <- subset_samples(insect3, Genus=="Peponapis")

insect_b <- subset_samples(insect3, Species=="Acalymma_vittatum")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")

insect_b <- subset_samples(insect3, Species=="Acalymma_trivittatum")
insect_b <- subset_samples(insect3, Species=="Diabrotica_undecimpunctata")
insect_b <- subset_samples(insect3, Spotted_Or_Striped=="Striped")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

pAlpha = plot_richness(insect_b,
                       shape = "Species",
                       color = "State2",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "Alpha diversity in bug guts")
pAlpha + geom_point(size = 1)

### More detailed alpha diversity

# Store as a new data variable
alphadt = data.table(pAlpha$data)
setnames(alphadt, "variable", "measure")
# Subset to just the Shannon index
alphadt <- alphadt[(measure == "Shannon")]
# Order by Days
alphadt <- alphadt[order(Species)]
# Define the plot
ggplot(data = alphadt, 
       mapping = aes(State2, value, color = State2)) + 
  geom_point(size = 5) + 
  geom_path() +
  geom_point(data = alphadt[(Genus == "Diabrotica")], 
             size = 8, alpha = 0.35) +
  #  facet_wrap(~State, ncol = 2, scales = "free_y") +
  ylab("Shannon Index") +
  ggtitle("Shannon index by state Peponapis (alpha diversity)")

###############################
#### Differential abundance ###
###############################

## To use the unsubsetted dataset
insect_b <-insect3

insect
insect_b <- subset_samples(insect3, Type=="beetle")
insect_b <- subset_samples(insect3, Genus=="Acalymma")
insect_b <- subset_samples(insect3, Genus=="Diabrotica")
insect_b <- subset_samples(insect3, Genus=="Peponapis")

insect_b <- subset_samples(insect3, Species=="Acalymma_vittatum")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")

insect_b <- subset_samples(insect3, Species=="Acalymma_trivittatum")
insect_b <- subset_samples(insect3, Species=="Diabrotica_undecimpunctata")
insect_b <- subset_samples(insect3, Spotted_Or_Striped=="Striped")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

ent.p.table <- mt(insect, "State", test = "f")
head(insect.fwer.table)
print(head(ent.p.table, 10))
names(ent.p.table)
et_Ftest <- ent.p.table[order(ent.p.table$adjp),]
write.table(et_Ftest, "Et_State_Beetle_Ftest.txt", sep="\t")



#### Differential abundance in DESeq2
# Convert to DESeq2's DESeqDataSet class
library("DESeq2")

packageVersion("DESeq2")

head(sample_data(insect3)$Region, 5)

hist(log10(apply(otu_table(insect3), 1, var)), xlab = "log10(variance)", main = "A large fraction of OTUs have very low variance")

diagdds = phyloseq_to_deseq2(insect3, ~ State)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(insect3)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.table(sigtab, "Peponapis_Diff_Abun.txt", sep="\t")

## Inpsect the OTUs that were significantly different between the two tissues. The following makes a nice ggplot2 summary of the results.

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
x
sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels = names(x))

p <- ggplot(sigtab, aes(x = Phylum, y = log2FoldChange, color = Phylum)) + geom_point(size = 6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
p

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
x
sigtab$Genus <- factor(as.character(sigtab$Genus), levels = names(x))

p <- ggplot(sigtab, aes(x = Genus, y = log2FoldChange, color = Genus)) + geom_point(size = 6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
p

# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
x
sigtab$Family <- factor(as.character(sigtab$Family), levels = names(x))

p <- ggplot(sigtab, aes(x = Family, y = log2FoldChange, color = Family)) + geom_point(size = 6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
p 

## Log transform data for variance stabilization
# pslog = transform_sample_counts(insect, function(x) log(1 + x))


## Distribution of reads
## THis is post filtering - should be very even (is)

sdt = data.table(as(sample_data(insect), "data.frame"),
                 TotalReads = sample_sums(insect), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
pSeqDepth + facet_wrap(~State)
pSeqDepth + 
  facet_grid(~State) + 
  ggtitle("Seq. Depth by Collection period")

tdt = data.table(tax_table(insect),
                 TotalCounts = taxa_sums(insect),
                 OTU = taxa_names(insect))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts")
# How many singletons?
tdt[(TotalCounts <= 0), .N]
tdt[(TotalCounts <= 1), .N]
# None. This data has been filtered already. 
# How many doubletons?
tdt[(TotalCounts <= 2), .N]
# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pCumSum
# Zoom-in
pCumSum + xlim(0, 100)

####
#### Taxa prevalence histogram, and fast_melt()
####

insect
insect_b <- subset_samples(insect, Type=="beetle")
insect_b <- subset_samples(insect, Genus=="Acalymma")
insect_b <- subset_samples(insect, Genus=="Diabrotica")
insect_b <- subset_samples(insect, Genus=="Peponapis")
insect_b <- subset_samples(insect, Species=="Acalymma_vittatum")
insect_b <- subset_samples(insect, Species=="Acalymma_trivittatum")
insect_b <- subset_samples(insect, Species=="Diabrotica_undecimpunctata")
insect_b <- subset_samples(insect, Spotted_Or_Striped=="Striped")
insect_b

insect_b <- subset_taxa(insect_b, Kingdom=="k__Bacteria")
insect_b <- subset_taxa(insect_b, Phylum=="p__Firmicutes")
insect_b <- subset_taxa(insect_b, Phylum=="p__Proteobacteria")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

source("../taxa_summary.R", local = TRUE)
mdt = fast_melt(insect_b)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Sort by total counts in ascending order 
prevdt <- prevdt[order(-prevdt$TotalCounts),]
prevdt
write.table(prevdt, "Top_bacteria_in_Striped_Cucumber_Beetles.txt", sep="\t")

# Sort by total counts in descending order 
prevdt[order(-prevdt$TotalCounts),]

# Sort by prevalence in ascending order 
prevdt[order(prevdt$Prevalence),]

# Sort by prevalence in descending order
prevdt[order(-prevdt$Prevalence),]


mdt = fast_melt(insect)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
ggplot(prevdt, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Histogram of Taxa Prevalence")
# How many singletons?
prevdt[(Prevalence <= 0), .N]
prevdt[(Prevalence <= 1), .N]
# None. This data has been filtered already. 
# How many doubletons?
prevdt[(Prevalence <= 2), .N]
# taxa cumulative sum
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pPrevCumSum

### Prevalence vs. Total Count Scatter plot

ggplot(prevdt, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

## Subset to the most abundant 9 phyla, and map these to color in a ggplot2 scatter plot.

addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
showPhyla = prevdt[, sum(TotalCounts), by = Phylum][order(-V1)][1:9]$Phylum
setkey(prevdt, Phylum)
ggplot(prevdt[showPhyla], 
       mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

## Alpha diversity in beetle guts
pAlpha = plot_richness(insect,
                       #               shape = "State",
                       color = "State",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "Alpha Diveristy, Moving Pictures Data")
pAlpha + geom_point(size = 1)

### More detailed alpha diversity

# Store as a new data variable
alphadt = data.table(pAlpha$data)
setnames(alphadt, "variable", "measure")
# Subset to just the Shannon index
alphadt <- alphadt[(measure == "Shannon")]
# Order by Days
alphadt <- alphadt[order(State)]
# Define the plot
ggplot(data = alphadt, 
       mapping = aes(State, value, color = State)) + 
  geom_point(size = 5) + 
  geom_path() +
  geom_point(data = alphadt[(Genus == "Acalymma")], 
             size = 8, alpha = 0.35) +
  facet_wrap(~measure, ncol = 2, scales = "free_y") +
  ylab("Shannon Index") +
  ggtitle("Shannon Index in Cucumber Beetles")

### Filter all taxa that occur in only one sample
# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.

keepTaxa = prevdt[(Prevalence >= 2 & TotalCounts > 3), TaxaID]
# Define new object with relative abundance
mpra = transform_sample_counts(insect, function(x) x / sum(x))
# Filter this new object
mpraf = prune_taxa(keepTaxa, mpra)
# Calculate distances
DistBC = phyloseq::distance(insect3, method = "bray")
DistUF = phyloseq::distance(insect3, method = "wUniFrac")

# Multidimensional Scaling (MDS aka PCoA)
## Perform MDS ("PCoA") Ordination 

ordBC = ordinate(insect3, method = "PCoA", distance = DistBC)
ordUF = ordinate(insect3, method = "PCoA", distance = DistUF)

## Scree plot

plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_scree(ordUF, "Scree Plot: Weighted UniFrac MDS")

# Plot a PCoA for each distance-and-ordination using the `plot_ordination()` function in phyloseq.

plot_ordination(insect3, ordBC, color = "Genus") + 
  #  geom_point(mapping = aes(shape = factor(Genus))) +
  geom_path() +
  ggtitle("PCoA: Bray-Curtis")

###

plot_ordination(mpraf, ordUF, color = "State") + 
  geom_point(mapping = aes(size = DaysSinceExperimentStart, 
                           shape = factor(Subject))) +
  ggtitle("PCoA: Weigthed Unifrac")

### Bray-Curtis

# Join sample data and ordination axes together in one data.table
ordBCdt = data.table(ordBC$vectors, keep.rownames = TRUE)
setnames(ordBCdt, "rn", "SampleID")
setkey(ordBCdt, SampleID)
setkey(sdt, SampleID)
ordBCdt <- ordBCdt[sdt]
setorder(ordBCdt, State)
# Axis 1
ggplot(ordBCdt, aes(State, Axis.1, color = factor(State))) +
  geom_point(size = 5) + 
  geom_path() +
  facet_wrap(~Genus) + ggtitle("Bray-Curtis PCoA Axis-1 vs. Day")
# Axis 2
ggplot(ordBCdt, aes(State, Axis.2, color = factor(State))) +
  geom_point(size = 5) + 
  geom_path() +
  facet_wrap(~Genus) + ggtitle("Bray-Curtis PCoA Axis-2 vs. Day")

## Subset based on Erwinia

sample_data(insect)$State <- factor(sample_data(insect)$State, levels = c("Cambridge_MA_Late", "Cambridge_MA_Mid", "Cambridge_MA_early"))
levels(sample_data(insect)$State)

# Family Enterobacteriaceae
insect_Et <- subset_taxa(insect, Family == "f__Enterobacteriaceae")
et_fam <- plot_bar(insect_Et, x="State",fill="Genus") + coord_flip() + 
  ggtitle("Enterobacteriaceae in Cucumber Beetle Guts")
xlab("Time Period of Collection")  
et_fam + theme(axis.title.y = element_text(face="bold", size=14),
               axis.text.y  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(lineheight=1.0, face="bold",size=10,)) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, siz=18))

# Genus Erwinia
insect_Et <- subset_taxa(insect, Genus == "g__Erwinia")
et_et <-plot_bar(insect_Et, x="State",fill="Genus") + coord_flip() + 
  ggtitle("Erwinia in Cucumber Beetle Guts")
xlab("Time Period of Collection")  
et_et + theme(axis.title.y = element_text(face="bold", size=14),
              axis.text.y  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(lineheight=1.0, face="bold",size=10,)) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x  = element_text(vjust=0.5, size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, siz=18))