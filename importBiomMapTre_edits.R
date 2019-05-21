library(biomformat)
library(vegan)
library(plotly)
library(phyloseq)
library(ggplot2)
library(genefilter)
library(PMA)
library(ade4)
library(ggrepel)
library(metagenomeSeq)
library(ape)
library(tidyverse)
library(purrr)
library(magrittr)
library(stringr)
library(stringi)
library(RColorBrewer)
library(viridis)
library(viridisLite)
library(plotly)
library(colormap)

theme_set(theme_bw())

setwd("~/Dropbox/Matrix_forsharing_16S_backup/Bee-beetle-only-beta")

## Pick global colors for all graphs

display.brewer.all() 
insects.12 <- brewer.pal(12, name = 'Paired')
# "#A6CEE3" Light Blue - states.7
# "#1F78B4" Dark Blue - colors.time.course
# "#B2DF8A" Light Green - colors.time.course
# "#33A02C" Dark Green - states.7
# "#FB9A99" Light Red - states.7
# "#E31A1C" Dark Red - states.7
# "#FDBF6F" Light Orange - bees.2
# "#FF7F00" Dark Orange - states.7
# "#CAB2D6" Light Purple - states.7
# "#6A3D9A" Dark Purple - states.7
# "#FFFF99" Yellow - colors.time.course
# "#B15928" Brown - bees.2

## For Acalymma, without massachusetts
states.7 <- c("#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#A6CEE3", "#6A3D9A")
print(states.7)

## For massachusetts time course
colors.time.course <- c("#A6CEE3", "#FFFF99", "#B2DF8A")
print(colors.time.course)

## For the bees
bees.2 <- c("#B15928", "#FDBF6F")
print(bees.2)

tree_colors <- c(
  "#B15928", "#FDBF6F",
  "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#A6CEE3", "#6A3D9A"
)

## Read in data
biom_file = read_biom("otu_table_filter_meta__Experiment_insect__.biom")
# map_file = ("New_mapping_file2.txt")
map_file = ("New_mapping_file_edits.txt")
trefile = ("rep_set.tre")
insect.tree <- read.tree(trefile)

biomot = import_biom(biom_file, parseFunction = parse_taxonomy_default)
bmsd = import_qiime_sample_data(map_file)
class(bmsd)
dim(bmsd)

#tax_table <- ("seqs_rep_set_tax_assignments.txt")
insect = merge_phyloseq(biomot, bmsd, insect.tree)

colnames(tax_table(insect)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(insect)

insect2 <- insect


## Prune the sparse taxa 
#insect2 <- prune_taxa(taxa_sums(insect2) > 0, insect2)
#insect2 <- prune_taxa(taxa_sums(insect2) > 1, insect2)
insect2 <- prune_taxa(taxa_sums(insect2) > 5, insect2) 
ntaxa(insect2)

### Alternate tidyverse solution to omitting otus

insect.filt <- insect2 %>%
  subset_taxa(
      Kingdom == "k__Bacteria" &
      Family  != "f__mitochondria" &
      Class   != "c__Chloroplast" &
      Genus   != "g__Wolbachia"
  )

##################
## Otus to omit ##
##################

# https://github.com/joey711/phyloseq/issues/652
badTaxa = read.table("wolbachia.chloroplast.mitochondria.otus", header=F) ## Wolbachia; should go back and look at distribution of htis later
badTaxa<-as.character(badTaxa[[1]])
allTaxa = taxa_names(insect2)
length(allTaxa)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
length(allTaxa)
insect_filt <- prune_taxa(allTaxa, insect2)
ntaxa(insect_filt)
insect2 <- prune_taxa(taxa_sums(insect_filt) > 10, insect_filt) 
ntaxa(insect2)

# Subset the dataset to just insects from Peponapis and Acalymma
acalymma.peponapis <- insect2 %>% 
  subset_samples(
    Genus=="Peponapis" | Genus=="Acalymma" &
    State_detail != "Cambridge_MA_Late-2"   # Omit this sample due to insufficient sequencing coverage
  ) %>%
  prune_taxa(taxa_sums(.) > 0, .) 
acalymma.peponapis


###################
### Rarefy data ###
###################

set.seed(28132)
sample_sums(acalymma.peponapis)
acalymma.peponapis.R = rarefy_even_depth(acalymma.peponapis, sample.size = 1000)
sample_sums(acalymma.peponapis.R)

acalymma.peponapis # phyloseq object filtered to remove non-bacterial reads
acalymma.peponapis.R #phyloseq object where wolbachia + chloroplasts have been removed, normalized to 1000 reads per sample

## Some summary statistics

# Number of distinct phyla and families
tax_table(acalymma.peponapis)[,"Phylum"] %>% unique %>% na.exclude %>% length
tax_table(acalymma.peponapis)[,"Family"] %>% unique %>% na.exclude %>% length

sum(sample_sums(acalymma.peponapis))

sample_variables(acalymma.peponapis)
levels(sample_data(acalymma.peponapis)$State2)

## Number of samples per variable
table(sample_data(acalymma.peponapis)$State)
table(sample_data(acalymma.peponapis)$State2)
table(sample_data(acalymma.peponapis)$Species)
table(sample_data(acalymma.peponapis)$Genus)
table(sample_data(acalymma.peponapis)$Region)
table(sample_data(insect2)$Spotted_Or_Striped)

######################
### Summarize taxa ###
######################

## acalymma.peponapis is UNRARIFIED data

## insectR is RAREFIED data use the rarefied data 
## for computing ecological indices and comparing richness, 
## alpha diversity, and things along those lines.

### Acalymma -----------------------------------------------------

acalymma <- subset_samples(acalymma.peponapis, Genus=="Acalymma")
acalymma <- prune_taxa(taxa_sums(insect2) > 5, insect2)

## Top OTUs in the entire dataset, nothing filtered out 
source("../taxa_summary.R", local = TRUE)
mdt = fast_melt(acalymma)

prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Sort by total counts in ascending order 
prevdt <- prevdt[order(-prevdt$TotalCounts),]
prevdt

write.table(prevdt, "Acalymma.otus.filtered.txt", sep="\t")

## Add a column that is the percentage of the total that each of the OTUs represents

prevdt.tbl<-as_tibble(prevdt)
str(prevdt.tbl)
prevdt.prop <- prevdt.tbl %>% mutate(
  Percentage =  TotalCounts / sum(TotalCounts) * 100
)

write.table(prevdt.prop, "Acalymma.filtered.otus.percent.txt", sep="\t")


### Peponapis -----------------------------------------------------

peponapis <- subset_samples(acalymma.peponapis, Genus=="Peponapis")
peponapis <- prune_taxa(taxa_sums(peponapis) > 5, peponapis)

## Top OTUs in the entire dataset, nothing filtered out 
source("../taxa_summary.R", local = TRUE)
mdt = fast_melt(peponapis)

prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Sort by total counts in ascending order 
prevdt <- prevdt[order(-prevdt$TotalCounts),]
prevdt

write.table(prevdt, "Peponapis.filtered.otus.txt", sep="\t")

## Add a column that is the percentage of the total that each of the OTUs represents

prevdt.tbl<-as_tibble(prevdt)
str(prevdt.tbl)
prevdt.prop <- prevdt.tbl %>% mutate(
  Percentage =  TotalCounts / sum(TotalCounts) * 100
)

write.table(prevdt.prop, "Peponapis.filtered.otus.percent.txt", sep="\t")

### Peponapis by state
peponapis.CA <- subset_samples(peponapis, State2=="California")
peponapis.CA <- prune_taxa(taxa_sums(peponapis.CA) > 5, peponapis)

mdt = fast_melt(peponapis.CA)

prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Sort by total counts in ascending order 
prevdt <- prevdt[order(-prevdt$TotalCounts),]
prevdt

write.table(prevdt, "Peponapis.CA.filtered.otus.txt", sep="\t")

peponapis.PA <- subset_samples(acalymma.peponapis, State2=="Pennsylvania")
peponapis.PA <- prune_taxa(taxa_sums(peponapis.PA) > 5, peponapis.PA)

mdt = fast_melt(peponapis.PA)

prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Sort by total counts in ascending order 
prevdt <- prevdt[order(-prevdt$TotalCounts),]
prevdt

write.table(prevdt, "Peponapis.PA.filtered.otus.txt", sep="\t")


