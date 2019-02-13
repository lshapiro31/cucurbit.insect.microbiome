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

setwd("/Users/lorishapiro/Dropbox/Matrix_forsharing_16S/Bee-beetle-only-beta")

biom_file = read_biom("otu_filter_meta__Type2_insect__.biom")
map_file = ("New_mapping_file2__Type2_insect__.txt")
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

## Prune the sparse taxa that show up only 2, 1 or zero times

#insect2 <- prune_taxa(taxa_sums(insect2) > 0, insect2)
#insect2 <- prune_taxa(taxa_sums(insect2) > 1, insect2)
insect2 <- prune_taxa(taxa_sums(insect2) > 2, insect2)
ntaxa(insect2)

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

######################
### Summarize taxa ###
######################

source("../taxa_summary.R", local = TRUE)
mdt = fast_melt(insect2)

prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Sort by total counts in ascending order 
prevdt <- prevdt[order(-prevdt$TotalCounts),]
prevdt

write.table(prevdt, "otus.txt", sep="\t")

###################
### Rarefy data ###
###################

set.seed(28132)
sample_sums(insect2)
insectR = rarefy_even_depth(insect2, sample.size = 1000)
sample_sums(insectR)

### phyloseq object insectR is now ready for use!

insect # phyloseq object with no filtering
insectR #phyloseq object where wolbachia + chloroplasts have been removed, normalized to 1000 reads per sample
