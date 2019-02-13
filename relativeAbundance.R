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

