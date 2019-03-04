#################################
### Acalymma - MA time course ###
#################################

#insect_b <- insectR
insect_b <- insect2
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")
insect_b <- prune_taxa(taxa_sums(insect_b) > 0, insect_b)

display.brewer.all() 
colors_vec <- brewer.pal(3, name = 'Paired')
print(colors_vec)

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
insect_b = prune_taxa(TopNOTUs, insect_b)

###  Change order of samples
# https://github.com/joey711/phyloseq/issues/291
sample_data(insect_b)$State2
sample_data(insect_b)$State2 <- factor(sample_data(insect_b)$State2, levels = c("Massachusetts-early", "Massachusetts-mid", "Massachusetts-late"))

# Transform to percentages of total available.
insect_b_tr = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))

title = "PCoA of Bray-Curtis distance\nof Acalymma beta-diversity over one growing\nseason in Massachusetts"
insectb_ord = ordinate(insect_b_tr, "PCoA", "bray")
p <- plot_ordination(insect_b_tr, insectb_ord, shape = "sample_Species", color = "State2")
p + geom_point(size = 6, alpha = 0.8) + ggtitle(title) +
  scale_color_manual(values=colors_vec) +
  scale_fill_manual(values =c("black", "black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank())

ggsave("Acalymma.timeCourse.PcOA.pdf", height=5, width=7)


#################
### Peponapis ###
#################

insect_b <- insectR
insect_b <- subset_samples(insect_b, Genus=="Peponapis")
insect_b <- subset_samples(insect_b, State2!="Guanajuato")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
insect_b = prune_taxa(TopNOTUs, insect_b)

# Transform to percentages of total available.
insect_b_tr = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))

title = "PCoA of Bray-Curtis distance\nPeponapis pruinosa"
insectb_ord = ordinate(insect_b, "PCoA", "bray")
p <- plot_ordination(insect_b, insectb_ord, shape = "sample_Species", color = "State2")
p + geom_point(size = 6, alpha = 0.8) + ggtitle(title) +
  scale_color_manual(values = c("orange4", "orange")) +
#  scale_color_manual(values=c("black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank())

ggsave("Peponapis.PcOA.pdf", height=5, width=7)


## Configuration 2
title = "PCoA of Bray-Curtis distance, everything"
insectR_ord = ordinate(insectR, "PCoA", "bray")
p = plot_ordination(insectR, insectR_ord, color = "sample_Species", shape = "State2")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
p + facet_wrap(~sample_Species)

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