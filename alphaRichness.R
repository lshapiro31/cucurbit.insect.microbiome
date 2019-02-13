### Alpha richness

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