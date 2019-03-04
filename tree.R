insect_b <- insect2
insect.striped <- subset_samples(insect_b, Spotted_Or_Striped=="Striped")  
insect.bee <- subset_samples(insect_b, insect_Genus=="Peponapis")
insect_b <- merge_phyloseq(insect.striped, insect.bee)
insect_b <- prune_taxa(taxa_sums(insect_b) > 0, insect_b)

#insect_b <- subset_samples(insect2, State2!="Guanajuato")
#insect_b <- subset_samples(insect_b, Genus=="Peponapis")
#insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
#insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:50])
insect_b = prune_taxa(TopNOTUs, insect_b)

insect_b <- transform_sample_counts(insect_b, function(x) x / sum(x))

bee.beetle.tree <- plot_tree(insect_b, 
                             label.tips="taxa_names", 
                             text.size = 2,
                             plot.margin = 0.2,
                             nodelabf=nodeplotblank, 
                             ladderize="left", 
                             size = "Abundance", 
                             shape="Phylum", 
                             justify = "yes please",  ## This option is what is making the labels disappear
                             color="sample_Species",
                             base.spacing=0.03,
                             #scale_colour_hue("#e5f5f9", "#99d8c9", "#2ca25f"),
                             title="Distribution of 30 most abundant OTUs in\nAcalymma and Peponapis") +
scale_size_continuous(range = c(1, 5))
bee.beetle.tree + 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=14),
        legend.title=element_blank())

ggsave("tree.comparison.pdf", height=10, width=10)


dev.off()
