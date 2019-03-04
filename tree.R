insect_b <- insect2
insect.striped <- subset_samples(insect_b, Spotted_Or_Striped=="Striped")  
insect.bee <- subset_samples(insect_b, insect_Genus=="Peponapis")
insect_b <- merge_phyloseq(insect.striped, insect.bee)
insect_b <- prune_taxa(taxa_sums(insect_b) > 0, insect_b)

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:50])
insect_b = prune_taxa(TopNOTUs, insect_b)

insect_b <- transform_sample_counts(insect_b, function(x) x / sum(x))

#Select custom colors
#display.brewer.all() 
#colors_vec <- brewer.pal(2, name = 'Greens')
#print(colors_vec)

tree_colors <- c(
  "#A1D99B", "#31A354", "orange2"
)
  
bee.beetle.tree <- plot_tree(insect_b, 
                             label.tips="taxa_names", 
                             text.size = 2,
                             plot.margin = 0.2,
                             nodelabf=nodeplotblank, 
                             ladderize="left", 
                             size = "Abundance", 
                             shape="Phylum", 
                             justify = "right",  ## This option is what is making the labels disappear
                             color="sample_Species",
                             base.spacing=0.03,
                             #scale_fill_manual(tree_colors)
                             #scale_colour_hue(values=tree_colors),
                             title="Comparison of the 30 most abundant OTUs in\nAcalymma and Peponapis") +
scale_size_continuous(range = c(1, 5))

bee.beetle.tree + 
  scale_color_manual(values=c(tree_colors)) +
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=16),
        legend.title=element_blank())

ggsave("tree.comparison.pdf", height=12, width=10)

dev.off()

### Same tree with tip labels; being cut off for some reason with justify option above 

plot_tree(insect_b, label.tips="taxa_names")

ggsave("tree.comparison.otus-labels.pdf", height=12, width=10)

dev.off()
