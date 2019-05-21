## Plot which insects have Wolbachia, and how much Wolbachia they have
## Use the original unfiltered original phyloseq object

insect.wolbachia <- insect %>%
  subset_samples(Genus=="Peponapis" | Genus=="Acalymma") %>%
  subset_taxa(Genus   == "g__Wolbachia") %>%
  prune_taxa(taxa_sums(.) > 1, .) %>%   
  tax_glom(taxrank = "Genus") %>% 
  prune_samples(sample_sums(.) >=10, .) 

wolbachia.plot <- plot_bar(insect.wolbachia, x="State_detail", fill="Genus") + 
  ggtitle("Distribution of Wolbachia in Acalymma and Peponapis (unrarefied)") +
  ylab("Wolbachia") +
  xlab("Individual Insects") +
  coord_flip()
wolbachia.plot + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
          axis.text.x  = element_text(face="bold", color = "black", hjust=0.5, size=10),
          axis.title.y = element_text(face="bold", color = "black", size=12),
          axis.title.x = element_text(face="bold", color = "black", size=12),
          axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
          legend.title = element_text(size=12, face="bold", hjust = 0.5),
          legend.text = element_text(size = 10, face = "bold"))

ggsave("wolbachia.barchart.pdf", height = 6, width = 8)

dev.off()
