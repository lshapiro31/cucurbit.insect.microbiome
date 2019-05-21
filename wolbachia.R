## Make Wolbachia only dataset with unfiltered original phyloseq object

insect.wolbachia <- insect %>%  # Original unfiltered phyloseq object that still contains wolbachia
  subset_samples(Genus=="Peponapis" | Genus=="Acalymma") %>%
  subset_taxa(Genus   == "g__Wolbachia") %>%
  prune_taxa(taxa_sums(.) > 3, .) %>%   
  tax_glom(taxrank = "Genus") %>%
  prune_samples(sample_sums(.) >=20, .)
insect.wolbachia


wolbachia.plot <- plot_bar(insect.wolbachia, x="State_detail", fill="Genus") + 
  ggtitle("Insects with Wolbachia (unrarefied)") +
  ylab("Total Abundance of Wolbachia") +
  xlab("Individual Insects") +
  coord_flip()
wolbachia.plot + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
          axis.text.x  = element_text(face="bold", color = "black", hjust=0.5, size=10),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face="bold", color = "black", size=12),
          axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, face = "bold")) +
  scale_fill_discrete(labels=c(" Wolbachia"))

ggsave("Manuscript/wolbachia.barchart.pdf", height = 6, width = 8)

dev.off()