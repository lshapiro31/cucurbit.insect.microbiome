## Tree files of otus

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
insect_b = prune_taxa(TopNOTUs, insect_b)

title="Distribution of 30 most abundant OTUs in Acalymma and Peponapis"
plot_tree(insect_b, nodelabf=nodeplotblank, ladderize="left", label.tips="taxa_names", shape="Phylum", color="sample_Species",title=title)

title="Distribution of 30 most abundant OTUs in Acalymma and Peponapis"
plot_tree(insect_b, nodelabf=nodeplotblank, ladderize="left", label.tips="taxa_names", shape="Phylum", color="State2",title=title)

dev.off()