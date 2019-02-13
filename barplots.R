### Barplots
insect_b <- insect2
insect_b <- subset_samples(insect2, State2!="Guanajuato")
insect_b <- subset_samples(insect_b, Genus=="Peponapis")
insect_b <- subset_taxa(insect_b, Family=="f__Enterobacteriaceae")
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia")
#insect_b <- subset_samples(insect2, State=="California")

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)


p <- plot_bar(insect_b, x="State_detail", fill="OTU")
p

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 20 OTUs in Peponapis"
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2
