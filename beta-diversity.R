####### Simple ordination of everything

## Configuration 1
title = "PCoA of Bray-Curtis distance, everything"
insectR_ord = ordinate(insectR, "PCoA", "bray")
p = plot_ordination(insect_R, insectR_ord, shape = "sample_Species", color = "State2")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
p + facet_wrap(~sample_Species)