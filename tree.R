insect_b <- insect2 %>%
  subset_samples(
#    Spotted_Or_Striped =="Striped" | insect_Genus == "Peponapis"  ## For paper, with Peponapis
    Spotted_Or_Striped =="Striped"    ## For poster; w/o Peponapis
  ) %>%
  prune_taxa(taxa_sums(.) > 0, .)


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

########## With all the states individually labeled

## Pick some colors!
display.brewer.all() 
tree_colors <- brewer.pal(7, name = 'Dark2')
print(tree_colors)

insect_b <- insectR %>% 
  subset_samples(
    Specialist=="Yes" &
      Genus=="Acalymma" & 
      Time_Course=="No"
  )


TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:100])
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
                             justify = "right",  ## This option is what is making the labels disappear
                             color="State2",
                             base.spacing=0.03,
                             #scale_fill_manual(tree_colors)
                             #scale_colour_hue(values=tree_colors),
                             title="Frequency of the 100 most abundant OTUs in Acalymma") +
  scale_size_continuous(range = c(1, 5))


bee.beetle.tree + 
  scale_color_manual(values=c(tree_colors)) +
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=16),
        legend.title=element_blank())

ggsave("tree.State2.comparison.pdf", height=12, width=10)

dev.off()


################# 
### Peponapis ###
#################

## Subset top 20 otus in Peponapis
# Start with Peponapis
# Proportional barchart

insect_b <- insect2 %>%
  subset_samples(
    Genus=="Peponapis" & State2!="Guanajuato"
  ) %>% 
  prune_samples(sample_sums(.) > 0, .) 
  
TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:50])
insect_b = prune_taxa(TopNOTUs, insect_b)

insect_b <- transform_sample_counts(insect_b, function(x) x / sum(x))

bee.tree <- plot_tree(insect_b, 
                             label.tips="taxa_names", 
                             text.size = 2,
                             plot.margin = 0.2,
                             nodelabf=nodeplotblank, 
                             ladderize="left", 
                             size = "Abundance", 
                             shape="Phylum", 
                             justify = "right",  ## This option is what is making the labels disappear
                             color="State2",
                             base.spacing=0.03,
                             #scale_fill_manual(tree_colors)
                             #scale_colour_hue(values=tree_colors),
                             title="Frequency of the 50 most abundant OTUs in Peponapis") +
  scale_size_continuous(range = c(1, 5))

bee.tree + 
  scale_color_manual(values = c("orange4", "orange")) +
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=16),
        legend.title=element_blank())

ggsave("tree.bee.State2.comparison.pdf", height=12, width=10)

dev.off()

############################ 
## Acalymma - Time Course ##
############################

display.brewer.all() 
colors_vec <- brewer.pal(3, name = 'Paired')
print(colors_vec)

insect_b <- insect2 %>%
  subset_samples(
    Genus=="Acalymma" &
      Time_Course=="Yes"
  )  %>%
  prune_samples(sample_sums(.) > 0, .)

###  Change order of samples
sample_data(insect_b)$State2
sample_data(insect_b)$State2 <- factor(sample_data(insect_b)$State2, levels = c("Massachusetts-early", "Massachusetts-mid", "Massachusetts-late"))

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)

insect_b <- transform_sample_counts(insect_b, function(x) x / sum(x))

bee.tree <- plot_tree(insect_b, 
                      label.tips="taxa_names", 
                      text.size = 2,
                      plot.margin = 0.2,
                      nodelabf=nodeplotblank, 
                      ladderize="left", 
                      size = "Abundance", 
                      shape="Phylum", 
                      justify = "right",  ## This option is what is making the labels disappear
                      color="State2",
                      base.spacing=0.03,
                      #scale_fill_manual(tree_colors)
                      #scale_colour_hue(values=tree_colors),
                      title="Frequency of the 20 most abundant OTUs in Acalymma vittatum\nthroughout one growing season in Massachusetts") +
  scale_size_continuous(range = c(1, 5))

bee.tree + 
  scale_color_manual(values = colors_vec) +
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=16),
        legend.title=element_blank())

ggsave("tree.timeCourse.comparison.pdf", height=12, width=10)

dev.off()
