############################################
### Visual comparison of shared and ########
### unique OTUs between bees and beetles ###
############################################

TopNOTUs = names(sort(taxa_sums(acalymma.peponapis), TRUE)[1:100])

acalymma.peponapis.tree <- acalymma.peponapis %>% 
 # subset_samples(Time_Course=="No") %>%
  prune_taxa(TopNOTUs, .) 

###  Change order of samples
sample_data(acalymma.peponapis.tree)$State2
sample_data(acalymma.peponapis.tree)$State2 <- factor(sample_data(acalymma.peponapis.tree)$State2, levels = c("California", "Pennsylvania", "Vermont", "Alabama", "Kentucky", "Iowa", "New_Mexico", "Arizona", "Guanajuato"))

plot_tree(acalymma.peponapis.tree, 
          nodelabf=nodeplotblank, 
          ladderize="left", 
          label.tips="taxa_names", 
          shape="Phylum", 
          method = "sampledodge",
          color="State2",
          size="Abundance",
          base.spacing=0.1,
          justify = TRUE,
          sizebase = 10,
          title=title)  +
  ggtitle("Distribution of ### most abundant OTUs in Acalymma and Peponapis") +
  scale_color_manual(values=c(tree_colors)) +
  scale_size_continuous(range = c(2, 8)) +
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=16),
        legend.title=element_blank())

ggsave("bee.beetle.tree.pdf", height = 10, width = 10)

#### ----------------------------------------------------------------------------------
# For the MA samples; have to change the mapping data in the biom file to get this to plot together
acalymma.timeCourse.tree <- acalymma.peponapis %>% 
  subset_samples(Time_Course=="Yes") %>%
  prune_taxa(TopNOTUs, .) 

###  Change order of samples
sample_data(acalymma.timeCourse.tree)$State
sample_data(acalymma.timeCourse.tree)$State <- factor(sample_data(acalymma.timeCourse.tree)$State, levels = c("Cambridge_MA_early", "Cambridge_MA_Mid", "Cambridge_MA_Late"))

plot_tree(acalymma.timeCourse.tree, 
          nodelabf=nodeplotblank, 
          ladderize="left", 
          label.tips="taxa_names", 
          shape="Phylum", 
          method = "sampledodge",
          color="State",
          size="Abundance",
          base.spacing=0.1,
          justify = TRUE,
          sizebase = 10,
          title=title)  +
  ggtitle("Distribution of ### most abundant OTUs in Acalymma and Peponapis") +
  scale_color_manual(values=c(colors.time.course)) +
  scale_size_continuous(range = c(2, 8)) +
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=16),
        legend.title=element_blank())

ggsave("beetle.timeCourse.pdf", height = 10, width = 10)

## Plot the tree without aligning labels to get the OTU labels to print
plot_tree(acalymma.peponapis.tree, 
          nodelabf=nodeplotblank, 
          ladderize="left", 
          label.tips="taxa_names")

ggsave("bee.beetle.tree.node.labels.pdf", height = 10, width = 10)

dev.off()

-------------------
