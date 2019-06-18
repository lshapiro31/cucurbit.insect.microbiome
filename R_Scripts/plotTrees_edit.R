############################################
### Visual comparison of shared and ########
### unique OTUs between bees and beetles ###
############################################

TopNOTUs.acalymma = names(sort(taxa_sums(acalymma), TRUE)[1:50])
TopNOTUs.peponapis = names(sort(taxa_sums(peponapis), TRUE)[1:40])
Top.60.OTUs <- c(TopNOTUs.acalymma, TopNOTUs.peponapis)

acalymma.peponapis.tree <- acalymma.peponapis %>% 
 # subset_samples(Time_Course=="No") %>%
  prune_taxa(Top.60.OTUs, .) 

###  Change order of samples

sample_data(acalymma.peponapis.tree)$State
sample_data(acalymma.peponapis.tree)$State <- factor(sample_data(acalymma.peponapis.tree)$State, levels = c("CA", "Dively_Farm", "UVM_Vermont", "Mobile_AL", "Kentucky", "Iowa", "New_Mexico", "Tucson_AZ", "Guanajuato", "Cambridge_MA_early", "Cambridge_MA_Mid", "Cambridge_MA_Late"))

#sample_data(acalymma.peponapis.tree)$State2
#sample_data(acalymma.peponapis.tree)$State2 <- factor(sample_data(acalymma.peponapis.tree)$State2, levels = c("California", "Pennsylvania", "Vermont", "Alabama", "Kentucky", "Iowa", "New_Mexico", "Arizona", "Guanajuato"))

plot_tree(acalymma.peponapis.tree, 
          nodelabf=nodeplotblank, 
          ladderize="left", 
          label.tips="taxa_names", 
 #         shape="Phylum", 
          method = "sampledodge",
          color="State",
          size="Abundance",
          base.spacing=0.1,
          justify = TRUE,
          sizebase = 10)  +
#  ggtitle("Comparison of most abundant OTUs in Acalymma and Peponapis") +
  scale_color_manual(values=c(tree_colors)) +
  scale_size_continuous(range = c(2, 8)) +
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5, size=16),
        legend.title=element_blank())
 

ggsave("Output/bee.beetle.tree.pdf", height = 10, width = 10)

#### ----------------------------------------------------------------------------------

## Plot the tree without aligning labels to get the OTU labels to print
plot_tree(acalymma.peponapis.tree, 
          nodelabf=nodeplotblank, 
          ladderize="left", 
          label.tips="taxa_names")

ggsave("bee.beetle.tree.node.labels.pdf", height = 10, width = 10)

dev.off()

#### ----------------------------------------------------------------------------------