##################
## Erwinia OTUS ##
##################

## The only three otus that are probably Et are: denovo577182, denovo590454, denovo634643
insect_b <- insect2 # phyloseq object unrarefied; filtered out wolbachia and chloroplasts
Et <- c("denovo577182", "denovo634643") # The two OTUs that may be Et
insect2 <- insect

insect.Et = prune_taxa(Et, insect_b)
insect.Et.merged = merge_taxa(insect.Et, taxa_names(insect.Et)[1:2])

#########################################################
### Detail of Et dynamics in vectors in Massachusetts ###
### throughout a single growing season ##################
#########################################################


prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(TopNOTUs, .) %>%
  tax_glom(., taxrank="Family") %>% 
  transform_sample_counts(., function(x) 100 * x/sum(x))

insect_b <- insect2 %>%
  subset_samples(Time_Course=="Yes") %>% ## Subset to beetles from the time course in Cambridge
  prune_samples(sample_sums(.) > 100, .) ## Restrict to samples that have at least 100 reads

insect.Et = prune_taxa(Et, insect_b)
insect.Et.merged = merge_taxa(insect.Et, taxa_names(insect.Et)[1:2])

sample_data(insect.Et.merged)$State2 <- factor(sample_data(insect_b)$State2, levels = c("Massachusetts-early", "Massachusetts-mid", "Massachusetts-late"))
levels(sample_data(insect.Et.merged)$State2)
#insect.Et.merged <- merge_samples(insect.Et.merged, "State2") ### Merge samples by category of state 

p2 <- plot_bar(insect.Et.merged, x="State2", fill="State2") + 
  ggtitle("Abundance of Erwinia tracheiphila:\nTotal unrarefied reads per sampling period") +
  coord_flip() +
  ylab("Abundance of Erwinia tracheiphila:\nTotal unrarefied reads per sampling period") +
 # xlab("Sampling Period") +
  scale_fill_manual(values=colors.time.course) 
p2 + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=14),
           axis.text.x  = element_text(face="bold", angle=360, color = "black", hjust=0.5, size=14),
           axis.title.y = element_text(face="bold", color = "white", size=12, hjust=0.5),
           axis.title.x = element_text(face="bold", color = "white", size=12),
           axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=12),
           legend.position="none")

ggsave("Erwinia.MA.Acalymma.pdf", height=3, width=6)

## All OTHER bacteria

Top.20.OTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(Top.20.OTUs, insect_b)
insect_b = merge_taxa(insect_b, Top.20.OTUs)

sample_data(insect_b)$State2 <- factor(sample_data(insect_b)$State2, levels = c("Massachusetts-early", "Massachusetts-mid", "Massachusetts-late"))
levels(sample_data(insect_b)$State2)

p2 <- plot_bar(insect_b, x="State2", fill="State2") + 
  ggtitle("Abundance of ALL bacteria in Acalymma:\nTotal unrarefied reads per sampling period") +
  coord_flip() +
  ylab("Abundance of all bacteria:\nTotal unrarefied reads per sampling period") +
  # xlab("Sampling Period") +
  scale_fill_manual(values=colors_vec) 
p2 + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
           axis.text.x  = element_text(face="bold", angle=360, color = "black", hjust=0.5, size=14),
           axis.title.y = element_text(face="bold", color = "white", size=12, hjust=0.5),
           axis.title.x = element_text(face="bold", color = "white", size=14),
           axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=12),
           legend.position="none")

ggsave("All.bacteria.MA.Acalymma.pdf", height=3, width=6)


## Count number of Et OTU reads --------------------------------
Top.20.OTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(Top.20.OTUs, insect_b)

table(sample_data(insect_b)$State) ## Top 20 OTUs
table(sample_data(insect.Et)$State) ## Just Et OTUs

insect_b_MA <- subset_samples(insect_b, State=="Cambridge_MA_Late")
print(insect_b_MA.abundance <- sample_sums(insect_b_MA))

##################################
## Detail of Et in ALL Acalymma ##
##################################

insect_b <- insect2 %>%
  subset_samples(
      Time_Course=="No" &
      Genus=="Acalymma"
    ) %>% 
  prune_samples(sample_sums(.) > 1, .) 

insect.Et = prune_taxa(Et, insect_b)
insect.Et.merged = merge_taxa(insect.Et, taxa_names(insect.Et)[1:2])

## Count number of Et OTU reads --------------------------------
Top.20.OTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(Top.20.OTUs, insect_b)

table(sample_data(insect_b)$State) ## Top 20 OTUs
table(sample_data(insect.Et)$State) ## Just Et OTUs

insect_b_MA <- subset_samples(insect.Et, State=="Tucson_AZ")
print(insect_b_MA.abundance <- sample_sums(insect_b_MA))

######## Bar charts ---------------------------

p2 <- plot_bar(insect.Et.merged, x="State2", fill="State2") + 
  ggtitle("Abundance of Erwinia tracheiphila:\nTotal unrarefied reads per sampling period") +
  coord_flip() +
  ylab("Abundance of Erwinia tracheiphila:\nTotal unrarefied reads per state") +
  # xlab("Sampling Period") +
  scale_fill_manual(values=colors_vec) 
p2 + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
           axis.text.x  = element_text(face="bold", angle=360, color = "black", hjust=0.5, size=14),
           axis.title.y = element_text(face="bold", color = "white", size=12, hjust=0.5),
           axis.title.x = element_text(face="bold", color = "white", size=14),
           axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=12),
           legend.position="none")

ggsave("Erwinia.all.Acalymma.pdf", height=3, width=6)

## All OTHER bacteria

TopNOTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect_b = prune_taxa(TopNOTUs, insect_b)
insect_b = merge_taxa(insect_b, TopNOTUs)

#sample_data(insect_b)$State2 <- factor(sample_data(insect_b)$State2, levels = c("Massachusetts-early", "Massachusetts-mid", "Massachusetts-late"))
#levels(sample_data(insect_b)$State2)
#insect.Et.merged <- merge_samples(insect.Et.merged, "State2") ### Merge samples by category of state 

p2 <- plot_bar(insect_b, x="State2", fill="State2") + 
  ggtitle("Abundance of ALL bacteria in Acalymma:\nTotal unrarefied reads per state") +
  coord_flip() +
  ylab("Abundance of all bacteria:\nTotal unrarefied reads per state") +
  # xlab("Sampling Period") +
  scale_fill_manual(values=colors_vec) 
p2 + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
           axis.text.x  = element_text(face="bold", angle=360, color = "black", hjust=0.5, size=14),
           axis.title.y = element_text(face="bold", color = "white", size=12, hjust=0.5),
           axis.title.x = element_text(face="bold", color = "white", size=14),
           axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=12),
           legend.position="none")

ggsave("All.bacteria.all.Acalymma.pdf", height=3, width=8)















#################################
## Graphs with unrarefied data ##
#################################

## Overall summary of all bacteria in Massachusetts Acalymma
insect_b <- insect2 %>%
  subset_samples(
    Genus=="Acalymma" &
      Time_Course=="Yes"
  ) %>%
  prune_samples(sample_sums(.) > 100, .) ## Restrict to samples that have at least 100 reads


Top.20.erwinia.OTUs = names(sort(taxa_sums(insect_b), TRUE)[1:20])
insect.20 = prune_taxa(Top.50.erwinia.OTUs, insect_b)

### Merge samples by category of state 
insectRm.State2 = merge_samples(insect.20, "State2")

# Repair the merged values associated with each surface after merge.
sample_data(insectRm.State2)$sample_Species <- levels(sample_data(insect_b)$sample_Species)
insectRm.Species.Family <- tax_glom(insectRm.State2, taxrank="Family")

# Transform to percentages of total available.
#insectRm.prop.Family = transform_sample_counts(insectRm.Species.Family, function(x) 100 * x/sum(x))

# 6x3 inches plot saved to pdf
title = "Community composition proportion in Acalymma in Massachusetts over a season (unrarefied)"
Striped.State <-plot_bar(insectRm.Species.Family, fill = "Family", title = title) + coord_flip() 
Striped.State + 
  theme(plot.title = element_text(lineheight=.8, face="bold",hjust=0.5))




###############################
## Graphs with rarefied data ##
###############################

insect_b <-insectR # Unrarefied data
insect_b <- subset_samples(insect_b, Time_Course=="Yes") ## Subset to beetles from the time course in Cambridge
insect_b <- subset_taxa(insect_b, Genus=="g__Erwinia") ## Subset to Erwinia genus
insect_b <- prune_taxa(taxa_sums(insect_b) > 0, insect_b) ## Restrict to taxa that are present at least once

ntaxa(insect_b)
sample_sums(insect_b)  # How many Erwinia reads are in each sample?
sample_variables(insect_b)

## This is less Erwinia than with unrarified data

title="Erwinia spp. in Acalymma"
p2 <- plot_bar(insect_b, x="State2", fill="OTU", title=title)
p2
## Very closely related to a plant associated Erwinia! Present in Western species as well

insect_prop = transform_sample_counts(insect_b, function(x) 100 * x/sum(x))
title = "Normalized Counts"

title="Normalized counts of top 20 OTUs in Acalymma"
p2 <- plot_bar(insect_prop, x="State_detail", fill="OTU", title=title)
p2


