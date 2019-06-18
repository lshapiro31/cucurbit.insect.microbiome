################# 
### Peponapis ###
#################

peponapis <- subset_samples(acalymma.peponapis, Genus=="Peponapis")
TopNOTUs = names(sort(taxa_sums(acalymma.peponapis), TRUE)[1:50]) # Top 50 overall
#TopNOTUs = names(sort(taxa_sums(peponapis), TRUE)[1:30]) # Top 50 in Peponapis

peponapis <- acalymma.peponapis %>%
  subset_samples(
      Genus=="Peponapis" & State2!="Guanajuato"
                 ) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(TopNOTUs, .) %>%
  tax_glom(., taxrank="Order") %>% 
  transform_sample_counts(., function(x) 100 * x/sum(x))
  
peponapis.bar <- plot_bar(peponapis, x="State_detail", fill="Order") + 
  ggtitle("Order taxonomy of top ### OTUS in in Peponapis (unrarefied)") +
  ylab("Proportion of Reads per Taxonomic Group") +
  xlab("Individual Bee") +
  coord_flip()
peponapis.bar + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
      axis.text.x  = element_text(face="bold", color = "black", hjust=0.5, size=10),
      axis.title.y = element_text(face="bold", color = "black", size=12),
      axis.title.x = element_text(face="bold", color = "black", size=12),
      axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
      legend.title = element_text(size=12, face="bold", hjust = 0.5),
      legend.text = element_text(size = 10, face = "bold"))

ggsave("Peponapis.proportional.barchart.pdf", height = 6, width = 8)

dev.off()

################################ 
## Acalymma - non Time Course ##
################################

## Subset top 30 otus in Acalymma

acalymma <- subset_samples(acalymma.peponapis, Genus=="Acalymma" & Time_Course=="No")
#TopNOTUs = names(sort(taxa_sums(acalymma), TRUE)[1:30])
TopNOTUs = names(sort(taxa_sums(acalymma.peponapis), TRUE)[1:50]) # Top 50 overall

acalymma.non.timecourse <- acalymma %>% 
  subset_samples(
    Genus=="Acalymma" &
    Time_Course=="No"
    ) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(TopNOTUs, .) %>%
  tax_glom(., taxrank="Order") %>% 
  transform_sample_counts(., function(x) 100 * x/sum(x))

acalymma.non.timecourse.bar <- plot_bar(acalymma.non.timecourse, x="State_detail", fill="Order") + 
  ggtitle("Order breakdown of top ### OTUS\nin Acalymma spp. (unrarefied)") +
  ylab("Proportion of Reads per Taxonomic Group") +
  xlab("Individual Striped Beetle") +
  coord_flip()
acalymma.non.timecourse.bar + theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
          axis.text.x  = element_text(face="bold", color = "black", hjust=0.5, size=10),
          axis.title.y = element_text(face="bold", color = "black", size=12),
          axis.title.x = element_text(face="bold", color = "black", size=12),
          axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
          legend.title = element_text(size=12, face="bold", hjust = 0.5),
          legend.text = element_text(size = 10, face = "bold"))

ggsave("Acalymma.proportional.barchart.pdf", height = 8, width = 12)
dev.off()

############################ 
## Acalymma - Time Course ##
############################

acalymma.timeCourse <- subset_samples(acalymma.peponapis, Genus=="Acalymma" & Time_Course == "Yes")
# TopNOTUs = names(sort(taxa_sums(acalymma), TRUE)[1:30]) # Top 30 in Acalymma

TopNOTUs = names(sort(taxa_sums(acalymma.peponapis), TRUE)[1:50]) # Top 50 overall

acalymma.timeCourse <- acalymma.peponapis %>% 
  subset_samples(
    Genus == "Acalymma" &
    Time_Course == "Yes"
  ) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(TopNOTUs, .) %>%
  tax_glom(., taxrank="Order") %>% 
  transform_sample_counts(., function(x) 100 * x/sum(x))

###  Change order of samples
# https://github.com/joey711/phyloseq/issues/291
sample_data(acalymma.timeCourse)$State_detail
sample_data(acalymma.timeCourse)$State_detail <- factor(sample_data(acalymma.timeCourse)$State_detail, levels = c("Cambridge_MA_Late-7",
                                                                                                                  "Cambridge_MA_Late-6",
                                                                                                                  "Cambridge_MA_Late-5",
                                                                                                                  "Cambridge_MA_Late-4",
                                                                                                                  "Cambridge_MA_Late-3",
                                                                                                                  "Cambridge_MA_Late-1",
                                                                                                                  "Cambridge_MA_Mid-1",
                                                                                                                  "Cambridge_MA_Mid-2",
                                                                                                                  "Cambridge_MA_Mid-3",
                                                                                                                  "Cambridge_MA_Mid-4",
                                                                                                                  "Cambridge_MA_Mid-5",
                                                                                                                  "Cambridge_MA_Mid-6",
                                                                                                                  "Cambridge_MA_early-1", 
                                                                                                                  "Cambridge_MA_early-2", 
                                                                                                                  "Cambridge_MA_early-3",
                                                                                                                  "Cambridge_MA_early-4",
                                                                                                                  "Cambridge_MA_early-5",
                                                                                                                  "Cambridge_MA_early-6",
                                                                                                                  "Cambridge_MA_early-7"))

acalymma.timeCourse.bar <- plot_bar(acalymma.timeCourse, x="State_detail", fill="Order") + 
  ggtitle("Order breakdown of top ### OTUS\nin Acalymma spp. in Massachussetts (unrarefied)") +
  ylab("Proportion of Reads per Taxonomic Group") +
  xlab("Individual Striped Beetle") +
  coord_flip() + 
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
                     axis.text.x  = element_text(face="bold", angle=0, color = "black", hjust=0.5, size=10),
                     axis.title.y = element_text(face="bold", color = "black", size=12),
                     axis.title.x = element_text(face="bold", color = "black", size=12),
                     axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
                     legend.title = element_text(size=12, face="bold", hjust = 0.5),
                     legend.text = element_text(size = 10, face = "bold"))
acalymma.timeCourse.bar

ggsave("Acalymma.timeCourse.proportional.barchart.pdf", height = 8, width = 12)

dev.off()

