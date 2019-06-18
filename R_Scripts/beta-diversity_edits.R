#####################################
### Acalymma - Non-MA time course ###
#####################################

acalymma.non.timecourse <- acalymma.peponapis %>% 
  subset_samples(
      Specialist=="Yes" &
      Genus=="Acalymma" & 
      Time_Course=="No"
  )

TopNOTUs = names(sort(taxa_sums(acalymma.non.timecourse), TRUE)[1:30])
acalymma.non.timecourse = prune_taxa(TopNOTUs, acalymma.non.timecourse)
acalymma.non.timecourse = transform_sample_counts(acalymma.non.timecourse, function(x) x / sum(x))

sample_data(acalymma.non.timecourse)$State2
sample_data(acalymma.non.timecourse)$State2 <- factor(sample_data(acalymma.non.timecourse)$State2, levels = c("California", "Pennsylvania", "Vermont", "Alabama", "Kentucky", "Iowa", "New_Mexico", "Arizona", "Guanajuato"))

# Calculate distances and ordinate 
acalymma_pcoa_bray = ordinate(acalymma.non.timecourse, method = "PCoA", distance = "bray") # Bray-curtis
acalymma_pcoa_wUF = ordinate(acalymma.non.timecourse, method = "PCoA", distance = "wUniFrac") # Weighted unifrac
acalymma_pcoa_UF = ordinate(acalymma.non.timecourse, method = "PCoA", distance = "UniFrac") # Original (unweighted unifrac)

# Plot scree
# plot_scree(acalymma_pcoa_UF, "Scree Plot: unWeighted UniFrac MDS")

# Plot acalymma_pcoa_bray, acalymma_pcoa_wUF, acalymma_pcoa_UF
acalymma.non.timecourse.ordination <- plot_ordination(acalymma.non.timecourse, acalymma_pcoa_bray, color = "State2") +
  geom_point(color = "grey40", size = 8) + 
  geom_point(aes(color = State2), size = 6) +
#  ggtitle("PCoA of Bray-Curtis distance\nof geographically distinct Acalymma populations") +
#  ggtitle("PCoA of Unifrac (weighted) distance\nof geographically distinct Acalymma populations") +
  ggtitle("PCoA of Unifrac (unweighted) distance\nof geographically distinct Acalymma populations") +
    scale_color_manual(values=states.7) +
#  scale_fill_manual(values =c("black", "black", "black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_text(face="bold", color = "black", size=12),
        axis.text.x  = element_text(face="bold", angle=360, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank())
acalymma.non.timecourse.ordination

#ggsave("Output/Acalymma.all.PcOA.BrayCurtis.pdf", height=9, width=11)
#ggsave("Output/Acalymma.all.PcOA.wUniFrac.pdf", height=9, width=11)
ggsave("Output/Acalymma.all.PcOA.UniFrac.pdf", height=9, width=11)

dev.off()

#################################
### Acalymma - MA time course ###
#################################

acalymma.timecourse <- acalymma.peponapis %>% 
  subset_samples(Time_Course=="Yes")

TopNOTUs = names(sort(taxa_sums(acalymma.timecourse), TRUE)[1:30])
acalymma.timecourse = prune_taxa(TopNOTUs, acalymma.timecourse)
acalymma.timecourse = transform_sample_counts(acalymma.timecourse, function(x) x / sum(x))

###  Change order of samples
# https://github.com/joey711/phyloseq/issues/291
sample_data(acalymma.timecourse)$State
sample_data(acalymma.timecourse)$State <- factor(sample_data(acalymma.timecourse)$State, levels = c("Cambridge_MA_early", "Cambridge_MA_Mid", "Cambridge_MA_Late"))

# Calculate distances and ordinate 
acalymma_season_pcoa_bray = ordinate(acalymma.timecourse, method = "PCoA", distance = "bray") # Bray-curtis
acalymma_season_pcoa_wUF = ordinate(acalymma.timecourse, method = "PCoA", distance = "wUniFrac") # Weighted unifrac
acalymma_season_pcoa_UF = ordinate(acalymma.timecourse, method = "PCoA", distance = "UniFrac") # Original (unweighted unifrac)

# Plot scree
plot_scree(acalymma_season_pcoa_bray, "Scree Plot: Bray-Curtis MDS")
plot_scree(acalymma_season_pcoa_wUF, "Scree Plot: Weighted UniFrac MDS")
plot_scree(acalymma_season_pcoa_UF, "Scree Plot: Weighted UniFrac MDS")

# Plot acalymma_season_pcoa_bray, acalymma_season_pcoa_wUF, acalymma_season_pcoa_UF
plot_ordination(acalymma.timecourse, acalymma_season_pcoa_UF, color = "State") +
#  stat_ellipse(type = "t", level=.9) +
  geom_point(color = "grey40", size = 8) + 
  geom_point(aes(color = State), size = 6) +
#  ggtitle("PCoA of Bray-Curtis distance\nof geographically distinct Acalymma populations") +
#  ggtitle("PCoA of Unifrac (weighted) distance\nof geographically distinct Acalymma populations") +
#  ggtitle("PCoA of Unifrac (unweighted) distance\nof geographically distinct Acalymma populations") +
  scale_color_manual(values=colors.time.course) +
  #  scale_fill_manual(values =c("black", "black", "black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_text(face="bold", color = "black", size=12),
        axis.text.x  = element_text(face="bold", angle=360, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank())

#ggsave("Output/Acalymma.timeCourse.PcOA.BrayCurtis.pdf", height=9, width=11)
#ggsave("Output/Acalymma.timeCourse.PcOA.wUniFrac.pdf", height=9, width=11)
ggsave("Output/Acalymma.timeCourse.PcOA.UniFrac.pdf", height=9, width=11)


dev.off()

#################
### Peponapis ###
#################

peponapis <- acalymma.peponapis %>%
  subset_samples(
      Genus=="Peponapis" &
      State2!="Guanajuato"
      ) %>%
  prune_taxa(taxa_sums(.) > 0, .) 

TopNOTUs = names(sort(taxa_sums(peponapis), TRUE)[1:30])
peponapis = prune_taxa(TopNOTUs, peponapis)
peponapis = transform_sample_counts(peponapis, function(x) 100 * x/sum(x))

# Calculate distances and ordinate 
peponapis_pcoa_bray = ordinate(peponapis, method = "PCoA", distance = "bray") # Bray-curtis
peponapis_pcoa_wUF = ordinate(peponapis, method = "PCoA", distance = "wUniFrac") # Weighted unifrac
peponapis_pcoa_UF = ordinate(peponapis, method = "PCoA", distance = "UniFrac") # Original (unweighted unifrac)

# Plot scree
plot_scree(peponapis_pcoa_bray, "Scree Plot: Bray-Curtis MDS")
plot_scree(peponapis_pcoa_wUF, "Scree Plot: Weighted UniFrac MDS")
plot_scree(peponapis_pcoa_UF, "Scree Plot: Weighted UniFrac MDS")

# Plot peponapis_pcoa_bray, peponapis_pcoa_wUF, peponapis_pcoa_UF
plot_ordination(peponapis, peponapis_pcoa_UF, color = "State2") +
  geom_point(color = "grey40", size = 8) + 
  geom_point(aes(color = State2), size = 6) +
#  ggtitle("PCoA of Bray-Curtis distance\nof geographically distinct Peponapis populations") +
#  ggtitle("PCoA of Unifrac (weighted) distance\nof geographically distinct Peponapis populations") +
  ggtitle("PCoA of Unifrac (unweighted) distance\nof geographically distinct Peponapis populations") +
  scale_color_manual(values = c(bees.2)) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank())

#ggsave("Output/Peponapis.BrayCurtis.pdf", height=5, width=7)
#ggsave("Output/Peponapis.PcOA.wUniFrac.pdf", height=5, width=7)
ggsave("Output/Peponapis.PcOA.UniFrac.pdf", height=5, width=7)


