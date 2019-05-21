#################
### Peponapis ###
#################

peponapis.R <- acalymma.peponapis.R %>%  ## Using rarefied data for diversity analysis
  subset_samples(
      Genus=="Peponapis" &
      State2!="Guanajuato"
    )

plot_richness(peponapis.R, "State2", measures=c("Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot(aes(x=State2, y=value, color=State2, fill=State2), alpha=0.9) +
  scale_fill_manual(values = c(bees.2)) +
  ggtitle("Geographic variation in Peponapis pruinosa alpha diversity") +
  scale_color_manual(values=c("black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title = element_blank(),
        legend.text = element_text(face="bold", color = "black", size=12),
        legend.position="right")

ggsave("Peponapis.alpha.pdf", height=5, width=10)

dev.off()

######################
### Acalymma - all ###
######################

acalymma.R <- acalymma.peponapis.R %>% 
  subset_samples(
    Specialist=="Yes" &
      Genus=="Acalymma" & 
      Time_Course=="No"
  )

## Alpha richness 
alpha_meas = c("Shannon", "Simpson")

sample_data(acalymma.R)$State2
sample_data(acalymma.R)$State2 <- factor(sample_data(acalymma.R)$State2, levels = c("California", "Pennsylvania", "Vermont", "Alabama", "Kentucky", "Iowa", "New_Mexico", "Arizona", "Guanajuato"))

plot_richness(acalymma.R, "State2", measures=c("Shannon", "Simpson", "InvSimpson"), title=title) +
  geom_boxplot(aes(x=State2, y=value, color=State2, fill=State2), alpha=0.9) +
  scale_fill_manual(values=states.7) +
  ggtitle("Geographic variation in Acalymma alpha diversity") +
  scale_color_manual(values =c("black", "black", "black", "black", "black", "black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank(),
        legend.text = element_text(face="bold", color = "black", size=12),
        legend.position="right")

ggsave("Acalymma.alpha.pdf", height=5, width=10)

dev.off()

#################################
### Acalymma - MA time course ###
#################################

acalymma.timecourse.R <- acalymma.peponapis.R %>%
  subset_samples(
    Genus == "Acalymma" & Time_Course == "Yes"
    )

###  Change order of samples
# https://github.com/joey711/phyloseq/issues/291
sample_data(acalymma.timecourse.R)$State
sample_data(acalymma.timecourse.R)$State <- factor(sample_data(acalymma.timecourse.R)$State, levels = c("Cambridge_MA_early", "Cambridge_MA_Mid", "Cambridge_MA_Late"))

alpha_meas = c("Shannon", "Simpson")
plot_richness(acalymma.timecourse.R, "State", measures=c("Shannon", "Simpson", "InvSimpson"), title=title) + 
  geom_boxplot(aes(x=State, y=value, color=State, fill=State), alpha=0.9) +
  scale_fill_manual(values=colors.time.course) +
  ggtitle("Variation in Acalymma spp. alpha diversity\nacross a single season in Massachusetts") +
  scale_color_manual(values =c("black", "black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title = element_blank(),
        legend.text = element_text(face="bold", color = "black", size=12),
        legend.position="right")

ggsave("Acalymma.timeCourse.alpha.pdf", height=5, width=10)

dev.off()

##############################
### Acalymma vs. Peponapis ###
##############################

plot_richness(acalymma.peponapis.R, "insect_Type", measures=c("Shannon", "Simpson", "InvSimpson"), title=title) +
  geom_boxplot(aes(x=insect_Type, y=value, color=insect_Type, fill=insect_Type), alpha=0.7) +
  ggtitle("Alpha diversity in Peponapis vs. Acalymma") +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank(),
        legend.position="none")

ggsave("Peponapis_v_Acalymma.alpha.pdf", height=5, width=10)

dev.off()

### ----------------------------------------------------------------------

data("soilrep")
plot_richness(soilrep, measures=c("InvSimpson", "Fisher"))
plot_richness(soilrep, "Treatment", "warmed", measures=c("Chao1", "ACE", "Simpson", "InvSimpson"))
data("GlobalPatterns")
plot_richness(GlobalPatterns, x="SampleType", measures=c("InvSimpson", "Simpson"))
plot_richness(GlobalPatterns, x="SampleType", measures=c("Chao1", "ACE", "Simpson", "InvSimpson"), nrow=3)
plot_richness(GlobalPatterns, x="SampleType", measures=c("Chao1", "ACE", "InvSimpson"), nrow=3, sortby = "Chao1")

