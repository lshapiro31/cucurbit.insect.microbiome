
#################
### Peponapis ###
#################

insect_b <- insectR
insect_b <- subset_samples(insect_b, Genus=="Peponapis")
insect_b <- subset_samples(insect_b, State2!="Guanajuato")

## Alpha richness Peponapis
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
title="Geographic variation in Peponapis pruinosa alpha diversity"
(p <- plot_richness(insect_b, "State2", measures=alpha_meas, title=title))
p + geom_boxplot(data=p$data, aes(x=State2, y=value, color=State2, fill=State2), alpha=0.7) +
  scale_fill_manual(values = c("orange4", "orange")) +
  scale_color_manual(values=c("black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank()) +
  theme(legend.position="none")

ggsave("Peponapis.alpha.pdf", height=5, width=10)

dev.off()

######################
### Acalymma - all ###
######################

insect_b <- insectR
insect_b <- subset_samples(insect_b, Specialist=="Yes")
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_samples(insect_b, Time_Course=="No")

display.brewer.all() 
colors_vec <- brewer.pal(7, name = 'BuG')
print(colors_vec)

## Alpha richness 
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
title="Geographic variation in Acalymma spp. alpha diversity"
(p <- plot_richness(insect_b, "State2", measures=alpha_meas, title=title))
p + geom_boxplot(data=p$data, aes(x=State2, y=value, color=State2, fill=State2), alpha=0.7) +
  scale_fill_manual(values=colors_vec) +
  scale_color_manual(values =c("black", "black", "black", "black", "black", "black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank()) +
  theme(legend.position="none")

ggsave("Acalymma.alpha.pdf", height=5, width=10)


dev.off()


#################################
### Acalymma - MA time course ###
#################################

insect_b <- insectR
insect_b <- subset_samples(insect_b, Genus=="Acalymma")
insect_b <- subset_samples(insect_b, Time_Course=="Yes")

display.brewer.all() 
colors_vec <- brewer.pal(3, name = 'Paired')
print(colors_vec)

###  Change order of samples
# https://github.com/joey711/phyloseq/issues/291
sample_data(insect_b)$State2
sample_data(insect_b)$State2 <- factor(sample_data(insect_b)$State2, levels = c("Massachusetts-early", "Massachusetts-mid", "Massachusetts-late"))

## Alpha richness 
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
title="Variation in Acalymma spp. alpha diversity\nacross a single season in Massachusetts"
(p <- plot_richness(insect_b, "State2", measures=alpha_meas, title=title))
p + geom_boxplot(data=p$data, aes(x=State2, y=value, color=State2, fill=State2), alpha=0.7) +
  scale_fill_manual(values=colors_vec) +
  scale_color_manual(values =c("black", "black", "black")) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank()) +
  theme(legend.position="none")

ggsave("Acalymma.timeCourse.alpha.pdf", height=5, width=10)


dev.off()

##############################
### Acalymma vs. Peponapis ###
##############################


insect_b <- insectR
insect_b <- subset_samples(insect_b, Specialist=="Yes")
insect_b <- subset_samples(insect_b, insect_Genus!="Peponapis_Maybe")

## Alpha richness Peponapis
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
title="Alpha diversity in Peponapis vs. Acalymma"
(p <- plot_richness(insect_b, "insect_Genus", measures=alpha_meas, title=title))
p + geom_boxplot(data=p$data, aes(x=insect_Genus, y=value, color=insect_Genus, fill=insect_Genus), alpha=0.7) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold", angle=90, color = "black", hjust=0.5, size=10),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title=element_blank()) +
  theme(legend.position="none")

ggsave("Peponapis_v_Acalymma.alpha.pdf", height=5, width=10)

dev.off()
