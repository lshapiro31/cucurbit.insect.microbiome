#################
### Peponapis ###
#################

peponapis.R <- acalymma.peponapis.R %>%  ## Using rarefied data for diversity analysis
  subset_samples(
      Genus=="Peponapis" &
      State2!="Guanajuato"
    )

peponapis.alpha <- plot_richness(peponapis.R, "State2", measures=c("Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot(aes(x=State2, y=value, color=State2, fill=State2), alpha=0.7) +
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
peponapis.alpha

ggsave("Peponapis.alpha.pdf", height=5, width=10)

dev.off()

## Test for statistical differences in alpha diversity 
peponapis.alpha.tbl <- peponapis.alpha$data %>% as_tibble(.)

peponapis.alpha.test.stat <- peponapis.alpha$data %>% 
  as_tibble(.) %>%
  filter(variable == "Simpson") %>%   ## Filter by results of "Shannon", "Simpson", "InvSimpson"
  select(State2, variable, value)

peponapis.alpha.test.stat

peponapis.alpha.test.stat$State2 <- as.factor(peponapis.alpha.test.stat$State2)
kruskal.test(value ~ State2, peponapis.alpha.test.stat)

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
sample_data(acalymma.R)$State2
sample_data(acalymma.R)$State2 <- factor(sample_data(acalymma.R)$State2, levels = c("California", "Pennsylvania", "Vermont", "Alabama", "Kentucky", "Iowa", "New_Mexico", "Arizona", "Guanajuato"))

acalymma.alpha <- plot_richness(acalymma.R, "State2", measures=c("Shannon", "Simpson", "InvSimpson"), title=title) +
  geom_boxplot(aes(x=State2, y=value, color=State2, fill=State2), alpha=0.7) +
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
acalymma.alpha

ggsave("Acalymma.alpha.pdf", height=5, width=10)

dev.off()

## Test for statistical differences in alpha diversity 
acalymma.alpha.tbl <- acalymma.alpha$data %>% as_tibble(.)

acalymma.alpha.test.stat <- acalymma.alpha$data %>% 
  as_tibble(.) %>%
  filter(variable == "Simpson") %>%   ## Filter by results of "Shannon", "Simpson", "InvSimpson"
  select(State2, Species, variable, value, Region)
acalymma.alpha.test.stat

acalymma.alpha.test.stat$Region <- as.factor(acalymma.alpha.test.stat$Region)
acalymma.alpha.test.stat$State2 <- as.factor(acalymma.alpha.test.stat$State2)
acalymma.alpha.test.stat$Species <- as.factor(acalymma.alpha.test.stat$Species)

kruskal.test(value ~ Region, acalymma.alpha.test.stat)
kruskal.test(value ~ State2, acalymma.alpha.test.stat)
kruskal.test(value ~ Species, acalymma.alpha.test.stat)

pairwise.wilcox.test(acalymma.alpha.test.stat$value, acalymma.alpha.test.stat$Region, p.adjust.method="fdr")
pairwise.wilcox.test(acalymma.alpha.test.stat$value, acalymma.alpha.test.stat$State2, p.adjust.method="fdr")

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

acalymma.time.course.alpha <- plot_richness(acalymma.timecourse.R, "State", measures=c("Shannon", "Simpson", "InvSimpson"), title=title) + 
  scale_color_manual(values =c("black", "black", "black")) +
  scale_fill_manual(values=colors.time.course) +
  ggtitle("Variation in Acalymma spp. alpha diversity\nacross a single season in Massachusetts") +
  geom_boxplot(aes(x=State, y=value, color=State, fill=State), alpha=0.7) +
  theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size=16),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.y  = element_text(face="bold", color = "black", vjust=0.5, size=10),
        legend.title = element_blank(),
        legend.text = element_text(face="bold", color = "black", size=12),
        legend.position="right")
acalymma.time.course.alpha

ggsave("Acalymma.timeCourse.alpha.pdf", height=5, width=10)

dev.off()

## Test for statistical differences in alpha diversity 
acalymma.time.course.alpha.tbl <- acalymma.time.course.alpha$data %>% as_tibble(.)

acalymma.time.course.alpha.stat <- acalymma.time.course.alpha$data %>% 
  as_tibble(.) %>%
  filter(variable == "Simpson") %>%   ## Filter by results of "Shannon", "Simpson", "InvSimpson"
  select(State, variable, value)
acalymma.time.course.alpha.stat

acalymma.time.course.alpha.stat$State <- as.factor(acalymma.time.course.alpha.stat$State)
kruskal.test(value ~ State, acalymma.time.course.alpha.stat)
pairwise.wilcox.test(acalymma.time.course.alpha.stat$value, acalymma.time.course.alpha.stat$State, p.adjust.method="fdr")

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
