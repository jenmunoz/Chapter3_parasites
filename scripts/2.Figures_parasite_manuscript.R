#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Figures                                                                       ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
### Last update: Nov 2022                                                     ###
################################################################################

# packages
install.packages(c("MASS","ggplot2","readr","dplyr","HDInterval","scales","e1071","tidyverse"))
install.packages(c("boot","patchwork","ggpubr", "grid", "gridExtra"))


# packages
library(MASS)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(HDInterval)
library(scales)
library(e1071)
library(tidyverse)
library(boot)
library(patchwork)
library(gridExtra)
library(ggpubr)
library(grid)
#install.packages("patchwork")
#library(Rmisc)


# Exploratory presence absence --------------------------------------------

ectos_df<-read.csv("data/5.ectos_pres_abs_df.csv") # data on presence absence
ectos_df<-ectos_df%>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations

ectos_df$sociality<-as.factor(ectos_df$sociality)

str(ectos_df)

ggplot(ectos_df, aes(x = sociality, y=ectoparasites_PA_max, group=sociality, color=sociality)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.02)+
  labs(title="a)P/A All ectoparasites",  y="Ocurrence P/A", x = "sociality")+
  ylim(0,1.5)+
  theme_classic(20)

ggsave("figures/1.P_A_all_parasites_sociality.pdf", width=10, height=10, units="in")


ggplot(ectos_df, aes(x = sociality, y=sample_size, group=sociality, color=sociality)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.02)+
  labs(title="a)P/A All ectoparasites",  y="Sample size", x = "sociality")+
  ylim(0,100)+
  theme_classic(20)

ggsave("figures/1.P_A_all_parasites_sample_size.pdf",
       width=15, height=15, units="in")

ggplot(ectos_df, aes(x = elevation_cat, y=ectoparasites_PA_max, group=elevation_cat, color=elevation_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.02)+
  labs(title="a)P/A All ectoparasites",  y="Ocurrence P/A", x = "elevation")+
  ylim(0,2)+
  theme_classic(20)

ggsave("figures/1.P_A_all_parasites_elevation.pdf",
       width=15, height=15, units="in")

ggplot(ectos_df, aes(x = foraging_cat, y=ectoparasites_PA_max, group=foraging_cat, color=foraging_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.02)+
  labs(title="a)P/A All ectoparasites",  y="Ocurrence P/A", x = "foraging")+
  ylim(0,2)+
  theme_classic(20)


ggsave("figures/1.P_A_all_parasites_foraging.pdf",
       width=15, height=15, units="in")




# Exploratory abundance ---------------------------------------------------

lice_df_abundance<-read.csv("data/5.lice_df_abundance_manu.csv")

lice_df_abundance$sociality<-as.factor(lice_df_abundance$sociality)


#Lice
# Exploring data patterns on lice abundance

mean_by_sociality<- lice_df_abundance %>% 
  group_by(sociality) %>% 
  summarise(mean_group=mean(total_lice , na.rm=TRUE))


hist(lice_df_abundance$total_lice,  breaks = 50)

ggplot(lice_df_abundance, aes(x = total_lice)) +
  geom_histogram(fill = "white", colour ="black") +
  facet_grid (sociality~.)+
  theme_classic(20)


ggplot(lice_df_abundance, aes(x = sociality, y=total_lice, group=sociality, color=sociality)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="b)Lice abundance",  y="Total abundance per individual", x = "sociality")+
  ylim(0,100)+
  theme_classic(20)

ggsave("figures/1.Abundance_Lice_sociality_simple.pdf",
       width=15, height=15, units="in")

ggplot(lice_df_abundance, aes(x = elevation_cat, y=total_lice, group=elevation_cat, color=elevation_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="b)Lice abundance",  y="Total abundance per individual", x = "elevation")+
  ylim(0,100)+
  theme_classic(20)

ggsave("figures/1.Abundance_Lice_elevation.pdf",
       width=15, height=15, units="in")


ggplot(lice_df_abundance, aes(x = foraging_cat, y=total_lice, group=foraging_cat, color=foraging_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.02)+
  labs(title="a)P/A All ectoparasites",  y="Ocurrence P/A", x = "foraging")+
  ylim(0,80)+
  theme_classic(20)


# Map smoke to fill, make the bars NOT stacked, and make them semitransparent

plot1<-ggplot(lice_df_abundance, aes(x = total_lice, fill= sociality)) +
  geom_histogram(position = "identity", 
                 alpha = 0.5,
                 na.rm=TRUE,
                 bins = 30)+
  geom_vline(data=mean_by_sociality, aes(xintercept=mean_group, color=sociality),linetype="dashed",size=1.5)+
  scale_color_manual(values=c("#999999", "turquoise4"))+
  scale_fill_manual(values=c("#999999", "turquoise4"))+
  labs(x="Total lice (per individual)", y = "Count")+
  theme_classic(20)+
  theme(axis.text=element_text(color="black"))


plot2<-ggplot(lice_df_abundance, aes(y=(total_lice), x = sociality, fill= sociality, color=sociality, alpha=sociality))+
  geom_boxplot(width=0.8,alpha=0.3)+
  geom_point(size=1, alpha=0.6,position = position_jitter(width = .2))+
  stat_summary(fun=mean, geom="crossbar", linetype = "dashed",size=0.1)+
  #geom_hline(linetype = "dashed", yintercept = overall_mean,size=1.5) + 
  geom_hline(data=mean_by_sociality, aes(yintercept=mean_group, color=sociality),linetype="dashed",size=1.5)+
  coord_flip()+
  scale_fill_manual(values=c("#999999", "turquoise4"))+
  scale_color_manual(values=c("#999999", "turquoise4"))+ ##E69F00
  scale_alpha_manual(values=c(0.3,0.4))+
  labs(title="a) Lice abundance (per individual)",y="Total lice", x = "sociality")+
  theme_classic()+
  theme_classic(20)

plot2<- plot2 + clean_theme()

figure<-ggarrange(plot2 + rremove("xlab"), plot1  + rremove("xlab"),
                  ncol=1, nrow=2,
                  align = "hv", 
                  common.legend = TRUE,legend = "right",
                  widths = c(5, 5), heights = c(3, 3),
                  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

lice_abundance_figure<-annotate_figure(figure,bottom = textGrob("Total lice", gp = gpar(cex = 1.3)))

ggsave("figures/1.Abundance_Lice_sociality.pdf",
       width=15, height=15, units="in")

# Mites

ggplot(mites_df_abundance, aes(x = sociality, y=total_mesostigmatidae, group=sociality, color=sociality)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="b)Mites abundance (Mesostigmatidae)",  y="Total abundance per individual", x = "sociality")+
  ylim(0,500)+
  theme_classic(20)

ggsave("figures/1.Abundance_mites_sociality_Mesostigmatidae.pdf",
       width=15, height=15, units="in")

ggplot(mites_df_abundance, aes(x = elevation_cat, y=total_mesostigmatidae, group=elevation_cat, color=elevation_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="b)mites abundance (Mesostigmatidae) ",  y="Total abundance per individual", x = "elevation")+
  ylim(0,100)+
  theme_classic(20)

ggsave("figures/1.Abundance_Lice_elevation_Mesostigmatidae.pdf",
       width=15, height=15, units="in")


ggplot(mites_df_abundance, aes(x = foraging_cat, y=total_mesostigmatidae, group=foraging_cat, color=foraging_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.02)+
  labs(title="a)mites abundance (Mesostigmatidae)",  y="total abundance per individual", x = "foraging")+
  ylim(0,80)+
  theme_classic(20)

ggsave("figures/1.Abundance_Lice_foraging_Mesostigmatidae.pdf",
       width=15, height=15, units="in")

# Mites means mesostigmatida

ggplot(mean_mites_abundance_mesostigmatidae, aes(x = sociality, y=mean_mites, group=sociality, color=sociality)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="b)Mites abundance (mean/ species) Mesostigmatidae",  y="Total abundance per species (mean)", x = "sociality")+
  ylim(0,75)+
  theme_classic(20)

ggsave("figures/1.Abundance_mites_sociality_Mesostigmatidae_means.pdf",
       width=15, height=15, units="in")


# Exploratory diversity ---------------------------------------------------


### lICE~Genera


### Exploratory plots
ggplot(lice_df_richness_summary, aes(x = richness)) +
  geom_histogram(fill = "white", colour ="black") +
  facet_grid (sociality~.)+
  theme_classic(20)

ggplot(lice_df_richness_genera, aes(x = sociality, y=richness, group=sociality, color=sociality)) +
  geom_boxplot() +
  labs(title="a) Lice diversity (total genera)",y="Total genera", x = "sociality")+
  theme_classic(20)

View(lice_df_richness_genera)
ggsave("figures/exploratory_lice_diversity_genera.pdf",
       width=15, height=15, units="in")

# Lice~species
###
lice_richness_manu_sp<-read.csv("data/5.lice_richness_sp_df_manu.csv")
lice_richness_manu_sp<-lice_richness_manu_sp %>% distinct( species_jetz,.keep_all = TRUE)

lice_richness_manu_sp$sociality<-as.factor(lice_richness_manu_sp$sociality)



ggplot(lice_richness_manu_sp, aes(x = sociality, y=richness_sp, group=sociality, color=sociality)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="b)Lice richness (sp)",  y="Total ectos morpho-species by bird species", x = "sociality")+
  ylim(0,10)+
  theme_classic(20)

ggsave("figures/1.Richness_lice_sociality_sp.pdf",
       width=15, height=15, units="in")


ggplot(lice_richness_manu_sp, aes(x = sociality, y=n_samples_lice, group=sociality, color=sociality)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="a)Richness_lice_sp",  y="Sample size", x = "sociality")+
  ylim(0,25)+
  theme_classic(20)

ggsave("figures/1.Richness_lice_sample_size.pdf",
       width=15, height=15, units="in")

ggplot(lice_richness_manu_sp, aes(x = elevation_cat, y=richness_sp, group=elevation_cat, color=elevation_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.01)+
  labs(title="a)Richness_lice",  y="Total ectos morpho-species by bird species", x = "elevation")+
  ylim(0,10)+
  theme_classic(20)

ggsave("figures/1.Richness_lice_elevation_sp.pdf",
       width=15, height=15, units="in")

ggplot(lice_richness_manu_sp, aes(x = foraging_cat, y=richness_sp, group=foraging_cat, color=foraging_cat)) +
  geom_point(alpha=0.5)+
  geom_boxplot(alpha = 0.10) +
  geom_jitter(height = 0.02)+
  labs(title="a)P/A All ectoparasites",  y="Total ectos morpho-species by bird specie", x = "foraging")+
  ylim(0,10)+
  theme_classic(20)


ggsave("figures/1.Richness_lice_foraging_sp.pdf",
       width=15, height=15, units="in")








ggplot(lice_df_richness_sp_summary, aes(x = richness_sp)) +
  geom_histogram(fill = "white", colour ="black") +
  facet_grid (sociality~.)+
  theme_classic(20)



# Correlation between sample size and diversity



# # Models prevalence -----------------------------------------------------


# Models occurrence -------------------------------------------------------


# Model abundance ---------------------------------------------------------

# PLOTTING TRANDOM EFFECTS an example 
plot.communityPGLMM(m_a_meso,predicted=TRUE)

a<-communityPGLMM.plot.re(x=m_a_meso, sp.var = "species_jetz", site.var = "total_mesostigmatidae",predicted=TRUE)

communityPGLMM.plot.re()
