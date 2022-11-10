#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Figures                                                                       ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
### Last update: Nov 2022                                                     ###
################################################################################

# packages
# packages
library(MASS)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(HDInterval)
library(scales)
library(e1071)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(boot)
library(patchwork)
library(gridExtra)
library(ggpubr)
library(grid)
#install.packages("patchwork")
#library(Rmisc)


# Exploratory presence absence --------------------------------------------


# Exploratory abundance ---------------------------------------------------

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
