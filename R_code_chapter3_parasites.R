###########R-code for: Chapter 3 parasites############################################################
# This contains the code to explore an analyse data for chapter 3 
#
#
#last modeified Oct 8 2021
###
# Loading packages --------------------------------------------------------
# libraries for easier manipulation of data
install.packages("tidyr") 
install.packages("tidyverse") 
install.packages("dplyr")
install.packages ("data.table")
install.packages ("extrafont")
installed.packages("lubridate")  #for dates
install.packages("car")
#Other libraries for data analyses
install.packages("vegan")
install.packages("ggplot2")
install.packages("devtools")
install.packages("lme4")
install.packages("knitr")
install.packages("ts")


library(vegan)
library(tidyverse)
library(tidyr)
library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(extrafont)
library(visreg)
library(lubridate)
#library(ts)


# Read Data --------------------------------------------------------------------

#Selecting species of interest

parasitescountf<-read.csv("FlockerGenus_FlockerSpecies_TotalSamples.csv",stringsAsFactors = FALSE, header=TRUE, strip.white=TRUE, na.strings = c("NA",""))
taxonomy<-read.csv("taxonomy.csv", stringsAsFactors = FALSE, header=TRUE, strip.white=TRUE, na.strings = c("NA",""))
parasitescountf_taxonomy<-right_join(parasitescountf, taxonomy, by = "species")
write.csv (parasitescountf_taxonomy, "parasites_birds_flocks_manu.csv")

#Data from Allen's lab
species_flock_priority<-read.csv("mixed_flocks_project_priority_list_1_output.csv", header=TRUE,strip.white=TRUE, na.strings = c("NA",""))
species_flock_interest<-read.csv("mixed_flocks_project_priority_list_2_output.csv", header=TRUE, strip.white=TRUE, na.strings = c("NA",""))
sociality_classification<-read.csv("sociality_classification.csv",header=TRUE, na.strings = c("NA","")) 
# notes on sociality 1= obligate species, 2= regular participants of flcoks, 3#occasional participants of flocsk, 4# accidental participant of flocks 5# solitary species
#sociality_binomial 0= solitary or accidentally  found in flcoks 1= joing mixed species flocks at some extent
#sociality conspecific refers to species that foraga in conspecific groups ( it is independet of mixed species flcoks participation, a species can be part of mixed psecies flcoks and conspecific flocks too or ant-followers flocks) 2= ant-followers 1=yes family groups 0=no

str(sociality_classification)
str(species_flock_priority_sociality)
str(species_flock_priority)   

#merge

species_flock_priority_sociality<-full_join(species_flock_priority,sociality_classification, by="species")

species_flock_priority_sociality %>% mutate_at(vars(Lice_total, Mites1_total,Mites2_total,Ticks_Total ), list(as.integer))

#species_flock_priority_sociality %>%mutate(total_ectoparasites=(Lice_total+Mites1_total+Mites2_total+Ticks_Total)) #create a new column with teh total of ectoparasites
#species_flock_priority_sociality_n_cero<-species_flock_priority_sociality %>%filter(total_ticks_mites_lice != "0")
species_flock_priority_sociality_upperlayers<-species_flock_priority_sociality %>%filter(foraging!= "on-near ground")
#species_flock_priority_sociality_n_sclerurus<-species_flock_priority_sociality %>%filter(species!= "Sclerurus caudacutus")
#If getting NaS
#which (is.na(species_flock_priority_sociality$foraging))

View(species_flock_priority_sociality)

write.csv(species_flock_priority_sociality,"2.species_flock_priority_sociality.csv")

#new collumns for species analyese
species_flock_priority_sociality$species_factor<-as.factor(species_flock_priority_sociality$species)

species_flock_priority_sociality_upperlayers$species_factor<-as.factor(species_flock_priority_sociality_upperlayers$species)

species_flock_priority_sociality$sociality_binomial<-as.factor(species_flock_priority_sociality$sociality_binomial)
species_flock_priority_sociality_upperlayers$sociality_binomial<-as.factor(species_flock_priority_sociality_upperlayers$sociality_binomial)


#plots


#general

ggplot(data=species_flock_priority_sociality, aes(x=sociality, y=total_ticks_mites_lice)) +
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=15,size=2, color="black", fill="black")+
  theme_classic()


ggplot(data=species_flock_priority_sociality_upperlayers, aes(x=sociality,y=total_ticks_mites_lice))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun.y=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

species_flock_priority_sociality_upperlayers%>% group_by(sociality) %>% 
  summarize(mean=mean(total_ticks_mites_lice, na.rm=TRUE),sd=sd(total_ticks_mites_lice, na.rm=TRUE))

species_flock_priority_sociality%>% group_by(sociality) %>% 
  summarize(mean=mean(total_ticks_mites_lice, na.rm=TRUE),sd=sd(total_ticks_mites_lice, na.rm=TRUE))


# Social others_ socilaity binomial
ggplot(data=species_flock_priority_sociality, aes(x=sociality_binomial , y=total_ticks_mites_lice)) +
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()


ggplot(data=species_flock_priority_sociality_upperlayers, aes(x=sociality_binomial,y=total_ticks_mites_lice))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun.y=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()


species_flock_priority_sociality_upperlayers%>% group_by(sociality_binomial) %>% 
  summarize(mean=mean(total_ticks_mites_lice, na.rm=TRUE))

species_flock_priority_sociality%>% group_by(sociality_binomial) %>% 
  summarize(mean=mean(total_ticks_mites_lice, na.rm=TRUE), sd=sd(total_ticks_mites_lice, na.rm=TRUE))
  
# This is the anlyses for lice only

###Lice
ggplot(data=species_flock_priority_sociality, aes(x=sociality, y=Lice_total))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

ggplot(data=species_flock_priority_sociality_upperlayers, aes(x=sociality, y=Lice_total))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun.y=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

species_flock_priority_sociality_upperlayers%>% group_by(sociality) %>% 
  summarize(mean=mean(Lice_total, na.rm=TRUE),sd=sd(Lice_total, na.rm=TRUE))

species_flock_priority_sociality%>% group_by(sociality) %>% 
  summarize(mean=mean(Lice_total, na.rm=TRUE),sd=sd(Lice_total, na.rm=TRUE))

# Social others_ socilaity binomial#lice

ggplot(data=species_flock_priority_sociality, aes(x=sociality_binomial , y=Lice_total)) +
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

ggplot(data=species_flock_priority_sociality_upperlayers, aes(x=sociality_binomial , y=Lice_total)) +
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()
 

species_flock_priority_sociality_upperlayers%>% group_by(sociality_binomial) %>% 
  summarize(mean=mean(Lice_total, na.rm=TRUE))

species_flock_priority_sociality%>% group_by(sociality_binomial) %>% 
  summarize(mean=mean(Lice_total, na.rm=TRUE))

str(species_flock_priority_sociality)
##Mites
ggplot(data=species_flock_priority_sociality, aes(x=sociality, y=Totals_mites))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()
ggplot(data=species_flock_priority_sociality_upperlayers, aes(x=sociality, y=Totals_mites))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

species_flock_priority_sociality_upperlayers%>% group_by(sociality) %>% 
  summarize(mean=mean(Totals_mites, na.rm=TRUE),sd=sd(Totals_mites, na.rm=TRUE))

species_flock_priority_sociality%>% group_by(sociality) %>% 
  summarize(mean=mean(Totals_mites, na.rm=TRUE),sd=sd(Totals_mites, na.rm=TRUE))

#Mites binomial

ggplot(data=species_flock_priority_sociality, aes(x=sociality_binomial, y=Totals_mites))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()
ggplot(data=species_flock_priority_sociality_upperlayers, aes(x=sociality_binomial, y=Totals_mites))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  stat_summary(fun=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

species_flock_priority_sociality_upperlayers%>% group_by(sociality_binomial) %>% 
  summarize(mean=mean(Totals_mites, na.rm=TRUE),sd=sd(Totals_mites, na.rm=TRUE))

species_flock_priority_sociality%>% group_by(sociality_binomial) %>% 
  summarize(mean=mean(Totals_mites, na.rm=TRUE),sd=sd(Totals_mites, na.rm=TRUE))

#Ticks
ggplot(data=species_flock_priority_sociality, aes(x=sociality, y=Ticks_Total ))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  scale_color_manual(values=c("#00798c", "#edae49", "#66a182","#d1495b","#8d96a3"))+
  theme_classic()


#Analyze variation by species

ggplot(data=species_flock_priority_sociality, aes(x=sociality, y=total_ticks_mites_lice))+
  geom_jitter(aes(colour=species),width = 0.3)+
  stat_summary(fun.y=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

ggplot(data=species_flock_priority_sociality, aes(x=species, y=total_ticks_mites_lice))+
  geom_jitter(aes(colour=sociality),width = 0.3)+
  theme_classic()

#without sclerurus

ggplot(data=species_flock_priority_sociality_n_sclerurus, aes(x=sociality,y=total_ticks_mites_lice))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  stat_summary(fun.y=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()


####Lets plot by speceis

ggplot(data=species_flock_priority_sociality, aes(x=species_factor, y=total_ticks_mites_lice)) +
  geom_jitter(aes(colour=sociality),width = 0.3)+
  theme_classic()

ggplot(data=species_flock_priority_sociality_upperlayers, aes(x=species_factor, y=total_ticks_mites_lice)) +
  geom_jitter(aes(colour=sociality),width = 0.3)+
  theme_classic()


####
# Ectoparadite diversity --------------------------------------------------
####

matrix_diversity<-read.csv("3_simulated_matrix_diversity_ectos.csv")
as.matrix.data.frame(matrix_diversity)
#lets play with Myrsidea

species_flock_priority_sociality_myrsidea<-species_flock_priority_sociality %>% filter(Lice_genus=="Myrsidea")

ggplot(data=species_flock_priority_sociality_myrsidea, aes(x=sociality, y=Lice_total))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  stat_summary(fun.y=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

#Mesostigmatid

species_flock_priority_sociality_mesostigmatid<-species_flock_priority_sociality %>% filter(Mites_group2=="Mesostigmatid")

ggplot(data=species_flock_priority_sociality_mesostigmatid, aes(x=sociality, y=Mites2_total))+
  geom_jitter(aes(colour=foraging),width = 0.3)+
  stat_summary(fun.y=mean, geom="point",shape=20,size=5, color="black", fill="black")+
  theme_classic()

#How similar is the composition of lies in the bird species, color code by sociality

#community by species matrix, where each species is a communit and each genus of parasite is a "taxonomic unit

example_NMDS=metaMDS(community_matrix,k=2,trymax=100)


NMDS

example_NMDS<-metaMDS(matrix_diversity, distance ="bray", k=2)

meta

orditorp(example_NMDS)




