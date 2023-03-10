#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Data formating adn compilations                                                                         ###
### R-code                                                                          ###
### Jenny Munoz      
###
### Last update: feb 2023                                                   ###
################################################################################

# Loading packages --------------------------------------------------------
# libraries for easier manipulation of data
#install.packages("tidyr") 
install.packages("tidyverse") 
install.packages("dplyr")
install.packages ("data.table")
install.packages ("extrafont")
installed.packages("lubridate")  #for dates
#data cleaning
install.packages ("janitor")

#Other libraries for data analyses and visualizations
install.packages("vegan")
install.packages("ggplot2")
install.packages("devtools")

install.packages("knitr")
install.packages("ts")

install.packages("RColorBrewer")
#for models
install.packages("car") #Anova command
install.packages("lattice") #preliminary plots
install.packages("lme4") #for glmer (generalized linear mixed models) 
install.packages("visreg")  #extract confidence intervals and trend lines from GLMMs
install.packages("lsmeans") #least squared means
install.packages("MuMIn") #pseudo R squared for GLMMs
install.packages("emmeans")

# Phylogenetic component

install.packages("ape")
#install.packages("here")
install.packages("phytools")
install.packages("tidyverse")
install.packages("metafor")
install.packages("phangorn") # to reconstruct a maximum clade credibility tree
install.packages("rr2")
install.packages ( "MCMCglmm")
install.packages("phyr")
install.packages("TreeTools")

library(TreeTools)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) # to be able to run bayesian inference on the PGLMM models more info here https://www.r-inla.org


library("devtools")
#devtools::install_github(repo = "https://github.com/hrue/r-inla", ref = "stable", subdir = "rinla", build = FALSE)

#Libraries for data
library(tidyverse)
#library(tidyr)
#library(plyr)
library(ggplot2)
library(dplyr)
library(data.table)
library(extrafont)
library(lubridate)
# data cleaning 

library(janitor)
# for data visualization
library(vegan)
library(ggplot2)
library(devtools)
library(knitr)
library(ts)
library(RColorBrewer)


#libraries for models and visualizations

library(lattice) #preliminary plots
library(car) #Anova command
library(lsmeans) #least squared means
library(lme4) #for glmer (generalized linear mixed models) 
library(visreg) #extract confidence intervals and trend lines from GLMMs
library(MuMIn) #pseudo R squared for GLMMs
library(emmeans)

#library(dplyr) 

#Phylogenetic component

library(ape)
#library(here)
library(phyr)
library(phytools)
library(metafor)
library (phangorn) # to reconstruct a maximum clade credibility tree
library (rr2)
library (MCMCglmm)

library(tidyverse)
library(skimr)
install.packages("TreeTools")
library(TreeTools)

#Libraries ofr plots
library(gridExtra)
library(ggpubr)
library(grid)


#[Overall] Part 0 Data summary and exploration -------------------------------------------------
# tHis part is to calculate some general descriptive data summaries 
# Presence absence
ectos_df<-read.csv("data/7.ectoparasite_df_presence_absence.csv") # data on presence absence
ectos_df<-ectos_df %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
ectos_df<-ectos_df%>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations

# we have presence absence for 62 host social species and for 27 non social species  in manu

ectos_df%>% group_by(Family,species_clean,species_jetz ) %>%  # we have abundance for 62 host social species and for 27 non social species  in manu
  
  
  as.data.frame(unique( ectos_df$Family)) %>% View()
ectos_df<-read.csv("data/5.ectos_pres_abs_df.csv") # data on presence absence
ectos_df<-ectos_df%>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations

View(ectos_df)
ectos_df  %>% group_by(sociality) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu

non_social<-ectos_df %>% filter( sociality=="0") # we have OCURRENCE  for 88 host social species and for 45 non social species  in manu
social<-ectos_df %>% filter( sociality=="1") # we have OCURRENCE for 88 host social species and for 45  non social species  in manu

mean(non_social$sample_size)
mean(social$sample_size)


# Create  summary tables of bird species sample size and family and type of parasite
names(ectos_df)
families<-as.data.frame(ectos_df%>% group_by(BLFamilyLatin) %>% summarize(n())) 
write.csv(families, "tables/1.Family_sample_size_manu.csv")
individuals <-as.data.frame(ectos_df%>% group_by(species_clean) %>% summarize(n())) 
write.csv(individuals, "tables/1.Individuals_sample_size_manu.csv")

ectos_df%>% group_by(BLFamilyLatin, species_clean) %>% summarize(total=sum(bill_tidy))



# Abundance 

# Lice

lice_df_abundance<-read.csv("data/7.lice_df_abundance.csv")
# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
lice_df_abundance<-lice_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

lice_df_abundance<-lice_df_abundance %>% distinct( species_jetz,.keep_all = TRUE)

lice_df_abundance %>% filter( sociality=="0") %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu
lice_df_abundance %>% filter( sociality=="1") %>%  View() # we have abundance for  host social species and for 45 non social species  in manu

lice_df_abundance  %>% group_by(sociality) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu

lice_df_abundance  %>% group_by(elevation_cat) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu


lice_df_abundance<-read.csv("data/5.lice_df_abundance_manu.csv") #  lice abundance manu only 

lice_df_abundance %>% group_by(sociality) %>% summarize(ave=sd(total_lice)) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu


# Mean abundance
mean_lice_abundance<-read.csv("data/5.lice_df_abundance_means.csv")
mean_lice_abundance<-mean_lice_abundance %>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations

mean_lice_abundance %>% filter( sociality=="1") # we have abundance for 62 host social species and for 27 non social species  in mau
mean_lice_abundance %>% filter( sociality=="0") # we have abundance for 62 host social species and for 27 non social species  in mau

mean_lice_abundance  %>% group_by(sociality) %>% summarize(ave=sd(mean_lice)) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu

#Mites

mites_df_abundance<-read.csv("data/7.mites_df_abundance.csv")
names(mites_df_abundance)
mites_df_abundance <-mites_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

mites_df_abundance<-mites_df_abundance %>% distinct( species_jetz,.keep_all = TRUE)

mites_df_abundance %>% filter( sociality=="0") %>% View()# we have abundance for  host social species and for  non social species  in manu
mites_df_abundance %>% filter( sociality=="1") %>%  View() # we have abundance for  host social species and for 45 non social species  in manu

mites_df_abundance  %>% group_by(sociality) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu
mites_df_abundance  %>% group_by(elevation_cat) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu

mites_df_abundance %>% group_by(sociality) %>% summarize(ave=sd(total_mites)) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu


###
# Diversity
##

lice_richness_manu_sp<-read.csv("data/5.lice_richness_sp_df_manu.csv")
lice_richness_manu_sp<-lice_richness_manu_sp %>% distinct( species_jetz,.keep_all = TRUE)

lice_richness_manu_sp %>% filter( sociality=="0") %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu
lice_richness_manu_sp %>% filter( sociality=="1") %>%  View() # we have abundance for  host social species and for 45 non social species  in manu

lice_richness_manu_sp %>% filter( sociality=="0") %>% summarize(total_samples=sum(n_samples_lice))# we have abundance for 62 host social species and for 27 non social species  in manu
lice_richness_manu_sp %>% filter( sociality=="1") %>% summarize(total_samples=sum(n_samples_lice))# we have abundance for 62 host social species and for 27 non social species  in manu

lice_richness_manu_sp %>% group_by(elevation_cat) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu
lice_richness_manu_sp %>% group_by(foraging_cat) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu

# [Presence_absence DATA] Part1_Ectoparasites_ presence__using binary data  1/0-----------------------------

###_###_###_###_##
# The data
###_###_###_###_##

# Gettting the data ready
ectoparasites_df<-read.csv("data/data_analyses/7.ectoparasite_df_presence_absence.csv") # data on presence absence
names( ectoparasites_df)
unique(ectoparasites_df$Mites)
View(phylogeny)
class(phylogeny)

# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
ectoparasites_df<-ectoparasites_df %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

# Re-strudture the data
# Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?

ectoparasites_df <-ectoparasites_df  %>% mutate_at("species_jetz", str_replace, " ", "_")
str( ectoparasites_df)

ectoparasites_df$elevation_cat<-as.factor(ectoparasites_df$elevation_cat)
ectoparasites_df$foraging_cat<-as.factor(ectoparasites_df$foraging_cat)
ectoparasites_df$species_jetz<-as.factor(ectoparasites_df$species_jetz)
ectoparasites_df$sociality<-as.numeric(ectoparasites_df$sociality)

###_###_###_###_##
#The models for presence absence
###_###_###_###_##

# all ectos together

# WARNING STILL NEED TO INCLUDE SAMPLE SIZE IN THE ANALYSES AS A RANDOM EFFECT but maybe not here because it is individual samples because 
#the higher the sample size the higher teh probability to find ectos
# this will be required in the diversity analyses
ectoparasites_df<- ectoparasites_df %>% mutate(ectoparasites_PA=Lice+Mites+Ticks)
ectoparasites_df$ectoparasites_PA[ectoparasites_df$ectoparasites_PA>=1]<-1   # convert the numerical values that we have without core to lowland iquitos
unique(ectoparasites_df$ectoparasites_PA)

ectos_pres_abs<-ectoparasites_df %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_PA_max=max(ectoparasites_PA), sample_size=(n())) 

species_atributes<-ectoparasites_df %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean,species_binomial,BLFamilyLatin )
species_attributes_distict<-distinct( species_atributes)

ectos_pres_abs_df<-right_join(species_attributes_distict, ectos_pres_abs, by="species_jetz")   %>% arrange(elevation_cat)

ectos_pres_abs_df<-ectos_pres_abs_df %>% distinct(species_jetz,ectoparasites_PA_max,.keep_all = TRUE) # Eliminates species taht ocurr at two elevatiosn and keep only one ( in alphabetical order of the stations, e.g for galbula low_elevations is kept and montane is deleted, but the samples size is mantained this is jut to be able to do the phylogenetic analyses properly)

View(ectos_pres_abs_df)

#write.csv(ectos_pres_abs_df, "data/5.ectos_pres_abs_df.csv")

names(ectoparasites_df)

# [Prevalence DATA ]Part1 Ectoparasite models_using prevalence data (Corrections after meeting with Jullie et al  using probabilities)_Ectoparasite models_  -----------------------------

###_###_###_###_##
# The data
###_###_###_###_##

#Getting the data ready
ectoparasites_df<-read.csv("data/data_analyses/7.ectoparasite_df_presence_absence.csv") # data on presence absence with all variables 
names( ectoparasites_df)
#unique(ectoparasites_df$Mites)
#View(phylogeny)
#class(phylogeny)

# including other data
#number of detections in flocks
elevation_midpoint<-read.csv("data/1.elevation_midpoint_manu_species.csv")
anti_join(ectoparasites_df,elevation_midpoint, by=c("species_clean"="species_taxonomy_SACC_2021")) # speceis that are in the ectoparasite list that do not have a matcj in b 
ectoparasites_df<-left_join(ectoparasites_df,elevation_midpoint,by=c("species_clean"="species_taxonomy_SACC_2021"))

# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
ectoparasites_df<-ectoparasites_df %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

###_###_###_###_##
# all ectos together  species level
###_###_###_###_##

# INCLUDE SAMPLE SIZE IN THE ANALYSES AS A RANDOM EFFECT but maybe not here because it is individual samples because 
#the higher the sample size the higher teh probability to find ectos

# Re-strudture the data
# Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?
ectoparasites_df <-ectoparasites_df  %>% mutate_at("species_jetz", str_replace, " ", "_")

# Sumarize the parasite prevalence 
ectoparasites_df<- ectoparasites_df %>% mutate(ectoparasites_PA=Lice+Mites+Ticks)
ectoparasites_df$ectoparasites_PA[ectoparasites_df$ectoparasites_PA>=1]<-1   # convert the numerical values that we have without core to lowland iquitos
unique(ectoparasites_df$ectoparasites_PA)

#str(ectoparasites_df$ectoparasites_PA)
ectoparasites_df %>% group_by(ectoparasites_PA) %>% 
summarise(n())

ectos_pres_abs<-ectoparasites_df %>% group_by(species_jetz ) %>% 
summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation_midpoint=max(elevation_midpoint)) %>% 
mutate(proportion_ectoparasites=ectoparasites_presence/sample_size)

species_atributes<-ectoparasites_df %>% select( species_jetz, species_clean,species_binomial,BLFamilyLatin, sociality,mass_tidy,elevation_cat, foraging_cat )
species_attributes_distict<-distinct(species_atributes)

ectos_pres_abs_df<-right_join(species_attributes_distict, ectos_pres_abs, by="species_jetz")   %>% arrange(elevation_cat)

ectos_pres_abs_df<-ectos_pres_abs_df %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE) # Eliminates species taht ocurr at two elevatiosn and keep only one ( in alphabetical order of the stations, e.g for galbula low_elevations is kept and montane is deleted, but the samples size is mantained this is jut to be able to do the phylogenetic analyses properly)

ectos_pres_abs_df<- ectos_pres_abs_df %>% rename(family=BLFamilyLatin)
#View(ectos_pres_abs_df)

#write.csv(ectos_pres_abs_df, "data/data_analyses/7.dff_ectos_pres_abs_species.csv")


###_###_###_###_##
# all ectos together  (individual level)
###_###_###_###_##

ectoparasites_df<-read.csv("data/data_analyses/7.ectoparasite_df_presence_absence.csv") # data on presence absence with all variables 

# including other data
#number of detections in flocks
elevation_midpoint<-read.csv("data/1.elevation_midpoint_manu_species.csv")
anti_join(ectoparasites_df,elevation_midpoint, by=c("species_clean"="species_taxonomy_SACC_2021")) # speceis that are in the ectoparasite list that do not have a matcj in b 
ectoparasites_df<-left_join(ectoparasites_df,elevation_midpoint,by=c("species_clean"="species_taxonomy_SACC_2021"))

# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
ectoparasites_df<-ectoparasites_df %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

# INCLUDE SAMPLE SIZE IN THE ANALYSES AS A RANDOM EFFECT but maybe not here because it is individual samples because 
#the higher the sample size the higher teh probability to find ectos

# Re-strudture the data
# Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?
ectoparasites_df <-ectoparasites_df  %>% mutate_at("species_jetz", str_replace, " ", "_")

# Summarize the parasite prevalence per individual
ectoparasites_df<- ectoparasites_df %>% mutate(ectoparasites_PA=Lice+Mites+Ticks)
ectoparasites_df$ectoparasites_PA[ectoparasites_df$ectoparasites_PA>=1]<-1   # convert the numerical values that we have without core to lowland iquitos
unique(ectoparasites_df$ectoparasites_PA)

# select relevant columns and delete repited rows ( surpridinly we have repited rows!!!!)

ectos_pres_abs_df_individuals<-ectoparasites_df %>% 
  select(sample,Sample_label,Sample_No,Full_Label,species_jetz, species_clean,species_binomial,BLFamilyLatin, sociality,mass_tidy,elevation_cat, foraging_cat,Powder.lvl,  elevation_midpoint, ectoparasites_PA, Lice, Mites, Ticks) %>% 
  rename(family=BLFamilyLatin)%>% 
  distinct(Full_Label,.keep_all = TRUE)

dim(ectos_pres_abs_df_individuals )


# Importnat to make sure  that there are not repited samples 

#write.csv(ectos_pres_abs_df_individuals, "data/data_analyses/7.dff_ectos_pres_abs_individuals.csv")

names(ectoparasites_df )
species_attributes_distict<-distinct(species_atributes)


# [DATA Abundance]Part2_Ectoparasite models_  -----------------------------
###_###_###_###_
# Lice
###_###_###_###_

#Response variable count 
#Random effects species ( we have multiple measurements per species), account to variation between species other than phylogenetic, maybe redundat with foraging high
#Random effect phylogeny
#Random effect elevation 

# Read the files
ectos_df<-read.csv("data/data_analyses/7.dff_ectos_pres_abs_individuals.csv")  # general samples list to include teh zeros 
unique (ectos_df$Lice)
  
# data on lice prevalence at the ind level
lice_df_abundance<-read.csv("data/data_analyses/7.lice_df_abundance.csv",na.strings=c("_","","<NA>")) %>% 
  filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos") %>% 
  select(total_lice,Female, Male,Nimph, Eggs, Full_Label, species_binomial) %>% 
  distinct(Full_Label,.keep_all = TRUE) 


# include the zeros 
ectos_dff<-full_join(ectos_df, lice_df_abundance, by="Full_Label")
ectos_dff$total_lice[is.na(ectos_dff$total_lice)] <- 0
ectos_dff$Female[is.na(ectos_dff$Female)] <- 0
ectos_dff$Male[is.na(ectos_dff$Male)] <- 0
ectos_dff$Nimph[is.na(ectos_dff$Nimph)] <- 0
ectos_dff$Eggs[is.na(ectos_dff$Eggs)] <- 0

# Exclude the zeros
#

###_###_###_###_
# Mites
###_###_###_###_

mites_df_abundance<-read.csv("data/data_analyses/7.mites_df_abundance.csv",na.strings=c("_","","<NA>")) %>% 
  filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos") %>% 
  select(total_mites,	total_mesostigmatidae,	total_no_feathers_mites, Full_Label, species_binomial) %>% 
  distinct(Full_Label,.keep_all = TRUE) 


View(ectos_dff)
ectos_df_abun<-left_join(ectos_dff, mites_df_abundance, by="Full_Label")
ectos_df_abun$total_mites[is.na(ectos_df_abun$total_mites)] <- 0
ectos_df_abun$total_mesostigmatidae[is.na(ectos_df_abun$total_mesostigmatidae)] <- 0
ectos_df_abun$total_no_feathers_mites[is.na(ectos_df_abun$total_no_feathers_mites)] <- 0

#### The data file # Yupi!! ( now just nees to include the diversity recalculated)
#write.csv(ectos_df_abun,"data/data_analyses/7.dff_all_ectos_prevalence_abundance_individual.csv")

###_###_###_###_
# Get the elevations for each sample
###_###_###_###_


# get the elevations on the different files at the individual sample level

# Extracting the elevation
captures<-read.csv("data/data_analyses/0_netting_manu_complete_for_ectoparasite_samples.csv") %>%
  select(ectoparasite_code,species,species_taxonomy_SACC_2021,jetz_species,elevation_all_years,elevation_net_elevation)
names(captures)

ecto_parameters<-read.csv("data/data_analyses/0_ectoparasites_samples_metadata.csv") %>%
  select(ectoparasite_code,Species,Event.Number,Station) %>% 
  rename(Species_parameters=Species)

  
samples_elevation<-right_join(captures,ecto_parameters, by="ectoparasite_code", multiple="all") %>% 
  distinct(Event.Number,.keep_all = TRUE)
#write.csv(samples_elevation, "data/data_analyses/0_dff_ectoparasite_samples_elevation_metadata.csv")

#Integrating  data with elevations 

ectos_samples_elevations<-read.csv("data/data_analyses/0_dff_ectoparasite_samples_elevation_metadata.csv")
ectoparasites_dff<-read.csv("data/data_analyses/7.dff_all_ectos_prevalence_abundance_individual.csv")

a<-(as.data.frame(ectoparasites_dff$Full_Label))%>% mutate(sample_code=ectoparasites_dff$Full_Label) %>% select(sample_code) %>% arrange(desc(sample_code))
b<-(as.data.frame(ectos_samples_elevations$Event.Number))%>% mutate(sample_code=ectos_samples_elevations$Event.Number) %>% select(sample_code) %>% arrange(desc(sample_code))

missing<-as.list(setdiff(a,b))

# For this samples  we do not have information  We needed to included manually 
# [1] "SP130013" "P121006"  "P120889"  "P120085"  "P115912" 


#SP130013  Chiroxiphia_boliviana
#P121006	Rupicola_peruvianus 1362 The correct name of this is SP121006 AND THE ELEVATION 1362
#"P120889" is Leptopogon superciliaris CORRRECT NAME SP120889  AND THE ELEVATION IS 1395
#P120085 Formicarius rufipectus   correct name  SP120085 elev 1362
#P115912 Leptopogon amaurocephalus elevation 381 


# we will add the elevation manually to this samples

ectoparasites_dff_elevations<-left_join(ectoparasites_dff,ectos_samples_elevations, by=c("Full_Label"="Event.Number"))
write.csv(ectoparasites_dff_elevations,"data/data_analyses/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv")

setdiff()




