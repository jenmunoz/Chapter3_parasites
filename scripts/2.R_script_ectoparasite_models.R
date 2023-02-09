#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models                                                                           ###
### R-code                                                                          ###
### Jenny Munoz      
###
### Last update: Oct 26 2022                                                     ###
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

#Libraries ofr plots
library(gridExtra)
library(ggpubr)
library(grid)

# #[Overall] Part 0 Data summary -------------------------------------------------
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



# [Presence_absence] Part1_Ectoparasite models_using binary data -----------------------------
# To underestand teh models and out puts better see this paper: https://www.iecolab.org/wp-content/uploads/2020/10/phyr_2020.pdf
# also see this blog with an example https://daijiang.github.io/phyr/articles/phyr_example_empirical.html
#All models fitted with pglmm() have class of communityPGLMM. Here is a list of functions that can be used to these models.

#Notes for this analyses we are trying to incorportae phylogeny and random effect but our response variable is binomial (0,1), or counts (abundance)
# Because of this we can not do a PGLS ( which only takes contonous data as response varaible)
# I run teh model with a gneral mixed effect model with out the phylogenetic correction first and then PGLMM and tehn McmcPGLM 
# The options are PGLMM AND MCMCpglmmm, there is also the opction binarypglmm but I am having troble underestanding the sintax
# MCMCpglmm uses bayesian approach 
# PGLMM looks straightforward to underetand but when using R2 to evaluta goodness of fit seem very low...
# prevalence is a 1 or 0 so we can use binomial # but elevation can not be continuos to be entered as a random effect logit 

#notes blog for PGLS when response variable is contonous blog from liam revel http://www.phytools.org/Cordoba2017/ex/4/PGLS.html
#Binarypglmm blog https://rdrr.io/cran/ape/man/binaryPGLMM.html
#MCMCglmm https://ourcodingclub.github.io/tutorials/mcmcglmm/

###_###_###_###_##
# The data
###_###_###_###_##

# Getttingthe data ready
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

species_atributes<-ectoparasites_df %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean)
species_attributes_distict<-distinct( species_atributes)

ectos_pres_abs_df<-right_join(species_attributes_distict, ectos_pres_abs, by="species_jetz")   %>% arrange(elevation_cat)

ectos_pres_abs_df<-ectos_pres_abs_df %>% distinct(species_jetz,ectoparasites_PA_max,.keep_all = TRUE) # Eliminates species taht ocurr at two elevatiosn and keep only one ( in alphabetical order of the stations, e.g for galbula low_elevations is kept and montane is deleted, but the samples size is mantained this is jut to be able to do the phylogenetic analyses properly)
  

View(ectos_pres_abs_df)

#write.csv(ectos_pres_abs_df, "data/5.ectos_pres_abs_df.csv")

names(ectoparasites_df)


# Reading the files 

ectos_df<-read.csv("data/5.ectos_pres_abs_df.csv") # data on presence absence
unique(ectos_df$species_jetz)
str(ectos_df) # should have same number of rows than teh phylogeny 
#phylogeny<- read.nexus("data/phylo_data/1_host_consensus_tree_Manuspecies.nex") 
phylogeny_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_ectos_pres_abs.nex")

tips<- as.data.frame(phylogeny_rooted$tip.label)
names<-as.data.frame(ectos_df$species_jetz)
anti_join(tips,names, by=c("phylogeny_rooted$tip.label"="ectos_df$species_jetz"))

# Using pglmmm
#I am not sue about including foraging cat since that some how is included in teh variation per species + 
#elevation_cat # only hs three categories so i am not sure I can used as a ranmod effect

ecto_PA_pglmm <-  phyr::pglmm(ectoparasites_PA_max ~ sociality+(1|species_jetz__)+(1|sample_size)+(1|foraging_cat)+(1|elevation_cat), 
                        data = ectos_df, 
                        family = "binomial",
                        cov_ranef = list(species_jetz= phylogeny_rooted), #class phylo
                        #bayes = TRUE,
                        REML = TRUE, 
                        verbose = TRUE,
                        s2.init = .25) # what is this last parameter for

summary(ecto_PA_pglmm)
predict(ecto_PA)
rr2::R2(ecto_PA_pglmm)
class(ecto_PA_pglmm )

# Using a glmm
ectos_PA_glmm<-  lme4::glmer(ectoparasites_PA_max ~ sociality+ (1|elevation_cat)+(1|species_jetz)+(1|sample_size),  #+1(1|foraging_cat)
                             data = ectos_df, 
                             family = "binomial",
                             verbose = TRUE)

summary(ecto_PA_glmm)
rr2::R2(ecto_PA_glmm)

#binary without random effects
ecto_PA_pglmm <- binaryPGLMM(ectoparasites_PA_max ~ sociality, 
                              data=ectos_df, 
                              phy=phylogeny_rooted)
                      
summary(ecto_PA_pglmm )
predict(ecto_PA)
rr2::R2(ecto_PA)

# Using MCMC without phylogeny 
MCMC<-MCMCglmm(ectoparasites_PA_max ~ sociality, random = ~foraging_cat+elevation_cat+species_jetz,
               family="categorical", data=ectos_df)

summary(MCMC)
plot(MCMC$Sol)
plot(MCMC$VCV)

# Using MCMC with Phylogeny ( this in not working yet)
# Some instructiosn in how to fit the data here https://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-2-multiple-measurements-model-mcmcglmm


V<-vcv(phylogeny_rooted)
V<-V/max(V)
detV<-exp(determinant(V)$modulus[1])
V <- V/detV^(1/n)

invV <- Matrix(solve(V),sparse=T)
ectos_df$species_jetz <- phylogeny_rooted$tip.label
rownames(invV) <- ectos_df$species_jetz

nitt <- 43000
thin <- 10
burnin <- 3000

prior <- list(R=list(V=0.5, fix=1), G=list(G1=list(V=0.5, nu=1, alpha.mu=0, alpha.V=0.5)))

ectos_pa_MCMC<-MCMCglmm(ectoparasites_PA_max ~ sociality, random=~species_jetz+elevation_cat, ginvers=list(species_jetz=invV),
                 data=ectos_df, slice=TRUE, nitt=nitt, thin=thin, burnin=burnin,
                 family="categorical", prior=prior, verbose=FALSE)

rr2::R2(ectos_pa_MCMC)

summary(ectos_pa_MCMC)


# GOODNESS OF FIT  
# Read this paper to underestand better https://academic.oup.com/sysbio/article/68/2/234/5098616?login=false
# An additional complication is created by model with random effects. Given that random effects are very flexible model components 
#(for example, nothing stops you from fitting a random effect for each observation in your dataset), a straight-up calculation of variance explains isn’t meaningful. T
#That said, methods that can produce a useful R2 metric in the complex situation have been developed. 
#The package rr2 is able to calculate several flavors of R2, and supports phyr’s pglmm model object. Let’s try it!
#
install.packages("rr2")
library(rr2)
rr2::R2(ecto_PA)
rr2::R2(ecto_PA_model_glmm)
rr2::R2(MCMC)




# Part1 Ectoparasite models_using prevalence data (Corrections after meeting with Jullie et al  using probabilities)_Ectoparasite models_  -----------------------------
#All models fitted with pglmm() have class of communityPGLMM. Here is a list of functions that can be used to these models.
# In this case we will analyses proportion data arising from counts (counts base data, See https://fukamilab.github.io/BIO202/04-B-binary-data.html#glm_for_proportional_data)

#Notes for this analyses we are trying to incorportae phylogeny and random effect but our response variable is binomial (0,1), or probaility ( 0 to 1) or counts (abundance)
# Because of this we can not do a PGLS ( which only takes contonous data as response varaible)
# I run teh model with a gneral mixed effect model with out the phylogenetic correction first glmm and then PGLMM and tehn McmcPGLM 
# The options are PGLMM AND MCMCpglmmm, there is also the opction binarypglmm but I am having troble underestanding the sintax
# MCMCpglmm uses bayesian approach 
# PGLMM looks straightforward to underetand but when using R2 to evaluta goodness of fit seem very low...
# probability of parasite ocurrence is bounded a 1 or 0 so we can use binomial # but elevation can not be continuos to be entered as a random effect logit 

#notes blog for PGLS when response variable is contonous blog from liam revel http://www.phytools.org/Cordoba2017/ex/4/PGLS.html
#Binarypglmm blog https://rdrr.io/cran/ape/man/binaryPGLMM.html
#MCMCglmm https://ourcodingclub.github.io/tutorials/mcmcglmm/


###_###_###_###_##
# The data
###_###_###_###_##

# Getting the data ready
ectoparasites_df<-read.csv("data/data_analyses/7.ectoparasite_df_presence_absence.csv") # data on presence absence with all variables 
names( ectoparasites_df)
unique(ectoparasites_df$Mites)
View(phylogeny)
class(phylogeny)

# including other data
#number of detections in flocks
elevation_midpoint<-read.csv("data/1.elevation_midpoint_manu_species.csv")
#anti_join(ectoparasites_df,elevation_midpoint, by=c("species_clean"="species_taxonomy_SACC_2021")) # speceis that are in the ectoparasite list that do not have a matcj in b 

ectoparasites_df<-left_join(ectoparasites_df,elevation_midpoint,by=c("species_clean"="species_taxonomy_SACC_2021"))

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
ectoparasites_df$elevation_midpoint<-as.numeric(ectoparasites_df$elevation_midpoint)

###_###_###_###_##
#The models for presence absence
###_###_###_###_##

# all ectos together

# INCLUDE SAMPLE SIZE IN THE ANALYSES AS A RANDOM EFFECT but maybe not here because it is individual samples because 
#the higher the sample size the higher teh probability to find ectos
# this will be required in the diversity analyses
ectoparasites_df<- ectoparasites_df %>% mutate(ectoparasites_PA=Lice+Mites+Ticks)
ectoparasites_df$ectoparasites_PA[ectoparasites_df$ectoparasites_PA>=1]<-1   # convert the numerical values that we have without core to lowland iquitos
unique(ectoparasites_df$ectoparasites_PA)

str(ectoparasites_df$ectoparasites_PA)
ectoparasites_df %>% group_by(ectoparasites_PA) %>% 
  summarise(n())

ectos_pres_abs<-ectoparasites_df %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation_midpoint=max(elevation_midpoint)) %>% 
  mutate(proportion_ectoparasites=ectoparasites_presence/sample_size)

species_atributes<-ectoparasites_df %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean, mass_tidy,Family )
species_attributes_distict<-distinct(species_atributes)

ectos_pres_abs_df<-right_join(species_attributes_distict, ectos_pres_abs, by="species_jetz")   %>% arrange(elevation_cat)

ectos_pres_abs_df<-ectos_pres_abs_df %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE) # Eliminates species taht ocurr at two elevatiosn and keep only one ( in alphabetical order of the stations, e.g for galbula low_elevations is kept and montane is deleted, but the samples size is mantained this is jut to be able to do the phylogenetic analyses properly)

View(ectos_pres_abs_df)

# write.csv(ectos_pres_abs_df, "data/data_analyses/7.dff_ectos_pres_abs.csv")

###_###_###
# [Prevalence_Presence_absence] THE MODELS
###_###_###

ectos_df<-read.csv("data/data_analyses/7.dff_ectos_pres_abs.csv")# data on prevalence FOR THE MODEL
# Converting variables as needed
ectos_df$sociality<-as.factor(ectos_df$sociality)

#as.data.frame(unique(ectos_df$species_jetz)) %>% View()
str(ectos_df) # should have same number of rows than teh phylogeny 
#phylogeny<- read.nexus("data/phylo_data/1_host_consensus_tree_Manuspecies.nex") 
phylogeny_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_ectos_pres_abs.nex")

is.ultrametric(phylogeny_rooted)
# Make sure the tips and the names on the file coincide
tips<- as.data.frame(phylogeny_rooted$tip.label)
names<-as.data.frame(ectos_df$species_jetz)
anti_join(tips,names, by=c("phylogeny_rooted$tip.label"="ectos_df$species_jetz"))

# Important!!!!!Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_rooted$tip.label,rownames(ectos_df)),]

#MODEL FITTING AND EVALUALLING MODEL FIT GRAOHICALLY

###_###_###_ Important for model fitting
#predict ( model) # plot predicted vlues and confidence limits on the logit scale, the points are not logit tranformed data, they are "working values to fit the model.....dificult to interpret 
#visreg(model) # Similarly  plot predicted values and confidence limits on the logit scale, dificult to interpret
#fitted (model) # willl give us the predicted values transformed to the original scale 
#visreg(model, scale="response") # willl give us the predicted values transformed to the original scale 
###_###_###_ ###_###_###_ ###_###_###_ 


# lETS CONSIDER A MIXED EFFECT MODEL first
# # Using a glmm
ectos_prevalence_glmm<-  lme4::glmer(proportion_ectoparasites ~ sociality+ sample_size+ (1|elevation_cat)+(1|species_jetz),  #+1(1|foraging_cat)
                                     data = ectos_df, 
                                     family = "binomial" (link = "logit"),
                                     verbose = TRUE)
summary(ectos_prevalence_glmm)
rr2::R2(ectos_prevalence_glmm)

log(1.21978)

# How well does it predict the data 
newdata <- data.frame(elevation_cat = ectos_df$elevation_cat,
                      species_jetz = ectos_df$species_jetz,
                      sociality = ectos_df$sociality,
                      sample_size = ectos_df$sample_size)

newdata$Predicted_Response <- predict(ectos_prevalence_glmm, newdata = newdata, type = "response")


# model  predicts overall patterns well?

# I can plot all predicted values or just the means 

# the mean predicted value
summary_predictions <- newdata %>%   # Specify data frame
  group_by(sociality)  %>%    # Specify group indicator
  dplyr::summarise(proportion_ectoparasites <- mean(Predicted_Response))

colnames(summary_predictions) <- c("sociality", "proportion_ectoparasites")

png("predicted", width = 2500, height = 3100, res = 300, units = "px")
ggplot(summary_predictions, aes(x = sociality, y = proportion_ectoparasites))+
  geom_point()+
  ylim(0,1)
dev.off()

# Using vis reg
visreg(ectos_prevalence_glmm, scale="response", ylim=c(0,1)) # If scale='response' for a glm, the inverse link function will be applied so that the model is plotted on the scale of the original response.


# all predicted values 
ggplot(newdata, aes(x = sociality, y =Predicted_Response))+
  geom_point()+
  ylim(0,1)
Predicted_Response


visreg(ectos_prevalence_glmm, scale="response", ylim=c(0,1)) # If scale='response' for a glm, the inverse link function will be applied so that the model is plotted on the scale of the original response.
# This is useful to obtain predicted values on the 0riginal scale 
emmeans(ectos_prevalence_glmm,"sociality", type="response")



# Using PGLMM [ Lets correct for phylogeny ]
#I am not sue about including foraging cat since that some how is included in teh variation per species + 
#elevation_cat # only hs three categories so i am not sure I can used as a ranmod effect
# we revomed (1|foraging_cat) because it was not significant in individual models 

unique(ectos_df$elevation_cat)
mean(ectos_df$proportion_ectoparasites) # mean prevalnece

lm(proportion_ectoparasites~1, data=ectos_df) # even the average is a linear regression

str( ectos_df)
ecto_prevalence_pglmm <-  phyr::pglmm(proportion_ectoparasites~sociality+sample_size+(1|species_jetz__)+(1|elevation_cat), #+elevation_midpoint
                              data = ectos_df, 
                              family = "binomial",
                              cov_ranef = list(species_jetz= phylogeny_rooted), #class phylo
                              #bayes = TRUE,
                              REML = TRUE, 
                              verbose = TRUE,
                              s2.init = .25) # what is this last parameter for
summary(ecto_prevalence_pglmm)
predict(ecto_prevalence_pglmm)
rr2::R2(ecto_prevalence_pglmm)
class(ecto_prevalence_pglmm )

# Explore if the random effects are important?
# One way to get an idea is to run a likelihood ratio test on the random effect. 
# This can be achieved by using the pglmm_profile_LRT() function, at least for binomial models.

LRTest <- sapply(1:3, FUN = function(x) phyr::pglmm_profile_LRT(ecto_prevalence_pglmm, re.number = x))
colnames(LRTest) <- names(ecto_prevalence_pglmm$ss)
t(LRTest)

#Pages Lambda to see if there is phylogenetic signal
#[ WARNING THIS CAN ONLY BE USED IN CONTINOUS DATA NO IN DISCRETE DATA (E.G NO IN COUNTS DATA), becase discrete data does not follow Brownian Motion mode of evo. One optionis to use kappa model but in non ultrametric trees (  with genetic info?)
#Phylogenetic signal is generally recognized to be the tendency of related species to resemble one another; and Blomberg et al.'s (2003) K and Pagel's (1999) λ are two quantitative measures of this pattern.
#The metrics are quite different from one another. λ is a scaling parameter for the correlations between species, relative to the correlation expected under Brownian evolution.
#K is a scaled ratio of the variance among species over the contrasts variance (the latter of which will be low if phylogenetic signal is high). 
#λ has a nice natural scale between zero (no correlation between species) and 1.0 (correlation between species equal to the Brownian expectation).
#λ itself is not a correlation, but a scaling factor for a correlation, so λ>1.0 is theoretically possible. However, depending on the structure of the tree, λ>>1.0 is usually not defined.
#K, a variance ratio, is rescaled by dividing by the Brownian motion expectation. This gives it the property of having an expected value of 1.0 under Brownian evolution, but K for empirical and simulated datasets can sometimes be >>1.0. 
#Since they measure the qualitative tendency towards similarity of relatives in entirely different ways, we have no real expectation that they will be numerically equivalent - except (by design) under Brownian motion evolution.

# Make sure the names  are in the same order in the phylogeny and in the traits

phylogeny_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_ectos_pres_abs.nex")
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_rooted$tip.label,rownames(ectos_df)),]

tips<- as.data.frame(phylogeny_rooted$tip.label)
names<-as.data.frame(ectos_df$species_jetz)

# caLCULATE Lambda option 1 
# https://www.zoology.ubc.ca/~schluter/R/Phylogenetic.html

mydat2 <- ectos_df[, c("proportion_ectoparasites", "mass_tidy")]
phytools::phylosig(phylogeny_rooted, as.matrix(mydat2)[,c("proportion_ectoparasites")], method = "lambda") # phylogenetic signal in trait "x"

#phytools::phylosig(phylogeny_rooted, as.matrix(mydat2)[,c("sociality")], method = "lambda") # phylogenetic signal in trait "x"

# Option 2 not used here but similar results
names_sp = match(ectos_df$species_jetz, phylogeny_rooted$tip.label)
z = nlme::gls(proportion_ectoparasites ~ sociality+ sample_size, data = ectos_df, correlation = corPagel(0.5,  phylogeny_rooted, fixed=FALSE))
summary(z)


###_###_###
# PLOTS PREVALENCE [Presence Absence]
###_###_###

# PLotting the model predictions 
#### How well does PGLMM it predict the data
newdata <- data.frame(elevation_cat = ectos_df$elevation_cat,
                      species_jetz = ectos_df$species_jetz,
                      sociality = ectos_df$sociality,
                      sample_size = ectos_df$sample_size)

predictions<- predict(ecto_prevalence_pglmm,newdata = newdata, type = "response" ) ##newdata = newdata
ectos_df_predicted <- cbind(ectos_df, predictions)

str(ectos_df_predicted)

# lets calculae the mean of the predited values on the untransformed scale
summary_model <- ectos_df_predicted %>% 
  group_by(sociality) %>%      
  dplyr::summarise(prevalence= mean(proportion_ectoparasites), sd=sd(proportion_ectoparasites), n = n()) %>% 
  mutate(se= sd/(sqrt(n)))

colnames(summary_model) <- c("sociality", "prevalence", "sd", "n", "se")

head(summary_model)

# make a plot of model predictions (that also shows data)
png("predicted3", width = 3000, height = 3000, res = 300, units = "px")
ggplot(data = ectos_df_predicted, aes(x = sociality, y = Y_hat ))+
  # geom_point(data = ectos_df, aes(x=sociality, y = proportion_ectoparasites),color="grey",size=2)+
  geom_jitter(data = ectos_df, aes(x=sociality, y = proportion_ectoparasites),color="grey",size=3,width = 0.06)+
  geom_segment(data = summary, aes(x = sociality, y = prevalence, xend = sociality, yend = prevalence+se, color="red"),show_guide = FALSE)+
  geom_segment(data = summary, aes(x = sociality, y = prevalence, xend = sociality, yend = prevalence-se, color="red"),show_guide = FALSE)+
  geom_point(data = summary, aes(x=sociality, y = prevalence), color="red", size=4,shape=19)+
  scale_y_continuous("Ectoparasites prevalence", limits = c(0,1)) +
  scale_x_discrete("Sociality")+
  geom_hline(yintercept = 0.775, linetype = "dashed")+
  theme_classic(40)
dev.off()


# Vis eg does not support PGLMM
#visreg(ecto_prevalence_pglmm, scale="response", ylim=c(0,1)) # If scale='response' for a glm, the inverse link function will be applied so that the model is plotted on the scale of the original response.
#emmeans(ecto_prevalence_glmm, type="response")

mean(ectos_df$proportion_ectoparasites) # mean prevalnece observed

# Plotting the traits things in the  phylogeny 
#Some examples fr plotting http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html
# Plot phylogenetic tree and  the trait continous trait ( prevalence)
#ContMap #Function plots a tree with a mapped continuous character. 
#The mapping is accomplished by estimating states at internal nodes using ML with fastAnc, and then interpolating the states along each edge using equation [2] of Felsenstein (1985).

#contMap(tree, x, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL,
#lims=NULL, outline=TRUE, sig=3, type="phylogram", direction="rightwards", 
#plot=TRUE, ...)

ectos_df<-read.csv("data/data_analyses/7.dff_ectos_pres_abs.csv")# data on prevalence # should have same number of rows than teh phylogeny 
phylogeny_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_ectos_pres_abs.nex")
list.names=setNames(ectos_df$proportion_ectoparasites, ectos_df$species_jetz)

# Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_rooted$tip.label,rownames(ectos_df)),]

object = contMap(phylogeny_rooted, list.names, direction = "leftwards", plot=FALSE)

png("figures/figures_manuscript/Fig2.Prevalence_phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))
plot(object_color, mar=c(5.1,1,1,1), ylim=c(1-0.09*(Ntip(object$tree)-1), Ntip(object$tree)), legend=FALSE,lwd=4,fsize=0.5)
#add.color.bar(5, object2, title = "Infection Prevalence", lims = object$lims, digits = 3, prompt=FALSE,x=0,y=1-0.08*(Ntip(object$tree)-1), lwd=4,fsize=1,subtitle="")
add.color.bar(5, object_color$cols, title = "", lims = object$lims, digits = 3, prompt=FALSE,x=10,y=-5, lwd=4,fsize=1,subtitle="Ectoparasites Prevalence")
dev.off()
title(xlab = "Time from root (Ma)")


####  Plot flcoking as a binomial trait on the phylogeny

phylogeny_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_ectos_pres_abs.nex")
ectos_df<-read.csv("data/data_analyses/7.dff_ectos_pres_abs.csv")# data on prevalence # should have same number of rows than teh phylogeny 
fmode<-as.factor(setNames(ectos_df$sociality,ectos_df$species_jetz))

# Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_rooted$tip.label,rownames(ectos_df)),]

View(fmode)

png("figures/figures_manuscript/Fig1.sociality_phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
dotTree(phylogeny_rooted,fmode,colors=setNames(c("red","black"),                                              c("1","0")),ftype="i",fsize=0.5, lwd=4)
dev.off()

# Combining both plots

fmode<-as.factor(setNames(ectos_df$sociality,ectos_df$species_jetz))

object = contMap(phylogeny_rooted, list.names, direction = "leftwards", plot=FALSE)
object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))

png("figures/figures_manuscript/Fig3.Sociality_and_prevalence_phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
plot(dotTree(phylogeny_rooted,fmode,colors=setNames(c("red","black"),c("1","0")),ftype="i",fsize=0.5, lwd=4),text(x=10,y=-5,"Mixed-species flocks",pos=1))

plot(object_color$tree,colors=object_color$cols,add=TRUE,ftype="off",lwd=5,fsize=0.5,
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)

add.color.bar(10, object_color$cols, title = "", lims = object$lims, digits = 3, prompt=FALSE,x=70,y=-5, lwd=4,fsize=1,subtitle="Ectoparasites Prevalence",pos=4)

dev.off()


#GADGETS 

# Color palette
require(RColorBrewer)
## Loading required package: RColorBrewer
ColorPalette <- brewer.pal(n = 9, name = "Blues")
# Plotting
plot(seedplantstree, type="p", use.edge.length = TRUE, label.offset=0.01,cex=0.6)
tiplabels(pch=21,bg=ColorPalette[maxH.cat],col="black",cex=1,adj=0.505)


###_





# [Abundance]Part2_Ectoparasite models_  -----------------------------
# Lice
#Response variable count 
#Random effects species ( we have multiple measurements per species), account to variation between species other than phylogenetic, maybe redundat with foraging high
#Random effect phylogeny
#Random effect elevation 

library(ape)
library(MCMCglmm)

# Read teh files
#lice_df_abundance<-read.csv("data/data_analyses/7.lice_df_abundance.csv")
# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
#lice_df_abundance<-lice_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
#unique (lice_df_abundance$species_jetz) 
#lice_df_abundance<-lice_df_abundance %>% distinct( species_jetz,.keep_all = TRUE)

#lice_df_abundance<-read.csv("data/5.lice_df_abundance_manu.csv") # data for manu only



#Add powder lever
#powder_bias<-read.csv("data/data_analyses/7.ectos_samples_powder_level.csv")
#lice_df_abundance_powder<-inner_join( lice_df_abundance, powder_bias, by="Full_Label") 
#lice_df_abundance_powder$powder_level<-as.factor(lice_df_abundance_powder$powder_level)

#write.csv(lice_df_abundance_powder,"data/data_analyses/7.dff_lice_abundance.csv")
str(lice_df_abundance_powder)

# Using pglmmm
lice_df_abundance<-read.csv("data/data_analyses/7.dff_lice_abundance.csv") 
unique(lice_df_abundance$powder_level)
phylo_lice_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")

# Make sure number in teh data and the phylogenty are consistent
lice_df_abundance <-lice_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
lice_df_abundance$elevation_cat<-as.factor(lice_df_abundance$elevation_cat)
lice_df_abundance$foraging_cat<-as.factor(lice_df_abundance$foraging_cat)
lice_df_abundance$species_jetz<-as.factor(lice_df_abundance$species_jetz)
lice_df_abundance$sociality<-as.factor(lice_df_abundance$sociality)

str(lice_df_abundance)

# Important! make sure data traist and phylogeny are in teh same order

rownames(lice_df_abundance) <- lice_df_abundance$species_jetz # first make it the row names 
lice_df_abundance<- lice_df_abundance[match(phylo_lice_rooted$tip.label,rownames(lice_df_abundance)),]

tips<- as.data.frame(phylo_lice_rooted$tip.label)
names<-as.data.frame(lice_df_abundance$species_jetz)

#I am not sue about including foraging cat since that some how is included in teh variation per species , 
#also sample size seems relevant but not sure how to included since the observations are individual.
# iTHINK THE BEST APPROACH IS TO calculate the mean abundance adn nclude sample size as a random effect

 
l_abun_pglmm <-  phyr::pglmm(total_lice ~ sociality+ (1|elevation_cat)+(1|species_jetz__)+ (1|powder_level), 
                              data = lice_df_abundance, 
                              family = "poisson",
                              cov_ranef = list(species_jetz= phylo_lice_rooted), #class phylo
                              #bayes = TRUE,
                              REML = TRUE,  # NOT SURE WHEN TO USE ML
                              verbose = TRUE,
                              s2.init = .25) # what is this last parameter for

summary(l_abun_pglmm)
predict(l_abun_pglmm)
rr2::R2(l_abun_pglmm)
R2.lik(l_abun_pglmm)
R2.pred(l_abun_pglmm)
R2.resid(l_abun_pglmm)
class(l_abun_pglmm) # i THINK r2 ONLY ALOW TO calculate goodness of fit in gaussian (community PGLMM)

# Using Glmm

#Similar results to above 
lice_abundance_glmm<-  lme4::glmer(total_lice ~ sociality+  (1|elevation_cat)+(1|species_jetz)+ (1|powder_level),  #+1(1|foraging_cat)
                                     data = lice_df_abundance, 
                                     family = "poisson",
                                     verbose = TRUE)

summary(lice_abundance_glmm)
rr2::R2(ectos_prevalence_glmm)


lice_abundance_glm<-glm(total_lice ~ sociality+ powder_level,  #+1(1|foraging_cat)
                                   data = lice_df_abundance, 
                                   family = "poisson")
                                   

#But we know that sample size might have an effect 
#Modeling the mean abundance to take care of the sample size
#### Sample size


mean_lice<-lice_df_abundance %>% group_by (species_jetz) %>% 
  summarize(mean_lice=mean(total_lice), sample_size=(n()))


View(mean_lice)
View(lice_df_abundance)
  
  #ectos_df  %>% group_by(sociality) %>% summarize(total_samples=n()) %>% View()# we have abundance for 62 host social species and for 27 non social species  in manu


species_atributes<-lice_df_abundance %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean) %>%  filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
species_attributes_distict<-distinct( species_atributes)

mean_lice_abundance<-right_join(species_attributes_distict, mean_lice, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 


###_###_####_###_
#write_csv(mean_lice_abundance,"data/data_analyses/7.dff_lice_abundance_means.csv")

#The data
mean_lice_abundance<-read.csv("data/data_analyses/7.dff_lice_abundance_means.csv")
mean_lice_abundance<-mean_lice_abundance %>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations


phylo_lice_rooted<-read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")


mean_lice_abundance <-mean_lice_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
mean_lice_abundance$elevation_cat<-as.factor(mean_lice_abundance$elevation_cat)
mean_lice_abundance$foraging_cat<-as.factor(mean_lice_abundance$foraging_cat)
mean_lice_abundance$species_jetz<-as.factor(mean_lice_abundance$species_jetz)
mean_lice_abundance$sociality<-as.factor(mean_lice_abundance$sociality)

str(mean_lice_abundance)
# The model 

hist(mean_lice_abundance$mean_lice)
unique(mean_lice_abundance$mean_lice)
names ( mean_lice_abundance)
l_abun_mean_pglmm <-  phyr::pglmm(mean_lice ~ sociality+sample_size+(1|elevation_cat)+(1|species_jetz__), 
                             data = mean_lice_abundance, 
                             family = "gaussian", #( maybe gaussian will be better in this case ?)
                             cov_ranef = list(species_jetz= phylo_lice_rooted), #class phylo
                             #bayes = TRUE,
                             REML= TRUE,  # NOT SURE WHEN TO USE ML
                             verbose = TRUE,
                             s2.init = .25) # what is this last parameter for

summary (l_abun_mean_pglmm)
rr2::R2(l_abun_mean_pglmm)

l_abun_mean_glmm <-  lme4::lmer(mean_lice ~ sociality+ (1|elevation_cat)+(1|species_jetz)+ (1|foraging_cat), 
                                  data = mean_lice_abundance)


#Mites: this is all the data not manu only

###_###_###_###_###_###_###_###
# Modeling no feather mites only (cause feathers mites can be misleding)
###_###_###_###_###_###_###_###

#mites_df_abundance<-read.csv("data/data_analyses/7.mites_df_abundance.csv") 
#powder_bias<-read.csv("data/data_analyses/7.ectos_samples_powder_level.csv")
#mites_df_abundance_powder<-inner_join(mites_df_abundance, powder_bias, by="Full_Label") 
#mites_df_abundance_powder$powder_level<-as.factor(mites_df_abundance_powder$powder_level)

#write.csv(mites_df_abundance_powder,"data/data_analyses/7.dff_mites_abundance.csv")

mites_df_abundance<-read.csv("data/data_analyses/7.dff_mites_abundance.csv") 
View(mites_df_abundance)
names(mites_df_abundance)
phylogeny_for_mites<- read.nexus("data/phylo_data/1_host_consensus_tree_mites.nex")
# Abundance is counts so we can use  a poisson but it is zero infladed (no opticon in PGLMM to take care of this)
#poisson error; fixed effect sociality=categories of variable of direct interest; random effect=foraging type
# species and elevation site has to be factors
#sociality 1, 0
#elevation as a site (as factor) several levels bamboo, lowlands, montane, high andes
mites_df_abundance <-mites_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

# Lets check for zero infalted data 
100*sum(mites_df_abundance$total_no_feathers_mites== 0)/nrow(mites_df_abundance)

# 21 % of our data is zeros( truth zeros)? i guess yes cause we collected the sample for teh individual

#### Abundance
# Modeling the individual abundances
# Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?
mites_df_abundance<-mites_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
mites_df_abundance$elevation_cat<-as.factor(mites_df_abundance$elevation_cat)
mites_df_abundance$foraging_cat<-as.factor(mites_df_abundance$foraging_cat)
mites_df_abundance$species_jetz<-as.factor(mites_df_abundance$species_jetz)
mites_df_abundance$sociality<-as.factor(mites_df_abundance$sociality)
mites_df_abundance$powder_level<-as.factor(mites_df_abundance$powder_level)

unique(mites_df_abundance$powder_level)

mean(mites_df_abundance$total_no_feathers_mites)
sd(mites_df_abundance$total_no_feathers_mites) # overdispersed variance> mean

# IMPORTANT !!!make sure traits data nad phylogeny are in teh same order
rownames(mites_df_abundance) <- mites_df_abundance$species_jetz # first make it the row names 
mites_df_abundance<- mites_df_abundance[match(phylogeny_mites$tip.label,rownames(mites_df_abundance)),]

tips<- as.data.frame(phylogeny_mites$tip.label)
names<-as.data.frame(mites_df_abundance$species_jetz)

# Modeling the data # I would prefer to use a zero inflated model however that is only aviallable in a gassioan approach bt that does no work with my model ( not sure why ye)

names( mites_df_abundance)
m_a_no_f<-  phyr::pglmm(total_no_feathers_mites~ sociality + (1|elevation_cat)  +(1|species_jetz__)+(1|powder_level), 
                        data = mites_df_abundance, 
                        family ="poisson", # use when bayes=true "zeroinflated.poisson",
                        cov_ranef = list(species_jetz=phylogeny_for_mites), #class phylo
                        #bayes = TRUE,
                        REML = TRUE, 
                        verbose = TRUE, 
                        s2.init = .25)
summary(m_a_no_f)


m_a_meso<-phyr::pglmm(total_mesostigmatidae~ sociality + (1|elevation_cat)+(1|species_jetz__)+(1|powder_level), 
                        data = mites_df_abundance, 
                        family ="poisson", # use when bayes=true "zeroinflated.poisson",
                        cov_ranef = list(species_jetz=phylogeny_for_mites), #class phylo
                        #bayes = TRUE,
                        REML = TRUE, 
                        verbose = TRUE, 
                        s2.init = .25)
summary( m_a_meso)

View(mites_df_abundance)
#### mean non feather mites 

mean_mites_mesostigmatidae<-mites_df_abundance %>% group_by (species_jetz) %>% 
  summarize(mean_mites=mean(total_mesostigmatidae)) 

species_atributes<-mites_df_abundance %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean)
species_attributes_distict<-distinct( species_atributes)

mean_mites_abundance_mesostigmatidae<-right_join(species_attributes_distict, mean_mites_mesostigmatidae, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 

mean_mites_abundance_mesostigmatidae<-mean_mites_abundance_mesostigmatidae  %>% mutate_at("species_jetz", str_replace, " ", "_")
mean_mites_abundance_mesostigmatidae$elevation_cat<-as.factor(mean_mites_abundance_mesostigmatidae$elevation_cat)
mean_mites_abundance_mesostigmatidae$foraging_cat<-as.factor(mean_mites_abundance_mesostigmatidae$foraging_cat)
mean_mites_abundance_mesostigmatidae$species_jetz<-as.factor(mean_mites_abundance_mesostigmatidae$species_jetz)
mean_mites_abundance_mesostigmatidae$sociality<-as.factor(mean_mites_abundance_mesostigmatidae$sociality)



mean_m_a_meso<-  phyr::pglmm(mean_mites~ sociality  + (1|foraging_cat)+(1|species_jetz__), 
                        data = mean_mites_abundance_mesostigmatidae, 
                        family ="gaussian", # use when bayes=true "zeroinflated.poisson",
                        cov_ranef = list(species_jetz=phylogeny_for_mites), #class phylo
                        #bayes = TRUE,
                        REML = TRUE, 
                        verbose = TRUE, 
                        s2.init = .25)
summary( mean_m_a_meso)

# [Diversity] Part 3 Ectoparasite models_ -----------------------------
###_###_###_###_
# at the species level 
###_###_###_###_

lice_df_diversity_species<-read.csv("data/7.lice_df_diversity_species.csv",header=TRUE,na.strings=c("_",""))
# Filtering only manu data  
lice_df_diversity_sp<-lice_df_diversity_species %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
# Formatting

# REMOVE EMPTY ROWS (for some samples we do not have Lice genus2 and tree but when summarizing the table the field were created)
lice_df_diversity_sp<-lice_df_diversity_sp %>% drop_na(valuecol)
View(lice_df_diversity_sp)


#richness by host species
lice_df_diversity_sp_sum<-lice_df_diversity_sp %>% group_by(species_jetz) %>% 
  summarise(richness_sp=n_distinct(valuecol))

species_atributes<-lice_df_diversity_species %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean) %>%  filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

species_attributes_distict<-distinct(species_atributes)

lice_df_richness_sp<-right_join(species_attributes_distict, lice_df_diversity_sp_sum, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 

# We are worried that diversity increases with sample size, lets check is that is the case
# summarize the number of samples per species  
n_samples_lice_sp<-lice_df_diversity_sp %>% group_by(species_jetz) %>% 
  summarise(n_samples_lice=n_distinct(Full_Label))

lice_richness_sp_df<-right_join(lice_df_richness_sp, n_samples_lice_sp, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 


#write.csv(lice_richness_sp_df, "data/5.lice_richness_sp_df_manu.csv")

####Correlation sample size and diversity #WARNING

plot(richness_sp~n_samples_lice, data=lice_richness_manu_sp)

model<-lm(richness_sp~n_samples_lice, data=lice_richness_manu_sp)
abline(model)

str(lice_diversity_sp_samplesize)

# the data

lice_richness_manu_sp<-read.csv("data/5.lice_richness_sp_df_manu.csv")
lice_richness_manu_sp<-lice_richness_manu_sp %>% distinct( species_jetz,.keep_all = TRUE)
unique(lice_richness_manu_sp$species_jetz)

phylo_lice_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex") # this is the same one that we use for lice abundance cause it has the same species

## Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?
lice_richness_manu_sp <-lice_richness_manu_sp %>% mutate_at("species_jetz", str_replace, " ", "_")
lice_richness_manu_sp$elevation_cat<-as.factor(lice_richness_manu_sp$elevation_cat)
lice_richness_manu_sp$foraging_cat<-as.factor(lice_richness_manu_sp$foraging_cat)
lice_richness_manu_sp$species_jetz<-as.factor(lice_richness_manu_sp$species_jetz)
lice_richness_manu_sp$sociality<-as.factor(lice_richness_manu_sp$sociality)


#The model pglmm

l_d_sp_pglmm<-  phyr::pglmm(richness_sp~ sociality + n_samples_lice+ (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                      data = lice_richness_manu_sp, 
                      family ="poisson", # use when bayes=true "zeroinflated.poisson",
                      cov_ranef = list(species_jetz=phylo_lice_rooted), #class phylo
                      #bayes = TRUE,
                      REML = TRUE, 
                      verbose = TRUE, 
                      s2.init = .25)
summary(l_d_sp_pglmm)
rr2::R2(l_d_sp_pglmm)

# SAMPLE SIZE AS a random effect.
l_d_sp_pglmm<-  phyr::pglmm(richness_sp~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__) + (1|n_samples_lice),
                            data = lice_richness_manu_sp, 
                            family ="poisson", # use when bayes=true "zeroinflated.poisson",
                      cov_ranef = list(species_jetz=phylo_lice_rooted), #class phylo
                      #bayes = TRUE,
                      REML = TRUE, 
                      verbose = TRUE, 
                      s2.init = .25)
summary(l_d_sp_pglmm)




l_d_sp_model_glmm<-  lme4::glmer(richness_sp~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz) +(1|n_samples_lice), 
                                 data = lice_diversity_sp_samplesize, 
                                 family=poisson(), # use when bayes=true "zeroinflated.poisson",
                                 #REML = TRUE, 
                                 verbose = TRUE)
summary(l_d_sp_model_glmm)



# [Infestation] Part 4 ----------------------------------------------------


plotTree.datamatrix # To plot two or more discrete characters withthe phylogeny see page 169 of liam's revel book, sociality and infested yes or no


#[Network] general calculation] ----------------------------------------




#[Network] Degree and ectoparasite abundance ----------------------------------------

#The data MEAN

degree_df<-read.csv("data/8.network_outputdegree_network_all_sp_manu_jetz_tax.csv")

mean_lice_abundance<-read.csv("data/5.lice_df_abundance_means.csv")
mean_lice_abundance<-mean_lice_abundance %>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations

mean_lice_abundance %>% filter( sociality=="1") # we have abundance for 62 host social species and for 27 non social species  in mau

phylo_lice_rooted<-read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")

# Reformating 
degree_df<-degree_df%>% mutate_at("species_jetz", str_replace, " ", "_")
degree_df$deg_binary <-as.numeric(degree_df$deg_binary )

#mean_lice_abundance <-mean_lice_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
mean_lice_abundance$elevation_cat<-as.factor(mean_lice_abundance$elevation_cat)
mean_lice_abundance$foraging_cat<-as.factor(mean_lice_abundance$foraging_cat)
mean_lice_abundance$species_jetz<-as.factor(mean_lice_abundance$species_jetz)
mean_lice_abundance$sociality<-as.factor(mean_lice_abundance$sociality)

sociality_continous_mean_lice<- inner_join(degree_df, mean_lice_abundance, by="species_jetz") 

#sociality_continous_mean_lice<-sociality_continous_mean_lice %>% rename("species_jetz"="species_taxonomy_SACC_2021")
sociality_continous_mean_lice$species_jetz<-as.factor(sociality_continous_mean_lice$species_jetz)

View(sociality_continous_mean_lice)
str(sociality_continous_mean_lice )
str( degree_df)

#The model  PGLS 9 SINCE THE MEAN IS  CONTIINOUS VARIBLE)??
l_abun_mean_degree<-  phyr::pglmm(mean_lice ~ deg_binary +(1|elevation_cat)+(1|species_jetz__)+ (1|foraging_cat),
                                  data = sociality_continous_mean_lice, 
                                  family = "poisson", 
                                  cov_ranef = list(species_jetz= phylo_lice_rooted), #class phylo
                                  #bayes = TRUE,
                                  REML= TRUE,  # NOT SURE WHEN TO USE ML
                                  verbose = TRUE,
                                  s2.init = .25) # what is this last parameter for
summary(l_abun_mean_degree)

rr2::R2(l_abun_mean_degree)


ggplot(sociality_continous_mean_lice, aes(x=deg_binary, y=mean_lice)) +
  geom_point()+
  geom_jitter(height = 0.01)+
  labs(title="b) Degree")+
  theme_classic(20)

# INDIVIDUAL abundance and binomial degree

degree_df<-read.csv("data/8.network_outputdegree_network_all_sp_manu_jetz_tax.csv")
lice_df_abundance<-read.csv("data/5.lice_df_abundance_manu.csv")

#lice_df_abundance<-read.csv("data/7.lice_df_abundance.csv")
# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
#lice_df_abundance<-lice_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
#write.csv(lice_df_abundance, "data/5.lice_df_abundance_manu.csv")

lice_df_abundance %>% filter(sociality=="1") # for 63 individuals
unique (lice_df_abundance$species_jetz)  # check that number coinide with the phylogeny

phylo_lice_rooted_abund_networks_manu<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abunmanu_networks.nex")
str(phylo_lice_rooted_abund_networks_manu)
phylo_lice_rooted_abund_networks_manu$tip.label
# Make sure number in teh data and the phylogenty are consistent

# Reformating 
#lice_df_abundance <-lice_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
lice_df_abundance$elevation_cat<-as.factor(lice_df_abundance$elevation_cat)
lice_df_abundance$foraging_cat<-as.factor(lice_df_abundance$foraging_cat)
lice_df_abundance$species_jetz<-as.factor(lice_df_abundance$species_jetz)
lice_df_abundance$sociality<-as.factor(lice_df_abundance$sociality)

degree_df<-degree_df%>% mutate_at("species_jetz", str_replace, " ", "_")
degree_df$deg_binary <-as.numeric(degree_df$deg_binary )

sociality_continous_abundance_lice<- inner_join(degree_df, lice_df_abundance, by="species_jetz") 
unique(sociality_continous_abundance_lice$species_jetz)

names(sociality_continous_abundance_lice)

l_abun_ind_degree<-  phyr::pglmm(total_lice ~deg_binary +(1|elevation_cat)+(1|species_jetz__)+ (1|foraging_cat),
                                  data = sociality_continous_abundance_lice,
                                  family = "poisson", 
                                  cov_ranef = list(species_jetz=phylo_lice_rooted_abund_networks_manu), #class phylo
                                  #bayes = TRUE,
                                  REML= FALSE,  # NOT SURE WHEN TO USE ML
                                  verbose = TRUE,
                                  s2.init = .25) # what is this last parameter for
rr2::R2(l_abun_ind_degree)

ggplot(sociality_continous_abundance_lice, aes(x=deg_binary, y=total_lice)) +
  geom_point()+
  geom_smooth(method = lm)
  geom_jitter(height = 0.01)+
  labs(title="b) Degree per ind ")+
  theme_classic(20)
  

# DEGREE weighted 

degree_w_df<-read.csv("data/8.network_degree_weighted_network_all_sp_manu_jetz_tax.csv")%>% mutate_at("species_jetz", str_replace, " ", "_")
str(degree_w_df)
degree_w_df$deg_weighted<-as.numeric(degree_w_df$deg_weighted)

mean_lice_abundance<-read.csv("data/5.lice_df_abundance_means.csv")

mean_lice_abundance<-mean_lice_abundance %>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations
#mean_lice_abundance %>% filter(sociality=="1") # we have 62 soacials speci

phylo_lice_rooted<-read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")

#mean_lice_abundance <-mean_lice_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
mean_lice_abundance$elevation_cat<-as.factor(mean_lice_abundance$elevation_cat)
mean_lice_abundance$foraging_cat<-as.factor(mean_lice_abundance$foraging_cat)
mean_lice_abundance$species_jetz<-as.factor(mean_lice_abundance$species_jetz)
mean_lice_abundance$sociality<-as.factor(mean_lice_abundance$sociality)

sociality_continous_w_mean_lice<- inner_join(degree_w_df, mean_lice_abundance, by="species_jetz") 
sociality_continous_w_mean_lice$species_jetz<-as.factor(sociality_continous_w_mean_lice$species_jetz)

#the model

l_abun_mean_w_degree<-  phyr::pglmm(mean_lice ~ deg_weighted +(1|elevation_cat)+(1|species_jetz__)+ (1|foraging_cat),
                                  data = sociality_continous_w_mean_lice, 
                                  family = "poisson", 
                                  cov_ranef = list(species_jetz= phylo_lice_rooted), #class phylo
                                  #bayes = TRUE,
                                  REML= TRUE,  # NOT SURE WHEN TO USE ML
                                  verbose = TRUE,
                                  s2.init = .25) # what is this last parameter for


rr2::R2(l_abun_mean_degree)

ggplot(sociality_continous_w_mean_lice, aes(x=deg_weighted, y=mean_lice)) +
  geom_point()+
  geom_jitter(height = 0.01)+
  labs(title="b) w_Degree")+
  theme_classic(20)



#[Network] Degree and ectoparasite diversity ----------------------------------------
degree_df<-read.csv("data/5.degree_network_all_sp_manu.csv")
degree_df <-degree_df%>% mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_")
degree_df$degree<-as.numeric(degree_df$degree)

lice_richness_manu_sp<-read.csv("data/5.lice_richness_sp_df_manu.csv")
lice_richness_manu_sp<-lice_richness_manu_sp %>% distinct( species_jetz,.keep_all = TRUE)
unique(lice_richness_manu_sp$species_jetz)

phylo_lice_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex") # this is the same one that we use for lice abundance cause it has the same species

#

lice_richness_manu_sp<-lice_richness_manu_sp %>% mutate_at("species_jetz", str_replace, " ", "_")
lice_richness_manu_sp$elevation_cat<-as.factor(lice_richness_manu_sp$elevation_cat)
lice_richness_manu_sp$foraging_cat<-as.factor(lice_richness_manu_sp$foraging_cat)
lice_richness_manu_sp$species_jetz<-as.factor(lice_richness_manu_sp$species_jetz)
lice_richness_manu_sp$sociality<-as.factor(lice_richness_manu_sp$sociality)

# joining files
sociality_continous_diversity_lice<- right_join(degree_df, lice_richness_manu_sp, by=c("species_taxonomy_SACC_2021"="species_jetz")) 
sociality_continous_diversity_lice<-sociality_continous_diversity_lice %>% rename("species_jetz"="species_taxonomy_SACC_2021")
sociality_continous_diversity_lice$species_jetz<-as.factor(sociality_continous_diversity_lice$species_jetz)

l_diver_mean_degree<-  phyr::pglmm( richness_sp ~ degree + n_samples_lice+(1|elevation_cat)+(1|species_jetz__)+ (1|foraging_cat),
                                  data = sociality_continous_diversity_lice, 
                                  family = "poisson", 
                                  cov_ranef = list(species_jetz= phylo_lice_rooted), #class phylo
                                  #bayes = TRUE,
                                  REML= TRUE,  # NOT SURE WHEN TO USE ML
                                  verbose = TRUE,
                                  s2.init = .25) # what is this last parameter for


ggplot(sociality_continous_diversity_lice, aes(x=degree, y=richness_sp )) +
  geom_point()+
  geom_jitter(height = 0.01)+
  geom_smooth()+
  labs(title="b) Degree")+
  theme_classic(20)

ggplot(sociality_continous_diversity_lice, aes(x=n_samples_lice, y=richness_sp )) +
  geom_point()+
  geom_jitter(height = 0.01)+
  geom_smooth(method = lm)+
  labs(title="b) Degree")+
  theme_classic(20)


# Crazy idea


##_###_###_###_###_###_###_###_###_###_###_###_###_###_
# Flocks
#Preparing flock data 
flocks<-read.csv ("data/0.flocks_manu_complete_18052022.csv")
flocks<-flocks %>% filter(database_decision=="include") 
#Look for duplicates there are 755 observations, 4 repited entrees
{flocks} %>%
  dplyr::group_by( flock_id, species_taxonomy_SACC_2021) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 
#flocks %>% group_by( flock_id, species_taxonomy_SACC_2021) %>%filter(n() > 1)

#eliminate duplicates 
flocks<-flocks %>% distinct( flock_id, species_taxonomy_SACC_2021, .keep_all = TRUE)
unique(flocks$flock_id)
# create list of flocking species
flocks_list<-flocks %>% 
  distinct( species_taxonomy_SACC_2021, keep_all=FALSE) %>% 
  add_column(flocking="yes") %>% 
  select(-keep_all)
View(flocks)

names(flocks)


# Number of individual per species

conspecific_groups_mean<-flocks %>% group_by(species_taxonomy_SACC_2021) %>% 
  summarize (mean_conspecific=mean(ind_per_sp)) 

conspecific_groups_max<-flocks %>% group_by(species_taxonomy_SACC_2021) %>% 
  summarize (max_conspecific=max(ind_per_sp)) 

conspecifics<-full_join (conspecific_groups_mean,conspecific_groups_max, by="species_taxonomy_SACC_2021")

as.data.frame(conspecifics)

View(conspecifics)

hist(conspecific_groups$mean_conspecific)
# make sure the formating is consistent

conspecifics<-conspecifics %>% mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_")
# Joining with jettx taxonomy 

jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv")%>% 
  mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_") %>% 
  select(species_taxonomy_SACC_2021,species,species_jetz,TipLabel) %>% 
  distinct(species_taxonomy_SACC_2021,.keep_all = TRUE )


conspecific_jetz<-left_join (conspecifics, jetz_taxonomy_manu_only, by="species_taxonomy_SACC_2021")
conspecific_jetz<-conspecific_jetz %>% mutate_at("species_jetz", str_replace, " ", "_")



View(conspecific_jetz)

lice_df_abundance_means<-read.csv("data/5.lice_df_abundance_means.csv")

lice_df_abundance_means_conspecifics<-inner_join (conspecific_jetz, lice_df_abundance_means, by="species_jetz")

names (lice_df_abundance_means_conspecifics)
ggplot(lice_df_abundance_means_conspecifics, aes(x = mean_lice, y=max_conspecific)) +
  geom_point(alpha=0.5)+
  geom_smooth(method=lm)
ggplot(lice_df_abundance_means_conspecifics, aes(x = mean_lice, y=mean_conspecific)) +
  geom_point(alpha=0.5)+
  geom_smooth(method=lm)


lice_df_abundance<-read.csv("data/5.lice_df_abundance_manu.csv")

lice_df_abundance_conspecifics<-inner_join (conspecific_jetz, lice_df_abundance, by="species_jetz")

ggplot(lice_df_abundance_conspecifics, aes(x =total_lice, y=mean_conspecific)) +
  geom_point(alpha=0.5)+
  geom_smooth(method=lm)


mites_df_abundance<-read.csv("data/7.mites_df_abundance.csv")
mites_df_abundance<-mites_df_abundance %>% mutate_at("species_jetz", str_replace, " ", "_")

mites_df_abundance_conspecifics<-inner_join (conspecific_jetz, mites_df_abundance, by="species_jetz")



ggplot(mites_df_abundance_conspecifics, aes(x =total_mites, y=max_conspecific)) +
  geom_point(alpha=0.5)+
  geom_smooth(method=lm)

write.csv(degree_manu_jetz, "data/8.network_outputdegree_network_all_sp_manu_jetz_tax.csv")
write.csv(degree_w_manu_jetz, "data/8.network_degree_weighted_network_all_sp_manu_jetz_tax.csv")



View(deg_binary)
