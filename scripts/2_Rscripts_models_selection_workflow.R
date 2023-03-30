#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models workflow                                                                       ###
### R-code                                                                          ###
### Jenny Munoz      
###
### Last update: March 2023                                                ###
################################################################################
# 0.Libraries  --------------------------------------------------------------

# libraries for easier manipulation of data
#install.packages("tidyr") 
install.packages("tidyverse") 
install.packages("dplyr")
install.packages ("data.table")
install.packages ("extrafont")
installed.packages("lubridate")  #for dates
#data cleaning
install.packages ("janitor")
install.packages("assertr")

#Other libraries for data analyses and visualizations
install.packages("vegan")
install.packages("ggplot2")
install.packages("devtools")
install.packages("knitr")
install.packages("ts")
install.packages("RColorBrewer")
install.packages("ggridges")
install.packages("ggtree")
install.packages("aplot")

#for models
install.packages("car") #Anova command
install.packages("lattice") #preliminary plots
install.packages("lme4") #for glmer (generalized linear mixed models) 
install.packages("visreg")  #extract confidence intervals and trend lines from GLMMs
install.packages("lsmeans") #least squared means
install.packages("MuMIn") #pseudo R squared for GLMMs
install.packages("emmeans")
#model assumptions
install.packages("DHARMa")

# bayes models 
install.packages('tidybayes')
install.packages ('bayesplot')
install.packages('brms') # bayesian approach to model phylogenetic data with repides observations
#install.packages('brmstools')
remotes::install_github("Pakillo/DHARMa.helpers")
devtools::install_github("mvuorre/brmstools")


# install.packages("remotes") # DHARMA FOR BRMS 

# Phylogenetic component
install.packages("ape")
#install.packages("here")
install.packages("phytools")
install.packages("tidyverse")
install.packages("metafor")
install.packages("phangorn") # to reconstruct a maximum clade credibility tree
install.packages("rr2")
install.packages ( "MCMCglmm")
#install.packages("phyr")
remotes::install_github("daijiang/phyr")
#options(repos = c(
# phyr = 'https://daijiang.r-universe.dev',
# CRAN = 'https://cloud.r-project.org'))
#install.packages('phyr')
install.packages("TreeTools")

# Libraries 
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
library(assertr)
# for data visualization
library(vegan)
library(ggplot2)
library(devtools)
library(knitr)
library(ts)
library(RColorBrewer)
library(ggridges)
library(ggtree)
library(aplot)

#libraries for models and visualizations
library(lattice) #preliminary plots
library(car) #Anova command
library(lsmeans) #least squared means
library(lme4) #for glmer (generalized linear mixed models) 
library(visreg) #extract confidence intervals and trend lines from GLMMs
library(MuMIn) #pseudo R squared for GLMMs
library(emmeans)
#mode assumption check
library (DHARMa)
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
library(TreeTools)
#Libraries for plots
library(gridExtra)
library(ggpubr)
library(grid)
#libaries for bayes 
library(bayesplot)
library(tidybayes)
library(brms) # bayesian approach to model phylogenetic data with repides observations
library(DHARMa.helpers)
library(brmstools)


# ##### 1.Data processing prevalence ectos ----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality ) %>% 
  na.omit() %>% 
  filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # Need to prunne the tree 


# Data structure
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$year_seasonality<-as.numeric(ectos_birds_dff$year_seasonality)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$ectoparasites_PA<-as.numeric(ectos_birds_dff$ectoparasites_PA)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one

names(ectos_birds_dff)
is.ultrametric(phylo)

# Make sure the tips and the names on the file coincide and formatting of name is consistent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_birds_dff$species_jetz)) %>% mutate(name=ectos_birds_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 
phy_cov<-ape::vcv(phylo, corr=TRUE)

# ##### 1.2.Model selection prevalence ectos --------------------------------------------------------
###_###_###

ecto_p_brms_bayes_no_int<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
                               (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                               (1|Powder.lvl)+
                               (1|species),
                             data=ectos_birds_dff,
                             family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                             data2 = list(phy_cov=phy_cov),
                             iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                             thin=2,
                             control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(ecto_p_brms_bayes_no_int, "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")
#ecto_p_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")

ecto_p_brms_bayes_sociality_interactions<-brms::brm(ectoparasites_PA~
                                                sociality+
                                                scale(elevation)+
                                                scale(year_seasonality)+
                                                sociality:scale(elevation)+
                                                sociality:scale(year_seasonality)+
                                                (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                (1|Powder.lvl)+
                                                (1|species),
                                              data=ectos_birds_dff,
                                              family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                              data2 = list(phy_cov=phy_cov),
                                              iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                              thin=2,
                                              control=list(adapt_delta=0.99, max_treedepth=12)) 
saveRDS(ecto_p_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions.RDS")
ecto_p_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions.RDS")

ecto_p_brms_bayes_all_interactions<-brms::brm(ectoparasites_PA~
                               sociality+
                               scale(elevation)+
                               scale(year_seasonality)+
                               sociality:scale(elevation)+
                               sociality:scale(year_seasonality)+
                               scale(elevation):scale(year_seasonality)+
                               sociality:scale(year_seasonality):scale(elevation)+
                               (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                               (1|Powder.lvl)+
                               (1|species),
                             data=ectos_birds_dff,
                             family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                             data2 = list(phy_cov=phy_cov),
                             iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                             thin=2,
                             control=list(adapt_delta=0.99, max_treedepth=12)) 

#saveRDS(ecto_p_brms_bayes_all_interactions, "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")
#ecto_p_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")

loo(ecto_p_brms_bayes_no_int,ecto_p_brms_bayes_sociality_interactions, ecto_p_brms_bayes_all_interactions,compare=TRUE)

summary (ecto_p_brms_bayes)
fixef() # to get more detailed values for estimates
coef(ecto_p_brms_bayes) # if you have group-level effects (hierarchical data)
bayes_R2(ecto_p_brms_bayes_all_interactions) # R2 0.1529

plot(ecto_p_brms_bayes) # paarameter distributio and convergence
mcmc_plot(ecto_p_brms_bayes_no_int) # Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict(ecto_p_brms_bayes)
ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 5), ylim = c(0,30))

# ##### 1.Data processing abundance lice ----------------------------------------------------

ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, total_lice,year_seasonality ) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens") 
#filter(total_lice<60)  # removing outliers 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

# Finding putliers
#ectos_birds_dff <-ectos_birds_dff  %>% mutate_at("species_jetz", str_replace, " ", "_")
#ectos_birds_dff$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$total_lice<-as.numeric(ectos_birds_dff$total_lice)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one
names(ectos_birds_dff)
is.ultrametric(phylo)

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_birds_dff$species_jetz)) %>% mutate(name=ectos_birds_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 

#Phylogenetic covariance matrix
phy_cov<-ape::vcv(phylo, corr=TRUE)

# ##### 2.2.Model selection abundance lice --------------------------------------------------------

zip_a_lice_brms_bayes_no_int<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                      (1|Powder.lvl)+
                                      (1|species),
                                    data=ectos_birds_dff,
                                    family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                    data2 = list(phy_cov=phy_cov),
                                    iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                    thin=2,
                                    control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(zip_a_lice_brms_bayes_no_int, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS")
zip_a_lice_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS")

zip_a_lice_brms_bayes_sociality_interactions<-brms::brm(total_lice~
                                                      sociality+
                                                      scale(elevation)+
                                                      scale(year_seasonality)+
                                                      sociality:scale(elevation)+
                                                      sociality:scale(year_seasonality)+
                                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                      (1|Powder.lvl)+
                                                      (1|species),
                                                    data=ectos_birds_dff,
                                                    family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                    data2 = list(phy_cov=phy_cov),
                                                    iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                    thin=2,
                                                    control=list(adapt_delta=0.99, max_treedepth=12)) 
saveRDS(zip_a_lice_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")
zip_a_lice_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms__LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")

zip_a_lice_brms_bayes_all_interactions<-brms::brm(total_lice~
                                                sociality+
                                                scale(elevation)+
                                                scale(year_seasonality)+
                                                sociality:scale(elevation)+
                                                sociality:scale(year_seasonality)+
                                                scale(elevation):scale(year_seasonality)+
                                                sociality:scale(year_seasonality):scale(elevation)+
                                                (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                (1|Powder.lvl)+
                                                (1|species),
                                              data=ectos_birds_dff,
                                              family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                              data2 = list(phy_cov=phy_cov),
                                              iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                              thin=2,
                                              control=list(adapt_delta=0.99, max_treedepth=12)) 

#saveRDS(ecto_p_brms_bayes_all_interactions, "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")
#ecto_p_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")


mcmc_plot(zip_a_lice_brms_bayes_no_int)

