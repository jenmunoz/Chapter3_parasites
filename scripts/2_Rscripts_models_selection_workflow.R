#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models selection workflow                                                                       ###
### R-code                                                                          ###
### Jenny Munoz      
###
### Last update: March 2023                                                ###
################################################################################
#As the proportion of zeros is quite high in the data, it is worthwhile to test also a zero-inflated negative-binomial model, 
#which is a mixture of two models - logistic regression to model the proportion of extra zero counts - negative-binomial model
#Notes check this blog for discussion in STAN BRMS situations: https://discourse.mc-stan.org
# check this blog for discussion in model selection https://discourse.mc-stan.org/t/model-selection-in-brms/30492
# about loo package interpretation https://mc-stan.org/loo/reference/loo-glossary.html
# some examples here : https://avehtari.github.io/modelselection/roaches.html

# Notes for model selection with BRMS, [Work in progress]
# definition [elpd_loo] The ELPD is the theoretical expected log pointwise predictive density for a new dataset 
# definition [elpd_loo SE] This standard error is a coarse description of our uncertainty about the predictive performance for unknown future data. 
# In well behaving cases p_loo < N and p_loo < p, where p is the total number of parameters in the model and n ( number of observations)
#elpd_loo differences of more than 4 magnitude are consider lareg enough to make a difference in teh model
# elpd_loo large compared to se_diff indicate changes in the model selection too
# notes more model overfiiting in BRMS 
# Look at the number of observations (226 is shown in the first line), 
# The effective number of parameters, p_loo ()  p_loo < N (number of observations) and p_loo < p (parameters of the model)
# and the number of high Pareto k’s ( All Pareto k’s are good, so there are no highly influential individual observations, This indicates that loo computation is reliable)

# also important to check teh priors prior_summary()
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

#MODEL COMPARISON
#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.
loo(m1,m2,m3, compare=TRUE)

#PLOTS
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

# model convergence 
png("data/data_analyses/models/model_plots/M.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int)
dev.off()

# model fit
png("data/data_analyses/models/model_plots/M.png",width = 3000, height = 3000, res = 300, units = "px")
pp_m<- brms::posterior_predict(MODEL)
ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
dev.off()

pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

color_scheme_set("blue")

estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions Ectos prevalence", subtitle ="Ectos prevalence with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50")+
  xlab("Estimate")

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions Ectos prevalence", subtitle ="Ectos prevalence with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig2.plot_model_parameters_intervals_ecto_PREVALENCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

png("data/data_analyses/model_selection/plots_model_selected/Fig2.plot_model_parameters_ecto_PREVALENCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

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

saveRDS(zip_a_lice_brms_bayes_all_interactions, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")
zip_a_lice_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")

loo(zip_a_lice_brms_bayes_all_interactions,zip_a_lice_brms_bayes_sociality_interactions,zip_a_lice_brms_bayes_no_int)



#MODEL COMPARISON
loo(m1,m2,m3, compare=TRUE)

#PLOTS
color_scheme_set("teal")

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

# model convergence 
png("data/data_analyses/model_selection/plots_model_selected/Fig2.plot_model_CONVERGENCE_intervals_ecto_PREVALENCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int)
dev.off()

# model fit
png("data/data_analyses/model_selection/plots_model_selected/Fig2.plot_model_CONVERGENCE_intervals_ecto_PREVALENCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
pp_m<- brms::posterior_predict(MODEL)
ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
dev.off()

pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(zip_a_lice_brms_bayes_no_int,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions LICE_ABUNDANCE", subtitle ="LICE_ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50")+
  xlab("Estimate")

estimates_plot_intervals<-mcmc_plot(zip_a_lice_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions LICE_ABUNDANCE", subtitle ="LICE_ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig2l.plot_model_parameters_intervals_LICE_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

png("data/data_analyses/model_selection/plots_model_selected/Fig2l.plot_model_parameters_LICE_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()



# ##### 3.1.Data processing abundance mites ----------------------------------------------------

ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality,total_mites, total_mesostigmatidae, total_no_feathers_mites,year_seasonality ) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens") 
#filter(total_lice<60)  # removing outliers 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$total_mites<-as.numeric(ectos_birds_dff$total_mites)
ectos_birds_dff$total_mesostigmatidae<-as.numeric(ectos_birds_dff$total_mesostigmatidae)
ectos_birds_dff$total_no_feathers_mites<-as.numeric(ectos_birds_dff$total_no_feathers_mites)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one

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

# ##### 2.2.Model selection abundance mites --------------------------------------------------------

loop<-loo(zip_a_nf_mites_brms_bayes_no_int)
zip_a_nf_mites_brms_bayes_no_int<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                          (1|Powder.lvl)+
                                          (1|species),
                                        data=ectos_birds_dff,
                                        family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                        data2 = list(phy_cov=phy_cov),
                                        iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                        thin=2,
                                        control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(zip_a_nf_mites_brms_bayes_no_int, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")
zip_a_nf_mites_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")

mcmc_plot(zip_a_nf_mites_brms_bayes_sociality_interactions)
zip_a_nf_mites_brms_bayes_sociality_interactions<-brms::brm(total_no_feathers_mites~
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
#saveRDS(zip_a_nf_mites_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")
zip_a_nf_mites_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")

zip_a_nf_mites_brms_bayes_all_interactions<-brms::brm(total_no_feathers_mites~
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

#saveRDS(zip_a_nf_mites_brms_bayes_all_interactions, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")
zip_a_nf_mites_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")

loo(zip_a_lice_brms_bayes_all_interactions,zip_a_lice_brms_bayes_sociality_interactions,zip_a_lice_brms_bayes_no_int)

mcmc_plot(zip_a_lice_brms_bayes_no_int)


# ##### 1.Data processing NETWORKS prevalence ectos ----------------------------------------------------

#seasonality_date<-read.csv("data/data_analyses/data_manuscript/7.dff_parasites_seasonality.csv")
#dff_ectos_network<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>%
#inner_join(read.csv("data/data_analyses/data_manuscript/7.dff_parasites_seasonality.csv"), by="Full_Label")
# write.csv(dff_ectos_network,"data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv" )


dff_ectos_network_individual_metrics<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
  select(elevation_extrapolated_date, species_jetz, Powder.lvl,foraging_cat, sociality,ectoparasites_PA, degree, w_degree, year_seasonality) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

unique(dff_ectos_network_individual_metrics$species_jetz) # this is teh total species that are in flocks taht we have samples for

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos social and non social so we need to trim it 

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(dff_ectos_network_individual_metrics$species_jetz)) %>% mutate(name=dff_ectos_network_individual_metrics$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 

# phylogenetic correlation structure, create a covariance matrix of species
phy_cov<-ape::vcv(phylo, corr=TRUE)

# data structure 

#dff_ectos_network_individual_metrics$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
dff_ectos_network_individual_metrics$foraging_cat<-as.factor(dff_ectos_network_individual_metrics$foraging_cat)
dff_ectos_network_individual_metrics$species_jetz<-as.factor(dff_ectos_network_individual_metrics$species_jetz)
dff_ectos_network_individual_metrics$elevation<-as.numeric(dff_ectos_network_individual_metrics$elevation)
dff_ectos_network_individual_metrics$degree<-as.numeric(dff_ectos_network_individual_metrics$degree)
dff_ectos_network_individual_metrics$w_degree<-as.numeric(dff_ectos_network_individual_metrics$w_degree)
dff_ectos_network_individual_metrics$year_seasonality<-as.numeric(dff_ectos_network_individual_metrics$year_seasonality)

#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
dff_ectos_network_individual_metrics$sociality<-as.factor(dff_ectos_network_individual_metrics$sociality)
dff_ectos_network_individual_metrics$Powder.lvl<-as.factor(dff_ectos_network_individual_metrics$Powder.lvl)
dff_ectos_network_individual_metrics$ectoparasites_PA<-as.numeric(dff_ectos_network_individual_metrics$ectoparasites_PA)
dff_ectos_network_individual_metrics$species<-dff_ectos_network_individual_metrics$species_jetz # create a column for the species effect different to the phylogenetic one

names(dff_ectos_network_individual_metrics)
is.ultrametric(phylo)
LICE 



# ##### 1.2.Model selection NETWORKS prevalence ectos --------------------------------------------------------
###_###_###
mcmc_plot(ecto_p_brms_bayes_no_int_degree)
ecto_p_brms_bayes_no_int_degree<-brms::brm(ectoparasites_PA~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                      (1|Powder.lvl)+
                                      (1|species),
                                    data=dff_ectos_network_individual_metrics,
                                    family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                    data2 = list(phy_cov=phy_cov),
                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                    thin=2,
                                    control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(ecto_p_brms_bayes_no_int_degree, "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE.RDS")
ecto_p_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE.RDS")

mcmc_plot(ecto_p_brms_bayes_sociality_interactions_degree)

ecto_p_brms_bayes_sociality_interactions_degree<-brms::brm(ectoparasites_PA~
                                                      scale(degree)+
                                                      scale(elevation)+
                                                      scale(year_seasonality)+
                                                      scale(degree):scale(elevation)+
                                                      scale(degree):scale(year_seasonality)+
                                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                      (1|Powder.lvl)+
                                                      (1|species),
                                                    data=dff_ectos_network_individual_metrics,
                                                    family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                    data2 = list(phy_cov=phy_cov),
                                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                    thin=2,
                                                    control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(ecto_p_brms_bayes_sociality_interactions_degree, "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")
ecto_p_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")

ecto_p_brms_bayes_all_interactions_degree<-brms::brm(ectoparasites_PA~
                                                scale(degree)+
                                                scale(elevation)+
                                                scale(year_seasonality)+
                                                scale(degree):scale(elevation)+
                                                scale(degree):scale(year_seasonality)+
                                                scale(elevation):scale(year_seasonality)+
                                                scale(degree):scale(year_seasonality):scale(elevation)+
                                                (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                (1|Powder.lvl)+
                                                (1|species),
                                              data=dff_ectos_network_individual_metrics,
                                              family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                              data2 = list(phy_cov=phy_cov),
                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                              thin=2,
                                              control=list(adapt_delta=0.99, max_treedepth=12)) 

saveRDS(ecto_p_brms_bayes_all_interactions_degree, "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree.RDS")
ecto_p_brms_bayes_all_interactions_degree<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree.RDS")

loo(ecto_p_brms_bayes_no_int_degree,ecto_p_brms_bayes_sociality_interactions_degree,ecto_p_brms_bayes_all_interactions_degree,compare=TRUE)

summary (ecto_p_brms_bayes)
fixef() # to get more detailed values for estimates
coef(ecto_p_brms_bayes) # if you have group-level effects (hierarchical data)
bayes_R2(ecto_p_brms_bayes_sociality_interactions) # R2 0.1529

plot(ecto_p_brms_bayes) # paarameter distributio and convergence
mcmc_plot(ecto_p_brms_bayes_sociality_interactions_degree) # Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict(ecto_p_brms_bayes)
ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 5), ylim = c(0,30))

# ##### 1.Data processing NETWORKS abundance lice ----------------------------------------------------


# Lice

dff_ectos_network_individual_metrics<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
  select(elevation_extrapolated_date, species_jetz, Powder.lvl,foraging_cat, sociality,total_lice, degree, w_degree, year_seasonality) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

unique(dff_ectos_network_individual_metrics$species_jetz) # this is teh total species that are in flocks taht we have samples for

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos social and non social so we need to trim it 

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(dff_ectos_network_individual_metrics$species_jetz)) %>% mutate(name=dff_ectos_network_individual_metrics$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 

# phylogenetic correlation structure, create a covariance matrix of species
phy_cov<-ape::vcv(phylo, corr=TRUE)

# data structure 

#dff_ectos_network_individual_metrics$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
dff_ectos_network_individual_metrics$foraging_cat<-as.factor(dff_ectos_network_individual_metrics$foraging_cat)
dff_ectos_network_individual_metrics$species_jetz<-as.factor(dff_ectos_network_individual_metrics$species_jetz)
dff_ectos_network_individual_metrics$elevation<-as.numeric(dff_ectos_network_individual_metrics$elevation)
dff_ectos_network_individual_metrics$degree<-as.numeric(dff_ectos_network_individual_metrics$degree)
dff_ectos_network_individual_metrics$w_degree<-as.numeric(dff_ectos_network_individual_metrics$w_degree)
dff_ectos_network_individual_metrics$year_seasonality<-as.numeric(dff_ectos_network_individual_metrics$year_seasonality)

#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
dff_ectos_network_individual_metrics$sociality<-as.factor(dff_ectos_network_individual_metrics$sociality)
dff_ectos_network_individual_metrics$Powder.lvl<-as.factor(dff_ectos_network_individual_metrics$Powder.lvl)
dff_ectos_network_individual_metrics$total_lice<-as.numeric(dff_ectos_network_individual_metrics$total_lice)
dff_ectos_network_individual_metrics$species<-dff_ectos_network_individual_metrics$species_jetz # create a column for the species effect different to the phylogenetic one
names(dff_ectos_network_individual_metrics)
is.ultrametric(phylo)


# ##### 2.2.Model selection NETWORKS abundance lice --------------------------------------------------------

zip_a_lice_brms_bayes_no_int_degree<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                          (1|Powder.lvl)+
                                          (1|species),
                                        data=dff_ectos_network_individual_metrics,
                                        family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                        data2 = list(phy_cov=phy_cov),
                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                        thin=2,
                                        control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(zip_a_lice_brms_bayes_no_int_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")
zip_a_lice_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")

zip_a_lice_brms_bayes_sociality_interactions_degree<-brms::brm(total_lice~
                                                          scale(degree)+
                                                          scale(elevation)+
                                                          scale(year_seasonality)+
                                                          scale(degree):scale(elevation)+
                                                          scale(degree):scale(year_seasonality)+
                                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                          (1|Powder.lvl)+
                                                          (1|species),
                                                        data=dff_ectos_network_individual_metrics,
                                                        family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                        data2 = list(phy_cov=phy_cov),
                                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                        thin=2,
                                                        control=list(adapt_delta=0.99, max_treedepth=12)) 
prior_summary(model)

#saveRDS(zip_a_lice_brms_bayes_sociality_interactions_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")
zip_a_lice_brms_bayes_sociality_interactions_degree<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms__LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")

zip_a_lice_brms_bayes_all_interactions_degree<-brms::brm(total_lice~
                                                    scale(degree)+
                                                    scale(elevation)+
                                                    scale(year_seasonality)+
                                                    scale(degree):scale(elevation)+
                                                    scale(degree):scale(year_seasonality)+
                                                    scale(elevation):scale(year_seasonality)+
                                                    scale(degree):scale(year_seasonality):scale(elevation)+
                                                    (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                    (1|Powder.lvl)+
                                                    (1|species),
                                                  data=dff_ectos_network_individual_metrics,
                                                  family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                  data2 = list(phy_cov=phy_cov),
                                                  iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                  thin=2,
                                                  control=list(adapt_delta=0.99, max_treedepth=12)) 

saveRDS(zip_a_lice_brms_bayes_all_interactions_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_DEGREE.RDS")
zip_a_lice_brms_bayes_all_interactions_degree<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_DEGREE.RDS")

loo(zip_a_lice_brms_bayes_sociality_interactions,zip_a_lice_brms_bayes_no_int)

# ##### 3.1.Data processing NETWORKS abundance mites ----------------------------------------------------

# Mites
dff_ectos_network_individual_metrics<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
  select(elevation_extrapolated_date, species_jetz, Powder.lvl,foraging_cat, sociality,total_mites, total_no_feathers_mites, degree, w_degree, year_seasonality) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites


unique(dff_ectos_network_individual_metrics$species_jetz) # this is teh total species that are in flocks taht we have samples for

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos social and non social so we need to trim it 


# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(dff_ectos_network_individual_metrics$species_jetz)) %>% mutate(name=dff_ectos_network_individual_metrics$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 

# phylogenetic correlation structure, create a covariance matrix of species
phy_cov<-ape::vcv(phylo, corr=TRUE)

# data structure 

#dff_ectos_network_individual_metrics$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
dff_ectos_network_individual_metrics$foraging_cat<-as.factor(dff_ectos_network_individual_metrics$foraging_cat)
dff_ectos_network_individual_metrics$species_jetz<-as.factor(dff_ectos_network_individual_metrics$species_jetz)
dff_ectos_network_individual_metrics$elevation<-as.numeric(dff_ectos_network_individual_metrics$elevation)
dff_ectos_network_individual_metrics$degree<-as.numeric(dff_ectos_network_individual_metrics$degree)
dff_ectos_network_individual_metrics$w_degree<-as.numeric(dff_ectos_network_individual_metrics$w_degree)
dff_ectos_network_individual_metrics$year_seasonality<-as.numeric(dff_ectos_network_individual_metrics$year_seasonality)

#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
dff_ectos_network_individual_metrics$sociality<-as.factor(dff_ectos_network_individual_metrics$sociality)
dff_ectos_network_individual_metrics$Powder.lvl<-as.factor(dff_ectos_network_individual_metrics$Powder.lvl)
#dff_ectos_network_individual_metrics$total_lice<-as.numeric(dff_ectos_network_individual_metrics$total_lice)
dff_ectos_network_individual_metrics$total_mites<-as.numeric(dff_ectos_network_individual_metrics$total_mites)
dff_ectos_network_individual_metrics$total_no_feathers_mites<-as.numeric(dff_ectos_network_individual_metrics$total_no_feathers_mites)
dff_ectos_network_individual_metrics$species<-dff_ectos_network_individual_metrics$species_jetz # create a column for the species effect different to the phylogenetic one
names(dff_ectos_network_individual_metrics)
is.ultrametric(phylo)

# ##### 2.2.Model selection NETWORKS abundance mites --------------------------------------------------------

zip_a_nf_mites_brms_bayes_no_int_degree<-brms::brm(total_no_feathers_mites~degree+ scale(elevation)+ scale(year_seasonality)+
                                              (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                              (1|Powder.lvl)+
                                              (1|species),
                                            data=dff_ectos_network_individual_metrics,
                                            family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                            data2 = list(phy_cov=phy_cov),
                                            iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                            thin=2,
                                            control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(zip_a_nf_mites_brms_bayes_no_int_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_DEGREE.RDS")
zip_a_nf_mites_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_DEGREE.RDS")

zip_a_nf_mites_brms_bayes_sociality_interactions_degree<-brms::brm(total_no_feathers_mites~
                                                              degree+
                                                              scale(elevation)+
                                                              scale(year_seasonality)+
                                                              sociality:scale(elevation)+
                                                              sociality:scale(year_seasonality)+
                                                              (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                              (1|Powder.lvl)+
                                                              (1|species),
                                                            data=dff_ectos_network_individual_metrics,
                                                            family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                            data2 = list(phy_cov=phy_cov),
                                                            iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                            thin=2,
                                                            control=list(adapt_delta=0.99, max_treedepth=12)) 
saveRDS(zip_a_nf_mites_brms_bayes_sociality_interactions_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")
zip_a_nf_mites_brms_bayes_sociality_interactions_degree<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")

zip_a_nf_mites_brms_bayes_all_interactions_degree<-brms::brm(total_lice~
                                                        degree+
                                                        scale(elevation)+
                                                        scale(year_seasonality)+
                                                        sociality:scale(elevation)+
                                                        sociality:scale(year_seasonality)+
                                                        scale(elevation):scale(year_seasonality)+
                                                        sociality:scale(year_seasonality):scale(elevation)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      control=list(adapt_delta=0.99, max_treedepth=12)) 

saveRDS(zip_a_nf_mites_brms_bayes_all_interactions_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions_DEGREE.RDS")
zip_a_nf_mites_brms_bayes_all_interactions_degree<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions_DEGREE.RDS")

loo(zip_a_lice_brms_bayes_all_interactions,zip_a_lice_brms_bayes_sociality_interactions,zip_a_lice_brms_bayes_no_int)
mcmc_plot(zip_a_nf_mites_brms_bayes_sociality_interactions)




