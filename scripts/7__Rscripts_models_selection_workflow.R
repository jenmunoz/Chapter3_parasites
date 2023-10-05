#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models selection workflow                                                                       ###
### R-code                                                                          ###
### Jenny Munoz      
### R version 4.2.2 (2022-10-31) $nickname [1] "Innocent and Trusting"
### Last update: June 14 2023                                                ###
################################################################################

# Underestanding parameters in brms models https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

# sd(RANDOM EFFECTS) is an estimate of the standard deviation across sociality slopes, thereby showing how much birds mean number of parasites differ in response to social context, elevation and ###
# Cross validation for model selection https://mc-stan.org/loo/articles/online-only/faq.html#modelselection

#Cross-validation is a family of techniques that try to estimate how well a model would predict previously unseen data by using fits of the model to a subset of the data to predict the rest of the data.
#Cross-validation can be used to:
#Asses the predictive performance of a single model
#Asses model misspecification or calibration of the predictive distribution of a single model
#Compare multiple models
#Select a single model from multiple candidates
#Combine the predictions of multiple models
#Even if the goal of the model is not to make predictions, a model which makes bad or badly calibrated predictions is less likely to provide useful insights to a phenomenon studied.

#Using cross-validation for a single model
#Two basic cases why to use cross-validation for one model are:
# We want to know how good predictions the model can make for future or otherwise unseen observations.
 # We want to know if the model describes the observed data well, but we are not going make any predictions for the future.

#How to use cross-validation for model selection?
#First avoid model selection by using the model which includes all predictors and includes all uncertain things. Then optimal thing is to integrate over all the uncertainties. 
#When including many components to a model, it is useful to think more carefully about the prior. For example, if there are many predictors, it is useful to use priors that
#a) state that only some of the effects are big, or b) many effects are big and correlating (it is not possible to have a large number of big independent effects Tosh et al. (2021)).
#If there is explicit utility or loss for observing future predictor values (e.g. medical tests) use decision theory.
#If there is implicit cost for bigger models (e.g. bigger model more difficult to explain or costs of feature measurements are unknown), 
#choose a smaller model which similar predictive performance as the biggest model. If there are only a small number of models, overfitting due to selection process is small. 
#If there are a large number of models, as for example often in variable selection, then the overfitting due to the selection process can be a problem (Piironen and Vehtari, 2017) and more elaborate approaches, such as projection predictive variable selection is recommended.
#If there is application specific utility or loss function, use that to assess practically relevant difference in predictive performance of two models.
#If there is no application specific utility or loss function, use log score, ie elpd. If elpd difference (elpd_diff in loo package) is less than 4, the difference is small (Sivula, Magnusson and Vehtari, 2020)).
#If elpd difference (elpd_diff in loo package) is larger than 4, then compare that difference to standard error of elpd_diff (provided e.g. by loo package) (Sivula, Magnusson and Vehtari, 2020). See also Section How to interpret in Standard error (SE) of elpd difference (elpd_diff)?.
#If there is a large number of models compared, there is possibility of overfitting in model selection.

#Can cross-validation be used for hierarchical / multilevel models?
#The short answer is “Yes”. Hierarchical model is useful, for example, if there are several subjects and for each subject several trials. As discussed in When is cross-validation valid?, it is useful to think of the prediction task or generalizability over different exchangeable entities. We can use different types of cross-validation to choose the focus. This means that also different forms of cross-validation are valid for hierarchical models

# workflow for model selection after reading some of the demo studies
#A case study Roaches cross-validation demo with “random effects” models Cross validation demo https://avehtari.github.io/modelselection/roaches.html

#Model selection note: 
# 1) Analyze the posterior, do they make sense
#2) Cross validation checking (use loo)
#2.1) if random effect present, high pareto K could indicate  misspecification but also flexibility of the model. importance sampling in PSIS-LOO can fail for “random effect” model. 
#2.2) Option 1, refi the model using reloo, but takes very long time
#2.3) Use k-fold and re-fit the model 10 times, each time leaving out 10 % of teh observations and compare the k-models with loo_compare()
#Use posterior predictive checking pp_check  to see if the model can predict the proportion of zeros well

# Prior selection https://paul-buerkner.github.io/brms/reference/set_prior.html
# very informative blog https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
#other interesting read with  many reference that i explored https://statmodeling.stat.columbia.edu/2018/04/03/justify-my-love/
#5 levels of priors
#Flat prior (not usually recommended);
#Super-vague but proper prior: normal(0, 1e6) (not usually recommended);
#Weakly informative prior, very weak: normal(0, 10);
#Generic weakly informative prior: normal(0, 1);
#Specific informative prior: normal(0.4, 0.2) or whatever. Sometimes this can be expressed as a scaling followed by a generic prior: theta = 0.4 + 0.2*z; z ~ normal(0, 1);

# 

#Avoiding model refits in leave-one-out cross-validation with moment matching https://cran.r-project.org/web/packages/loo/vignettes/loo2-moment-matching.html
# also important to check teh priors prior_summary()
# 0.Libraries  --------------------------------------------------------------
R.Version()
# libraries for easier manipulation of data
install.packages("pacman")
library("pacman")
install.packages('BiocManager')
library("BiocManager")



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
install.packages("rstan")
#install.packages('brmstools')
remotes::install_github("Pakillo/DHARMa.helpers")
devtools::install_github("mvuorre/brmstools")
install.packages('rstanarm')
install.packages('loo')
devtools::install_github("paul-buerkner/brms") # see if this helps when using reloo=TRUE to remove pareto
#install.packages('brms') # bayesian approach to model phylogenetic data with repides observations


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
library(rstan)
library(rstanarm)
library(loo)

# # ## 0. Models selected -------------------------------------------------

#Models

#Infection
selected_ecto_infection_brms_bayes_no_int

#Infection networks
selected_ecto_p_brms_bayes_no_int_degree_prior

#Abundance Lice
selected_zinb_a_lice_brms_bayes_no_int_priors

#Abundance Lice networks
selected_zinb_a_lice_brms_bayes_no_int_degree_prior

#Abundance non-feather mites 
selected_zinb_a_nf_mites_brms_bayes_no_int_prior

#Abundance non-feather mites networks
selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior

#Abundance all mites 
selected_zinb_a_all_mites_brms_bayes_no_int_prior

#Prevalence
ecto_p_brms_bayes_no_int_species_priors_zobi

#Prevalence networks
ecto_p_brms_bayes_no_int_species_priors_degree_zobi

#Lice richness
selected_poisson_lice_diversity_sociality_no_int_priors

#Lice richness networks 
selected_poisson_lice_diversity_degree_no_int_priors
# ### 0. Data summaries  -----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
select(general_diversity, family,elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality, mass_tidy_species, mass_ind_tidy,mass_ind_comp ) %>% 
  na.omit()
View(ectos_birds_dff)

#total_lice, total_mites, total_no_feathers_mites, total_ticks

names(ectos_birds_dff)

#View(ectos_birds_dff)
dim(ectos_birds_dff)
unique(ectos_birds_dff$species_jetz)
unique(ectos_birds_dff$family)

# sociality
ectos_birds_dff %>% 
  group_by(sociality,species_jetz) %>% 
  summarise(n=n(), P=(n()/871*100) )

a<-ectos_birds_dff %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), sociality=max(sociality)) 

a %>% group_by(sociality) %>% 
  summarise(n=n(), P=(n()/871*100) )

ectos_birds_dff<-ectos_pres_abs %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE)%>% 
  filter(sample_size>4) 

pp_check()
R.Version()


# overall prevalence 
ectos_birds_dff %>% 
  group_by(ectoparasites_PA) %>% 
summarise(n=n(), P=(n()/871*100) )

# Prevalence by family 
prevalence_step1<-ectos_birds_dff %>% 
  group_by( family,ectoparasites_PA) %>% 
  summarise(n=n())
View(prevalence_step1)

# this calculation is not perfect because if one family only have zeros in the samples that will be counted as positive but since we are filtering families with more than 5 samples we shoudl be ok
# Following the paper (Jovani and Tella 2006) we modeled only families with more that 10 samples, we did not use the same aproach for species because it will exclude most of our non social species.
View(prevalence_step1)
prevalence<-prevalence_step1 %>% summarise(family_name=first(family), total_samples=sum(n), total_positive=last(n), prevalence=((total_positive/total_samples)*100)) %>% 
  filter(total_samples>4)

write.csv(prevalence, "data/data_manuscript/7_prevalence_ectoparasite_family_level.csv")

# dIVERSITY AND COINFECTIONS

ectos_birds_dff %>% 
  group_by(general_diversity) %>% 
  summarise(n=n(),P=(n()/871*100))

# DIVERSITY BY GROUP

# mites
ectos_birds_dff<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(general_diversity, family,elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality, mass_tidy_species, mass_ind_tidy,mass_ind_comp, total_mites, total_no_feathers_mites ) %>% 
  na.omit() 

ectos_birds_dff%>% summarize(max=max(total_mites), mean=mean(total_mites), sd=sd(total_mites), min=min(total_mites))

ectos_birds_dff%>% summarize(max=max(total_no_feathers_mites), mean=mean(total_no_feathers_mites), sd=sd(total_no_feathers_mites), min=min(total_no_feathers_mites))


ectos_birds_dff %>% 
  filter(total_mites!=0) %>% 
  summarise(samples_mites=sum(!is.na(total_mites)), (p=samples_mites/831*100))


# lice
ectos_birds_dff<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(general_diversity, family,elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality, mass_tidy_species, mass_ind_tidy,mass_ind_comp, total_lice ) %>% 
  na.omit() %>% 
  View()

ectos_birds_dff%>% summarize(max=max(total_lice), mean=mean(total_lice), sd=sd(total_lice), min=min(total_lice))


ectos_birds_dff %>% 
  filter(total_lice!=0) %>% 
  summarise(samples_lice=sum(!is.na(total_lice)), (p=samples_lice/831*100), max=max(total_lice))
  
# ticks

ectos_birds_dff<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(general_diversity, family,elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality, mass_tidy_species, mass_ind_tidy,mass_ind_comp, Ticks ) %>% 
  na.omit()

names(ectos_birds_dff)

ectos_birds_dff %>% 
  filter(Ticks!=0) %>% 
  summarise(samples_ticks=sum(!is.na(Ticks)), (p=samples_ticks/831*100))


# ##### 1.Data processing prevalence ectos ----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, sociality_groups, year_seasonality, mass_tidy_species, mass_ind_tidy,mass_ind_comp ) %>% 
  na.omit()
#ectos_birds_dff <- get.complete.cases(ectos_birds_dff) # mke sure we get all complete cases 



phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # Need to prunne the tree 

# Data structure
# species=variable accounts for any specific effect that would be independent of the phylogenetic relationship between species
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$year_seasonality<-as.numeric(ectos_birds_dff$year_seasonality)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$sociality_groups<-as.factor(ectos_birds_dff$sociality_groups) # including ant followers as social

ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$ectoparasites_PA<-as.numeric(ectos_birds_dff$ectoparasites_PA)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one
ectos_birds_dff$species<-as.factor(ectos_birds_dff$species)
ectos_birds_dff$mass_tidy_species<-as.numeric(ectos_birds_dff$mass_tidy_species)
ectos_birds_dff$mass_ind_tidy<-as.numeric(ectos_birds_dff$mass_ind_tidy)
ectos_birds_dff$mass_ind_comp<-as.numeric(ectos_birds_dff$mass_ind_comp)

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

#This may seem fairly informative, fixed effects are centers and scales your variables before fitting your model, so these are not as informative as they would be if the data were left  on their input scales!
  
#fixed_prior<- we will keep the dafault flat priors

ecto_p_brms_bayes_no_int<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
                               (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                               (1|Powder.lvl)+
                               (1|species),
                             data=ectos_birds_dff,
                             #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                             family= bernoulli(), # bernoulli() uses the (link = "logit")
                             data2 = list(phy_cov=phy_cov),
                             #prior = c(random_prior,intercept_prior),
                             iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                             thin=2,
                             control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(ecto_p_brms_bayes_no_int, "data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")
ecto_p_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")

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
                                              #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                              #prior = c(random_prior,intercept_prior),
                                              family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                              data2 = list(phy_cov=phy_cov),
                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                              thin=2,
                                              control=list(adapt_delta=0.99, max_treedepth=14)) 

#saveRDS(ecto_p_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions.RDS")
ecto_p_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions.RDS")

loo(ecto_p_brms_bayes_all_interactions)

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
                             iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                             thin=2,
                             control=list(adapt_delta=0.99, max_treedepth=14)) 

#saveRDS(ecto_p_brms_bayes_all_interactions, "data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")
ecto_p_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")

# PRIOrS SPECIFICATON 

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


# half student only allows positive values

# I am not sure why it does not like the intercept wider prior so I used the default
ecto_p_brms_bayes_no_int_prior<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
                                             (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                             (1|Powder.lvl)+
                                             (1|species),
                                           data=ectos_birds_dff,
                                           save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
                                           family= bernoulli(), # bernoulli() uses the (link = "logit")
                                           data2 = list(phy_cov=phy_cov),
                                           prior=c(prior_predictors,prior_random),
                                           iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                           thin=2,
                                           control=list(adapt_delta=0.99, max_treedepth=14))

#saveRDS(ecto_p_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions_priors_default_intercept.RDS")
ecto_p_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions_priors_default_intercept.RDS")


#saveRDS(ecto_p_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1P.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_priors.RDS")

ecto_p_brms_bayes_sociality_interactions_priors<-brms::brm(ectoparasites_PA~
                                                      sociality+
                                                      scale(elevation)+
                                                      scale(year_seasonality)+
                                                      sociality:scale(elevation)+
                                                      sociality:scale(year_seasonality)+
                                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                      (1|Powder.lvl)+
                                                      (1|species),
                                                    data=ectos_birds_dff,
                                                    save_pars = save_pars(all=  TRUE),
                                                    prior=c(prior_predictors,prior_random),
                                                    family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                    data2 = list(phy_cov=phy_cov),
                                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                    thin=2,
                                                    control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(ecto_p_brms_bayes_sociality_interactions_priors,"data/data_analyses/model_selection/M1P.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_priors.RDS")
ecto_p_brms_bayes_sociality_interactions_priors<-readRDS("data/data_analyses/model_selection/M1P.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_priors.RDS")

ecto_p_brms_bayes_all_interactions_priors<-brms::brm(ectoparasites_PA~
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
                                              prior=c(prior_predictors,prior_random),
                                              save_pars = save_pars(all=  TRUE),
                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                              thin=2,
                                              control=list(adapt_delta=0.99, max_treedepth=14)) 

#saveRDS(ecto_p_brms_bayes_all_interactions_priors,"data/data_analyses/model_selection/M1P.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_priors.RDS")
ecto_p_brms_bayes_all_interactions_priors<-readRDS("data/data_analyses/model_selection/M1P.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_priors.RDS")

#Other models removing individual variables 



ecto_p_brms_bayes_no_int_no_species<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                 (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                 (1|Powder.lvl),
                                               data=ectos_birds_dff,
                                               #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                               family= bernoulli(), # bernoulli() uses the (link = "logit")
                                               data2 = list(phy_cov=phy_cov),
                                               #prior = c(random_prior,intercept_prior),
                                               iter=4000, warmup=2000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                               thin=2,
                                               control=list(adapt_delta=0.99, max_treedepth=14)) 


ecto_p_brms_bayes_no_int_no_phylo<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
                                               (1|Powder.lvl)+
                                               (1|species),
                                             data=ectos_birds_dff,
                                             #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                             family= bernoulli(), # bernoulli() uses the (link = "logit")
                                             #prior = c(random_prior,intercept_prior),
                                             iter=4000, warmup=2000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                             thin=2,
                                             control=list(adapt_delta=0.99, max_treedepth=14)) 


mcmc_plot(ecto_p_brms_bayes_no_int_no_species)
loo(ecto_p_brms_bayes_no_int_nopowder,ecto_p_brms_bayes_no_int,ecto_p_brms_bayes_no_int_no_species,ecto_p_brms_bayes_no_int_no_phylo, compare=TRUE )

ecto_p_brms_bayes_no_int_prior_mass<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+scale(mass_tidy)+
                                                 (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                 (1|Powder.lvl)+
                                                 (1|species),
                                               data=ectos_birds_dff,
                                               save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
                                               family= bernoulli(), # bernoulli() uses the (link = "logit")
                                               data2 = list(phy_cov=phy_cov),
                                               prior=c(prior_predictors,prior_random),
                                               iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                               thin=2,
                                               control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(ecto_p_brms_bayes_no_int_prior_mass, "data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions_priors_default_intercept_mass.RDS")



###_###_###_##
#MODEL COMPARISON
###_###_###_##

#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.
bayes_R2(ecto_p_brms_bayes_no_int_prior)
bayes_R2(ecto_p_brms_bayes_sociality_interactions_priors)
bayes_R2(ecto_p_brms_bayes_all_interactions_priors)

# use loo cross validation 
loo(ecto_p_brms_bayes_no_int_prior, ecto_p_brms_bayes_sociality_interactions_priors, ecto_p_brms_bayes_all_interactions_priors,compare=TRUE)
loo(ecto_p_brms_bayes_no_int_prior,ecto_p_brms_bayes_no_int_prior_mass,compare=TRUE)
bayes_R2(ecto_p_brms_bayes_no_int_prior)
# R=eld_diff<4 so we keep the simplest model !
loo_compare(waic(ecto_p_brms_bayes_no_int_prior), waic(ecto_p_brms_bayes_sociality_interactions_priors)) # interesting warning

# use k-fold-cv validation instead because in models with random efects loo tends to fail # but this is taking forevwer so i WILL RUNT IT AT UBC TOMORROW
k_ecto_p_brms_no_int<-kfold(ecto_p_brms_bayes_no_int_prior, K=10)
k_ecto_p_brms_sociality_int<-kfold(ecto_p_brms_bayes_sociality_interactions_priors, K=10)
loo_compare(k_ecto_p_brms_no_int, k_ecto_p_brms_sociality_int)

# Model posterior predictive checks 
#The idea of posterior predictive checks is to compare our observed data to replicated data from the model. 
#If our model is a good fit, we should be able to use it to generate a dataset that resembles the observed data.
#gettigsamples from the posterior predictive distribution:
pp_check(ecto_p_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100) 

#or dependeing or the data use
#pp_check(ecto_p_brms_bayes_no_int_prior, type = "stat", stat = 'median', nsamples = 100)
#pp_m<- brms::posterior_predict(ecto_p_brms_bayes_no_int_prior)
#ppc_rootogram(y=ecto_p_brms_bayes_no_int_prior$data$ectoparasites_PA, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))

###_###_###_##
#PLOTS SELECTED MODEL
###_###_###_##


conditional_effects(ecto_p_brms_bayes_no_int_prior)
marginal_effects(ecto_p_brms_bayes_no_int_prior)
plot( conditional_effects(ecto_p_brms_bayes_no_int_prior), 
    points = TRUE, 
    point_args = list(width = .05, shape = 1))
  
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("blue")

# model convergence 
#png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
png("figures/figures_manuscript/models_selected_figures/best_models/Fig1P_ecto_INFECTION_brms_bayes_no_int_prior_CONVERGENCE.png",width = 4000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int_prior)
dev.off()

# model fit
#png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
png("figures/figures_manuscript/models_selected_figures/best_models/Fig1P_ecto_INFECTION_brms_bayes_no_int_prior_FIT.png",width = 4000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100) 
dev.off()

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

color_scheme_set("blue")

estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS INFECTION ")+
  theme_classic(30)+
  xlim(-4,4)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

#png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_ESTIMATES.png",width = 3000, height = 3000, res = 300, units = "px")
png("figures/figures_manuscript/models_selected_figures/best_models/Fig1P_ecto_INFECTION_brms_bayes_no_int_prior_ESTIMATES.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_prior_mass,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS INFECTION ")+
  theme_classic(30)+
  xlim(-4,4)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig1P_ecto_INFECTION_brms_bayes_no_int_prior_ESTIMATES_INTERVALS.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



# PREVALENCE ####### 1.3 ***Selected*** model prevalence ectos INFECTION (included mass) ----------------------------

# PRIORS SPECIFICATON 

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

# EXCLUDED POWDER LEVEL BECAUSE IT DOES NOT SIGNIFICANTLY IMPROVE MODEL FIT # alternatively use the model ecto_p_brms_bayes_no_int_prior

#[SELECTED MODEL]
# including the body mass of the individual scaled and ant-followers as social species 
selected_ecto_infection_brms_bayes_no_int<-brms::brm(ectoparasites_PA~sociality_groups+
                                                       scale(elevation)+ scale(year_seasonality)+
                                                       scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                       (1|gr(species_jetz, cov = phy_cov))+ (1|species), #(1|Powder.lvl) # excluded poweder level cause it does not significantly improve model fit
                                                     data=ectos_birds_dff,
                                                     save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
                                                     family= bernoulli(), # bernoulli() uses the (link = "logit")
                                                     data2 = list(phy_cov=phy_cov),
                                                     prior =c(prior_predictors,prior_random,prior_intercept),
                                                     iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                     thin=2,
                                                     control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(selected_ecto_infection_brms_bayes_no_int,"results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED_antfollowers_included.RDS")
selected_ecto_infection_brms_bayes_no_int<-readRDS("results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED_antfollowers_included.RDS")

coef(selected_ecto_infection_brms_bayes_no_int)


# including the body mass of the individual scaled
# selected_ecto_infection_brms_bayes_no_int<-brms::brm(ectoparasites_PA~sociality+
#                                                     scale(elevation)+ scale(year_seasonality)+
#                                                      scale(mass_tidy_species)+(mass_ind_comp)+
#                                                      (1|gr(species_jetz, cov = phy_cov))+ (1|species), #(1|Powder.lvl) # excluded poweder level cause it does not significantly improve model fit
#                                              data=ectos_birds_dff,
#                                              save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
#                                              family= bernoulli(), # bernoulli() uses the (link = "logit")
#                                              data2 = list(phy_cov=phy_cov),
#                                              prior =c(prior_predictors,prior_random,prior_intercept),
#                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                              thin=2,
#                                              control=list(adapt_delta=0.99, max_treedepth=14)) 
# saveRDS(selected_ecto_infection_brms_bayes_no_int,"results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED.RDS")
# selected_ecto_infection_brms_bayes_no_int<-readRDS("results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED.RDS")



# selected_ecto_infection_brms_bayes_no_int<-brms::brm(ectoparasites_PA~sociality+
#                                                        scale(elevation)+ scale(year_seasonality)+
#                                                        scale(mass_tidy_species)+scale(mass_ind_comp)+
#                                                        (1|gr(species_jetz, cov = phy_cov))+ (1|species), #(1|Powder.lvl) # excluded poweder level cause it does not significantly improve model fit
#                                                      data=ectos_birds_dff,
#                                                      save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
#                                                      family= bernoulli(), # bernoulli() uses the (link = "logit")
#                                                      data2 = list(phy_cov=phy_cov),
#                                                      prior =c(prior_predictors,prior_random,prior_intercept),
#                                                      iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                                      thin=2,
#                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
# 
# 
# saveRDS(selected_ecto_infection_brms_bayes_no_int,"results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED_ind_mass_scaled.RDS")
# selected_ecto_infection_brms_bayes_no_int<-readRDS("results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED_ind_mass_scaled.RDS")


###_###_###_##
#PLOTS SELECTED MODEL
###_###_###_##

conditional_effects(selected_ecto_infection_brms_bayes_no_int)
marginal_effects(selected_ecto_infection_brms_bayes_no_int)
plot( conditional_effects(selected_ecto_infection_brms_bayes_no_int), 
      points = TRUE, 
      point_args = list(width = .05, shape = 1))

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("blue")

# model convergence 
#png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
png("results/selected_models_figures/1_Fig1P_BEST_ecto_INFECTION_brms_bayes_no_int_prior_CONVERGENCE.png",width = 4000, height = 3000, res = 300, units = "px")
plot(selected_ecto_infection_brms_bayes_no_int)
dev.off()

# model fit
#png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
png("results/selected_models_figures/1_Fig1P_BEST_ecto_INFECTION_brms_bayes_no_int_prior_FIT.png",width = 4000, height = 3000, res = 300, units = "px")
pp_check(selected_ecto_infection_brms_bayes_no_int, type = "dens_overlay", ndraws = 100) 
dev.off()

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

color_scheme_set("blue")

estimates_plot<-mcmc_plot(selected_ecto_infection_brms_bayes_no_int,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS INFECTION ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/1_Fig1P_BEST_ecto_INFECTION_brms_bayes_no_int_prior_ESTIMATES.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_ecto_infection_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS INFECTION ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/1_Fig1P_BEST_ecto_INFECTION_brms_bayes_no_int_prior_ESTIMATES_INTERVALS_scaled_ind_mass.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

bayes_R2(selected_ecto_infection_brms_bayes_no_int)

# half student only allows positive values

# ######1.Data processing and  model selection for prevalence species level (zero_one_beta) --------
# Our prevalence data has an extra complication and is thatit has ones while beta does not allow u to model proportion that reach one
# I found this solutions and decided to implement 1- response ( but phi does not converge), as that will just model the probability of not getting infecteda https://github.com/paul-buerkner/brms/issues/942

ectos_birds_df<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
                           select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality,sociality_groups, year_seasonality, mass_tidy_species ) %>% 
                           na.omit() %>%  filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites
#note sociality groups include ant followers as mixed species flockers 

str(ectos_birds_df)

ectos_pres_abs<-ectos_birds_df %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation=mean(elevation), sociality=max(sociality_groups), mass=max(mass_tidy_species), foraging_cat=first(foraging_cat))%>% 
  mutate(proportion_ectoparasites=ectoparasites_presence/sample_size) %>% 
  na.omit() 
View(ectos_pres_abs)

ectos_birds_dff<-ectos_pres_abs %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE)%>% 
  filter(sample_size>4) 

unique(ectos_birds_dff$proportion_ectoparasites)
View(ectos_birds_dff)

#ectos_birds_dff <- get.complete.cases(ectos_birds_dff) # mke sure we get all complete cases 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # Need to prunne the tree 

# Data structure
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$proportion_ectoparasites<-as.numeric(ectos_birds_dff$proportion_ectoparasites)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)

#ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality) # including ant followers as social

str(ectos_birds_dff)
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

#PRIORS
prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl work, cause our predictors are scaled could use a wider prior if needed normal(0,10)
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
prior_phi<-prior("normal(0, 1)", class ="phi")


#residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default

# modeling the opossitive of prevalence of getting infected, meaninsn modeling the probability of not getting infected


# The models 

prior_summary(ecto_p_brms_bayes_no_int_species)

#Since we have zeros and ones we will use a zero-one inflated distribution
ecto_p_brms_bayes_no_int_species_priors_zobi<-brms::brm((proportion_ectoparasites)~sociality+ scale(elevation)+scale(sample_size)+scale(mass)+
                                                                   (1|gr(species_jetz, cov = phy_cov))+
                                                                   (1|species),
                                                                 data=ectos_birds_dff,
                                                                 #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                                                 family= zero_one_inflated_beta(), # # this family is for proportions that include zeros and ones                                         
                                                                 data2 = list(phy_cov=phy_cov),
                                                                 save_pars = save_pars(all = TRUE),
                                                                 prior = c(prior_predictors,prior_random,prior_intercept,prior_phi),
                                                                 iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                 thin=2,
                                                                 control=list(adapt_delta=0.99, max_treedepth=14))


saveRDS(ecto_p_brms_bayes_no_int_species_priors_zobi, "results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_zobi_antbirds_included.RDS")
ecto_p_brms_bayes_no_int_species_priors_zobi<-readRDS("results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_zobi_antbirds_included.RDS")

loo(ecto_p_brms_bayes_no_int_species_priors_zobi)

#PLOTS
plot(conditional_effects(ecto_p_brms_bayes_no_int_species_priors_zib, dpar="mu"), 
     points = TRUE, 
     point_args = list(width = .05, shape = 1))

color_scheme_set("viridisE") 

bayes_R2(ecto_p_brms_bayes_no_int_species_priors)

#model convergence 
png("results/selected_models_figures/2_1Fig2Ps.plot_model_CONVERGENCE_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES_zobi.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int_species_priors_zobi)
dev.off()

# model fit
png("results/selected_models_figures/2_1Fig2Ps.plot_model_FIT_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES_zobi.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_species_priors_zobi, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_zobi,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/2_1Fig2Ps..plot_model_parameters_PREVALENCE SPECIES_brms_bayes_no_int_zobi.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_zobi2,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/2_1Fig2Ps..plot_model_parameters_intervals_PREVALENCE SPECIES_brms_bayes_no_int_zobi2.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

###
# INCLUDING DEGREE
###

ectos_birds_df<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality, mass_tidy_species,degree_species,w_degree_species ) %>% 
  na.omit() %>%  filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites
##

ectos_pres_abs<-ectos_birds_df %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation=mean(elevation), sociality=max(sociality), mass=max(mass_tidy_species), degree_species=max(degree_species ), w_degree_species=max(w_degree_species))%>% 
  mutate(proportion_ectoparasites=ectoparasites_presence/sample_size) %>% 
  na.omit() 

ectos_birds_dff<-ectos_pres_abs %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE)%>% 
  filter(sample_size>4) 

unique(ectos_birds_dff$proportion_ectoparasites)
View(ectos_birds_dff)

#ectos_birds_dff <- get.complete.cases(ectos_birds_dff) # mke sure we get all complete cases 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # Need to prunne the tree 

# Data structure
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$proportion_ectoparasites<-as.numeric(ectos_birds_dff$proportion_ectoparasites)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one
ectos_birds_dff$species<-as.factor(ectos_birds_dff$species) # create a column for the species effect different to the phylogenetic one

ectos_birds_dff$degree_species<-as.numeric (ectos_birds_dff$degree_species) # create a column for the species effect different to the phylogenetic one
ectos_birds_dff$w_degree_species<-as.numeric (ectos_birds_dff$w_degree_species) # create a column for the species effect different to the phylogenetic one

ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
str(ectos_birds_dff)
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

#PRIORS
prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl work, cause our predictors are scaled could use a wider prior if needed normal(0,10)
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
prior_phi<-prior("normal(0, 1)", class ="phi")

View(ectos_birds_dff)

#ecto_p_brms_bayes_no_int_species_priors_degree_zib<-brms::brm((1-proportion_ectoparasites)~scale(degree_species)+ scale(elevation)+scale(sample_size)+scale(mass)+
#                                                                 (1|gr(species_jetz, cov = phy_cov))+
#                                                                (1|species),
 #                                                              data=ectos_birds_dff,
 #                                                              #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
#                                                             family= zero_inflated_beta(), # # this family is for proportions that include zeros and ones                                         
#                                                            data2 = list(phy_cov=phy_cov),
#                                                           save_pars = save_pars(all = TRUE),
#                                                               prior = c(prior_predictors,prior_random,prior_intercept,prior_phi),
#                                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                                               thin=2,
#                                                              control=list(adapt_delta=0.99, max_treedepth=14))

ecto_p_brms_bayes_no_int_species_priors_degree_zobi<-brms::brm((proportion_ectoparasites)~scale(degree_species)+ scale(elevation)+scale(sample_size)+scale(mass)+
                                              (1|gr(species_jetz, cov = phy_cov))+
                                              (1|species),
                                            data=ectos_birds_dff,
                                            #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                            family= zero_one_inflated_beta(), # # this family is for proportions that include zeros and ones                                         
                                            data2 = list(phy_cov=phy_cov),
                                            save_pars = save_pars(all = TRUE),
                                            prior = c(prior_predictors,prior_random,prior_intercept,prior_phi),
                                            iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                            thin=2,
                                            control=list(adapt_delta=0.99, max_treedepth=14))


saveRDS(ecto_p_brms_bayes_no_int_species_priors_degree_zobi, "results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_DEGREE_zobi.RDS")
ecto_p_brms_bayes_no_int_species_priors_degree_zobi<-readRDS("results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_DEGREE_zobi.RDS")

# Plots

color_scheme_set("viridisE") 

#model convergence 
png("results/selected_models_figures/2_2Fig2Ps.plot_model_CONVERGENCE_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES_zobi_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int_species_priors_degree_zobi)
dev.off()

# model fit
png("results/selected_models_figures/2_2Fig2Ps.plot_model_FIT_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES_zobi_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_species_priors_degree_zobi, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()

png("results/selected_models_figures/2_2Fig2Ps.plot_model_FIT_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES_zib_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_species_priors_degree_zib, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_degree_zobi,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~DEGREE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/2_2_Fig2Ps.plot_model_parameters_PREVALENCE SPECIES_brms_bayes_no_int_zobi_DEGREE.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_degree_zobi,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~DEGREE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/2_2Fig2Ps..plot_model_parameters_intervals_PREVALENCE SPECIES_brms_bayes_no_int_zobi_DEGREE.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

###
### Strenght 
###

ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi<-brms::brm((proportion_ectoparasites)~scale(w_degree_species)+ scale(elevation)+scale(sample_size)+scale(mass)+
                                                                 (1|gr(species_jetz, cov = phy_cov))+
                                                                 (1|species),
                                                               data=ectos_birds_dff,
                                                               #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                                               family= zero_one_inflated_beta(), # # this family is for proportions that include zeros and ones                                         
                                                               data2 = list(phy_cov=phy_cov),
                                                               save_pars = save_pars(all = TRUE),
                                                               prior = c(prior_predictors,prior_random,prior_intercept,prior_phi),
                                                               iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                               thin=2,
                                                               control=list(adapt_delta=0.99, max_treedepth=14))


saveRDS(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi, "results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_STRENGHT_zobi.RDS")
ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi<-readRDS("results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_STRENGHT_zobi.RDS")

# Plots

color_scheme_set("viridisE") 

#model convergence 
png("results/selected_models_figures/2_3Fig2Ps.plot_model_CONVERGENCE_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES_zobi_STRENGHT.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi)
dev.off()

# model fit
png("results/selected_models_figures/2_3Fig2Ps.plot_model_FIT_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES_zobi_STRENGHT.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~STRENGHT")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/2_3Fig2Ps.plot_model_parameters_PREVALENCE SPECIES_brms_bayes_no_int_zobi_STRENGHT.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~STRENGHT")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/2_3Fig2Ps..plot_model_parameters_intervals_PREVALENCE SPECIES_brms_bayes_no_int_zobi_STRENGHT.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()




# ##### 2.Data processing abundance lice ----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, sociality_groups, total_lice,year_seasonality, mass_tidy_species, mass_ind_tidy,mass_ind_comp) %>% 
  filter(species_jetz!="Premnoplex_brunnescens") %>% 
  na.omit()


#%>% filter(total_lice<60)  # removing outliers 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

# Finding putliers
#ectos_birds_dff <-ectos_birds_dff  %>% mutate_at("species_jetz", str_replace, " ", "_")
#ectos_birds_dff$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$sociality_groups<-as.factor(ectos_birds_dff$sociality_groups) # including ant followers as social

ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$total_lice<-as.numeric(ectos_birds_dff$total_lice)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one
ectos_birds_dff$species<-as.factor(ectos_birds_dff$species)
ectos_birds_dff$mass_tidy_species<-as.numeric(ectos_birds_dff$mass_tidy_species)
ectos_birds_dff$mass_ind_tidy<-as.numeric(ectos_birds_dff$mass_ind_tidy)
ectos_birds_dff$mass_ind_comp<-as.numeric(ectos_birds_dff$mass_ind_comp)

names(ectos_birds_dff)
is.ultrametric(phylo)
str(ectos_birds_dff)
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
prior_summary(zinb_a_lice_brms_bayes_no_int)

zinb_a_lice_brms_bayes_no_int<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                      (1|Powder.lvl)+
                                      (1|species),
                                    data=ectos_birds_dff,
                                    family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                    data2 = list(phy_cov=phy_cov),
                                    #prior = c(random_prior,intercept_prior, residual_prior,residual_prior_2),
                                    #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                    thin=2,
                                    control=list(adapt_delta=0.99, max_treedepth=14)) 



#saveRDS(zinb_a_lice_brms_bayes_no_int, "data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS")
zinb_a_lice_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS")

zinb_a_lice_brms_bayes_sociality_interactions<-brms::brm(total_lice~
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
                                                    #prior = c(random_prior,intercept_prior, residual_prior,residual_prior_2),
                                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                    thin=2,
                                                    control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(zip_a_lice_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/M1L.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")
zinb_a_lice_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/M1L.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")

zinb_a_lice_brms_bayes_all_interactions<-brms::brm(total_lice~
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
                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                              thin=2,
                                              control=list(adapt_delta=0.99, max_treedepth=14)) 

#saveRDS(zinb_a_lice_brms_bayes_all_interactions, "data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")
zinb_a_lice_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")

# PRIORS

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

prior_summary(zinb_a_lice_brms_bayes_no_int_priors)

zinb_a_lice_brms_bayes_no_int_priors<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                           (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                           (1|Powder.lvl)+
                                           (1|species),
                                         data=ectos_birds_dff,
                                         family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                         data2 = list(phy_cov=phy_cov),
                                         prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                         #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                         iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                         thin=2,
                                         control=list(adapt_delta=0.99, max_treedepth=14))

#saveRDS(zinb_a_lice_brms_bayes_no_int_priors, "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")
zinb_a_lice_brms_bayes_no_int_priors<-readRDS( "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")

loo(zinb_a_lice_brms_bayes_no_int)
mcmc_plot(zinb_a_lice_brms_bayes_sociality_interactions_priors)
mcmc_plot(zinb_a_lice_brms_bayes_sociality_interactions)


zinb_a_lice_brms_bayes_sociality_interactions_priors<-brms::brm(total_lice~
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
                                                         prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                         iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                         thin=2,
                                                         control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(zinb_a_lice_brms_bayes_sociality_interactions_priors, "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_sociality_interactions_priors.RDS")
zinb_a_lice_brms_bayes_sociality_interactions_priors<-readRDS( "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_sociality_interactions_priors.RDS")


zinb_a_lice_brms_bayes_all_interactions_priors<-brms::brm(total_lice~
                                                     sociality+
                                                     scale(elevation)+
                                                     scale(year_seasonality)+
                                                     sociality:scale(elevation)+
                                                     sociality:scale(year_seasonality)+
                                                     scale(elevation):scale(year_seasonality)+
                                                     sociality:scale(year_seasonality):scale(elevation)+
                                                     (1|gr(species_jetz, cov = phy_cov))+ 
                                                     (1|Powder.lvl)+
                                                     (1|species),                                                       
                                                     prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                     data=ectos_birds_dff,
                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

#saveRDS(zinb_a_lice_brms_bayes_all_interactions_priors, "data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_priors.RDS")
zinb_a_lice_brms_bayes_all_interactions_priors<-readRDS("data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_priors.RDS")

###_###_###_##
#MODEL COMPARISON
###_###_###_##

#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.
bayes_R2(zinb_a_lice_brms_bayes_no_int_priors)

bayes_R2(zinb_a_lice_brms_bayes_sociality_interactions_priors)
bayes_R2(zinb_a_lice_brms_bayes_all_interactions_priors)

# use loo cross validation [ In this case we are not interested inteh interaaction between the different factors because the sociality of a speceis does not change with seasonality or elevation, cause those are attributes to the individuals does not change ]
loo(zinb_a_lice_brms_bayes_no_int_priors) 


loo(zinb_a_lice_brms_bayes_no_int, zinb_a_lice_brms_bayes_sociality_interactions, zinb_a_lice_brms_bayes_all_interactions,compare=TRUE)

loo(zinb_a_lice_brms_bayes_no_int_priors, zinb_a_lice_brms_bayes_sociality_interactions_priors, zinb_a_lice_brms_bayes_all_interactions_priors,compare=TRUE)
loo(zinb_a_lice_brms_bayes_no_int_priors, zinb_a_lice_brms_bayes_sociality_interactions_priors,compare=TRUE)
# R=eld_diff<4 so we keep the simplest model !
loo_compare(waic(zinb_a_lice_brms_bayes_no_int_priors), waic(zinb_a_lice_brms_bayes_sociality_interactions_priors)) # interesting warning
mcmc_plot(zinb_a_lice_brms_bayes_no_int)
mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors)
# use k-fold-cv validation instead because in models with random efects loo tends to fail # but this is taking forevwer so i WILL RUNT IT AT UBC TOMORROW
k_lice_a_brms_no_int_prior<-kfold(zinb_a_lice_brms_bayes_no_int_priors, K=10)
k_ecto_p_brms_sociality_int<-kfold(ecto_p_brms_bayes_sociality_interactions_priors, K=10)
loo_compare(k_ecto_p_brms_no_int, k_ecto_p_brms_sociality_int)

# Model posterior predictive checks 
#The idea of posterior predictive checks is to compare our observed data to replicated data from the model. 
#If our model is a good fit, we should be able to use it to generate a dataset that resembles the observed data.
#gettigsamples from the posterior predictive distribution:

#or dependeing or the data use
pp_check(zinb_a_lice_brms_bayes_no_int_priors, type = "stat", stat = 'median', nsamples = 100)
pp_m<- brms::posterior_predict(zinb_a_lice_brms_bayes_no_int_priors)
ppc_rootogram(y=zinb_a_lice_brms_bayes_no_int_priors$data$total_lice, pp_m[1:200, ])  +   
coord_cartesian(xlim = c(0, 100), ylim = c(0,30))

###_###_###_##
#PLOTS
###_###_###_##

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2(zinb_a_lice_brms_bayes_no_int_priors) # R2 0.1529

color_scheme_set("teal")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig1L.ZINB_a_lice_brms_bayes_no_int_priors_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_lice_brms_bayes_no_int_priors)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig1L.ZINB_a_lice_brms_bayes_no_int_priors_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_lice_brms_bayes_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

# Choose depending on the model type
#png("figures/figures_manuscript/mT.png",width = 3000, height = 3000, res = 300, units = "px")
#pp_m<- brms::posterior_predict(MODEL)
#ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
#  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
#dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

color_scheme_set("teal")

estimates_plot<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions ZINB LICE ABUNDANCE", subtitle ="ZINB LICE ABUNDANCE with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1L.ZINB_LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions ZINB LICE ABUNDANCE", subtitle ="ZINB LICE ABUNDANCE with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1L.ZINB_LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES_INTERVALS.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


#MODEL COMPARISSOM warrning 20 pareto k are high!!!! what to do?
loo(ZIP_a_lice_brms_bayes_no_int_priors, zinb_a_lice_brms_bayes_no_int_priors)
looni<-loo(zip_a_lice_brms_bayes_no_int)
loonip<-loo(zinb_a_lice_brms_bayes_no_int_priors)

loosi<-loo(zip_a_lice_brms_bayes_sociality_interactions)
#looni1<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, moment_match = TRUE)
#looni2<-loo_moment_match(zinb_a_nf_mites_brms_bayes_no_int_outliers, loo=loo1,k_threshold = 0.7)
looni3<-loo(zinb_a_lice_brms_bayes_no_int, reloo = TRUE) # exclude the observation with high pareto to allow me to compare with other models 
looni4<-reloo(zinb_a_lice_brms_bayes_no_int,loo=looni, chains=4) # ### THIS ONE WORKS!!! abit faster that the one on top # and actually removes all pareto >0.7
loosi4<-reloo(zinb_a_lice_brms_bayes_sociality_interactions, loo=loosi, chains=4)

loo_compare(loon, loosi)



# #### ### #### ## 2 ** Model abundance excluding zeros Lice -------------------------------
ectos_birds_dff<-read.csv("data/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, total_lice,year_seasonality ) %>% 
   filter(species_jetz!="Premnoplex_brunnescens") %>% 
  filter (total_lice!=0) %>% 
  na.omit()
  
#%>% filter(total_lice<60)  # removing outliers 

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

# MODEL 
#prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("cauchy(0,2.5)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled

#prior_intercept # we will just keep the default
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

bayes_R2(NZ_a_lice_brms_bayes_no_int)
prior_summary(NZ_a_lice_brms_bayes_no_int)

NZ_a_lice_brms_bayes_no_int_prior<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                           (1|gr(species_jetz, cov = phy_cov))+ 
                                           (1|Powder.lvl)+
                                           (1|species),
                                         data=ectos_birds_dff,
                                         family=poisson(),  #zero_inflated_negbinomial()
                                         data2 = list(phy_cov=phy_cov),
                                         prior = c(prior_predictors, prior_random),
                                         save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
                                         iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                         thin=2,
                                         control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(NZ_a_lice_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1P_NZ.model_LICE_ABUNDANCE_brms_no_interactions_priors.RDS")
NZ_a_lice_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1P_NZ.model_LICE_ABUNDANCE_brms_no_interactions_priors.RDS")

prior_summary(NZ_a_lice_brms_bayes_no_int_prior_nb)
NZ_a_lice_brms_bayes_no_int_prior_nb<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                               (1|gr(species_jetz, cov = phy_cov))+ 
                                               (1|Powder.lvl)+
                                               (1|species),
                                             data=ectos_birds_dff,
                                             family=negbinomial(),  #zero_inflated_negbinomial()
                                             data2 = list(phy_cov=phy_cov),
                                             prior = c(prior_predictors, prior_random),
                                             save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
                                             iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                             thin=2,
                                             control=list(adapt_delta=0.99, max_treedepth=14)) 

mcmc_plo

#saveRDS(NZ_a_lice_brms_bayes_no_int_prior_nb, "data/data_analyses/model_selection/M1L_NZ_NB.model_LICE_ABUNDANCE_brms_no_interactions_priors.RDS")
NZ_a_lice_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1l_NZ_NB.model_LICE_ABUNDANCE_brms_no_interactions_priors.RDS")

bayes_R2(NZ_a_lice_brms_bayes_no_int_prior_nb)
# Also see withnormal priors and cauchy 
bayes_R2(NZ_a_lice_brms_bayes_no_int_prior_nb_normal) 

# Some model exploration 
class(NZ_a_lice_brms_bayes_no_int_prior)
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(NZ_a_lice_brms_bayes_no_int_prior_nb, integer = TRUE)
plot(simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations

#  model comparison and overfitting 
bayes_R2(NZ_a_lice_brms_bayes_no_int_prior)
bayes_R2(NZ_a_lice_brms_bayes_no_int_prior_nb)

loo_nz_lice<-loo(NZ_a_lice_brms_bayes_no_int_prior, moment_match = TRUE)
k_ecto_NZ_lice_brms_no_int_prior<-kfold(NZ_a_lice_brms_bayes_no_int_prior, K=10)
saveRDS(k_ecto_NZ_lice_brms_no_int_prior, "data/data_analyses/model_selection/k_fold/K_fold_M1P_NZ.model_LICE_ABUNDANCE_brms_no_interactions_priors.RDS")

loo(NZ_a_lice_brms_bayes_no_int_prior_nb)
loo_nz_lice_nb<-loo(NZ_a_lice_brms_bayes_no_int_prior_nb, moment_match = TRUE)
k_ecto_NZ_lice_brms_no_int_prior_nb<-kfold(NZ_a_lice_brms_bayes_no_int_prior_nb, K=10)
saveRDS(k_ecto_NZ_lice_brms_no_int_prior_nb, "data/data_analyses/model_selection/k_fold/K_fold_M1P_NZ.model_LICE_ABUNDANCE_brms_no_interactions_priors_nb.RDS")

loo_compare(loo_nz_lice,loo_nz_lice_nb)
loo_compare(k_ecto_NZ_lice_brms_no_int_prior_nb,k_ecto_NZ_lice_brms_no_int_prior)

mcmc_plot(NZ_a_lice_brms_bayes_no_int_prior_nb)
mcmc_plot(NZ_a_lice_brms_bayes_no_int_prior)

# plots 

color_scheme_set("red") 

# poisson 

#model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig1_ALNZ_NB.plot_model_CONVERGENCE_LICE_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
plot(NZ_a_lice_brms_bayes_no_int_prior)
dev.off()

# model fit
# model fit
png("figures/figures_manuscript/models_selected_figures/Fig1_ALNZ_NB.plot_modell_FIT__ecto_LICE_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(NZ_a_lice_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(NZ_a_lice_brms_bayes_no_int_prior_nb,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions NZ_LICE ABUNDANCE", subtitle ="NZ_LICE_ABUNDANCE_brms_bayes_no_int_prior with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1_ALNZ_NB.plot_model_parameters_LICE_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(NZ_a_lice_brms_bayes_no_int_prior_nb,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions NZ_LICE ABUNDANCE", subtitle ="NZ_LICE_ABUNDANCE with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1_ALNZ_NB.plot_model_parameters_intervals_LICE_ABUNDANCE_brms_bayes_no_int_.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



###_###_###_###_###_###_###_###
####### fITTING THE MODEL WITH PGLMM gave up on this approach!!
###_###_###_###_###_###_###_###


# zero inflated?
#b) model pglmm
###_###_###

NZ_a_lice_pglmm<-phyr::pglmm(total_lice~sociality+ scale(elevation)+ 
                               scale(year_seasonality)+
                               (1|species_jetz__)+
                                (1|Powder.lvl),
                             data=ectos_birds_dff,
                             family="poisson",  #zero_inflated_negbinomial()
                             cov_ranef = list(species_jetz= phylo), #class phylo
                             add.obs.re = TRUE,
                             REML = TRUE, 
                             verbose = TRUE,
                             s2.init = .25) # what is this last parameter for


histogram(ectos_birds_dff$total_nf_mites) # id some outliers 

summary(NZ_a_lice_pglmm)
rr2::R2(NZ_a_lice_pglmm)
fixef(nf_mites_a_pglmm)
predict(ecto_a_pglmm)

# Assumptions check

simulationOutput<- DHARMa::simulateResiduals(fittedModel=NZ_a_lice_pglmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput)
testUniformity(simulationOutput) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput) ## tests if there are more zeros in the data than expected from the simulations




# ##### 2.3 Selected model Abundance Lice(ZIP vs ZINB) ------------------------------------

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

ZIP_a_lice_brms_bayes_no_int_priors<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                  (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                  (1|Powder.lvl)+
                                                  (1|species),
                                                data=ectos_birds_dff,
                                                family=zero_inflated_poisson(),  #zero_inflated_negbinomial()
                                                data2 = list(phy_cov=phy_cov),
                                                prior = c(prior_predictors,prior_random,prior_intercept,residual_prior2),
                                                #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                                iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(ZIP_a_lice_brms_bayes_no_int_priors, "data/data_analyses/model_selection/1.ZIP_model_ABUNDANCE_LICE_brms_multiple_obs_all_interactions_priors.RDS")
ZIP_a_lice_brms_bayes_no_int_priors<-readRDS()

##_###_###
##_##_Besty_##_##
##_###_###
zinb_a_lice_brms_bayes_no_int_priors<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                  (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                  (1|Powder.lvl)+
                                                  (1|species),
                                                data=ectos_birds_dff,
                                                family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                data2 = list(phy_cov=phy_cov),
                                                prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                                iter=10000, warmup=5000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(zinb_a_lice_brms_bayes_no_int_priors, "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")
zinb_a_lice_brms_bayes_no_int_priors<-readRDS( "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")


names(ectos_birds_dff)

selected_zinb_a_lice_brms_bayes_no_int_priors_mass<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                                scale(mass_tidy_species)+mass_ind_comp+
                                                  (1|gr(species_jetz, cov = phy_cov))+  (1|species)+
                                                  (1|Powder.lvl)+
                                                data=ectos_birds_dff,
                                                family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                data2 = list(phy_cov=phy_cov),
                                                prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                                iter=5000, warmup=2500, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                control=list(adapt_delta=0.99, max_treedepth=14))
saveRDS(zinb_a_lice_brms_bayes_no_int_priors_mass, "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_mass.RDS")


zinb_a_lice_brms_bayes_no_int_priors_mass_ind<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+scale(mass_individual_day)+
                                                       (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                       (1|Powder.lvl)+
                                                       (1|species),
                                                     data=ectos_birds_dff,
                                                     family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                     data2 = list(phy_cov=phy_cov),
                                                     prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                     #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                                     iter=5000, warmup=2500, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                     thin=2,
                                                     control=list(adapt_delta=0.99, max_treedepth=14))



mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors_mass)

#explore model structure
simulate_residuals <- dh_check_brms(ZIP_a_lice_brms_bayes_no_int_priors, integer = TRUE)
plot(simulate_residuals, form = dff_ectos_network_individual_metrics$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations

#MODEL COMPARISSOM w
loo(ZIP_a_lice_brms_bayes_no_int_priors, zinb_a_lice_brms_bayes_no_int_priors, compare = TRUE)
loonip_nb<-loo(zinb_a_lice_brms_bayes_no_int_priors)
loonip_p<-loo(ZIP_a_lice_brms_bayes_no_int_priors)
loo_compare(loonip_nb, loonip_p)

#arrning 20 pareto k are high!!!! what to do? 
# Use k-fold cross validation instead
loo(zinb_a_lice_brms_bayes_no_int_priors_mass,zinb_a_lice_brms_bayes_no_int_priors, compare=TRUE)


k_ZIP_a_lice_no_int_prior<-kfold(ZIP_a_lice_brms_bayes_no_int_priors, K=10)
saveRDS(k_ZIP_a_lice_no_int_prior, "data/data_analyses/model_selection/k_fold/K_fold_1.ZIP_model_ABUNDANCE_LICE_brms_multiple_obs_all_interactions_priors_poisson.RDS")

k_zinb_a_lice_no_int<-kfold(zinb_a_lice_brms_bayes_no_int_priors, K=10)
saveRDS(k_zinb_a_lice_no_int, "data/data_analyses/model_selection/k_fold/K_fold_1_zinb_model_ABUNDANCE_LICE_brms_multiple_obs_all_interactions_priors_zinb.RDS")

k_zinb_a_lice_no_int_mass<-kfold(zinb_a_lice_brms_bayes_no_int_priors_mass, K=10)
saveRDS(k_zinb_a_lice_no_int_mass, "data/data_analyses/model_selection/k_fold/K_fold_1_zinb_model_ABUNDANCE_LICE_brms_multiple_obs_all_interactions_priors_zinb_poisson.RDS")

loo_compare(k_zinb_a_lice_no_int_mass, k_zinb_a_lice_no_int) # compare using elpd_diff


#looni1<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, moment_match = TRUE)
#looni2<-loo_moment_match(zinb_a_nf_mites_brms_bayes_no_int_outliers, loo=loo1,k_threshold = 0.7)
#looni3<-loo(zinb_a_lice_brms_bayes_no_int, reloo = TRUE) # exclude the observation with high pareto to allow me to compare with other models 
#looni4<-reloo(zinb_a_lice_brms_bayes_no_int,loo=looni, chains=4) # ### THIS ONE WORKS!!! abit faster that the one on top # and actually removes all pareto >0.7
#loosi4<-reloo(zinb_a_lice_brms_bayes_sociality_interactions, loo=loosi, chains=4)


zinb_a_lice_brms_bayes_no_int_priors_nopowder<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                  (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                  (1|species),
                                                data=ectos_birds_dff,
                                                family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                data2 = list(phy_cov=phy_cov),
                                                prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                save_pars = save_pars(all=  TRUE), # if i need to use moment match but makes the model heavier
                                                iter=4000, warmup=2000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                control=list(adapt_delta=0.99, max_treedepth=14))
saveRDS(zinb_a_lice_brms_bayes_no_int_priors_nopowder, "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_nopowder.RDS")
mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors_nopowder)

loo(zinb_a_lice_brms_bayes_no_int_priors_nopowder,zinb_a_lice_brms_bayes_no_int_priors, compare=TRUE)
###_###_###_##
#PLOTS
###_###_###_##

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("teal")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/best_models/Fig3L_BEST_zinb_a_lice_brms_bayes_no_int_priors_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_lice_brms_bayes_no_int_priors)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig3L_BEST_zinb_a_lice_brms_bayes_no_int_priors_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_lice_brms_bayes_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

# Choose depending on the model type
#png("figures/figures_manuscript/mT.png",width = 3000, height = 3000, res = 300, units = "px")
#pp_m<- brms::posterior_predict(MODEL)
#ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
#  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
#dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

  
color_scheme_set("teal")

estimates_plot<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ZINB LICE ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig3L_BEST_ZINB_LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ZINB LICE ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig3L_BEST_ZINB_LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES_INTERVALS.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# model check  example

yrepnzb <- posterior_predict(brm_glmznb)
(prop_zero_test4 <- ppc_stat(y=roaches$y, yrepnzb, stat=function(y) mean(y==0)))
(max_test_nb <- pp_check(stan_glmnb, plotfun = "stat", stat = "max"))



# ABUNDANCE ###### 2.4 ***Selected*** model abundance lice included body mass ---------------

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default so no need to speciied
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

##_###_###
##_##_Besty_##_##
##_###_###

selected_zinb_a_lice_brms_bayes_no_int_priors<-brms::brm(total_lice~sociality_groups+ scale(elevation)+ scale(year_seasonality)+
                                                            scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                            (1|gr(species_jetz, cov = phy_cov))+ (1|species)+
                                                            (1|Powder.lvl),
                                                          data=ectos_birds_dff,
                                                          family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                          data2 = list(phy_cov=phy_cov),
                                                          prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                          save_pars = save_pars(all=  TRUE),# if i need to use moment match but makes the model heavier
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          thin=2,
                                                          control=list(adapt_delta=0.99, max_treedepth=14))


saveRDS(selected_zinb_a_lice_brms_bayes_no_int_priors, "results/selected_models/3_M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_ind_mass_scaled_SELECTED_antbirds_included.RDS")


# selected_zinb_a_lice_brms_bayes_no_int_priors<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
#                                                            scale(mass_tidy_species)+(mass_ind_comp)+
#                                                        (1|gr(species_jetz, cov = phy_cov))+ (1|species)+
#                                                        (1|Powder.lvl),
#                                                      data=ectos_birds_dff,
#                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
#                                                      data2 = list(phy_cov=phy_cov),
#                                                      prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
#                                                      save_pars = save_pars(all=  TRUE),# if i need to use moment match but makes the model heavier
#                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                                      thin=2,
#                                                      control=list(adapt_delta=0.99, max_treedepth=14))
# 
# 
# saveRDS(selected_zinb_a_lice_brms_bayes_no_int_priors, "results/selected_models/3_M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_SELECTED.RDS")
# selected_zinb_a_lice_brms_bayes_no_int_priors<-readRDS( "results/selected_models/3_M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_SELECTED.RDS")
# 
# #_##_Besty_##_##
# ##_###_###
# 
# selected_zinb_a_lice_brms_bayes_no_int_priors2<-brms::brm(total_lice~sociality+ scale(elevation)+ scale(year_seasonality)+
#                                                            scale(mass_tidy_species)+scale(mass_ind_comp)+
#                                                            (1|gr(species_jetz, cov = phy_cov))+ (1|species)+
#                                                            (1|Powder.lvl),
#                                                          data=ectos_birds_dff,
#                                                          family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
#                                                          data2 = list(phy_cov=phy_cov),
#                                                          prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
#                                                          save_pars = save_pars(all=  TRUE),# if i need to use moment match but makes the model heavier
#                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                                          thin=2,
#                                                          control=list(adapt_delta=0.99, max_treedepth=14))
# 
# 
# saveRDS(selected_zinb_a_lice_brms_bayes_no_int_priors2, "results/selected_models/3_M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_ind_mass_scaled_SELECTED.RDS")
# 
# loo(selected_zinb_a_lice_brms_bayes_no_int_priors2,selected_zinb_a_lice_brms_bayes_no_int_priors, compare=TRUE, moment_match = TRUE)
# ###_###_###_##
#PLOTS
###_###_###_##

marginal_effects(selected_zinb_a_lice_brms_bayes_no_int_priors)
plot( conditional_effects(selected_zinb_a_lice_brms_bayes_no_int_priors), 
      points = TRUE, 
      point_args = list(width = .05, shape = 1))

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("teal")

# model convergence 
png("results/selected_models_figures/3_Fig3L_BEST_zinb_a_lice_brms_bayes_no_int_priors_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_zinb_a_lice_brms_bayes_no_int_priors)
dev.off()

# model fit
png("results/selected_models_figures/3_Fig3L_BEST_zinb_a_lice_brms_bayes_no_int_priors_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_zinb_a_lice_brms_bayes_no_int_priors2, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

# Choose depending on the model type
#png("figures/figures_manuscript/mT.png",width = 3000, height = 3000, res = 300, units = "px")
#pp_m<- brms::posterior_predict(MODEL)
#ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
#  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
#dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


color_scheme_set("teal")

estimates_plot<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ZINB LICE ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("results/selected_models_figures/3_Fig3L_BEST_ZINB_LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ZINB LICE ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("results/selected_models_figures/3_Fig3L_BEST_ZINB_LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES_INTERVALS.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# model check  example

yrepnzb <- posterior_predict(brm_glmznb)
(prop_zero_test4 <- ppc_stat(y=roaches$y, yrepnzb, stat=function(y) mean(y==0)))
(max_test_nb <- pp_check(stan_glmnb, plotfun = "stat", stat = "max"))


# ##### 3.1.Data processing abundance mites ----------------------------------------------------

ectos_birds_dff<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality,sociality_groups,total_mites, total_mesostigmatidae, total_no_feathers_mites,year_seasonality, mass_tidy_species, mass_ind_tidy,mass_ind_comp ) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens") %>% filter(total_no_feathers_mites<61)# removing outliers 

View(ectos_birds_dff)

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$sociality_groups<-as.factor(ectos_birds_dff$sociality_groups)

ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$total_mites<-as.numeric(ectos_birds_dff$total_mites)
ectos_birds_dff$total_mesostigmatidae<-as.numeric(ectos_birds_dff$total_mesostigmatidae)
ectos_birds_dff$total_no_feathers_mites<-as.numeric(ectos_birds_dff$total_no_feathers_mites)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one
ectos_birds_dff$species<-as.factor(ectos_birds_dff$species)
ectos_birds_dff$mass_tidy_species<-as.numeric(ectos_birds_dff$mass_tidy_species)
ectos_birds_dff$mass_ind_tidy<-as.numeric(ectos_birds_dff$mass_ind_tidy)
ectos_birds_dff$mass_ind_comp<-as.numeric(ectos_birds_dff$mass_ind_comp)

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

# ##### 3.2.Model selection abundance mites --------------------------------------------------------
prior_summary(zinb_a_nf_mites_brms_bayes_no_int)
zinb_a_nf_mites_brms_bayes_no_int<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                     (1|Powder.lvl)+
                                                     (1|species),
                                                   data=ectos_birds_dff,
                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   #prior=c(random_prior, intercept_prior, residual_prior, residual_prior_2), # I am keping the flat priors for the fixed effects
                                                   #save_pars = save_pars(all=  TRUE),
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

loo(zinb_a_nf_mites_brms_bayes_no_int)

#saveRDS(zinb_a_nf_mites_brms_bayes_no_int, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")
zinb_a_nf_mites_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")

prior_summary(zinb_a_nf_mites_brms_bayes_sociality_interactions)

zinb_a_nf_mites_brms_bayes_sociality_interactions<-brms::brm(total_no_feathers_mites~
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
                                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                        thin=2,
                                                        #prior=c(random_prior, intercept_prior, residual_prior, residual_prior_2), # I am keping the flat priors for the fixed effects
                                                        #save_pars = save_pars(all=  TRUE),
                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(zinb_a_nf_mites_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")
zinb_a_nf_mites_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")

zinb_a_nf_mites_brms_bayes_all_interactions<-brms::brm(total_no_feathers_mites~
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
                                                  iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                  thin=2,
                                                  control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(zinb_a_nf_mites_brms_bayes_all_interactions, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")
zinb_a_nf_mites_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")

#Priors

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int)
mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior)

View(ectos_birds_dff)


zinb_a_nf_mites_brms_bayes_no_int_prior<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                               (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                               (1|Powder.lvl)+
                                               (1|species),
                                             data=ectos_birds_dff,
                                             family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                             data2 = list(phy_cov=phy_cov),
                                             iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                             thin=2,
                                             prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                             #save_pars = save_pars(all=  TRUE),
                                             control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zinb_a_nf_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")
zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")

zinb_a_nf_mites_brms_bayes_sociality_interactions_prior<-brms::brm(total_no_feathers_mites~
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
                                                             iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                             thin=2,
                                                             prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                             #save_pars = save_pars(all=  TRUE),
                                                             control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(zinb_a_nf_mites_brms_bayes_sociality_interactions_prior, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions_prior.RDS")
zinb_a_nf_mites_brms_bayes_sociality_interactions_prior<-readRDS ("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions_prior.RDS")




zinb_a_nf_mites_brms_bayes_all_interactions_prior<-brms::brm(total_no_feathers_mites~
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
                                                      prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
zinb_a_nf_mites_brms_bayes_all_interactions_prior=zip_a_nf_mites_brms_bayes_all_interactions_prior

#saveRDS(zinb_a_nf_mites_brms_bayes_all_interactions_prior, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions_prior.RDS")
zinb_a_nf_mites_brms_bayes_all_interactions_prior<-readRDS( "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions_prior.RDS")


bayes_R2(zinb_a_nf_mites_brms_bayes_no_int_prior_hurdle)
loo(zinb_a_nf_mites_brms_bayes_no_int_prior_hurdle, zinb_a_nf_mites_brms_bayes_no_int_prior, compare = TRUE)
pp_check(zinb_a_nf_mites_brms_bayes_no_int_prior_hurdle, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior_hurdle)


prior_summary(zinb_a_nf_mites_brms_bayes_no_int_prior)
loo(zinb_a_nf_mites_brms_bayes_no_int_prior,zinb_a_nf_mites_brms_bayes_sociality_interactions_prior,zinb_a_nf_mites_brms_bayes_all_interactions_prior)

loo1<-loo(zinb_a_nf_mites_brms_bayes_no_int_prior)
loo2<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, moment_match = TRUE)
loo3<-loo_moment_match(zinb_a_nf_mites_brms_bayes_no_int_prior, loo=loo1,k_threshold = 0.7)
loo4<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, reloo = TRUE) # exclude the observation with high pareto to allow me to compare with other models 
loo5<-reloo(zinb_a_nf_mites_brms_bayes_no_int_prior,loo=loo1, chains=2) # ### THIS ONE WORKS abit faster that the one on top # actually removes all pareto >0.7

###_###_###_##
#MODEL COMPARISON
###_###_###_##

#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.
bayes_R2(zinb_a_nf_mites_brms_bayes_no_int_prior)
bayes_R2(zinb_a_nf_mites_brms_bayes_sociality_interactions_prior)
bayes_R2(zip_a_nf_mites_brms_bayes_all_interactions_prior)

# use loo cross validation [ In this case we are not interested inteh interaaction between the different factors because the sociality of a speceis does not change with seasonality or elevation, cause those are attributes to the individuals does not change ]
# R= the negaive binomial with non interactions is better that the hrdle model and models with ineteractions
loo(zinb_a_nf_mites_brms_bayes_no_int_prior, zinb_a_nf_mites_brms_bayes_sociality_interactions_prior, zip_a_nf_mites_brms_bayes_all_interactions_prior,compare=TRUE)
loo(zinb_a_nf_mites_brms_bayes_no_int, zinb_a_nf_mites_brms_bayes_sociality_interactions, zip_a_nf_mites_brms_bayes_all_interactions,compare=TRUE)
loo(zinb_a_nf_mites_brms_bayes_no_int_prior, zinb_a_nf_mites_brms_bayes_sociality_interactions_prior,compare=TRUE)
# R=eld_diff<4 so we keep the simplest model !
loo_compare(waic(zinb_a_lice_brms_bayes_no_int_priors), waic(zinb_a_lice_brms_bayes_sociality_interactions_priors)) # interesting warning
mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior)
# use k-fold-cv validation instead because in models with random efects loo tends to fail # but this is taking forevwer so i WILL RUNT IT AT UBC TOMORROW
k_ecto_p_brms_no_int<-kfold(zinb_a_nf_mites_brms_bayes_no_int_prior, K=10)
k_mites_p_brms_no_int_zip<-kfold(ecto_p_brms_bayes_sociality_interactions_priors, K=10)
loo_compare(k_ecto_p_brms_no_int, k_ecto_p_brms_sociality_int)

# Model posterior predictive checks 
#The idea of posterior predictive checks is to compare our observed data to replicated data from the model. 
#If our model is a good fit, we should be able to use it to generate a dataset that resembles the observed data.
#gettigsamples from the posterior predictive distribution:
pp_check(zinb_a_nf_mites_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100) +xlim(0,10)

#or dependeing or the data use
#pp_check(ecto_p_brms_bayes_no_int_prior, type = "stat", stat = 'median', nsamples = 100)
#pp_m<- brms::posterior_predict(ecto_p_brms_bayes_no_int_prior)
#ppc_rootogram(y=ecto_p_brms_bayes_no_int_prior$data$ectoparasites_PA, pp_m[1:200, ])  +   
coord_cartesian(xlim = c(0, 100), ylim = c(0,30))

###_###_###_##
#PLOTS
###_###_###_##

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2(zinb_a_nf_mites_brms_bayes_no_int_prior) # R2 0.1529

color_scheme_set("green")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig1M.ZINB_ABUNDANCE_MITES_brms_bayes_no_int_priors_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_nf_mites_brms_bayes_no_int_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig1M.ZINB_ABUNDANCE_MITES_brms_bayes_no_int_priors_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_nf_mites_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

# Choose depending on the model type
#png("figures/figures_manuscript/mT.png",width = 3000, height = 3000, res = 300, units = "px")
#pp_m<- brms::posterior_predict(MODEL)
#ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
#  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
#dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

color_scheme_set("green")

estimates_plot<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions ABUNDANCE_MITES", subtitle ="ABUNDANCE_MITES with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1M.ZINB_ABUNDANCE_MITES_brms_bayes_no_int_prior_ESTIMATES.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions ABUNDANCE_MITES", subtitle ="ABUNDANCE_MITES with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1M.ZINB_ABUNDANCE_MITES_brms_bayes_no_int_prior_ESTIMATES_INTERVALS.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()




# #### ### #### ## 3 ** Model abundance excluding zeros Mites -------------------------------
ectos_birds_dff<-read.csv("data/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality,total_mites, total_mesostigmatidae, total_no_feathers_mites,year_seasonality ) %>% 
  filter(species_jetz!="Premnoplex_brunnescens") %>% 
  filter(total_no_feathers_mites>0) %>% 
  na.omit() 
  
#%>% filter(total_lice<60)# removing outliers 
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

# The model 

#prior_predictors<-prior("student_t(3,0,1)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


NZ_a_mites_brms_bayes_no_int_prior<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                               (1|gr(species_jetz, cov = phy_cov))+
                                                (1|Powder.lvl)+
                                               (1|species),
                                             data=ectos_birds_dff,
                                             family=poisson(),  #zero_inflated_negbinomial()
                                             data2 = list(phy_cov=phy_cov),
                                             prior = c(prior_predictors, prior_random),
                                             save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
                                             iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                             thin=2,
                                             control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(NZ_a_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1M_NZ.model_MITES_ABUNDANCE_brms_no_interactions_priors.RDS")
NZ_a_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1M_NZ.model_MITES_ABUNDANCE_brms_no_interactions_priors.RDS")

NZ_a_mites_brms_bayes_no_int_prior_nb<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                  (1|gr(species_jetz, cov = phy_cov))+ 
                                                  (1|Powder.lvl)+
                                                  (1|species),
                                                data=ectos_birds_dff,
                                                family=negbinomial(),  #zero_inflated_negbinomial()
                                                data2 = list(phy_cov=phy_cov),
                                                prior = c(prior_predictors, prior_random),
                                                save_pars = save_pars(all=  TRUE), #if i need to use moment match but makes the model heavier
                                                iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(NZ_a_mites_brms_bayes_no_int_prior_nb, "data/data_analyses/model_selection/M1M_NZ.model_MITES_ABUNDANCE_brms_no_interactions_priors_nb.RDS")
NZ_a_mites_brms_bayes_no_int_prior_nb<-readRDS("data/data_analyses/model_selection/M1M_NZ.model_MITES_ABUNDANCE_brms_no_interactions_priors_nb.RDS")




# Some model exploration 
class(NZ_a_mites_brms_bayes_no_int_prior)

loo(NZ_a_mites_brms_bayes_no_int_prior)
loo_nz_mite<-loo(NZ_a_mites_brms_bayes_no_int_prior, moment_match = TRUE)
k_ecto_NZ_mites_brms_no_int_prior<-kfold(NZ_a_mites_brms_bayes_no_int_prior, K=10)
saveRDS(k_ecto_NZ_mites_brms_no_int_prior, "data/data_analyses/model_selection/k_fold/K_fold_M1P_NZ.model_MITES_ABUNDANCE_brms_no_interactions_priors.RDS")

loo(NZ_a_mites_brms_bayes_no_int_prior_nb)
loo_nz_mite_nb<-loo(NZ_a_mites_brms_bayes_no_int_prior_nb, moment_match = TRUE)
k_ecto_NZ_mites_brms_no_int_prior_nb<-kfold(NZ_a_mites_brms_bayes_no_int_prior_nb, K=10)
saveRDS(k_ecto_NZ_mites_brms_no_int_prior_nb, "data/data_analyses/model_selection/k_fold/K_fold_M1P_NZ_NB.model_MITES_ABUNDANCE_brms_no_interactions_priors_nb.RDS")

loo_compare(k_ecto_NZ_mites_brms_no_int_prior,k_ecto_NZ_mites_brms_no_int_prior_nb)
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(NZ_a_mites_brms_bayes_no_int_prior_nb, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations

# plots 
color_scheme_set("pink") 

bayes_R2(NZ_a_mites_brms_bayes_no_int_prior_nb)

#model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig1_AMNZ_NB.plot_model_CONVERGENCE_intervals_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
plot(NZ_a_mites_brms_bayes_no_int_prior_nb)
dev.off()

# model fit
# model fit
png("figures/figures_manuscript/models_selected_figures/Fig1_AMNZ_NB.plot_modell_FIT_intervals_ecto_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(NZ_a_mites_brms_bayes_no_int_prior_nb, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

pp_check(NZ_a_mites_brms_bayes_no_int_prior_nb, ndraws = 100)+ xlim(0, 10)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(NZ_a_mites_brms_bayes_no_int_prior_nb, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

mcmc_plot(NZ_a_mites_brms_bayes_no_int_prior_nb)

estimates_plot<-mcmc_plot(NZ_a_mites_brms_bayes_no_int_prior_nb,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions NZ_MITES ABUNDANCE", subtitle ="NZ_MITES_ABUNDANCE_brms_bayes_no_int_prior with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1_AMNZ_NB.plot_model_parameters_MITES_ABUNDANCE_brms_bayes_no_int_nb.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(NZ_a_mites_brms_bayes_no_int_prior_nb,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions NZ_MITES ABUNDANCE", subtitle ="NZ_MITES_ABUNDANCE with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1_AMNZ_NB_plot_model_parameters_intervals_MITES_ABUNDANCE_brms_bayes_no_int_nb.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# ##### 3.3 Selection model Abundance MITES (ZIP and ZINB) ------------------------------------
#Priors

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

prior_summary(ZIP_a_nf_mites_brms_bayes_no_int_prior)

View(ectos_birds_dff)
# Remembet that we removed one outlier with >80 
zinb_a_nf_mites_brms_bayes_no_int_prior<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                     (1|Powder.lvl)+
                                                     (1|species),
                                                   data=ectos_birds_dff,
                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                   #save_pars = save_pars(all=  TRUE),
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 


saveRDS(zinb_a_nf_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")
zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")

zinb_a_nf_mites_brms_bayes_no_int_prior_mass<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+scale(mass_tidy)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                     (1|Powder.lvl)+
                                                     (1|species),
                                                   data=ectos_birds_dff,
                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                   #save_pars = save_pars(all=  TRUE),
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 


saveRDS(zinb_a_nf_mites_brms_bayes_no_int_prior_mass, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_mass.RDS")

mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior_mass)
# 28 DIVERGENT TRANSICTIONS !!! WARNING # ithout outliersun wneed to ren
ZIP_a_nf_mites_brms_bayes_no_int_prior<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                     (1|Powder.lvl)+
                                                     (1|species),
                                                   data=ectos_birds_dff,
                                                   family=zero_inflated_poisson(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   prior = c(prior_predictors,prior_random,residual_prior2),
                                                   #save_pars = save_pars(all=  TRUE),
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(ZIP_a_nf_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1MNF.model_prevalence_ZIP_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")
ZIP_a_nf_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_ZIP_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")




# Excluding powder level 

zinb_a_nf_mites_brms_bayes_no_int_prior_nopowder<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                     (1|species),
                                                   data=ectos_birds_dff,
                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=4000, warmup=2000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                   #save_pars = save_pars(all=  TRUE),
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zinb_a_nf_mites_brms_bayes_no_int_prior_nopowder, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_nopowder.RDS")
zinb_a_nf_mites_brms_bayes_no_int_prior_nopowder<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_nopowder.RDS")



loo(zinb_a_nf_mites_brms_bayes_no_int_prior,zinb_a_nf_mites_brms_bayes_no_int_prior_nopowder, compare=TRUE)

###_###_###_##
#MODEL COMPARISON
###_###_###_##

#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.
bayes_R2()
bayes_R2(ZIP_a_nf_mites_brms_bayes_no_int_prior)

# use loo cross validation [ In this case we are not interested inteh interaaction between the different factors because the sociality of a speceis does not change with seasonality or elevation, cause those are attributes to the individuals does not change ]
loo(zinb_a_nf_mites_brms_bayes_no_int_prior,ZIP_a_nf_mites_brms_bayes_no_int_prior, compare=TRUE)
loozip_mites<-loo(zinb_a_nf_mites_brms_bayes_no_int_prior, moment_match=TRUE)
loozinb_mites<-loo(ZIP_a_nf_mites_brms_bayes_no_int_prior, moment_match=TRUE)

loo_compare(loo_pm,loo_nb)

# use k-fold-cv validation instead because in models with random efFects loo tends to fail # but this is taking forevwer so i WILL RUNT IT AT UBC TOMORROW


k_mites_brms_bayes_no_int_prior_zinb<-kfold(zinb_a_nf_mites_brms_bayes_no_int_prior, K=10)
saveRDS(k_mites_brms_bayes_no_int_prior_zinb, "data/data_analyses/model_selection/k_fold/K_fold_M1MNF_model_MITES_ABUNDANCE_no_interactions_priors_zinb.RDS")
k_mites_brms_bayes_no_int_prior_zinb<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M1MNF_model_MITES_ABUNDANCE_no_interactions_priors_zinb.RDS")

k_mites_brms_bayes_no_int_prior_p<-kfold(ZIP_a_nf_mites_brms_bayes_no_int_prior, K=10)
saveRDS(k_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/k_fold/K_fold_M1MNF_model_MITES_ABUNDANCE_no_interactions_priors_poisson.RDS")
k_mites_brms_bayes_no_int_prior_p<-readRDS()

k_mites_brms_bayes_no_int_prior_mass<-kfold(zinb_a_nf_mites_brms_bayes_no_int_prior_mass, K=10)
saveRDS(k_mites_brms_bayes_no_int_prior_mass, "data/data_analyses/model_selection/k_fold/K_fold_M1MNF_model_MITES_ABUNDANCE_no_interactions_priors_zinb_mass.RDS")
k_mites_brms_bayes_no_int_prior_mass<-readRDS()

loo_compare(k_ecto_p_brms_no_int, k_ecto_p_brms_sociality_int)

# plots 
color_scheme_set("green") 

# poisson 
#model convergence 
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2M_BEST_ZINB_plot_model_CONVERGENCE_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_nf_mites_brms_bayes_no_int_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2M_BEST_ZINB_plot_modell_FIT_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_nf_mites_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES

mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior_outlier)


estimates_plot<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE  ")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2M_BEST_plot_model_parameters_ZINB_MITES ABUNDANCE_brms_bayes_social_int.png",width =4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2M_BEST_plot_model_parameters_intervals_ZINB_MITES_ABUNDANCE_brms_bayes_social_int.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# Choose depending on the model type
#png("figures/figures_manuscript/mT.png",width = 3000, height = 3000, res = 300, units = "px")
#pp_m<- brms::posterior_predict(MODEL)
#ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
#  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
#dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 



# ABUNDANCE ####### 3.4 ***Selected*** model abundance mites included body mass ---------------
#Priors

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

#View(ectos_birds_dff) # Remember that we removed one outliers with >80 

selected_zinb_a_nf_mites_brms_bayes_no_int_prior<-brms::brm(total_no_feathers_mites~sociality_groups+ scale(elevation)+ scale(year_seasonality)+
                                                              scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  (1|species)+
                                                     (1|Powder.lvl),
                                                   data=ectos_birds_dff,
                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=10000, warmup=5000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                   save_pars = save_pars(all=  TRUE),
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, "results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_SELECTED_antfollowers_included.RDS")
selected_zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_SELECTED_antfollowers_included.RDS")


# selected_zinb_a_nf_mites_brms_bayes_no_int_prior2<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
#                                                                scale(mass_tidy_species)+scale(mass_ind_comp)+
#                                                                (1|gr(species_jetz, cov = phy_cov))+  (1|species)+
#                                                                (1|Powder.lvl),
#                                                              data=ectos_birds_dff,
#                                                              family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
#                                                              data2 = list(phy_cov=phy_cov),
#                                                              iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                                              thin=2,
#                                                              prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
#                                                              save_pars = save_pars(all=  TRUE),
#                                                              control=list(adapt_delta=0.99, max_treedepth=14)) 
# 
# saveRDS(selected_zinb_a_nf_mites_brms_bayes_no_int_prior2, "results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_ind_mass_scaled_SELECTED.RDS")
# selected_zinb_a_nf_mites_brms_bayes_no_int_prior2<-readRDS("results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_ind_mass_scaled_SELECTED.RDS")
#         
# bayes_R2(selected_zinb_a_nf_mites_brms_bayes_no_int_prior)
# loo(selected_zinb_a_nf_mites_brms_bayes_no_int_prior,selected_zinb_a_nf_mites_brms_bayes_no_int_prior2, compare = TRUE)
#         
        
  # plots 
color_scheme_set("green") 

# poisson 
#model convergence 
png("results/selected_models_figures/3_Fig2M_BEST_ZINB_plot_model_CONVERGENCE_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_zinb_a_nf_mites_brms_bayes_no_int_prior)
dev.off()

# model fit
png("results/selected_models_figures/3_Fig2M_BEST_ZINB_plot_modell_FIT_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_prior2, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES


estimates_plot<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE  ")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("results/selected_models_figures/3_Fig2M_BEST_plot_model_parameters_ZINB_MITES ABUNDANCE_brms_bayes_no_int.png",width =4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_prior2,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("results/selected_models_figures/3_Fig2M_BEST_plot_model_parameters_intervals_ZINB_MITES_ABUNDANCE_brms_bayes_no_int_scaled_individial_mass.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

### ALL MITES INCLUDED

#Priors

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

View(ectos_birds_dff) # Remember that we removed one outliers with >80 

selected_zinb_a_all_mites_brms_bayes_no_int_prior<-brms::brm(total_mites~sociality_groups+ scale(elevation)+ scale(year_seasonality)+
                                                              scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                              (1|gr(species_jetz, cov = phy_cov))+  (1|species)+
                                                              (1|Powder.lvl),
                                                            data=ectos_birds_dff,
                                                            family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                            data2 = list(phy_cov=phy_cov),
                                                            iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                            thin=2,
                                                            prior = c(prior_predictors,prior_random,prior_intercept,residual_prior,residual_prior2),
                                                            save_pars = save_pars(all=  TRUE),
                                                            control=list(adapt_delta=0.99, max_treedepth=14)) 


saveRDS(selected_zinb_a_all_mites_brms_bayes_no_int_prior, "results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_ALL_MITES_phylo_multiple_obs_no_interactions_prior_SELECTED_antfollowers_included.RDS")
selected_zinb_a_all_mites_brms_bayes_no_int_prior<-readRDS("results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_ALL_MITES_phylo_multiple_obs_no_interactions_prior_SELECTED_antfollowers_included.RDS")



# ##### 4.Data processing NETWORKS prevalence ectos ----------------------------------------------------

#seasonality_date<-read.csv("data/data_manuscript/7.dff_parasites_seasonality.csv")
#dff_ectos_network<-read.csv("data/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>%
#inner_join(read.csv("data/data_manuscript/7.dff_parasites_seasonality.csv"), by="Full_Label")
# write.csv(dff_ectos_network,"data/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv" )

dff_ectos_network_individual_metrics<-read.csv("data/data_manuscript/3_dff_all_ectos_network_metrics_individuals_FILE_TIDY.csv",na.strings =c("","NA"))%>% 
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality,ectoparasites_PA, degree, w_degree, year_seasonality, mass_tidy_species, mass_ind_tidy, mass_ind_comp) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

names(dff_ectos_network_individual_metrics)

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


dff_ectos_network_individual_metrics$species<-as.factor(dff_ectos_network_individual_metrics$species)
dff_ectos_network_individual_metrics$mass_tidy_species<-as.numeric(dff_ectos_network_individual_metrics$mass_tidy_species)
dff_ectos_network_individual_metrics$mass_ind_tidy<-as.numeric(dff_ectos_network_individual_metrics$mass_ind_tidy)
dff_ectos_network_individual_metrics$mass_ind_comp<-as.numeric(dff_ectos_network_individual_metrics$mass_ind_comp)

names(dff_ectos_network_individual_metrics)
is.ultrametric(phylo)


# ##### 4.1.Model selection NETWORKS prevalence ectos DEGREE AND STRENGHT --------------------------------------------------------
###_###_###
#2PND stands for PREVALENCE NETWORK DEGREE
prior_summary(ecto_p_brms_bayes_no_int_degree)
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
                                    save_pars = save_pars(all=  TRUE),
                                    control=list(adapt_delta=0.99, max_treedepth=14)) 
#saveRDS(ecto_p_brms_bayes_no_int_degree, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE.RDS")
ecto_p_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE.RDS")


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
                                                    save_pars = save_pars(all=  TRUE),
                                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                    thin=2,
                                                    control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(ecto_p_brms_bayes_sociality_interactions_degree, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")
ecto_p_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")

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
                                              save_pars = save_pars(all=  TRUE),
                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                              thin=2,
                                              control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(ecto_p_brms_bayes_all_interactions_degree, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree.RDS")
ecto_p_brms_bayes_all_interactions_degree<-readRDS( "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree.RDS")

# PRIORS
#2PND stands for PREVALENCE NETWORK DEGREE

prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_######_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_######_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
# 
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_##_###_###_###_###_###_###_#

ecto_p_brms_bayes_no_int_degree_prior<-brms::brm(ectoparasites_PA~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                             (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                             (1|Powder.lvl)+
                                             (1|species),
                                           data=dff_ectos_network_individual_metrics,
                                           family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                           data2 = list(phy_cov=phy_cov),
                                           iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                           thin=2,
                                           save_pars = save_pars(all=  TRUE),
                                           prior = c(prior_predictors,prior_random, prior_intercept),
                                           control=list(adapt_delta=0.999, max_treedepth=14)) 
#saveRDS(ecto_p_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")
ecto_p_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_####_###_###_###_###_###_###_##_###
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_######_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###

ecto_p_brms_bayes_sociality_interactions_degree_prior<-brms::brm(ectoparasites_PA~
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
                                                           save_pars = save_pars(all=  TRUE),
                                                           prior = c(prior_predictors,prior_random, prior_intercept),
                                                           iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                           thin=2,
                                                           control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(ecto_p_brms_bayes_sociality_interactions_degree_prior, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE_prior.RDS")
#ecto_p_brms_bayes_sociality_interactions_degree_prior<-readRDS( "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE_prior.RDS")

ecto_p_brms_bayes_all_interactions_degree_prior<-brms::brm(ectoparasites_PA~
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
                                                       prior = c(prior_predictors,prior_random, prior_intercept),
                                                     data=dff_ectos_network_individual_metrics,
                                                     family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                     data2 = list(phy_cov=phy_cov),
                                                     save_pars = save_pars(all=  TRUE),
                                                     iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                     thin=2,
                                                     control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(ecto_p_brms_bayes_all_interactions_degree_prior, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree_prior.RDS")
#ecto_p_brms_bayes_all_interactions_degree_prior<-readRDS( "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree_prior.RDS")

ecto_p_brms_bayes_no_int_degree_prior_mass<-brms::brm(ectoparasites_PA~scale(degree)+ scale(elevation)+ scale(year_seasonality)+scale(mass_tidy)+
                                                  (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                  (1|Powder.lvl)+
                                                  (1|species),
                                                data=dff_ectos_network_individual_metrics,
                                                family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                data2 = list(phy_cov=phy_cov),
                                                iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                save_pars = save_pars(all=  TRUE),
                                                prior = c(prior_predictors,prior_random, prior_intercept),
                                                control=list(adapt_delta=0.999, max_treedepth=14)) 

saveRDS(ecto_p_brms_bayes_no_int_degree_prior_mass, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_mass.RDS")

#p
loo(ecto_p_brms_bayes_no_int_degree,ecto_p_brms_bayes_sociality_interactions_degree,ecto_p_brms_bayes_all_interactions_degree,compare=TRUE)

LOO(brms.0, brms.0.RE, brms.noRE, brms.speciesRE, brms.phyloRE, brms.full,brms.reduced.noRE,brms.reduced.speciesRE,brms.reduced.phyloRE,brms.reduced.bothRE,reloo=TRUE)


bayes_R2(ecto_p_brms_bayes_no_int_degree_prior) 
bayes_R2(ecto_p_brms_bayes_sociality_interactions_degree_prior) 
bayes_R2(ecto_p_brms_bayes_all_interactions_degree_prior)

# use loo cross validation 
loo(ecto_p_brms_bayes_no_int_degree_prior, ecto_p_brms_bayes_sociality_interactions_degree_prior,ecto_p_brms_bayes_all_interactions_degree_prior,  moment_match = TRUE,compare=TRUE)
loo(ecto_p_brms_bayes_no_int_degree_prior, moment_match = TRUE)

# There is not large difference between the three models so we keep the simplest model 

# use k-fold-cv validation [since While importance sampling in PSIS-LOO can fail for “random effect” model]
k_zinb_a_lice_brms_bayes_no_int<-kfold(ecto_p_brms_bayes_no_int_degree_prior, K=10)
k_zinb_a_lice_brms_bayes_sociality<-kfold(ecto_p_brms_bayes_sociality_interactions_degree_prior, K=10)
k_zinb_a_lice_brms_bayes_sociality<-kfold(ecto_p_brms_bayes_all_interactions_degree_prior, K=10)

loo_compare(k_zinb_a_lice_brms_bayes_no_int, k_zinb_a_lice_brms_bayes_sociality_interactions)


###_####_###_###_###
#SELECTED MODEL the BESTY !
###_####_###_###_###

loo(ecto_p_brms_bayes_no_int_degree_nopowder,ecto_p_brms_bayes_no_int_degree_prior,compare=TRUE)

#2PND stands for PREVALENCE NETWORK DEGREE
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

mcmc_plot(ecto_p_brms_bayes_no_int_degree)

loo(ecto_p_brms_bayes_no_int_degree_prior,ecto_p_brms_bayes_no_int_degree_prior_nopowder,compare=TRUE)


selected_ecto_p_brms_bayes_no_int_degree_prior<-brms::brm(ectoparasites_PA~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                            (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                            (1|species),
                                                          data=dff_ectos_network_individual_metrics,
                                                          family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          thin=2,
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random, prior_intercept),
                                                          control=list(adapt_delta=0.999, max_treedepth=14))


saveRDS(selected_ecto_p_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree_prior_SELECTED.RDS")
#selected_ecto_p_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree_prior_SELECTED.RDS")


# ######## 4.0. Selected model NETWORKS PREVALENCE degree  (ZINB)  -------

#PLOTS
color_scheme_set("blue")

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

# model convergence 
png("figures/figures_manuscript/models_selected_figures//best_models/Fig2PND_BEST_Bernu_model_CONVERGENCE_intervals_ecto_INFECTION_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_ecto_p_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2PND_BEST_Bernu_plot_model_FIT_intervals_ecto_INFECTION_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_ecto_p_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION~DEGREE")+
  theme_minimal(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2PND_BEST_Bernu_plot_model_parameters_INFECTION_brms_bayes_no_int_DEGREE.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION ~DEGREE ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2PND_BEST_Bernu_plot_model_parameters_intervals_INFECTION_brms_bayes_no_int_degree.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


# ######## 4.0. SelecTING model NETWORKS PREVALENCE STRENGHT  (ZINB)  -------
###_###_###_###_###_##_###_###_###_###_
#### Prevalence VS strenght
###_###_###_###_###_##_###_###_###_###_

prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?



names(dff_ectos_network_individual_metrics)
ectos_infection_brms_bayes_no_int_strength_prior<-brms::brm(ectoparasites_PA~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                (1|Powder.lvl)+
                                                                (1|species),
                                                              data=dff_ectos_network_individual_metrics,
                                                              family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                              data2 = list(phy_cov=phy_cov),
                                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                              thin=2,
                                                              save_pars = save_pars(all=  TRUE),
                                                              prior = c(prior_predictors,prior_random,prior_intercept),
                                                              control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(ectos_infection_brms_bayes_no_int_strength_prior, "data/data_analyses/model_selection/M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_all_interactions_STRENGTH_prior.RDS")
ectos_infection_brms_bayes_no_int_strength_prior<-readRDS( "data/data_analyses/model_selection/M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_all_interactions_STRENGTH_prior.RDS")

# SElected model does not include powder level cause it does not significantly improve model fit
selected_ectos_infection_brms_bayes_no_int_strength_prior<-brms::brm(ectoparasites_PA~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                (1|species),
                                                              data=dff_ectos_network_individual_metrics,
                                                              family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                              data2 = list(phy_cov=phy_cov),
                                                              iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                              thin=2,
                                                              save_pars = save_pars(all=  TRUE),
                                                              prior = c(prior_predictors,prior_random,prior_intercept),
                                                              control=list(adapt_delta=0.99, max_treedepth=14))
saveRDS(selected_ectos_infection_brms_bayes_no_int_strength_prior, "data/data_analyses/model_selection/M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_all_interactions_STRENGTH_prior_SELECTED.RDS")
selected_ectos_infection_brms_bayes_no_int_strength_prior<-readRDS( "data/data_analyses/model_selection/M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_all_interactions_STRENGTH_prior_SELECTED.RDS")


loo(ectos_infection_brms_bayes_no_int_strength_prior,selected_ectos_infection_brms_bayes_no_int_strength_prior, compare=TRUE, moment_match = TRUE)
mcmc_plot(ectos_infection_brms_bayes_no_int_strength_prior_nopowder)

bayes_R2(ectos_infection_brms_bayes_no_int_strength_prior)
bayes_R2(selected_ectos_infection_brms_bayes_no_int_strength_prior)

#PLOTS
color_scheme_set("blue")

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

# model convergence 
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2PNS_BEST_Bernu_model_CONVERGENCE_intervals_ecto_INFECTION_brms_bayes_no_int_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_ectos_infection_brms_bayes_no_int_strength_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2PNS_BEST_Bernu_plot_model_FIT_intervals_ecto_INFECTION_brms_bayes_no_int_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_ectos_infection_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 100)+ xlim(0,5)
dev.off()


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
estimates_plot<-mcmc_plot(selected_ectos_infection_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION~STRENGTH")+
  theme_minimal(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2PNS_BEST_Bernu_plot_model_parameters_INFECTION_brms_bayes_no_int_STRENGTH.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_ectos_infection_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION ~STRENGTH ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2PNS_BEST_Bernu_plot_model_parameters_intervals_INFECTION_brms_bayes_no_int_STRENGTH.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


mcmc_plot(zinb_prevalence_brms_bayes_no_int_W_degree_prior)
msms_zinb_prevalence_brms_bayes_no_int_W_degree_prior()

# 4.2 ###### **** Selected**** model Infection Networks Degree and strength ----------

###_###_###_###_###_##_###_###_###_###_
#### Prevalence VS Degree
###_###_###_###_###_##_###_###_###_###_

prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

mcmc_plot(ecto_p_brms_bayes_no_int_degree)

names(dff_ectos_network_individual_metrics)

selected_ecto_p_brms_bayes_no_int_degree_prior<-brms::brm(ectoparasites_PA~scale(degree)+ scale(elevation)+ scale(year_seasonality)+ 
                                                            scale(mass_tidy_species)+(mass_ind_comp)+
                                                            (1|gr(species_jetz, cov = phy_cov))+(1|species), #(1|Powder.lvl)
                                                          data=dff_ectos_network_individual_metrics,
                                                          family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          thin=2,
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random, prior_intercept),
                                                          control=list(adapt_delta=0.999, max_treedepth=14))


saveRDS(selected_ecto_p_brms_bayes_no_int_degree_prior, "results/selected_models/4_1_M2PND_model_INFECTION_b_brms_phylo_multiple_obs_all_interactions_degree_prior_SELECTED.RDS")
#selected_ecto_p_brms_bayes_no_int_degree_prior<-readRDS("results/selected_models/4_1_M2PND_model_prevalence_b_brms_phylo_multiple_obs_all_interactions_degree_prior_SELECTED.RDS")

selected_ecto_p_brms_bayes_no_int_degree_prior2<-brms::brm(ectoparasites_PA~scale(degree)+ scale(elevation)+ scale(year_seasonality)+ 
                                                            scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                            (1|gr(species_jetz, cov = phy_cov))+(1|species), #(1|Powder.lvl)
                                                          data=dff_ectos_network_individual_metrics,
                                                          family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          thin=2,
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random, prior_intercept),
                                                          control=list(adapt_delta=0.999, max_treedepth=14))


saveRDS(selected_ecto_p_brms_bayes_no_int_degree_prior2, "results/selected_models/4_1_M2PND_model_INFECTION_b_brms_phylo_multiple_obs_no_interactions_degree_prior_SELECTED_ind_mass_scaled.RDS")


#PLOTS
color_scheme_set("blue")

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

# model convergence 
png("results/selected_models_figures/4_1_Fig4PND_BEST_Bernu_model_CONVERGENCE_intervals_ecto_INFECTION_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_ecto_p_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("results/selected_models_figures/4_1_Fig4PND_BEST_Bernu_plot_model_FIT_intervals_ecto_INFECTION_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_ecto_p_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION~DEGREE")+
  theme_minimal(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("results/selected_models_figures/4_1_Fig4PND_BEST_Bernu_plot_model_parameters_INFECTION_brms_bayes_no_int_DEGREE.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION ~DEGREE ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_1_Fig4PND_BEST_Bernu_plot_model_parameters_intervals_INFECTION_brms_bayes_no_int_degree.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


###_###_###_###_###_##_###_###_###_###_
#### Prevalence VS strenght
###_###_###_###_###_##_###_###_###_###_

prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

selected_ecto_p_brms_bayes_no_int_strength_prior<-brms::brm(ectoparasites_PA~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+ 
                                                            scale(mass_tidy_species)+(mass_ind_comp)+
                                                            (1|gr(species_jetz, cov = phy_cov))+(1|species), #(1|Powder.lvl)
                                                          data=dff_ectos_network_individual_metrics,
                                                          family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          thin=2,
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random, prior_intercept),
                                                          control=list(adapt_delta=0.999, max_treedepth=14))

saveRDS(selected_ecto_p_brms_bayes_no_int_strength_prior, "results/selected_models/4_1_2_M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior_SELECTED.RDS")
#selected_ectos_infection_brms_bayes_no_int_strength_prior<-readRDS( "results/selected_models/4_1_2_M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_all_interactions_STRENGTH_prior_SELECTED.RDS")


selected_ecto_p_brms_bayes_no_int_strength_prior2<-brms::brm(ectoparasites_PA~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+ 
                                                              scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                              (1|gr(species_jetz, cov = phy_cov))+(1|species), #(1|Powder.lvl)
                                                            data=dff_ectos_network_individual_metrics,
                                                            family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
                                                            data2 = list(phy_cov=phy_cov),
                                                            iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                            thin=2,
                                                            save_pars = save_pars(all=  TRUE),
                                                            prior = c(prior_predictors,prior_random, prior_intercept),
                                                            control=list(adapt_delta=0.999, max_treedepth=14))

saveRDS(selected_ecto_p_brms_bayes_no_int_strength_prior2, "results/selected_models/4_1_2_M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior_SELECTED_ind_mass_scaled.RDS")

bayes_R2(ectos_infection_brms_bayes_no_int_strength_prior)
bayes_R2(selected_ectos_infection_brms_bayes_no_int_strength_prior)

#PLOTS
color_scheme_set("blue")

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

# model convergence 
png("results/selected_models_figures/4_1_2_Fig4PNS_BEST_Bernu_model_CONVERGENCE_intervals_ecto_INFECTION_brms_bayes_no_int_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_ecto_p_brms_bayes_no_int_strength_prior)
dev.off()

# model fit
png("results/selected_models_figures/4_1_2_Fig4PNS_BEST_Bernu_plot_model_FIT_intervals_ecto_INFECTION_brms_bayes_no_int_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_ecto_p_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 100)+ xlim(0,5)
dev.off()


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
estimates_plot<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION~STRENGTH")+
  theme_minimal(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_1_2_Fig4PNS_BEST_Bernu_plot_model_parameters_INFECTION_brms_bayes_no_int_STRENGTH.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_strength_prior2,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION ~STRENGTH ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_1_2_Fig4PNS_BEST_Bernu_plot_model_parameters_intervals_INFECTION_brms_bayes_no_int_STRENGTH_ind_mass_scaled.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


# # ##### 4.1.1 Model selection NETWORKS prevalence ectos species  --------

dff_ectos_network_individual_metrics<-read.csv("data/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
  select(elevation_extrapolated_date, species_jetz, Powder.lvl,foraging_cat, sociality,ectoparasites_PA, degree, w_degree, degree_species, degree_w_species, year_seasonality) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

ectos_prevalence_neworks<-dff_ectos_network_individual_metrics %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation=mean(elevation), degree_species_level=max(degree_species), degree_w_species_level=max(degree_w_species)) %>% 
  mutate(proportion_ectoparasites=ectoparasites_presence/sample_size) %>% 
  na.omit() 

ectos_birds_dff<-ectos_prevalence_neworks %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE)%>% 
  filter(sample_size>4)  %>% 
  na.omit()    # Eliminates species taht ocurr at two elevatiosn and keep only one ( in alphabetical order of the stations, e.g for galbula low_elevations is kept and montane is deleted, but the samples size is mantained this is jut to be able to do the phylogenetic analyses properly)

#ectos_birds_dff <- get.complete.cases(ectos_birds_dff) # mke sure we get all complete cases 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # Need to prunne the tree 


# Data structure
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$degree_w_species_level<-as.numeric(ectos_birds_dff$degree_w_species_level)
ectos_birds_dff$degree_species_level<-as.numeric(ectos_birds_dff$degree_species_level)
ectos_birds_dff$proportion_ectoparasites<-as.numeric(ectos_birds_dff$proportion_ectoparasites)
ectos_birds_dff$species<- ectos_birds_dff$species_jetz # create a column for the species effect different to the phylogenetic one

str(ectos_birds_dff)
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


#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


ecto_prevalence_degree_brms_bayes_no_int_species_priors<-brms::brm(proportion_ectoparasites~scale(degree_species_level)+ scale(elevation)+scale(sample_size)+
                                                     (1|gr(species_jetz, cov = phy_cov))+
                                                     (1|species),
                                                   data=ectos_birds_dff,
                                                   #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                                   family= zero_one_inflated_beta(),                                           
                                                   prior = c(prior_predictors,prior_random,prior_intercept),
                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(ecto_prevalence_degree_brms_bayes_no_int_species_priors, "data/data_analyses/model_selection/M2PND_species.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")
ecto_prevalence_degree_brms_bayes_no_int_species_priors<-readRDS("data/data_analyses/model_selection/M2PND_species.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")


bayes_R2(ecto_prevalence_degree_brms_bayes_no_int_species_priors)
mcmc_plot(ecto_prevalence_degree_brms_bayes_no_int_species_priors)
pp_check(ecto_prevalence_degree_brms_bayes_no_int_species_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)

color_scheme_set("purple")




# ##### 5.1.Data processing NETWORKS abundance lice ----------------------------------------------------

dff_ectos_network_individual_metrics<-read.csv("data/data_manuscript/3_dff_all_ectos_network_metrics_individuals_FILE_TIDY.csv",na.strings =c("","NA"))%>% 
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality,total_lice, degree, w_degree, year_seasonality,mass_tidy_species, mass_ind_tidy, mass_ind_comp) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

names(dff_ectos_network_individual_metrics)
# Lice

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
dff_ectos_network_individual_metrics$species<-dff_ectos_network_individual_metrics$species_jetz # create a column for the species effect different to the phylogenetic one

dff_ectos_network_individual_metrics$species<-as.factor(dff_ectos_network_individual_metrics$species)
dff_ectos_network_individual_metrics$mass_tidy_species<-as.numeric(dff_ectos_network_individual_metrics$mass_tidy_species)
dff_ectos_network_individual_metrics$mass_ind_tidy<-as.numeric(dff_ectos_network_individual_metrics$mass_ind_tidy)
dff_ectos_network_individual_metrics$mass_ind_comp<-as.numeric(dff_ectos_network_individual_metrics$mass_ind_comp)

names(dff_ectos_network_individual_metrics)
is.ultrametric(phylo)


# ##### 5.2.Model selection NETWORKS abundance lice --------------------------------------------------------
bayes_R2(zinb_a_lice_brms_bayes_no_int_degree)
prior_summary(zinb_a_lice_brms_bayes_no_int_degree)

zinb_a_lice_brms_bayes_no_int_degree<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                          (1|Powder.lvl)+
                                          (1|species),
                                        data=dff_ectos_network_individual_metrics,
                                        family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                        data2 = list(phy_cov=phy_cov),
                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                        thin=2,
                                        save_pars = save_pars(all=  TRUE),
                                        control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zinb_a_lice_brms_bayes_no_int_degree, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")
zinb_a_lice_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")

prior_summary(zinb_a_lice_brms_bayes_sociality_interactions_degree)
zinb_a_lice_brms_bayes_sociality_interactions_degree<-brms::brm(total_lice~
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
                                                        save_pars = save_pars(all=  TRUE),
                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 
prior_summary(model)

saveRDS(zip_a_lice_brms_bayes_sociality_interactions_degree, "data/data_analyses/model_selection/M2LND.model_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")
zip_a_lice_brms_bayes_sociality_interactions_degree<-readRDS( "data/data_analyses/model_selection/M2LND.model_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions_DEGREE.RDS")

zinb_a_lice_brms_bayes_all_interactions_degree<-brms::brm(total_lice~
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
                                                  save_pars = save_pars(all=  TRUE),
                                                  thin=2,
                                                  control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zinb_a_lice_brms_bayes_all_interactions_degree, "data/data_analyses/model_selection/M2LND.model_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_DEGREE.RDS")
zinb_a_lice_brms_bayes_all_interactions_degree<-readRDS( "data/data_analyses/model_selection/M2LND.model_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_DEGREE.RDS")

#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

loo(zinb_a_lice_brms_bayes_no_int_degree_prior,zinb_a_lice_brms_bayes_sociality_interactions_degree_prior, compare=TRUE)

mcmc_plot(zinb_a_lice_brms_bayes_no_int_degree_prior)
mcmc_plot(zinb_a_lice_brms_bayes_sociality_interactions_degree_prior)
prior_summary(zinb_a_lice_brms_bayes_no_int_degree_prior1)

###_###_###_###_###_###_###_###_###_###_
#THE BESTY  LICE DEGREE###
###_###_###_###_###_###_###_###_###_###_
zinb_a_lice_brms_bayes_no_int_degree_prior_nopowder<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_no_int_degree_prior_nopowder, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_nopowder.RDS")
zinb_a_lice_brms_bayes_no_int_degree_prior_nopowder<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_nopowder.RDS")


zinb_a_lice_brms_bayes_no_int_degree_prior<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                  (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                  (1|Powder.lvl)+
                                                  (1|species),
                                                data=dff_ectos_network_individual_metrics,
                                                family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                data2 = list(phy_cov=phy_cov),
                                                iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                save_pars = save_pars(all=  TRUE),
                                                prior = c(prior_predictors,prior_random),
                                                control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")
zinb_a_lice_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")

# comparing with a poisson distribution
zip_a_lice_brms_bayes_no_int_degree_prior_p<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_poisson(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(zip_a_lice_brms_bayes_no_int_degree_prior_p, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_p_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_poisson.RDS")
zip_a_lice_brms_bayes_no_int_degree_prior_p<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_p_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_poisson.RDS")


prior_summary(zinb_a_lice_brms_bayes_sociality_interactions_degree_prior)
zinb_a_lice_brms_bayes_sociality_interactions_degree_prior<-brms::brm(total_lice~
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
                                                                save_pars = save_pars(all=  TRUE),
                                                                prior = c(prior_predictors,prior_random),
                                                                control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_sociality_interactions_degree_prior, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE_prior.RDS")
zinb_a_lice_brms_bayes_sociality_interactions_degree_prior<-readRDS( "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_sociality_interactions_DEGREE_prior.RDS")

zinb_a_lice_brms_bayes_all_interactions_degree_prior<-brms::brm(total_lice~
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
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random, prior_intercept),
                                                          thin=2,
                                                          control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zinb_a_lice_brms_bayes_all_interactions_degree_prior, "data/data_analyses/model_selection/M2LND.model_LICE ABUNDANCE_b_brms_phylo_multiple_obs_all_interactions_degree_prior.RDS")
zinb_a_lice_brms_bayes_all_interactions_degree_prior<-readRDS( "data/data_analyses/model_selection/M2LND.model_LICE_ABUNDANCE_b_brms_phylo_multiple_obs_all_interactions_degree_prior.RDS")

# MODEL COMPARISON 

loo(zinb_a_lice_brms_bayes_no_int_degree_prior, zinb_a_lice_brms_bayes_sociality_interactions_degree_prior,compare=TRUE)
loo(zinb_a_lice_brms_bayes_no_int_degree, zinb_a_lice_brms_bayes_sociality_interactions_degree_prior,compare=TRUE)
loo(zinb_a_lice_brms_bayes_no_int_degree_prior, zinb_a_lice_brms_bayes_sociality_interactions_degree_prior, zip_a_nf_mites_brms_bayes_all_interactions_degree_prior,compare=TRUE)
loo(zinb_a_lice_brms_bayes_no_int_degree, zinb_a_lice_brms_bayes_sociality_interactions_degree,compare=TRUE)
# R=eld_diff<4 so we keep the simplest model !
loo_compare(waic(zinb_a_lice_brms_bayes_no_int_priors), waic(zinb_a_lice_brms_bayes_sociality_interactions_priors)) # interesting warning
mcmc_plot(zinb_a_lice_brms_bayes_no_int)
mcmc_plot(zinb_a_lice_brms_bayes_no_int_priors)

# Model posterior predictive checks 
#The idea of posterior predictive checks is to compare our observed data to replicated data from the model. 
#If our model is a good fit, we should be able to use it to generate a dataset that resembles the observed data.
#gettigsamples from the posterior predictive distribution:
pp_check(zinb_a_nf_mites_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100) 

#or dependeing or the data use
#pp_check(ecto_p_brms_bayes_no_int_prior, type = "stat", stat = 'median', nsamples = 100)
#pp_m<- brms::posterior_predict(ecto_p_brms_bayes_no_int_prior)
#ppc_rootogram(y=ecto_p_brms_bayes_no_int_prior$data$ectoparasites_PA, pp_m[1:200, ])  +   
coord_cartesian(xlim = c(0, 100), ylim = c(0,30))

###_###_###_##
#PLOTS
###_###_###_##

conditional_effects(zinb_a_lice_brms_bayes_no_int_degree_prior)
marginal_effects()
plot(
  conditional_effects(zinb_a_lice_brms_bayes_no_int_degree_prior_p, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1)
)


summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("purple")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig2LND.zinb_ABUNDANCE_LICE_brms_bayes_int_sociality_priors_CONVERGENCE_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_lice_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig2LND.zinb_ABUNDANCE_LICE_brms_bayes_int_sociality_priors_FIT_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_lice_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES

estimates_plot<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE DEGREE with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2LND.ZINB_ABUNDANCE_LICE_brms_bayes_ESTIMATES_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="ZINB Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE NETWORKS with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2LND.ZINB_LICE_ABUNDANCE_brms_bayes_INTERVALS_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

#### strenght

#PRIORS
#prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

prior_summary(zinb_a_lice_brms_bayes_no_int_W_degree_prior)

zinb_a_lice_brms_bayes_no_int_W_degree_prior<-brms::brm(total_lice~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_no_int_W_degree_prior, "data/data_analyses/model_selection/M2LNWD.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_W_DEGREE_prior.RDS")
zinb_a_lice_brms_bayes_no_int_W_degree_prior<-readRDS("data/data_analyses/model_selection/M2LNWD.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_W_DEGREE_prior.RDS")

loo(zinb_a_lice_brms_bayes_no_int_W_degree_prior, moment_match = TRUE)
mcmc_plot(zinb_a_lice_brms_bayes_no_int_W_degree_prior)
pp_check(zinb_a_lice_brms_bayes_no_int_W_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)

zinb_a_lice_brms_bayes_no_int_wdegree_prior<-brms::brm(total_lice~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                          (1|Powder.lvl)+
                                                          (1|species),
                                                        data=dff_ectos_network_individual_metrics,
                                                        family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                        data2 = list(phy_cov=phy_cov),
                                                        iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                        thin=2,
                                                        save_pars = save_pars(all=  TRUE),
                                                        prior = c(prior_predictors,prior_random,prior_intercept),
                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 

bayes_R2(zinb_a_lice_brms_bayes_no_int_wdegree_prior)
mcmc_plot(zinb_a_lice_brms_bayes_no_int_wdegree_prior)
pp_check(zinb_a_lice_brms_bayes_no_int_wdegree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)


# #### ### #### ## 5.3 ** Model abundance NETWORKS excluding zeros Lice -------------------------------

dff_ectos_network_individual_metrics<-read.csv("data/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
  select(elevation_extrapolated_date, species_jetz, Powder.lvl,foraging_cat, sociality,total_lice, degree, w_degree, year_seasonality) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  filter (total_lice!=0) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

names(dff_ectos_network_individual_metrics)
# Lice

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

# MODEL 

#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
prior_summary(NZ_a_lice_brms_bayes_no_degree_int_prior)
NZ_a_lice_brms_bayes_no_degree_int_prior_nb<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                      (1|Powder.lvl)+
                                                      (1|species),
                                                    data=dff_ectos_network_individual_metrics,
                                                    family=negbinomial(),  #zero_inflated_negbinomial()
                                                    data2 = list(phy_cov=phy_cov),
                                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                    thin=2,
                                                    save_pars = save_pars(all=  TRUE),
                                                    prior = c(prior_predictors,prior_random),
                                                    control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(NZ_a_lice_brms_bayes_no_degree_int_prior_nb, "data/data_analyses/model_selection/M2LND_NZ.model_lICE_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")
NZ_a_lice_brms_bayes_no_degree_int_prior_nb<-readRDS("data/data_analyses/model_selection/M2LND_NZ.model_lICE_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")

NZ_a_lice_brms_bayes_no_degree_int_prior_p<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=poisson(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14))


NZ_a_lice_brms_bayes_no_degree_int_prior_truncated<-brms::brm(total_lice|trunc(lb=1)~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                 (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                 (1|Powder.lvl)+
                                                                 (1|species),
                                                               data=dff_ectos_network_individual_metrics,
                                                               family=negbinomial(),  #zero_inflated_negbinomial()
                                                               data2 = list(phy_cov=phy_cov),
                                                               iter=500, warmup=200, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                               thin=2,
                                                               save_pars = save_pars(all=  TRUE),
                                                               control=list(adapt_delta=0.99, max_treedepth=14)) 

# Evaludating model fit and convergence

simulate_residuals <- dh_check_brms(NZ_a_lice_brms_bayes_no_degree_int_prior_nb, integer = TRUE)
plot(simulate_residuals, form = dff_ectos_network_individual_metrics$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations

# Model comparison 
loo(NZ_a_lice_brms_bayes_no_degree_int_prior_nb, moment_match=TRUE)
k_ecto_NZ_lice_brms_no_int_prior_nb<-kfold(NZ_a_lice_brms_bayes_no_degree_int_prior, K=10)
saveRDS(k_ecto_NZ_lice_brms_no_int_prior_nb, "data/data_analyses/model_selection/k_fold/K_fold_ALNZ_DEGREE_model_lice_abundance_brms_no_interactions_priors_nbinomial.RDS")

loo(NZ_a_lice_brms_bayes_no_degree_int_prior_p,moment_match=TRUE)
k_ecto_NZ_lice_brms_no_int_prior_p<-kfold(NZ_a_lice_brms_bayes_no_degree_int_prior_p, K=10)
saveRDS(k_ecto_NZ_lice_brms_no_int_prior_nb, "data/data_analyses/model_selection/k_fold/K_fold_ALNZ_DEGREE_model_lice_abundance_brms_no_interactions_priors_poisson.RDS")

loo(NZ_a_lice_brms_bayes_no_degree_int_prior_nb,NZ_a_lice_brms_bayes_no_degree_int_prior_p, compare=TRUE)

# The negaive binomial performs much better than the poission ! model selected
loo_compare(k_ecto_NZ_lice_brms_no_int_prior_nb,k_ecto_NZ_lice_brms_no_int_prior_p)

# plots 
color_scheme_set("red") 

# poisson 
#model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig2_LNDNZ_NB_plot_model_CONVERGENCE_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_nb.png",width = 3000, height = 3000, res = 300, units = "px")
plot(NZ_a_lice_brms_bayes_no_degree_int_prior_nb)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig2_LNDNZ_NB_plot_modell_FIT__ecto_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_nb.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(NZ_a_lice_brms_bayes_no_degree_int_prior_nb, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(NZ_a_lice_brms_bayes_no_degree_int_prior_nb,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions NZ_LICE ABUNDANCE_DEGREE", subtitle ="NZ_LICE_ABUNDANCE_DEGREE with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2_LNDNZ_NB_plot_model_parameters_LICE_ABUNDANCE_brms_bayes_no_int_degree.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(NZ_a_lice_brms_bayes_no_degree_int_prior_nb,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions NZ_LICE ABUNDANCE", subtitle ="NZ_LICE_ABUNDANCE with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2_LNDNZ_NB_plot_model_parameters_intervals_LICE_ABUNDANCE_brms_bayes_no_int_degree.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# ##### 5.5. Selecting model NETWORKS LICE ABUNDANCE DEGREE (ZIP VS ZINB) ------------------------------------------

#prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

prior_summary(zinb_a_lice_brms_bayes_no_int_degree_prior)

###_###_###_###_###_###_###_###_###_###_
#THE BESTY  LICE DEGREE###
###_###_###_###_###_###_###_###_###_###_
zinb_a_lice_brms_bayes_no_int_degree_prior<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")
zinb_a_lice_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")
###_###_###_###_###_###_###_###_###_###_


# comparing with a poisson distribution
zip_a_lice_brms_bayes_no_int_degree_prior_p<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                         (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                         (1|Powder.lvl)+
                                                         (1|species),
                                                       data=dff_ectos_network_individual_metrics,
                                                       family=zero_inflated_poisson(),  #zero_inflated_negbinomial()
                                                       data2 = list(phy_cov=phy_cov),
                                                       iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                       thin=2,
                                                       save_pars = save_pars(all=  TRUE),
                                                       prior = c(prior_predictors,prior_random),
                                                       control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(zip_a_lice_brms_bayes_no_int_degree_prior_p, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_p_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_poisson.RDS")
zip_a_lice_brms_bayes_no_int_degree_prior_p<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_p_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_poisson.RDS")

###_###_### EXCLUDING POWDER 
zinb_a_lice_brms_bayes_no_int_degree_prior_no_powder<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_no_int_degree_prior_no_powder, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_nopower.RDS")
zinb_a_lice_brms_bayes_no_int_degree_prior_no_powder<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_nopower.RDS")

# including mass 
zinb_a_lice_brms_bayes_no_int_degree_prior_mass<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+scale(mass_tidy)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_no_int_degree_prior_mass, "data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_mass.RDS")
###_###_###_###_###_###_###_###_###_###_
mcmc_plot(zinb_a_lice_brms_bayes_no_int_degree_prior_no_powder)

bayes_R2(zinb_a_lice_brms_bayes_no_int_degree_prior_mass)
bayes_R2(zinb_a_lice_brms_bayes_no_int_degree_prior)


loo(zinb_a_lice_brms_bayes_no_int_degree_prior_mass,zinb_a_lice_brms_bayes_no_int_degree_prior, compare=TRUE)

###_###_###_###_###_###_###_###_###_###_
###_###_###_##
#MODEL COMPARISON
###_###_###_##

#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.
bayes_R2(zinb_a_lice_brms_bayes_no_int_degree_prior)
bayes_R2(zip_a_lice_brms_bayes_no_int_degree_prior_p)

# use loo cross validation [ In this case we are not interested inteh interaaction between the different factors because the sociality of a speceis does not change with seasonality or elevation, cause those are attributes to the individuals does not change ]
loo(zinb_a_lice_brms_bayes_no_int_degree_prior_p,zinb_a_lice_brms_bayes_no_int_degree_prior, compare=TRUE)
loo_p<-loo(zinb_a_lice_brms_bayes_no_int_degree_prior_p, moment_match=TRUE)
loo_nb<-loo(zinb_a_lice_brms_bayes_no_int_degree_prior, moment_match=TRUE)

loo_compare(loo_pm,loo_nb)

# use k-fold-cv validation instead because in models with random efFects loo tends to fail # but this is taking forevwer so i WILL RUNT IT AT UBC TOMORROW
k_lice_brms_bayes_no_int_degree_prior<-kfold(zinb_a_lice_brms_bayes_no_int_degree_prior, K=10)
saveRDS(k_lice_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2LND.model_LICE_ABUNDANCE_DEGREE_no_interactions_priors_nb.RDS")
k_lice_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2LND.model_LICE_ABUNDANCE_DEGREE_no_interactions_priors_nb.RDS")

k_lice_brms_bayes_no_int_degree_prior_p<-kfold(zip_a_lice_brms_bayes_no_int_degree_prior_p, K=10)
saveRDS(k_lice_brms_bayes_no_int_degree_prior_p,"data/data_analyses/model_selection/k_fold/K_fold_M2LND.model_LICE_ABUNDANCE_DEGREE_no_interactions_priors_poisson.RDS")
k_lice_brms_bayes_no_int_degree_prior_p<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2LND.model_LICE_ABUNDANCE_DEGREE_no_interactions_priors_poisson.RDS")

loo_compare(k_lice_brms_bayes_no_int_degree_prior, k_lice_brms_bayes_no_int_degree_prior_p)

###The best model is the zero inflated negatve binomialr 

# plots 
color_scheme_set("teal") 

# poisson 
#model convergence 
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LNDNZ_BEST_plot_model_CONVERGENCE_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_BEST_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_lice_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LNDNZ_BEST_plot_modell_FIT_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_BEST_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_lice_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE DEGREE ")+
  xlim(-5,5)+
   theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LND_BEST_plot_model_parameters_LICE ABUNDANCE_brms_bayes_social_int_DEGREE_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE DEGREE ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LND_BEST_plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_social_int_DEGREE_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



# ##### 5.6. Selected model NETWORKS LICE ABUNDANCE STRENGHT (ZIP VS ZINB) ------------------------------------------

prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

prior_summary(zinb_a_lice_brms_bayes_no_int_degree_prior)

###_###_###_###_###_###_###_###_###_###_
#THE BESTY  LICE DEGREE###
###_###_###_###_###_###_###_###_###_###_
zinb_a_lice_brms_bayes_no_int_strength_prior<-brms::brm(total_lice~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(zinb_a_lice_brms_bayes_no_int_strength_prior, "data/data_analyses/model_selection/M2LNS.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior.RDS")
zinb_a_lice_brms_bayes_no_int_strength_prior<-readRDS("data/data_analyses/model_selection/M2LNS.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior.RDS")
###_###_###_###_###_###_###_###_###_###_


# comparing with a poisson distribution
zip_a_lice_brms_bayes_no_int_strength_prior<-brms::brm(total_lice~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                         (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                         (1|Powder.lvl)+
                                                         (1|species),
                                                       data=dff_ectos_network_individual_metrics,
                                                       family=zero_inflated_poisson(),  #zero_inflated_negbinomial()
                                                       data2 = list(phy_cov=phy_cov),
                                                       iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                       thin=2,
                                                       save_pars = save_pars(all=  TRUE),
                                                       prior = c(prior_predictors,prior_random),
                                                       control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(zip_a_lice_brms_bayes_no_int_strength_prior, "data/data_analyses/model_selection/M2LNS.model_lICE_ABUNDANCE_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior_poisson.RDS")
zip_a_lice_brms_bayes_no_int_strength_prior<-readRDS("data/data_analyses/model_selection/M2LNS.model_lICE_ABUNDANCE_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior_poisson.RDS")

###_###_###_###_###_###_###_###_###_###_
###_###_###_##
#MODEL COMPARISON
###_###_###_##

#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.
bayes_R2(zinb_a_lice_brms_bayes_no_int_degree_prior)
bayes_R2(zip_a_lice_brms_bayes_no_int_degree_prior_p)


loo_compare(loo_pm,loo_nb)

# use k-fold-cv validation instead because in models with random efFects loo tends to fail # but this is taking forevwer so i WILL RUNT IT AT UBC TOMORROW
k_lice_brms_bayes_no_int_strength_prior<-kfold(zinb_a_lice_brms_bayes_no_int_strength_prior, K=10)
saveRDS(k_lice_brms_bayes_no_int_strength_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2LNS.model_LICE_ABUNDANCE_STRENGTH_no_interactions_priors_nb.RDS")
k_lice_brms_bayes_no_int_strength_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2LNs.model_LICE_ABUNDANCE_STRENGTH_no_interactions_priors_nb.RDS")

k_lice_brms_bayes_no_int_strength_prior<-kfold(zip_a_lice_brms_bayes_no_int_strength_prior, K=10)
saveRDS(k_lice_brms_bayes_no_int_strength_prior,"data/data_analyses/model_selection/k_fold/K_fold_M2LNS.model_LICE_ABUNDANCE_STRENGTH_no_interactions_priors_poisson.RDS")
k_lice_brms_bayes_no_int_strength_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2LNS.model_LICE_ABUNDANCE_STRENGTH_no_interactions_priors_poisson.RDS")

loo_compare(k_lice_brms_bayes_no_int_degree_prior, k_lice_brms_bayes_no_int_degree_prior_p)

###The best model is the zero inflated negatve binomialr 

# plots 
color_scheme_set("teal") 

# poisson 
#model convergence 
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LNS_BEST_plot_model_CONVERGENCE_LICE_ABUNDANCE_brms_bayes_no_int_STRENGT_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_lice_brms_bayes_no_int_strength_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LNS_BEST_plot_modell_FIT_LICE_ABUNDANCE_brms_bayes_no_int_STRENGTH_BEST_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_lice_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE STRENGTH ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LS_BEST_plot_model_parameters_LICE ABUNDANCE_brms_bayes_social_int_STRENGTH_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE STRENGTH ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LNS_BEST_plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_social_int_STRENGTH_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# NETWORK ABUNDANCE  5.7 **** Selected**** model abundance lice NETWORKS degree and strength --------

###_###_###_###_###_###_###_###_###_###_
#LICE ABUNDANCE DEGREE###
###_###_###_###_###_###_###_###_###_###_

#prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


selected_zinb_a_lice_brms_bayes_no_int_degree_prior<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                 scale(mass_tidy_species)+(mass_ind_comp)+
                                                        (1|gr(species_jetz, cov = phy_cov))+(1|species)+
                                                        (1|Powder.lvl),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(selected_zinb_a_lice_brms_bayes_no_int_degree_prior, "results/selected_models/4_3_M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")
selected_zinb_a_lice_brms_bayes_no_int_degree_prior<-readRDS("results/selected_model/4_3_M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")

selected_zinb_a_lice_brms_bayes_no_int_degree_prior2<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                 scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                                 (1|gr(species_jetz, cov = phy_cov))+(1|species)+
                                                                 (1|Powder.lvl),
                                                               data=dff_ectos_network_individual_metrics,
                                                               family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                               data2 = list(phy_cov=phy_cov),
                                                               iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                               thin=2,
                                                               save_pars = save_pars(all=  TRUE),
                                                               prior = c(prior_predictors,prior_random),
                                                               control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(selected_zinb_a_lice_brms_bayes_no_int_degree_prior2, "results/selected_models/4_3_M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_ind_mass_scaled.RDS")

###The best model is the zero inflated negatve binomialr 

# plots 
color_scheme_set("teal") 
# poisson 
#model convergence 
png("results/selected_models_figures/4_3_Fig4LNDNZ_BEST_plot_model_CONVERGENCE_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_BEST_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_zinb_a_lice_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("results/selected_models_figures/4_3_Fig4LNDNZ_BEST_plot_modell_FIT_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_BEST_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_zinb_a_lice_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE DEGREE ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_Fig4LND_BEST_plot_model_parameters_LICE ABUNDANCE_brms_bayes_social_int_DEGREE_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE DEGREE ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_Fig4LND_BEST_plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_social_int_DEGREE_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


###_###_###_###_###_###_###_###_###_###_
#LICE ABUNDANCE Strenght ###
###_###_###_###_###_###_###_###_###_###_

prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

prior_summary(zinb_a_lice_brms_bayes_no_int_degree_prior)

selected_zinb_a_lice_brms_bayes_no_int_strength_prior<-brms::brm(total_lice~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                   scale(mass_tidy_species)+(mass_ind_comp)+
                                                                   (1|gr(species_jetz, cov = phy_cov))+  
                                                                   (1|Powder.lvl)+
                                                                   (1|species),
                                                        data=dff_ectos_network_individual_metrics,
                                                        family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                        data2 = list(phy_cov=phy_cov),
                                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                        thin=2,
                                                        save_pars = save_pars(all=  TRUE),
                                                        prior = c(prior_predictors,prior_random),
                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(selected_zinb_a_lice_brms_bayes_no_int_strength_prior, "results/selected_models/4_3_2_M2LNS.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior.RDS")
selected_zinb_a_lice_brms_bayes_no_int_strength_prior<-readRDS("results/selected_model/4_3_2_M2LNS.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior.RDS")

selected_zinb_a_lice_brms_bayes_no_int_strength_prior2<-brms::brm(total_lice~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                    scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                                   (1|gr(species_jetz, cov = phy_cov))+ 
                                                                   (1|Powder.lvl)+
                                                                   (1|species),
                                                                 data=dff_ectos_network_individual_metrics,
                                                                 family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                                 data2 = list(phy_cov=phy_cov),
                                                                 iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                 thin=2,
                                                                 save_pars = save_pars(all=  TRUE),
                                                                 prior = c(prior_predictors,prior_random),
                                                                 control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(selected_zinb_a_lice_brms_bayes_no_int_strength_prior2, "results/selected_models/4_3_2_M2LNS.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior_ind_mass_scaled.RDS")
###_###_###_###_###_###_###_###_###_###_


###The best model is the zero inflated negatve binomialr 

# plots 
color_scheme_set("teal") 

# poisson 
#model convergence 
png("results/selected_models_figures/4_3_2Fig4LNDNZ_BEST_plot_model_CONVERGENCE_LICE_ABUNDANCE_brms_bayes_no_int_STRENGTH_BEST_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_zinb_a_lice_brms_bayes_no_int_strength_prior)
dev.off()

# model fit
png("results/selected_models_figures/4_3_2Fig4LNDNZ_BEST_plot_modell_FIT_LICE_ABUNDANCE_brms_bayes_no_int_STRENGTH_BEST_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_zinb_a_lice_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE STRENGTH ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_2Fig4LND_BEST_plot_model_parameters_LICE ABUNDANCE_brms_bayes_social_int_STRENGTH_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE STRENGTH ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_2Fig4LND_BEST_plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_social_int_STRENGTH_ZINB.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# ##### 6.1.Data processing NETWORKS abundance mites ----------------------------------------------------

# Mites
dff_ectos_network_individual_metrics<-read.csv("data/data_manuscript/3_dff_all_ectos_network_metrics_individuals_FILE_TIDY.csv",na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality,total_mites, total_no_feathers_mites, degree, w_degree, year_seasonality,mass_tidy_species, mass_ind_tidy, mass_ind_comp) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")%>% filter(total_no_feathers_mites<60)

View(dff_ectos_network_individual_metrics)
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
dff_ectos_network_individual_metrics$species<-as.factor(dff_ectos_network_individual_metrics$species)
dff_ectos_network_individual_metrics$mass_tidy_species<-as.numeric(dff_ectos_network_individual_metrics$mass_tidy_species)
dff_ectos_network_individual_metrics$mass_ind_tidy<-as.numeric(dff_ectos_network_individual_metrics$mass_ind_tidy)
dff_ectos_network_individual_metrics$mass_ind_comp<-as.numeric(dff_ectos_network_individual_metrics$mass_ind_comp)

names(dff_ectos_network_individual_metrics)
is.ultrametric(phylo)

# ##### 6.2.Model selection NETWORKS abundance mites --------------------------------------------------------
prior_summary(zip_a_nf_mites_brms_bayes_no_int_degree)
zinb_a_nf_mites_brms_bayes_no_int_degree<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                              (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                              (1|Powder.lvl)+
                                              (1|species),
                                            data=dff_ectos_network_individual_metrics,
                                            family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                            data2 = list(phy_cov=phy_cov),
                                            iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                            save_pars = save_pars(all=  TRUE),
                                            thin=2,
                                            control=list(adapt_delta=0.99, max_treedepth=14)) 


mcmc_plot(zip2_a_nf_mites_brms_bayes_no_int_degree2)
bayes_R2(zip_a_nf_mites_brms_bayes_no_int_degree)
bayes_R2(zip2_a_nf_mites_brms_bayes_no_int_degree2)

loo(zip_a_nf_mites_brms_bayes_no_int_degree,zip2_a_nf_mites_brms_bayes_no_int_degree2)

bayes_R2()

loo(zip_a_nf_mites_brms_bayes_no_int_degree)
#saveRDS(zip_a_nf_mites_brms_bayes_no_int_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_DEGREE.RDS")
zip_a_nf_mites_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_DEGREE.RDS")

zinb_a_nf_mites_brms_bayes_sociality_interactions_degree<-brms::brm(total_no_feathers_mites~
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


loo(zip_a_lice_brms_bayes_all_interactions,zip_a_lice_brms_bayes_sociality_interactions,zip_a_lice_brms_bayes_no_int)
mcmc_plot(zip_a_nf_mites_brms_bayes_sociality_interactions)


#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                     (1|Powder.lvl)+
                                                     (1|species),
                                                   data=dff_ectos_network_individual_metrics,
                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   save_pars = save_pars(all=  TRUE),
                                                   prior = c(prior_predictors,prior_random),
                                                   thin=2,
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

#saveRDS(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_ZINB_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")
zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-readRDS( "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_ZINB_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")

ZIP_a_nf_mites_brms_bayes_no_int_degree_prior<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                            (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                            (1|Powder.lvl)+
                                                            (1|species),
                                                          data=dff_ectos_network_individual_metrics,
                                                          family=zero_inflated_poisson(), 
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random),
                                                          thin=2,
                                                          control=list(adapt_delta=0.99, max_treedepth=14)) 

#saveRDS(ZIP_a_nf_mites_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_ZIP_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")
ZIP_a_nf_mites_brms_bayes_no_int_degree_prior<-readRDS( "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_ZIP_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")


### MODEL COMPARISON 

loo(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, moment_match=TRUE)
loo(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, ZIP_a_nf_mites_brms_bayes_no_int_degree_prior, compare=TRUE)

k_zinb_nf_mites_brms_no_int_degree_prior<-kfold(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, K=10)
saveRDS(k_zinb_nf_mites_brms_no_int_degree_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_nb.RDS")
k_zinb_nf_mites_brms_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_nb.RDS")

k_zip_nf_mites_brms_no_int_degree_prior<-kfold(zip_a_nf_mites_brms_bayes_no_int_degree_prior, K=10)
saveRDS(k_zip_nf_mites_brms_no_int_degree_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_poisson.RDS")
k_zip_nf_mites_brms_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_poisson.RDS")


###_###_###_##
#PLOTS
###_###_###_##

conditional_effects(zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
marginal_effects()
plot(
  conditional_effects(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1)
)


summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("purple")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig2MND.zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_CONVERGENCE_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig2MND.zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_FIT_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES
estimates_plot<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions MITES ABUNDANCE DEGREE", subtitle ="MITES ABUNDANCE DEGREE with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2MND.ZINB_ABUNDANCE_MITES_brms_bayes_ESTIMATES_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="ZINB Posterior distributions MITES ABUNDANCE DEGREE", subtitle ="MITES ABUNDANCE NETWORKS with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2MND.ZINB_MITES_ABUNDANCE_brms_bayes_INTERVALS_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()




###_###_###_###_###_###_###_###_###_###_###_###_

# ##### 6.3 Model selection NETWORKS ABUNDANCE MITES DEGREE (ZIP&ZINB) --------------------------------

#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                            (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                            (1|Powder.lvl)+
                                                            (1|species),
                                                          data=dff_ectos_network_individual_metrics,
                                                          family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random),
                                                          thin=2,
                                                          control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")
zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-readRDS( "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")

ZIP_a_nf_mites_brms_bayes_no_int_degree_prior<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                           (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                           (1|Powder.lvl)+
                                                           (1|species),
                                                         data=dff_ectos_network_individual_metrics,
                                                         family=zero_inflated_poisson(), 
                                                         data2 = list(phy_cov=phy_cov),
                                                         iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                         save_pars = save_pars(all=  TRUE),
                                                         prior = c(prior_predictors,prior_random),
                                                         thin=2,
                                                         control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(ZIP_a_nf_mites_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_zip_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")
ZIP_a_nf_mites_brms_bayes_no_int_degree_prior<-readRDS( "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_zip_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")


### MODEL COMPARISON 

loo(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, moment_match=TRUE)
loo(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, ZIP_a_nf_mites_brms_bayes_no_int_degree_prior, compare=TRUE)

k_zinb_nf_mites_brms_no_int_degree_prior<-kfold(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, K=10)
saveRDS(k_zinb_nf_mites_brms_no_int_degree_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_nb.RDS")
k_zinb_nf_mites_brms_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_nb.RDS")

k_zip_nf_mites_brms_no_int_degree_prior<-kfold(ZIP_a_nf_mites_brms_bayes_no_int_degree_prior, K=10)
saveRDS(k_zip_nf_mites_brms_no_int_degree_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_poisson.RDS")
k_zip_nf_mites_brms_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2MND.model_MITES_ABUNDANCE_DEGREE_no_interactions_priors_poisson.RDS")

loo_compare(k_zinb_nf_mites_brms_no_int_degree_prior,k_zip_nf_mites_brms_no_int_degree_prior)
#BESTY # Best model for mites degree 


###_###_###_##
#PLOTS
###_###_###_##

conditional_effects(zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
marginal_effects()
plot(
  conditional_effects(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1)
)


summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("green")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2MND_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_CONVERGENCE_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2MND_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_FIT_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES

estimates_plot<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE DEGREE")+
  xlim(-8,8)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures//best_models/Fig2MND_BEST_zinb_ABUNDANCE_MITES_brms_bayes_ESTIMATES_int_DEGREE.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE DEGREE")+
  xlim(-8,8)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures//best_models/Fig2MND_BEST_zinb_MITES_ABUNDANCE_brms_bayes_INTERVALS_no_int_DEGREE.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# ##### 6.4 Model selection NETWORKS ABUNDANCE MITES STRENGHT (ZIP&ZINB) --------------------------------


#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions 
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default so no need to speciied
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default so no need to speciied


prior_summary(zinb_a_nf_mites_brms_bayes_no_int_strength_prior)

mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_strength_prior)
zinb_a_nf_mites_brms_bayes_no_int_strength_prior<-brms::brm(total_no_feathers_mites~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                            (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                            (1|Powder.lvl)+
                                                            (1|species),
                                                          data=dff_ectos_network_individual_metrics,
                                                          family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random),
                                                          thin=2,
                                                          control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zinb_a_nf_mites_brms_bayes_no_int_strength_prior, "data/data_analyses/model_selection/M2MNS_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_STRENGHT_prior.RDS")
zinb_a_nf_mites_brms_bayes_no_int_strength_prior<-readRDS( "data/data_analyses/model_selection/M2MNS_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_STRENGHT_prior.RDS")

ZIP_a_nf_mites_brms_bayes_no_int_strength_prior<-brms::brm(total_no_feathers_mites~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                           (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                           (1|Powder.lvl)+
                                                           (1|species),
                                                         data=dff_ectos_network_individual_metrics,
                                                         family=zero_inflated_poisson(), 
                                                         data2 = list(phy_cov=phy_cov),
                                                         iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                         save_pars = save_pars(all=  TRUE),
                                                         prior = c(prior_predictors,prior_random),
                                                         thin=2,
                                                         control=list(adapt_delta=0.99, max_treedepth=14)) 


#saveRDS(ZIP_a_nf_mites_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2MNS_model_ZIP_MITES_ABUNDANCE_STRENGTH_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior.RDS")
ZIP_a_nf_mites_brms_bayes_no_int_strength_prior<-readRDS( "data/data_analyses/model_selection/M2MNS_model_ZIP_MITES_ABUNDANCE_STRENGTH_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior.RDS")


### MODEL COMPARISON 

loo(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, moment_match=TRUE)
loo(ZIP_a_nf_mites_brms_bayes_no_int_strength_prior, zinb_a_nf_mites_brms_bayes_no_int_strength_prior, compare=TRUE)

k_zinb_nf_mites_strength_prior<-kfold(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, K=10)
saveRDS(k_zinb_nf_mites_strength_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2MNS.model_ZINB_MITES_ABUNDANCE_STRENGTH_no_interactions_priors.RDS")
k_zinb_nf_mites_strength_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2MNS.model_ZINB_MITES_ABUNDANCE_STRENGTH_no_interactions_priors.RDS")

k_zip_nf_mites_strength_prior<-kfold(ZIP_a_nf_mites_brms_bayes_no_int_strength_prior, K=10)
saveRDS(k_zip_nf_mites_strength_prior, "data/data_analyses/model_selection/k_fold/K_fold_M2MNS.model_ZIP_MITES_ABUNDANCE_STRENGTH_no_interactions_priors.RDS")
k_zip_nf_mites_strength_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M2MNS.model_ZIP_MITES_ABUNDANCE_STRENGTH_no_interactions_priors.RDS")

loo_compare(k_zinb_nf_mites_strength_prior,k_zip_nf_mites_strength_prior)
#BESTY # Best model for mites degree 

###_###_###_##
#PLOTS
###_###_###_##

conditional_effects(zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
marginal_effects()
plot(
  conditional_effects(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1)
)


summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("green")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_CONVERGENCE_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_a_nf_mites_brms_bayes_no_int_strength_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/best_models/Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_FIT_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(zinb_a_nf_mites_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES

estimates_plot<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~STRENGTH")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_ESTIMATES_int_STRENGTH.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~STRENGTH")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2MNS_BEST_zinb_MITES_ABUNDANCE_brms_bayes_INTERVALS_no_int_STRENGTH.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


# #### ### #### ## 6.5 ** Model networks abundance excluding zeros Mites -------------------------------

dff_ectos_network_individual_metrics<-read.csv("data/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
  select(elevation_extrapolated_date, species_jetz, Powder.lvl,foraging_cat, sociality,total_no_feathers_mites, degree, w_degree, year_seasonality) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  filter (total_no_feathers_mites!=0) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

names(dff_ectos_network_individual_metrics)
# Lice

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
dff_ectos_network_individual_metrics$total_no_feathers_mites<-as.numeric(dff_ectos_network_individual_metrics$total_no_feathers_mites)
dff_ectos_network_individual_metrics$species<-dff_ectos_network_individual_metrics$species_jetz # create a column for the species effect different to the phylogenetic one
names(dff_ectos_network_individual_metrics)
is.ultrametric(phylo)

# MODEL 

#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

NZ_a_mites_brms_bayes_no_degree_int_prior_nb<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                      (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                      (1|Powder.lvl)+
                                                      (1|species),
                                                    data=dff_ectos_network_individual_metrics,
                                                    family=negbinomial(),  #zero_inflated_negbinomial()
                                                    data2 = list(phy_cov=phy_cov),
                                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                    thin=2,
                                                    save_pars = save_pars(all=  TRUE),
                                                    prior = c(prior_predictors,prior_random),
                                                    control=list(adapt_delta=0.99, max_treedepth=14)) 
saveRDS(NZ_a_mites_brms_bayes_no_degree_int_prior_nb, "data/data_analyses/model_selection/M2MND_NZ.model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_nb.RDS")
NZ_a_lice_brms_bayes_no_degree_int_prior_nb<-readRDS("data/data_analyses/model_selection/M2MND_NZ.model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_nb.RDS")


# extremily slow!! truncating the data will be teh way of doig this 
NZ_a_mites_brms_bayes_no_degree_int_prior_truncated<-brms::brm(total_no_feathers_mites|trunc(lb=1)~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                          (1|Powder.lvl)+
                                                          (1|species),
                                                        data=dff_ectos_network_individual_metrics,
                                                        family=negbinomial(),  #zero_inflated_negbinomial()
                                                        data2 = list(phy_cov=phy_cov),
                                                        iter=500, warmup=200, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                        thin=2,
                                                        save_pars = save_pars(all=  TRUE),
                                                        prior = c(prior_predictors,prior_random),
                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 

NZ_a_mites_brms_bayes_no_degree_int_prior_p<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                        (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                        (1|Powder.lvl)+
                                                        (1|species),
                                                      data=dff_ectos_network_individual_metrics,
                                                      family=poisson(),  #zero_inflated_negbinomial()
                                                      data2 = list(phy_cov=phy_cov),
                                                      iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                      thin=2,
                                                      save_pars = save_pars(all=  TRUE),
                                                      prior = c(prior_predictors,prior_random),
                                                      control=list(adapt_delta=0.99, max_treedepth=14))

# Evaludating model fit and convergence

simulate_residuals <- dh_check_brms(NZ_a_mites_brms_bayes_no_degree_int_prior_n, integer = TRUE)
plot(simulate_residuals, form = dff_ectos_network_individual_metrics$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations

# Model comparison 
loo_nb<-(loo(NZ_a_mites_brms_bayes_no_degree_int_prior_nb, moment_match=TRUE))
k_ecto_NZ_MITES_brms_no_int_prior_nb<-kfold(NZ_a_mites_brms_bayes_no_degree_int_prior_nb, K=10)
#saveRDS(NZ_a_mites_brms_bayes_no_degree_int_prior_nb, "data/data_analyses/model_selection/k_fold/K_fold_AMNZ_DEGREE_model_mites_abundance_brms_no_interactions_priors_nbinomial.RDS")

loo_p<-(loo(NZ_a_mites_brms_bayes_no_degree_int_prior_p,moment_match=TRUE))
k_ecto_NZ_MITES_brms_no_int_prior_p<-kfold(loo_p, K=10)
#saveRDS(k_ecto_NZ_lice_brms_no_int_prior_nb, "data/data_analyses/model_selection/k_fold/K_fold_AMNZ_DEGREE_model_mites_abundance_brms_no_interactions_priors_poisson.RDS")

loo(NZ_a_lice_brms_bayes_no_degree_int_prior_nb,NZ_a_lice_brms_bayes_no_degree_int_prior_p, compare=TRUE)

# The negaive binomial performs much better than the poission ! model selected
loo_compare(k_ecto_NZ_MITES_brms_no_int_prior_nb,k_ecto_NZ_MITES_brms_no_int_prior_p)

# plots 
color_scheme_set("red") 

conditional_effects(model)
marginal_effects()
plot(
  conditional_effects(fit, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1)
)
# poisson 
#model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig2_MNDNZ_NB_plot_model_MITES_ABUNDANCE_CONVERGENCEbrms_bayes_no_int_DEGREE_nb.png",width = 3000, height = 3000, res = 300, units = "px")
plot(NZ_a_mites_brms_bayes_no_degree_int_prior_nb)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig2_MNDNZ_NB_plot_modell_FIT_MITES_ABUNDANCE_brms_bayes_no_int_DEGREE_nb.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(NZ_a_mites_brms_bayes_no_degree_int_prior_nb, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(NZ_a_mites_brms_bayes_no_degree_int_prior_nb,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions NZ_LICE ABUNDANCE_DEGREE", subtitle ="NZ_LICE_ABUNDANCE_DEGREE with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2_MNDNZ_NB_plot_model_parameters_MITES_ABUNDANCE_brms_bayes_no_int_degree.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(NZ_a_mites_brms_bayes_no_degree_int_prior_nb,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions NZ_LICE ABUNDANCE", subtitle ="NZ_LICE_ABUNDANCE with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2_MNDNZ_NB_plot_model_parameters_intervals_MITES_ABUNDANCE_brms_bayes_no_int_degree.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()
















# NETWORK ABUNDANCE # #### 6.7****Selected ****model NETWORKS and mites abundance mass---------

###_###_###_###_###_###_###_###_###_###_
#MITES ABUNDANCE DEGREE###
###_###_###_###_###_###_###_###_###_###_
#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions 
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?
residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default so no need to speciied
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default so no need to speciied


selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                               scale(mass_tidy_species)+(mass_ind_comp)+
                                                            (1|gr(species_jetz, cov = phy_cov))+
                                                            (1|Powder.lvl)+
                                                            (1|species),
                                                          data=dff_ectos_network_individual_metrics,
                                                          family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                          data2 = list(phy_cov=phy_cov),
                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                          save_pars = save_pars(all=  TRUE),
                                                          prior = c(prior_predictors,prior_random),
                                                          thin=2,
                                                          control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, "results/selected_models/4_3_3_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")
selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-readRDS( "results/selected_models/4_3_3_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")

selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior2<-brms::brm(total_no_feathers_mites~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                     scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                                     (1|gr(species_jetz, cov = phy_cov))+
                                                                     (1|Powder.lvl)+
                                                                     (1|species),
                                                                   data=dff_ectos_network_individual_metrics,
                                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                                   data2 = list(phy_cov=phy_cov),
                                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                   save_pars = save_pars(all=  TRUE),
                                                                   prior = c(prior_predictors,prior_random),
                                                                   thin=2,
                                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior2, "results/selected_models/4_3_3_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior_individual_mass_scaled.RDS")

###_###_###_##
#PLOTS
###_###_###_##
marginal_effects()
plot(
  conditional_effects(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1)
)
bayes_R2() # R2 0.1529
color_scheme_set("green")

# model convergence 
png("results/selected_models_figures/4_3_3_Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_CONVERGENCE_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("results/selected_models_figures/4_3_3_Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_FIT_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES

estimates_plot<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~DEGREE")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_3_Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_ESTIMATES_int_DEGREE.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior2,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~DEGREE")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_3_Fig2MNS_BEST_zinb_MITES_ABUNDANCE_brms_bayes_INTERVALS_no_int_DEGREE_ind_mass_scaled.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


###_###_###_###_###_###_###_###_###_###_
#MITES ABUNDANCE Strength###
###_###_###_###_###_###_###_###_###_###_
#PRIORS
prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
#prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

selected_zinb_a_nf_mites_brms_bayes_no_int_strength_prior<-brms::brm(total_no_feathers_mites~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                     scale(mass_tidy_species)+(mass_ind_comp)+
                                                                     (1|gr(species_jetz, cov = phy_cov))+
                                                                     (1|Powder.lvl)+
                                                                     (1|species),
                                                                   data=dff_ectos_network_individual_metrics,
                                                                   family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                                   data2 = list(phy_cov=phy_cov),
                                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                   save_pars = save_pars(all=  TRUE),
                                                                   prior = c(prior_predictors,prior_random),
                                                                   thin=2,
                                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(selected_zinb_a_nf_mites_brms_bayes_no_int_strength_prior, "results/selected_models/4_3_4_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_strength_prior.RDS")
selected_zinb_a_nf_mites_brms_bayes_no_int_strength_prior<-readRDS( "results/selected_models/4_3_4_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_strength_prior.RDS")

selected_zinb_a_nf_mites_brms_bayes_no_int_strength_prior2<-brms::brm(total_no_feathers_mites~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                                       scale(mass_tidy_species)+scale(mass_ind_comp)+
                                                                       (1|gr(species_jetz, cov = phy_cov))+
                                                                       (1|Powder.lvl)+
                                                                       (1|species),
                                                                     data=dff_ectos_network_individual_metrics,
                                                                     family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                                     data2 = list(phy_cov=phy_cov),
                                                                     iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                     save_pars = save_pars(all=  TRUE),
                                                                     prior = c(prior_predictors,prior_random),
                                                                     thin=2,
                                                                     control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(selected_zinb_a_nf_mites_brms_bayes_no_int_strength_prior2, "results/selected_models/4_3_4_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_strength_prior_ind_mass_scaled.RDS")
selected_zinb_a_nf_mites_brms_bayes_no_int_strength_prior<-readRDS( "results/selected_models/4_3_4_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_strength_prior.RDS")


###_###_###_##
#PLOTS
###_###_###_##
marginal_effects()
plot(
  conditional_effects(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1)
)


summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("green")

# model convergence 
png("results/selected_models_figures/4_3_4_Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_CONVERGENCE_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
png("results/selected_models_figures/4_3_4_Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_int_sociality_priors_FIT_STRENGTH.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES

estimates_plot<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~STRENGTH")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_4_Fig2MNS_BEST_zinb_ABUNDANCE_MITES_brms_bayes_ESTIMATES_int_STRENGTH.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior2,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~STRENGTH")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("results/selected_models_figures/4_3_4_Fig2MNS_BEST_zinb_MITES_ABUNDANCE_brms_bayes_INTERVALS_no_int_STRENGTH_ind_mass_scaled.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()




# # ##### 7 Modeling LICE DIVERSITY -------------------------------------------

#The lice diversity file
#lice_diversity_raw<-read.csv("data/data_raw/7.ectoparasite_raw_lice_diversity.csv",na.strings =c("","NA"))

#lice_diversity_raw<-read.csv("data/data_raw/7.ectoparasite_raw_lice_diversity.csv",na.strings =c("","NA","-")) %>% 
 # filter(Collection=="JEJ")

#lice_diversity<-lice_diversity_raw %>% 
#  group_by(species_binomial) %>%
#  summarise(cumulative_richness=n_distinct(na.omit(lice_genus_tidy)), sample_size=n_distinct(Full_Label), host_family=first(host_family)) 
  

# the general parameters 
#ectos_birds_df<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
#select(elevation_midpoint,species_binomial,species_jetz,foraging_cat, sociality, mass_tidy_species) %>% 
#na.omit() %>% 
#group_by(species_jetz) %>% 
#summarise(elevation_midpoint=first(elevation_midpoint),
#         species_binomial=first(species_binomial),
#         foraging_cat=first(foraging_cat), 
#          sociality=first(sociality), 
#         mass_tidy_species=first(mass_tidy_species),
#         total_sample_size=n())


# dff_lice_diversity_general<-right_join(lice_diversity,ectos_birds_df, by="species_binomial" )
  
#write.csv(dff_lice_diversity_general, "data/data_manuscript/3_dff_ectos_diversity_species_prevalence_abundance_diversity_elevation_mass_FILE_TIDY.csv")

####
#THE DATA
####

dff_lice_diversity<-read.csv("data/data_manuscript/3_dff_ectos_diversity_species_prevalence_abundance_diversity_elevation_mass_FILE_TIDY.csv",na.strings =c("","NA"))%>% filter(total_sample_size>9)
dff_lice_diversity$cumulative_richness[is.na(dff_lice_diversity$cumulative_richness)] = 0

unique(dff_lice_diversity$cumulative_richness)
unique(dff_lice_diversity$species_jetz)


dff_lice_diversity %>% group_by(sociality) %>% 
  summarize(n())

# the phylo data
phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos social and non social so we need to trim it 

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(dff_lice_diversity$species_jetz)) %>% mutate(name=dff_lice_diversity$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 

# phylogenetic correlation structure, create a covariance matrix of species
phy_cov<-ape::vcv(phylo, corr=TRUE)

# the data structure

#dff_ectos_network_individual_metrics$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
dff_lice_diversity$cumulative_richness<-as.numeric(dff_lice_diversity$cumulative_richness)
dff_lice_diversity$sociality<-as.factor(dff_lice_diversity$sociality)
dff_lice_diversity$sociality_groups<-as.factor(dff_lice_diversity$sociality_groups)

dff_lice_diversity$elevation_midpoint<-as.numeric(dff_lice_diversity$elevation_midpoint)
dff_lice_diversity$mass_tidy_species<-as.numeric(dff_lice_diversity$mass_tidy_species)
dff_lice_diversity$total_sample_size<-as.numeric(dff_lice_diversity$total_sample_size)


dff_lice_diversity$foraging_cat<-as.factor(dff_lice_diversity$foraging_cat)
dff_lice_diversity$species_jetz<-as.factor(dff_lice_diversity$species_jetz)
dff_lice_diversity$species<-as.factor(dff_lice_diversity$species_jetz) # create a column for the species effect different to the phylogenetic one

unique(dff_lice_diversity$cumulative_richness)

#View(dff_lice_diversity)

# model 
#prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
#prior_random<- prior("normal(0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 

prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

prior_summary(selected_poisson_lice_diversity_sociality_no_int_priors_trunc)
unique(dff_lice_diversity$cumulative_richness)


#### non trncuating the distribution 


selected_poisson_lice_diversity_sociality_no_int_priors<-brms::brm(cumulative_richness~sociality_groups+
                                                                     scale(elevation_midpoint)+scale(total_sample_size)+
                                                                     scale(mass_tidy_species)+
                                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                     (1|species),
                                                                   data=dff_lice_diversity,
                                                                   family=poisson(),  #zero_inflated_negbinomial()
                                                                   data2 = list(phy_cov=phy_cov),
                                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                   thin=2,
                                                                   prior = c(prior_predictors,prior_random),
                                                                   save_pars = save_pars(all=  TRUE),
                                                                   control=list(adapt_delta=0.999, max_treedepth=14)) 

saveRDS(selected_poisson_lice_diversity_sociality_no_int_priors, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_NO_truncated_antfollowers_included.RDS")


# selected_poisson_lice_diversity_sociality_no_int_priors_trunc<-brms::brm(cumulative_richness|trunc(lb=0, ub=3)~sociality+
#                                                                            scale(elevation_midpoint)+scale(mass_tidy_species)+
#                                                                            scale(total_sample_size)+
#                                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
#                                                                     (1|species),
#                                                                   data=dff_lice_diversity,
#                                                                   family=poisson("log"), # "sqrt" #zero_inflated_negbinomial()
#                                                                   data2 = list(phy_cov=phy_cov),
#                                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                                                   thin=2,
#                                                                   prior = c(prior_predictors,prior_random),
#                                                                   save_pars = save_pars(all=  TRUE),
#                                                                   control=list(adapt_delta=0.999, max_treedepth=14)) 
# 
# 
# saveRDS(selected_poisson_lice_diversity_sociality_no_int_priors_trunc, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_trunc.RDS")
# selected_poisson_lice_diversity_sociality_no_int_priors_trunc<-readRDS("results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_trunc.RDS")

# similar analyses but with sample size 10 or more

# sd(scale(dff_lice_diversity$mass_tidy_species))
# sd()
# selected_poisson_lice_diversity_sociality_no_int_priors_trunc_sample10<-brms::brm(cumulative_richness|trunc(lb=0, ub=3)~sociality+
#                                                                            scale(elevation_midpoint)+scale(mass_tidy_species)+
#                                                                            scale(total_sample_size)+
#                                                                            (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
#                                                                            (1|species),
#                                                                          data=dff_lice_diversity,
#                                                                          family=poisson("log"), # "sqrt" #zero_inflated_negbinomial()
#                                                                          data2 = list(phy_cov=phy_cov),
#                                                                          iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
#                                                                          thin=2,
#                                                                          prior = c(prior_predictors,prior_random),
#                                                                          save_pars = save_pars(all=  TRUE),
#                                                                          control=list(adapt_delta=0.999, max_treedepth=14)) 
# 
# 
# saveRDS(selected_poisson_lice_diversity_sociality_no_int_priors_trunc_sample10, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_trunc_10_samples.RDS")
# selected_poisson_lice_diversity_sociality_no_int_priors_trunc<-readRDS("results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_trunc_10_samples.RDS")




simulate_residuals <- dh_check_brms(selected_poisson_lice_diversity_sociality_no_int_priors_trunc2, integer = TRUE)
plot( simulate_residuals, form = dff_lice_diversity$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations


# plots

### Our data is not zero inflated so we used poisson
# plots 
color_scheme_set("red") 

# poisson 


png("trunc.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_poisson_lice_diversity_sociality_no_int_priors_trunc, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()
png("trunc2.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

#model convergence 
png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_CONVERGENCE_LICE_DIVERSITY_poisson_trunc.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_poisson_lice_diversity_sociality_no_int_priors_trunc)
dev.off()

# model fit
png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_FIT_LICE_DIVERSITY_poisson_trunc.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_poisson_lice_diversity_sociality_no_int_priors_trunc, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES

estimates_plot<-mcmc_plot(selected_poisson_lice_diversity_sociality_no_int_priors_trunc,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="LICE DIVERSITY ")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_parameters_LICE_DIVERSITY_poisson_trunc.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_sociality_no_int_priors_trunc,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY ")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_parameters_intervals_LICE_DIVERSITY_poisson_trunc.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


# # #### Data processing LICE DIVERSITY  NETWORKS and models-------------------------------------------

# Getting the data ready

#The lice diversity file

#lice_diversity_raw<-read.csv("data/data_raw/7.ectoparasite_raw_lice_diversity.csv",na.strings =c("","NA")) %>% 
#  filter(Collection=="JEJ")

# lice_diversity<-lice_diversity_raw %>% 
#  group_by(species_binomial) %>%
# summarise(cumulative_richness=n_distinct(na.omit(lice_genus_tidy)), sample_size=n_distinct(Full_Label), host_family=first(host_family)) 
 
# the other parameters that we need for the model 

# ectos_birds_df<-read.csv("data/data_manuscript/3_dff_all_ectos_network_metrics_individuals_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
#  select(elevation_midpoint,species_binomial,species_jetz,foraging_cat, sociality, mass_tidy_species,degree_species, w_degree_species) %>% 
# na.omit() %>% 
# group_by(species_jetz) %>% 
# summarise(elevation_midpoint=first(elevation_midpoint),
#          species_binomial=first(species_binomial),
#          total_sample_size=n(),
#           foraging_cat=first(foraging_cat), 
#           sociality=first(sociality),
#           mass_tidy_species=first(mass_tidy_species),
#           degree_species=first(degree_species),
#           degree_w_species=first(w_degree_species))

#dff_lice_diversity_networks<-right_join(lice_diversity,ectos_birds_df, by="species_binomial" ) 

#write.csv(dff_lice_diversity_networks, "data/data_manuscript/3_dff_ectos_diversity_network_metrics_species_FILE_TIDY.csv")

### ### ###
# the data 
### ### ###

dff_lice_diversity_networks<-read.csv("data/data_manuscript/3_dff_ectos_diversity_network_metrics_species_FILE_TIDY.csv") %>% filter(total_sample_size>4)%>% 
  select(-sample_size)
dff_lice_diversity_networks$cumulative_richness[is.na(dff_lice_diversity_networks$cumulative_richness)] = 0
#str(dff_lice_diversity_networks)
phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos social and non social so we need to trim it 

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(dff_lice_diversity_networks$species_jetz)) %>% mutate(name=dff_lice_diversity_networks$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 

# phylogenetic correlation structure, create a covariance matrix of species
phy_cov<-ape::vcv(phylo, corr=TRUE)

names(dff_lice_diversity_networks)
# data structure 

dff_lice_diversity_networks$cumulative_richness<-as.numeric(dff_lice_diversity_networks$cumulative_richness)
dff_lice_diversity_networks$sociality<-as.factor(dff_lice_diversity_networks$sociality)
#dff_lice_diversity_networks$degree_species<-as.numeric(dff_lice_diversity_networks$degree_species)
#dff_lice_diversity_networks$degree_w_species<-as.numeric(dff_lice_diversity_networks$degree_w_species)

dff_lice_diversity_networks$elevation_midpoint<-as.numeric(dff_lice_diversity_networks$elevation_midpoint)
dff_lice_diversity_networks$mass_tidy_species<-as.numeric(dff_lice_diversity_networks$mass_tidy_species)
dff_lice_diversity_networks$total_sample_size<-as.numeric(dff_lice_diversity_networks$total_sample_size)

dff_lice_diversity_networks$foraging_cat<-as.factor(dff_lice_diversity_networks$foraging_cat)
dff_lice_diversity_networks$species_jetz<-as.factor(dff_lice_diversity_networks$species_jetz)
dff_lice_diversity_networks$species<-as.factor(dff_lice_diversity_networks$species_jetz) # create a column for the species effect different to the phylogenetic one


#model

#prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_random<- prior("normal(0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 

#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?



selected_poisson_lice_diversity_degree_no_int_priors<-brms::brm(cumulative_richness~scale(degree_species)+
                                                                     scale(elevation_midpoint)+
                                                                     scale(total_sample_size)+
                                                                     scale(mass_tidy_species)+
                                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                     (1|species),
                                                                   data=dff_lice_diversity_networks,
                                                                   family=poisson(),  #zero_inflated_negbinomial()
                                                                   data2 = list(phy_cov=phy_cov),
                                                                   iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                   thin=2,
                                                                   prior = c(prior_predictors,prior_random),
                                                                   save_pars = save_pars(all=  TRUE),
                                                                   control=list(adapt_delta=0.999, max_treedepth=14)) 




saveRDS(selected_poisson_lice_diversity_degree_no_int_priors, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_degree.RDS")


# same bu only with more tahn 10 samples no truncated 

dff_lice_diversity_networks<-read.csv("data/data_manuscript/3_dff_ectos_diversity_network_metrics_species_FILE_TIDY.csv") %>% filter(total_sample_size>10)%>% 
  select(-sample_size)
dff_lice_diversity_networks$cumulative_richness[is.na(dff_lice_diversity_networks$cumulative_richness)] = 0
#str(dff_lice_diversity_networks)
phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos social and non social so we need to trim it 

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(dff_lice_diversity_networks$species_jetz)) %>% mutate(name=dff_lice_diversity_networks$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 

# phylogenetic correlation structure, create a covariance matrix of species
phy_cov<-ape::vcv(phylo, corr=TRUE)

names(dff_lice_diversity_networks)
# data structure 

dff_lice_diversity_networks$cumulative_richness<-as.numeric(dff_lice_diversity_networks$cumulative_richness)
dff_lice_diversity_networks$sociality<-as.factor(dff_lice_diversity_networks$sociality)
#dff_lice_diversity_networks$degree_species<-as.numeric(dff_lice_diversity_networks$degree_species)
#dff_lice_diversity_networks$degree_w_species<-as.numeric(dff_lice_diversity_networks$degree_w_species)

dff_lice_diversity_networks$elevation_midpoint<-as.numeric(dff_lice_diversity_networks$elevation_midpoint)
dff_lice_diversity_networks$mass_tidy_species<-as.numeric(dff_lice_diversity_networks$mass_tidy_species)
dff_lice_diversity_networks$total_sample_size<-as.numeric(dff_lice_diversity_networks$total_sample_size)

dff_lice_diversity_networks$foraging_cat<-as.factor(dff_lice_diversity_networks$foraging_cat)
dff_lice_diversity_networks$species_jetz<-as.factor(dff_lice_diversity_networks$species_jetz)
dff_lice_diversity_networks$species<-as.factor(dff_lice_diversity_networks$species_jetz) # create a column for the species effect different to the phylogenetic one


#model

#prior_predictors<-prior("normal(0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_predictors<-prior("student_t(3,0,10)", class ="b") # This prior is generating divergent transitions os i move to amore weakly informative parameter
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
#prior_random<- prior("normal(0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 

#prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

selected_poisson_lice_diversity_degree_no_int_priors_10<-brms::brm(cumulative_richness~scale(degree_species)+
                                                                  scale(elevation_midpoint)+scale(mass_tidy_species)+ scale(total_sample_size)+
                                                                  (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                  (1|species),
                                                                data=dff_lice_diversity_networks,
                                                                family=poisson(),  #zero_inflated_negbinomial()
                                                                data2 = list(phy_cov=phy_cov),
                                                                iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                thin=2,
                                                                prior = c(prior_predictors,prior_random),
                                                                save_pars = save_pars(all=  TRUE),
                                                                control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(selected_poisson_lice_diversity_degree_no_int_priors_10, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_degree_10.RDS")

# plots 

# poisson 
#model convergence 
png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_CONVERGENCE_LICE_DIVERSITY_DEGREE_poisson_trunc.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_poisson_lice_diversity_degree_no_int_priors_10)
dev.off()

# model fit
png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_FIT_LICE_DIVERSITY_DEGREE_poisson.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_poisson_lice_diversity_degree_no_int_priors_10, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_poisson_lice_diversity_degree_no_int_priors_10,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="LICE DIVERSITY ~DEGREE")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_parameters_LICE_DIVERSITY_DEGREE_poisson_trunc.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_degree_no_int_priors_10,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY~DEGREE ")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_parameters_intervals_LICE_DIVERSITY_DEGREE_poisson.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

### TRUNCATED


#truncated is not working well several 50 divergent transitions
selected_poisson_lice_diversity_degree_no_int_priors_trunc<-brms::brm(cumulative_richness|trunc(lb=0,ub=3)~scale(degree_species)+
                                                                          scale(elevation_midpoint)+scale(mass_tidy_species)+
                                                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                           scale(total_sample_size)+
                                                                          (1|species),
                                                                        data=dff_lice_diversity_networks,
                                                                        family=poisson(),  #zero_inflated_negbinomial()
                                                                        data2 = list(phy_cov=phy_cov),
                                                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                        thin=2,
                                                                        prior = c(prior_predictors,prior_random),
                                                                        save_pars = save_pars(all=  TRUE),
                                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 



#saveRDS(selected_poisson_lice_diversity_degree_no_int_priors_trunc, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_degree.RDS")

selected_poisson_lice_diversity_degree_no_int_priors_trunc<-readRDS("results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_degree.RDS")

# poisson 
#model convergence 
png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_CONVERGENCE_LICE_DIVERSITY_DEGREE_poisson_trunc.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_poisson_lice_diversity_degree_no_int_priors_trunc)
dev.off()

# model fit
png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_FIT_LICE_DIVERSITY_DEGREE_poisson.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_poisson_lice_diversity_degree_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_poisson_lice_diversity_degree_no_int_priors_trunc,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="LICE DIVERSITY ~DEGREE")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_parameters_LICE_DIVERSITY_DEGREE_poisson_trunc.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_degree_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY~DEGREE ")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DLD_BEST_plot_model_parameters_intervals_LICE_DIVERSITY_DEGREE_poisson.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()
### w_degree
selected_poisson_lice_diversity_w_degree_no_int_priors<-brms::brm(cumulative_richness~scale(degree_w_species)+
                                                                          scale(elevation_midpoint)+ scale(total_sample_size)+scale(mass_tidy_species)+
                                                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                          (1|species),
                                                                        data=dff_lice_diversity_networks,
                                                                        family=poisson(),  #zero_inflated_negbinomial()
                                                                        data2 = list(phy_cov=phy_cov),
                                                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                        thin=2,
                                                                        save_pars = save_pars(all=  TRUE),
                                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 


saveRDS(selected_poisson_lice_diversity_w_degree_no_int_priors, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_w_degree.RDS")
selected_poisson_lice_diversity_w_degree_no_int_priors<-readRDS( "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_w_degree.RDS")


selected_poisson_lice_diversity_w_degree_no_int_priors_trunc<-brms::brm(cumulative_richness|trunc(lb=0,ub=3)~scale(degree_w_species)+
                                                                          scale(elevation_midpoint)+ scale(total_sample_size)+scale(mass_tidy_species)+
                                                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                                          (1|species),
                                                                        data=dff_lice_diversity_networks,
                                                                        family=poisson(),  #zero_inflated_negbinomial()
                                                                        data2 = list(phy_cov=phy_cov),
                                                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                                        thin=2,
                                                                        save_pars = save_pars(all=  TRUE),
                                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 



saveRDS(selected_poisson_lice_diversity_w_degree_no_int_priors_trunc, "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_w_degree_trunc.RDS")



# poisson 
#model convergence 
png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_CONVERGENCE_LICE_DIVERSITY_poisson_W_DEGREE_trunc.png",width = 3000, height = 3000, res = 300, units = "px")
plot(selected_poisson_lice_diversity_w_degree_no_int_priors_trunc)
dev.off()

# model fit
png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_FIT_LICE_DIVERSITY_poisson_W_DEGREE_trunc.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(selected_poisson_lice_diversity_w_degree_no_int_priors_trunc, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_poisson_lice_diversity_sociality_no_int_priors_trunc,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="LICE DIVERSITY ~W DEGREE")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_parameters_LICE_DIVERSITY_poisson_WDEGREE_trunc.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_w_degree_no_int_priors_trunc,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY ~W_degree")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

png("results/selected_models_figures/5_Fig_DL_BEST_plot_model_parameters_intervals_LICE_DIVERSITY_poisson_W_DEGREE_trunc.png",width = 4000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


bplot<-plot(
  conditional_effects(selected_poisson_lice_diversity_degree_no_int_priors_trunc, dpar = "mu"), 
  points = TRUE, 
  point_args = list(width = .05, shape = 1))



# #### Extrapolated richness values  --------------------------------------------

library(vegan)

data(dune)
View(dune)
data(dune.env)
View(dune.env)
pool <- with(dune.env, specpool(dune, Management))
pool
op <- par(mfrow=c(1,2))
boxplot(specnumber(dune) ~ Management, data = dune.env,
        col = "hotpink", border = "cyan3")
boxplot(specnumber(dune)/specpool2vect(pool) ~ Management,
        data = dune.env, col = "hotpink", border = "cyan3")
par(op)
data(BCI)
View(BCI)
## Accumulation model
pool <- poolaccum(BCI)
summary(pool, display = "chao")
plot(pool)
## Quantitative model
estimateR(BCI[1:5,])



specpool(BCI[1:5,])


#### LOAD PREVIOUSLY SAVED MODELS #####
###_###_###_###_###_###_###_###_###_###_###_###_

#models

#Prevalence

ecto_p_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")
ecto_p_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions_priors_default_intercept.RDS")

ecto_p_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE.RDS")
ecto_p_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")


#lICE
zinb_a_lice_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS") # lice abundance default priors             !!!! warrning 20 pareto k are high!!!! what to do?
zinb_a_lice_brms_bayes_all_interactions_priors<-readRDS("data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_priors.RDS")# Lice abundance weakly informative priors !!!! warrning 20 pareto k are high!!!! what to do?

zinb_a_lice_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")
zinb_a_lice_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")


#MITES
zinb_a_nf_mites_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")# Mites abundance default priors 
zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/1NFM.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")# Mites abundance weakly informative priors 


# Mites abundance weakly by degree informative priors 


# NOTES Checking model overfitting and missbehavior and model comparison  ----------------------------------------------

# they promote the use of  WAIC (Widely applicable information criterion), and LOO (leave-one-out cross-validation). These use calculations of the log-likelihood across the entire posterior.
#For predictive models, a common method for picking models is to test the predictive accuracy of a model on new data, or some held-out portion of  the data you already have (but not the portion you used to build the model).
#If your model does well predicting a random datapoint excluded from your model when you build it, it is likely not overfit. However, re-running your model  many times, dropping a different data point every time would take a lot of time making it unpractical for complex models.
#The "loo" package in R employs a algorithm to approximate what the performance of your model would be if you performed full leave-one-out cross-validation.
#If you are interested in how the approximation is calculated,  and how to determine if you can reasonably make this approximation  for you model, you can read the original paper by Aki Vehtari which  explains WAIC and LOO: https://arxiv.org/abs/1507.04544 and check out the example on CRAN: https://cran.r-project.org/web/packages/loo/index.html


#As the proportion of zeros is quite high in the data, it is worthwhile to test also a zero-inflated negative-binomial model, 
#which is a mixture of two models - logistic regression to model the proportion of extra zero counts - negative-binomial model
# notes from Max farell https://github.com/maxfarrell/qcbs_stan_workshop/blob/master/QCBS_stan.Rmd
# some step by step in brms and bayesian https://ourcodingclub.github.io/tutorials/brms/
#Notes check this blog for discussion in STAN BRMS situations: https://discourse.mc-stan.org
# check this blog for discussion in model selection https://discourse.mc-stan.org/t/model-selection-in-brms/30492
# about loo package interpretation https://mc-stan.org/loo/reference/loo-glossary.html
# Using the loo package https://mc-stan.org/loo/articles/loo2-example.html
# some examples here : https://avehtari.github.io/modelselection/roaches.html
# Infereing from P_loo https://discourse.mc-stan.org/t/a-quick-note-what-i-infer-from-p-loo-and-pareto-k-values/3446
# Avoidng model refit with moment match when pareko is high https://mc-stan.org/loo/articles/loo2-moment-matching.html
#how to use the loo package to carry out Pareto smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO) for purposes of model checking and model comparison. https://mc-stan.org/loo/articles/loo2-example.html

# Notes for model selection with BRMS, [Work in progress]
#loo gives us warnings about the Pareto diagnostics, which indicate that for some observations the leave-one-out posteriors are different enough from the full posterior that importance-sampling is not able to correct the difference. We can see more details by printing the loo object.
# definition [elpd_loo] The ELPD is the theoretical expected log pointwise predictive density for a new dataset 
# definition [elpd_loo SE] This standard error is a coarse description of our uncertainty about the predictive performance for unknown future data. 
# In well behaving cases p_loo < N and p_loo < p, where p is the total number of parameters in the model and n ( number of observations)
#elpd_loo differences of more than 4 magnitude are consider lareg enough to make a difference in teh model
# elpd_loo large compared to se_diff indicate changes in the model selection too
# notes more model overfiiting in BRMS 
# Look at the number of observations (226 is shown in the first line), 
# The effective number of parameters, p_loo ()  p_loo < N (number of observations) and p_loo < p (parameters of the model)
# and the number of high Pareto k’s ( All Pareto k’s are good, so there are no highly influential individual observations, This indicates that loo computation is reliable)

#To improve the accuracy of the loo() result above, we could perform leave-one-out cross-validation by explicitly leaving out single observations 
#and refitting the model using MCMC repeatedly. However, the Pareto 𝑘diagnostics indicate that there are 19 observations which are problematic. 
#This would require 19 model refits which may require a lot of computation time.
loo(zip_a_nf_mites_brms_bayes_no_int,reloo = TRUE) # refits the model leaving the observation with hgh pareto out , it is an slow alternative 
#Instead of refitting with MCMC, we can perform a faster moment matching correction to the importance sampling for the problematic observations. 
loo_2<-loo(zip_a_nf_mites_brms_bayes_no_int, moment_match = TRUE)
plot(loo_2)




#Steps to check model misspecification 
# i) Calculate loo 
loo1<-loo(ecto_p_brms_bayes_no_int,save_psis = TRUE)
#ii) Plot pareto k diagnostics
plot(loo1)
#iii) leave-one-out cross-validation marginal posterior predictive checks Gabry et al (2018) in combination with bayesplot()
# For a good model, the distribution of LOO-PIT values should be uniform
#LOO-PIT values for our model (thick curve) is compared to many independently generated samples (each the same size as our dataset) 
yrep <- posterior_predict(ecto_p_brms_bayes_no_int)

ppc_loo_pit_overlay(
  y = ectos_birds_dff$ectoparasites_PA,
  yrep = yrep,
  lw = weights(loo1$psis_object)
)

#What to do if I have many high Pareto 𝑘?

#If 𝑘̂<0.5 k  then the corresponding component of elpd_loo is estimated with high accuracy. If 0.5<𝑘̂<0.7the accuracy is lower, but still OK.
#If 𝑘̂>0.7, then importance sampling is not able to provide useful estimate for that component/observation. Pareto-𝑘̂is also useful as a measure of influence of an observation. 
#Highly influential observations have high 𝑘̂values. Very high 𝑘̂values often indicate model misspecification, outliers or mistakes in data processing. 
#If p_loo <𝑝 and the number of parameters 𝑝 is relatively large compared to the number of observations (e.g., 𝑝>𝑁/5), it is likely that the model is so flexible or the population prior so weak that it’s difficult to predict the left out observation (even for the true model). This happens, for example, in the simulated 8 schools (Vehtari, Gelman and Gabry, 2017), random effect models with a few observations per random effect, and Gaussian processes and spatial models with short correlation lengths.

# Solution 1 

#One option is to re-fit the model without these problematic observations, 
#and directly calculate the loo statistic directly for them.

#Exclude those with high values 
loo2<-loo(MODEL, moment_match = TRUE)


#loo2 <- loo(model, k_threshold=0.7)

# i THINK THIS IS MORE OF OUR CASES
#If p_loo < p and the number of parameters p is relatively large compared to the number of observations (e.g., p>N/5),
#it is likely that the model is so flexible or the population prior so weak that it’s difficult to predict the left out observation (even for the true model). This happens, for example, in the simulated 8 schools (in VGG2017), random effect models with a few observations per random effect, and Gaussian processes and spatial models with short correlation lengths.


# Using ELPD for model selection 
# ELPD is good for model comparison as it measures the goodness of the whole predictive distribution
#As quick rule: If elpd difference (elpd_diff in loo package) is less than 4, the difference is small (Sivula, Magnusson and Vehtari, 2020). 
#If elpd difference (elpd_diff in loo package) is larger than 4, then compare that difference to standard error of elpd_diff (provided e.g. by loo package) (Sivula, Magnusson and Vehtari, 2020).
#The value for deciding what is small or large can be based on connection to Pseudo-BMA+-weights (Yao et al., 2018). See also How to interpret in Standard error (SE) of elpd difference (elpd_diff)?. https://mc-stan.org/loo/articles/online-only/faq.html#se_diff
#SE tends to be underestimated especially if the number of observations is small or the models are badly misspecified.
# When the difference (elpd_diff) is larger than 4, the number of observations is larger than 100 and the model is not badly misspecified then normal approximation and SE are quite reliable description of the uncertainty in the difference. 
#Differences smaller than 4 are small and then the models have very similar predictive performance and it doesn’t matter if the normal approximation fails or SE is underestimated (Sivula, Magnusson and Vehtari, 2020).



#Why is scale not the same than  standardize
#install.packages("arm")
#library(arm)
#mean(standardize(dff_ectos_network_individual_metrics$degree))


# Fancy plots  ------------------------------------------------------------
#From https://www.andrewheiss.com/blog/2021/11/08/beta-regression-guide/#b-beta-regression-bayesian-style

  ame_fancy_zi_polyarchy_quota <- fancy_model %>% 
  emmeans(~ quota + polyarchy,
          at = list(polyarchy = seq(0, 1, by = 0.1)),
          epred = TRUE,
          re_formula = NULL) %>% 
  gather_emmeans_draws()
ggplot(ame_fancy_zi_polyarchy_quota,
       aes(x = polyarchy, y = .value, color = quota, fill = quota)) +
  stat_lineribbon(aes(fill_ramp = stat(level))) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  scale_fill_ramp_discrete(range = c(0.2, 0.7)) +
  facet_wrap(vars(quota), ncol = 2,
             labeller = labeller(quota = c(`TRUE` = "Quota",
                                           `FALSE` = "No Quota"))) +
  labs(x = "Polyarchy (democracy)",
       y = "Predicted proportion of women MPs",
       fill = "Quota", color = "Quota",
       fill_ramp = "Credible interval") +
  theme_clean() +
  

# Fancy plots  ------------------------------------------------------------

# combine the summaries for two models
bind_rows(
  # m11.10
  fixef(m11.10) %>% 
    data.frame() %>% 
    rownames_to_column("param") %>% 
    mutate(fit = "m11.10"),
  
  # m11.11
  fixef(m11.11) %>% 
    data.frame() %>% 
    rownames_to_column("param") %>% 
    mutate(fit = "m11.11")
) %>% 
  # wrangle
  filter(param != "Intercept") %>% 
  mutate(param = factor(param,
                        levels = str_c("X", 1:30),
                        labels = str_c("X[", 1:30, "]"))) %>% 
  
  # plot!
  ggplot(aes(x = param, y = Estimate, ymin = Q2.5, ymax = Q97.5, color = fit)) +
  geom_hline(yintercept = 0, linetype = 2, size = 1/4) +
  geom_linerange(position = position_dodge(width = 1/3), size = 1/4, key_glyph = "path") +
  geom_point(position = position_dodge(width = 1/3), size = 1/2) +
  scale_color_manual(NULL,
                     values = c("grey50", "black"),
                     guide = guide_legend(override.aes = list(angle = 90))) +
  scale_x_discrete(NULL, labels = ggplot2:::parse_safe) +
  ylab("marginal posterior") +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 0))


# plot conditional effects 

###_###_###_##
#PLOTS SELECTED MODEL
###_###_###_##


conditional_effects(ecto_p_brms_bayes_no_int_prior)
marginal_effects(selected_ecto_p_brms_bayes_no_int)
plot(conditional_effects(zinb_a_nf_mites_brms_bayes_no_int_prior), 
      points = TRUE, 
      point_args = list(width = .05, shape = 1))

