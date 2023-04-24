#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models selection workflow                                                                       ###
### R-code                                                                          ###
### Jenny Munoz      
###
### Last update: March 2023                                                ###
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


# ##### 1.Data processing prevalence ectos ----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality ) %>% 
  na.omit() %>% 
  filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

ectos_birds_dff <- get.complete.cases(ectos_birds_dff) # mke sure we get all complete cases 

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

# THIS IS INCLUDING THE WEAKLY INFORMATIVE SPECIFIED FOR THE INTERCEPT BUT SOME HOW IT DOES NOT LIKE IT

#ecto_p_brms_bayes_no_int_prior2<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
                                             (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                             (1|Powder.lvl)+
                                             (1|species),
                                           data=ectos_birds_dff,
                                           #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                           family= bernoulli(), # bernoulli() uses the (link = "logit")
                                           data2 = list(phy_cov=phy_cov),
                                           prior=c(prior_predictors,prior_random,prior_intercept),
                                           iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                           thin=2,
                                           control=list(adapt_delta=0.99, max_treedepth=14))

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
loo(ecto_p_brms_bayes_no_int_prior, ecto_p_brms_bayes_sociality_interactions_priors,compare=TRUE)
loo(ecto_p_brms_bayes_no_int_prior, ecto_p_brms_bayes_all_interactions_priors,compare=TRUE)
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
#PLOTS
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
png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int_prior)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100) 
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

color_scheme_set("blue")

estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions Ectos prevalence", subtitle ="Ectos prevalence with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_ESTIMATES.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions Ectos prevalence", subtitle ="Ectos prevalence with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1P.ecto_p_brms_bayes_no_int_prior_ESTIMATES_INTERVALS.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


# 1.Data processing and  model selection for prevalence species level (zero_one_beta) --------
ectos_birds_df<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality ) %>% 
  na.omit() %>% 
  filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites

ectos_pres_abs<-ectos_birds_df %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation=mean(elevation), sociality=max(sociality)) %>% 
  mutate(proportion_ectoparasites=ectoparasites_presence/sample_size) %>% 
  na.omit() 

ectos_birds_dff<-ectos_pres_abs %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE)%>% 
  filter(sample_size>4)     # Eliminates species taht ocurr at two elevatiosn and keep only one ( in alphabetical order of the stations, e.g for galbula low_elevations is kept and montane is deleted, but the samples size is mantained this is jut to be able to do the phylogenetic analyses properly)


#ectos_birds_dff <- get.complete.cases(ectos_birds_dff) # mke sure we get all complete cases 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # Need to prunne the tree 


# Data structure
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
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


# The models 
ecto_p_brms_bayes_no_int_species<-brms::brm(proportion_ectoparasites~sociality+ scale(elevation)+scale(sample_size)+
                                      (1|gr(species_jetz, cov = phy_cov))+
                                      (1|species),
                                    data=ectos_birds_dff,
                                    #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                    family= zero_one_inflated_beta(), # for proportions
                                    data2 = list(phy_cov=phy_cov),
                                    iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                    thin=2,
                                    control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(ecto_p_brms_bayes_no_int_species, "data/data_analyses/model_selection/P1s.model_prevalence_brms_phylo_SPECIES_no_interactions.RDS")
ecto_p_brms_bayes_no_int_species<-readRDS("data/data_analyses/model_selection/P1s.model_prevalence_brms_phylo_SPECIES_no_interactions.RDS")

prior_summary(ecto_p_brms_bayes_no_int_species)
bayes_R2(ecto_p_brms_bayes_no_int_species)
mcmc_plot(ecto_p_brms_bayes_no_int_species_priors)


#PRIORS
prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl work, cause our predictors are scaled could use a wider prior if needed normal(0,10)
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


ecto_p_brms_bayes_no_int_species_priors<-brms::brm(proportion_ectoparasites~sociality+ scale(elevation)+
                                              (1|gr(species_jetz, cov = phy_cov))+
                                              (1|species),
                                            data=ectos_birds_dff,
                                            #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                            family= zero_one_inflated_beta(), # # this family is for proportions that include zeros and ones                                         
                                            data2 = list(phy_cov=phy_cov),
                                            save_pars = save_pars(all = TRUE),
                                            prior = c(prior_predictors,prior_random,prior_intercept),
                                            iter=12000, warmup=5000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                            thin=2,
                                            control=list(adapt_delta=0.99, max_treedepth=14)) 


saveRDS(ecto_p_brms_bayes_no_int_species_priors, "data/data_analyses/model_selection/P1s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors.RDS")
ecto_p_brms_bayes_no_int_species_priors<-readRDS("data/data_analyses/model_selection/P1s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors.RDS")

# model comparison

View(ectos_birds_dff)

loo(ecto_p_brms_bayes_no_int_species_priors, moment_match = TRUE)

#TONS OF PAREKO >0.7 MAYBE TRY A DIFFERENT DISTRIBUTION !!!

#PLOTS
plot(conditional_effects(ecto_p_brms_bayes_no_int_species_priors, dpar="mu"), 
      points = TRUE, 
      point_args = list(width = .05, shape = 1))

color_scheme_set("orange") 

bayes_R2(ecto_p_brms_bayes_no_int_species_priors)

#model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig1Ps.plot_model_CONVERGENCE_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int_species_priors)
dev.off()

# model fit
# model fit
png("figures/figures_manuscript/models_selected_figures/Fig1Ps.plot_model_FIT_intervals_ecto_PREVALENCE_brms_bayes_no_int_SPECIES.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_species_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 10)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions PREVALENCE SPECIES", subtitle ="PREVALENCE SPECIES with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1Ps..plot_model_parameters_PREVALENCE SPECIES_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions PREVALENCE SPECIES", subtitle ="PREVALENCE SPECIES with medians and 95% intervals")+
  theme_classic(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig1Ps..plot_model_parameters_intervals_PREVALENCE SPECIES_brms_bayes_no_int_degree.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



# ##### 2.Data processing abundance lice ----------------------------------------------------

ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, total_lice,year_seasonality ) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens") 
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
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
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
                                                iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                thin=2,
                                                control=list(adapt_delta=0.99, max_treedepth=14))

#saveRDS(zinb_a_lice_brms_bayes_no_int_priors, "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")
zinb_a_lice_brms_bayes_no_int_priors<-readRDS( "data/data_analyses/model_selection/M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")

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

k_ZIP_a_lice_no_int_prior<-kfold(ZIP_a_lice_brms_bayes_no_int_priors, K=10)
saveRDS(k_ZIP_a_lice_no_int_prior, "data/data_analyses/model_selection/k_fold/K_fold_1.ZIP_model_ABUNDANCE_LICE_brms_multiple_obs_all_interactions_priors_poisson.RDS")

k_zinb_a_lice_no_int<-kfold(zinb_a_lice_brms_bayes_no_int_priors, K=10)
saveRDS(k_zinb_a_lice_no_int, "data/data_analyses/model_selection/k_fold/K_fold_1_zinb_model_ABUNDANCE_LICE_brms_multiple_obs_all_interactions_priors_poisson.RDS")

loo_compare(k_ZIP_a_lice_no_int_prior, k_zinb_a_lice_no_int) # compare using elpd_diff


#looni1<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, moment_match = TRUE)
#looni2<-loo_moment_match(zinb_a_nf_mites_brms_bayes_no_int_outliers, loo=loo1,k_threshold = 0.7)
looni3<-loo(zinb_a_lice_brms_bayes_no_int, reloo = TRUE) # exclude the observation with high pareto to allow me to compare with other models 
looni4<-reloo(zinb_a_lice_brms_bayes_no_int,loo=looni, chains=4) # ### THIS ONE WORKS!!! abit faster that the one on top # and actually removes all pareto >0.7
loosi4<-reloo(zinb_a_lice_brms_bayes_sociality_interactions, loo=loosi, chains=4)

###_###_###_##
#PLOTS
###_###_###_##

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

color_scheme_set("blue")

# model convergence 
png("figures/figures_manuscript/models_selected_figures/ZIP_Fig3L.zinb_a_lice_brms_bayes_no_int_priors_CONVERGENCE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ZIP_a_lice_brms_bayes_no_int_priors)
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/ZIP_Fig3L.zinb_a_lice_brms_bayes_no_int_priors_FIT.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ZIP_a_lice_brms_bayes_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
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

estimates_plot<-mcmc_plot(ZIP_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions ZIP LICE ABUNDANCE", subtitle ="ZIP LICE ABUNDANCE with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/ZIP_Fig3L.LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ZIP_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions ZIP LICE ABUNDANCE", subtitle ="ZIP LICE ABUNDANCE with medians and 95% intervals")+
  theme_classic(30)+
  xlim(-3,3)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/ZIP_Fig3L.LICE_ABUNDANCE_brms_bayes_no_int_prior_ESTIMATES_INTERVALS.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# model check  example

yrepnzb <- posterior_predict(brm_glmznb)
(prop_zero_test4 <- ppc_stat(y=roaches$y, yrepnzb, stat=function(y) mean(y==0)))
(max_test_nb <- pp_check(stan_glmnb, plotfun = "stat", stat = "max"))


# ##### 3.1.Data processing abundance mites ----------------------------------------------------

ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality,total_mites, total_mesostigmatidae, total_no_feathers_mites,year_seasonality ) %>% 
  na.omit() %>% filter(species_jetz!="Premnoplex_brunnescens") 
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


zinb_a_nf_mites_brms_bayes_no_int_prior_hurdle<-brms::brm(total_no_feathers_mites~sociality+ scale(elevation)+ scale(year_seasonality)+
                                                     (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                     (1|Powder.lvl)+
                                                     (1|species),
                                                   data=ectos_birds_dff,
                                                   family=hurdle_poisson(),  #zero_inflated_negbinomial()
                                                   data2 = list(phy_cov=phy_cov),
                                                   iter=1000, warmup=500, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                   thin=2,
                                                   prior = c(prior_predictors,prior_random,prior_intercept),
                                                   #save_pars = save_pars(all=  TRUE),
                                                   control=list(adapt_delta=0.99, max_treedepth=14)) 

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
k_ecto_p_brms_no_int<-kfold(ecto_p_brms_bayes_no_int_prior, K=10)
k_ecto_p_brms_sociality_int<-kfold(ecto_p_brms_bayes_sociality_interactions_priors, K=10)
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
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
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

# ##### 3.3 Selected model Abundance MITES (ZIP and ZINB) ------------------------------------
#Priors

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default

prior_summary(ZIP_a_nf_mites_brms_bayes_no_int_prior)

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

#saveRDS(zinb_a_nf_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")
zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")

# 28 DIVERGENT TRANSICTIONS !!! WARNING 
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
k_mites_brms_bayes_no_int_prior<-kfold(zinb_a_nf_mites_brms_bayes_no_int_prior, K=10)
saveRDS(k_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/k_fold/K_fold_M1MNF_model_MITES_ABUNDANCE_no_interactions_priors_nb.RDS")
k_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/k_fold/K_fold_M1MNF_model_MITES_ABUNDANCE_no_interactions_priors_nb.RDS")

k_mites_brms_bayes_no_int_prior_p<-kfold(ZIP_a_nf_mites_brms_bayes_no_int_prior, K=10)
saveRDS(k_mites_brms_bayes_no_int_prior, "data/data_analyses/model_selection/k_fold/K_fold_M1MNF_model_MITES_ABUNDANCE_no_interactions_priors_poisson.RDS")
k_mites_brms_bayes_no_int_prior_p<-readRDS()

loo_compare(k_ecto_p_brms_no_int, k_ecto_p_brms_sociality_int)

# plots 
color_scheme_set("red") 

# poisson 
#model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig2LNDNZ_BEST_plot_model_CONVERGENCE_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_BEST.png",width = 3000, height = 3000, res = 300, units = "px")
plot()
dev.off()

# model fit
png("figures/figures_manuscript/models_selected_figures/Fig2LNDNZ_BEST_plot_modell_FIT__ecto_LICE_ABUNDANCE_brms_bayes_no_int_DEGREE_BEST.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(MODEL, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(ZIP_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE DEGREE with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2LND_BEST_plot_model_parameters_LICE ABUNDANCE_brms_bayes_social_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(MODEL,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE DEGREE with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2LND_BEST_plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_social_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
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



# ##### 4.Data processing NETWORKS prevalence ectos ----------------------------------------------------

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


# ##### 4.1.Model selection NETWORKS prevalence ectos --------------------------------------------------------
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
# THE BESTY 
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

#PLOTS
color_scheme_set("purple")

summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)
bayes_R2() # R2 0.1529

# model convergence 
png("figures/figures_manuscript/models_selected_figures/Fig2PND.plot_model_CONVERGENCE_intervals_ecto_PREVALENCE_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
# model fit
png("figures/figures_manuscript/models_selected_figures/Fig2PND.plot_model_FIT_intervals_ecto_PREVALENCE_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 10)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions PREVALENCE DEGREE", subtitle ="PREVALENCE DEGREE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2PND.plot_model_parameters_PREVALENCE_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="dens_overlay") +
  labs(title="Posterior distributions PREVALENCE DEGREE", subtitle ="PREVALENCE DEGREE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/Fig2PND.plot_model_parameters_intervals_PREVALENCE_brms_bayes_no_int_degree.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

plot.O <- stanplot(brms.full)
my_O_title <- expression(paste(bold(" Overall Infected"),bold(" posterior estimates")))
plot.O + ggtitle(my_O_title) + theme(plot.title = element_text(family = "Arial", color="black",size=16)) + theme(text = element_text(family="Arial",size=14))

#### strenght
names(dff_ectos_network_individual_metrics)
zinb_prevalence_brms_bayes_no_int_W_degree_prior<-brms::brm(ectoparasites_PA~scale(w_degree)+ scale(elevation)+ scale(year_seasonality)+
                                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                                          (1|Powder.lvl)+
                                                          (1|species),
                                                        data=dff_ectos_network_individual_metrics,
                                                        family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                                                        data2 = list(phy_cov=phy_cov),
                                                        iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                        thin=2,
                                                        save_pars = save_pars(all=  TRUE),
                                                        prior = c(prior_predictors,prior_random, prior_intercept),
                                                        control=list(adapt_delta=0.99, max_treedepth=14)) 




mcmc_plot(zinb_prevalence_brms_bayes_no_int_W_degree_prior)
msms_zinb_prevalence_brms_bayes_no_int_W_degree_prior()
# # ##### 4.1.1 Model selection NETWORKS prevalence ectos species  --------

dff_ectos_network_individual_metrics<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
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

dff_ectos_network_individual_metrics<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
  select(elevation_extrapolated_date, species_jetz, Powder.lvl,foraging_cat, sociality,total_lice, degree, w_degree, year_seasonality) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
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

dff_ectos_network_individual_metrics<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
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

# ##### 5.5. Selected model NETWORKS LICE ABUNDANCE degree (ZIP VS ZINB) ------------------------------------------

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
  labs(title="Posterior distributions ZINBLICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE DEGREE with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LND_BEST_plot_model_parameters_LICE ABUNDANCE_brms_bayes_social_int_DEGREE_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions ZINB LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE DEGREE with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures/best_models/Fig2LND_BEST_plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_social_int_DEGREE_ZINB.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


# ##### 6.1.Data processing NETWORKS abundance mites ----------------------------------------------------

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

estimates_plot_intervals<-mcmc_plot(ZIP_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
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

#saveRDS(zinb_a_nf_mites_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")
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

#saveRDS(ZIP_a_nf_mites_brms_bayes_no_int_degree_prior, "data/data_analyses/model_selection/M2MND_model_MITES_ABUNDANCE_zip_brms_phylo_multiple_obs_no_interactions_degree_prior.RDS")
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
  labs(title="Posterior distributions ZINB MITES ABUNDANCE DEGREE", subtitle ="MITES ABUNDANCE DEGREE with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures//best_models/Fig2MND_BEST_zinb_ABUNDANCE_MITES_brms_bayes_ESTIMATES_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions ZINB MITES ABUNDANCE DEGREE", subtitle ="MITES ABUNDANCE NETWORKS with medians and 95% intervals")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

png("figures/figures_manuscript/models_selected_figures//best_models/Fig2MND_BEST_zinb_MITES_ABUNDANCE_brms_bayes_INTERVALS_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# #### ### #### ## 6.4 ** Model abundance excluding zeros Mites -------------------------------

dff_ectos_network_individual_metrics<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_network_metrics_individuals_FILE.csv",na.strings =c("","NA"))%>% 
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


# Checking model overfitting and missbehavior and model comparison  ----------------------------------------------

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

