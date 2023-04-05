#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models selection workflow                                                                       ###
### R-code                                                                          ###
### Jenny Munoz      
###
### Last update: March 2023                                                ###
################################################################################

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
library(rstan)
library(rstanarm)
library(loo)


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
#saveRDS(ecto_p_brms_bayes_no_int, "data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")
ecto_p_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")

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

#saveRDS(ecto_p_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions.RDS")
ecto_p_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_sociality_interactions.RDS")

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

#saveRDS(ecto_p_brms_bayes_all_interactions, "data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")
ecto_p_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions.RDS")

# PRIOrS SPECIFICATON 

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?


# half student only allows positive values
#prior_summary(ecto_p_brms_bayes_no_int)
ecto_p_brms_bayes_no_int_prior<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
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

#saveRDS(ecto_p_brms_bayes_no_int_prior, "data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_priors.RDS")


ecto_p_brms_bayes_no_int_prior2<-brms::brm(ectoparasites_PA~sociality+ scale(elevation)+ scale(year_seasonality)+
                                             (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                             (1|Powder.lvl)+
                                             (1|species),
                                           data=ectos_birds_dff,
                                           #save_pars = save_pars(all=  TRUE), if i need to use moment match but makes the model heavier
                                           family= bernoulli(), # bernoulli() uses the (link = "logit")
                                           data2 = list(phy_cov=phy_cov),
                                           prior=c(prior_predictors,prior_random),
                                           iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                           thin=2,
                                           control=list(adapt_delta=0.99, max_treedepth=14))

saveRDS(ecto_p_brms_bayes_no_int_prior2, "data/data_analyses/model_selection/P1.model_prevalence_b_brms_phylo_multiple_obs_all_interactions_priors_wo_i.RDS")

loo(ecto_p_brms_bayes_no_int_prior2,ecto_p_brms_bayes_no_int_prior, compare=TRUE)


#MODEL COMPARISON
#paper https://arxiv.org/abs/2008.10296v3 and forum https://discourse.mc-stan.org/t/relationship-between-elpd-differences-and-likelihood-ratios/23855
#we can compute the probability that one model has a better predictive performance than the other.
#The oft-cited rule of thumb is that Bayesian elpd_loo differences less than |4| are small.

bayes_R2(ecto_p_brms_bayes_no_int) 
bayes_R2(ecto_p_brms_bayes_sociality_interactions) 

# use loo cross validation 
loo(ecto_p_brms_bayes_sociality_interactions, ecto_p_brms_bayes_no_int, compare=TRUE)
# eld_diff<0.4 so we keep the simplest model 

# use k-fold-cv validation
k_ecto_p_brms_bayes_no_int<-kfold(ecto_p_brms_bayes_no_int, K=10)
k_ecto_p_brms_bayes_sociality_interactions<-kfold(ecto_p_brms_bayes_sociality_interactions, K=10)

loo_compare(k_ecto_p_brms_bayes_no_int, k_ecto_p_brms_bayes_sociality_interactions)

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

png("data/data_analyses/model_selection/plots_model_selected/Fig2.plot_model_parameters_ecto_PREVALENCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions Ectos prevalence", subtitle ="Ectos prevalence with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig2.plot_model_parameters_intervals_ecto_PREVALENCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



# ##### 1.Data processing abundance lice ----------------------------------------------------

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



#saveRDS(zinb_a_lice_brms_bayes_no_int, "data/data_analyses/model_selection/L1.model_prevalence_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS")
zinb_a_lice_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/l1.model_prevalence_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS")

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
#saveRDS(zip_a_lice_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/L1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")
zinb_a_lice_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/L1.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")

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

#saveRDS(zinb_a_lice_brms_bayes_all_interactions, "data/data_analyses/model_selection/1L.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")
zinb_a_lice_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/1L.model_prevalence_zip_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")

# PRIORS

prior_predictors<-prior("student_t(3,0,10)", class ="b") # Mean of 0 shoudl works, cause our predictors are scaled
prior_random<- prior("student_t(3,0,10)", class="sd",lb=0) # half student allows to only incorporate positive values 
prior_intercept<-prior("student_t(3,0,10)", class="Intercept")  # I am not sure what are good priors for an intercept shoudl I ALSO include negative values?

residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0) # this is the default
residual_prior2<-prior(beta(1,1), class = "zi",lb=0,ub=1) # this is teh default


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
#saveRDS(zinb_a_lice_brms_bayes_no_int_priors, "data/data_analyses/model_selection/1L.model_zip_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")
zinb_a_lice_brms_bayes_no_int_priors<-readRDS( "data/data_analyses/model_selection/1L.model_zip_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors.RDS")


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
#saveRDS(zinb_a_lice_brms_bayes_sociality_interactions_priors, "data/data_analyses/model_selection/1L.model_zip_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_sociality_interactions_priors.RDS")
zinb_a_lice_brms_bayes_sociality_interactions_priors<-readRDS( "data/data_analyses/model_selection/1L.model_zip_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_sociality_interactions_priors.RDS")

#MODEL COMPARISSOM
looni<-loo(zip_a_lice_brms_bayes_no_int)
loonip<-loo(zinb_a_lice_brms_bayes_no_int_priors)

loosi<-loo(zip_a_lice_brms_bayes_sociality_interactions)
#looni1<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, moment_match = TRUE)
#looni2<-loo_moment_match(zinb_a_nf_mites_brms_bayes_no_int_outliers, loo=loo1,k_threshold = 0.7)
looni3<-loo(zinb_a_lice_brms_bayes_no_int, reloo = TRUE) # exclude the observation with high pareto to allow me to compare with other models 
looni4<-reloo(zinb_a_lice_brms_bayes_no_int,loo=looni, chains=4) # ### THIS ONE WORKS!!! abit faster that the one on top # and actually removes all pareto >0.7
loosi4<-reloo(zinb_a_lice_brms_bayes_sociality_interactions, loo=loosi, chains=4)

loo_compare(loon, loosi)

#MODEL COMPARISON

# model check  example

yrepnzb <- posterior_predict(brm_glmznb)
(prop_zero_test4 <- ppc_stat(y=roaches$y, yrepnzb, stat=function(y) mean(y==0)))
(max_test_nb <- pp_check(stan_glmnb, plotfun = "stat", stat = "max"))


bayes_R2(zinb_a_lice_brms_bayes_no_int) 
bayes_R2(zinb_a_lice_brms_bayes_sociality_interactions) 

# use loo cross validation 
loo(zinb_a_lice_brms_bayes_no_int, zinb_a_lice_brms_bayes_sociality_interactions, compare=TRUE)
# eld_diff  so we keep the simplest model 

# use k-fold-cv validation [since While importance sampling in PSIS-LOO can fail for “random effect” model]
k_zinb_a_lice_brms_bayes_no_int<-kfold(ecto_p_brms_bayes_no_int, K=10)
k_zinb_a_lice_brms_bayes_sociality_interactions<-kfold(ecto_p_brms_bayes_sociality_interactions, K=10)

loo_compare(k_zinb_a_lice_brms_bayes_no_int, k_zinb_a_lice_brms_bayes_sociality_interactions)

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

png("data/data_analyses/model_selection/plots_model_selected/Fig2l.plot_model_parameters_LICE_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zip_a_lice_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions LICE_ABUNDANCE", subtitle ="LICE_ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig2l.plot_model_parameters_intervals_LICE_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()


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

# ##### 2.2.Model selection abundance mites --------------------------------------------------------

random_prior<-prior(student_t(3, 0, 10), class = "sd",lb=0)
intercept_prior<-prior(normal( 0, 10), class = "Intercept",lb=0)
residual_prior<-prior(gamma(0.01, 0.01), class = "shape",lb=0)
residual_prior_2<-prior(beta(1,1), class = "zi",lb=0,ub=1)


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

#saveRDS(zip_a_nf_mites_brms_bayes_no_int, "data/data_analyses/model_selection/1NFM.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")
zip_a_nf_mites_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/1NFM.model_prevalence_zip_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")

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
#saveRDS(zip_a_nf_mites_brms_bayes_sociality_interactions, "data/data_analyses/model_selection/1NFM.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")
zip_a_nf_mites_brms_bayes_sociality_interactions<-readRDS( "data/data_analyses/model_selection/1NFM.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_sociality_interactions.RDS")

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
                                                  iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                                  thin=2,
                                                  control=list(adapt_delta=0.99, max_treedepth=14)) 

saveRDS(zip_a_nf_mites_brms_bayes_all_interactions, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")
zip_a_nf_mites_brms_bayes_all_interactions<-readRDS( "data/data_analyses/model_selection/1.model_prevalence_zip_brms_nf_MITES_ABUNDANCE_phylo_multiple_obs_all_interactions.RDS")

prior_summary(zinb_a_nf_mites_brms_bayes_no_int)

loo1<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers)
loo2<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, moment_match = TRUE)
loo3<-loo_moment_match(zinb_a_nf_mites_brms_bayes_no_int_outliers, loo=loo1,k_threshold = 0.7)
loo4<-loo(zinb_a_nf_mites_brms_bayes_no_int_outliers, reloo = TRUE) # exclude the observation with high pareto to allow me to compare with other models 
loo5<-reloo(zinb_a_nf_mites_brms_bayes_no_int_outliers,loo=loo1, chains=1) # ### THIS ONE WORKS abit faster that the one on top # actually removes all pareto >0.7

#Model selection

## use loo cross validation 
loo(zip_a_nf_mites_brms_bayes_sociality_interactions,zip_a_nf_mites_brms_bayes_no_int)
# eld_diff  0.1 so we keep the simplest model Which is the one without interactions

# use k-fold-cv validation [since While importance sampling in PSIS-LOO can fail for “random effect” model]
k_zinb_a_lice_brms_bayes_no_int<-kfold(ecto_p_brms_bayes_no_int, K=10)
k_zinb_a_lice_brms_bayes_sociality_interactions<-kfold(ecto_p_brms_bayes_sociality_interactions, K=10)

loo_compare(k_zinb_a_lice_brms_bayes_no_int, k_zinb_a_lice_brms_bayes_sociality_interactions)

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(zip_a_nf_mites_brms_bayes_no_int,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions nf_MITES_ABUNDANCE", subtitle ="MITES_ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig2m.plot_model_parameters_nf_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zip_a_nf_mites_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions nf_MITES_ABUNDANCE", subtitle ="MITES_ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig2m.plot_model_parameters_intervals_nf_MITES_ABUNDANCE_brms_bayes_no_int.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



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

# model check 
loo_lice_d<-loo(zip_a_lice_brms_bayes_sociality_interactions_degree,save_psis = TRUE)

plot(loo_lice_d)

yrep <- posterior_predict(zip_a_lice_brms_bayes_sociality_interactions_degree)

ppc_loo_pit_overlay(
  y = dff_ectos_network_individual_metrics$total_lice,
  yrep = yrep,
  lw = weights(loo_lice_d$psis_object))

rstan::get_num_upars(zip_a_lice_brms_bayes_no_int_degree$fit)
get_num_upars(zip_a_lice_brms_bayes_no_int_degree)

install.packages("rstan")
library(rstan)
loo(zip_a_lice_brms_bayes_no_int_degree,zip_a_lice_brms_bayes_sociality_interactions_degree, compare=TRUE)
zip_a_lice_brms_bayes_no_int_degree<-brms::brm(total_lice~scale(degree)+ scale(elevation)+ scale(year_seasonality)+
                                          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
                                          (1|Powder.lvl)+
                                          (1|species),
                                        data=dff_ectos_network_individual_metrics,
                                        family=zero_inflated_poisson(),  #zero_inflated_negbinomial()
                                        data2 = list(phy_cov=phy_cov),
                                        iter=8000, warmup=4000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
                                        thin=2,
                                        control=list(adapt_delta=0.99, max_treedepth=12)) 

prior_summary(zip_a_lice_brms_bayes_no_int_degree)
#saveRDS(zinb_a_lice_brms_bayes_no_int_degree, "data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")
zinb_a_lice_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/1.model_prevalence_zip_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")

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

loo(zip_a_lice_brms_bayes_no_int_degree,zip_a_lice_brms_bayes_sociality_interactions_degree)
loo(zip_a_lice_brms_bayes_no_int_degree,zip_a_lice_brms_bayes_sociality_interactions_degree,zip_a_lice_brms_bayes_all_interactions_degree)


loo(zip_a_lice_brms_bayes_no_int_degree, moment_match=TRUE)
loo_moment_match(zip_a_lice_brms_bayes_no_int_degree)

# all interactions
estimates_plot<-mcmc_plot(zip_a_lice_brms_bayes_all_interactions_degree,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig3l.plot_model_parameters_LICE ABUNDANCE_brms_bayes_all_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zip_a_lice_brms_bayes_all_interactions_degree,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE NETWORKS with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig3l.plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_all_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# non interactions 

estimates_plot<-mcmc_plot(zip_a_lice_brms_bayes_sociality_interactions_degree,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig3l.plot_model_parameters_LICE ABUNDANCE_brms_bayes_sociality_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zip_a_lice_brms_bayes_sociality_interactions_degree,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE NETWORKS with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig3l.plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_sociality_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()

# social interactions

# non interactions 

estimates_plot<-mcmc_plot(zip_a_lice_brms_bayes_no_int_degree,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig3l.plot_model_parameters_LICE ABUNDANCE_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

estimates_plot_intervals<-mcmc_plot(zip_a_lice_brms_bayes_no_int_degree,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions LICE ABUNDANCE DEGREE", subtitle ="LICE ABUNDANCE NETWORKS with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/model_selection/plots_model_selected/Fig3l.plot_model_parameters_intervals_LICE ABUNDANCE_brms_bayes_no_int_DEGREE.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



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




