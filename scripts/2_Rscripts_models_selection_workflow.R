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


ecto_p_brms_bayes_all_interactions<-brms::brm(ectoparasites_PA~
                               sociality+
                               scale(elevation)+
                               scale(year_seasonality)+
                               socality:scale(elevation)+
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

loo_compare(b_ecto_p_brms_bayes,ecto_p_brms_bayes)

m_no_interactions<-loo(b_ecto_p_brms_bayes)
m_interactions<-loo(m_no_interactions,m_interactions)
loo_compare(m_no_interactions,m_interactions)
saveRDS(ecto_p_brms_bayes, "data/data_analyses/models/seasonality/1.model_prevalence_brms_phylo_multiple_obs_031623_seasonality_inter_elev_season.RDS")
saveRDS(b_ecto_p_brms_bayes, "data/data_analyses/models/seasonality/1.model_prevalence_brms_phylo_multiple_obs_031623_seasonality.RDS")
b_ecto_p_brms_bayes<-readRDS("data/data_analyses/models/seasonality/1.model_prevalence_brms_phylo_multiple_obs_031623_seasonality.RDS")

#saveRDS(ecto_p_brms_bayes, "data/data_analyses/models/1.model_prevalence_brms_phylo_multiple_obs_031623.RDS")
#ecto_p_brms_bayes<-readRDS("data/data_analyses/models/1.model_prevalence_brms_phylo_multiple_obs_031623.RDS")
# Summarize the model
summary (ecto_p_brms_bayes)
levels(ectos_birds_dff$sociality)

fixef() # to get more detailed values for estimates
coef(ecto_p_brms_bayes) # if you have group-level effects (hierarchical data)

# interpret the model 
#To test whether all regression coefficients are different from zero, we can look at the Credible Intervals that are listed in the summary output or we can visually represent them in density plots.
#If we do so, we clearly see that zero is not included in any of the density plots, meaning that we can be reasonably certain the regression coefficients are different from zero.
#INTERPRETATION:In the model, the parameter for Sociality means the expected difference between non_social(0) and social (1) with all other covariates held constant. we clearly see that zero is included in the density plot for sociality so there is not effect of sociality??
bayes_R2(ecto_p_brms_bayes) # R2 0.1529
plot(ecto_p_brms_bayes) # paarameter distributio and convergence

mcmc_plot(ecto_p_brms_bayes) # Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
launch_shinystan()
pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict(ecto_p_brms_bayes)
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  + 
  coord_cartesian(xlim = c(0, 10))

ppc_rootogram(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 5), ylim = c(0,30))

ppc_stat(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m, stat = "prop_zero")

#Assumptions check model pglmm_ bayes
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(ecto_p_brms_bayes, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations

# assumptiosn look good!   



# Plots BRMS 
color_scheme_set("blue")

# convergence 
png("data/data_analyses/models/model_plots/1.parameters_distribution_convergence_plot_model_PREVALENCE_brms_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
plot(ecto_p_brms_bayes)
dev.off()

# model fit
png("data/data_analyses/models/model_plots/1.model_fit_PREVALENCE_LICE_brms_phylo_multiple_obs_032123.png.png",width = 3000, height = 3000, res = 300, units = "px")
pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
dev.off()


# Posterior distributions

color_scheme_set("green")
#brmstools::forest(ecto_p_brms_bayes, level = 0.95, show_data = TRUE)

estimates_plot<-mcmc_plot(ecto_p_brms_bayes ,prob=0.90, prob_outer=0.95,
                          variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                          type="areas") +
  labs(title="Posterior distributions ectoparasite prevalence", subtitle ="with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes ,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                                    type="intervals") +
  labs(title="Posterior distribution ectoparasite prevalence", subtitle ="with means and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/models/model_plots/1.parameters_plot_model_PREVALENCE_brms_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()


