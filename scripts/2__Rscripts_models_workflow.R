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
install.packages('brmstools')
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
install.packages("phyr")
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
library(brms) # bayesian approach to model phylogenetic data with repides observations
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
library(brms)
library(DHARMa.helpers)
library(brmstools)

# # 1.Data import-----------------------------------------------------------------
#DATA
# This dataset contains all parasites samples (887) after removing exact duplicated rows, for which we have an assigned elevation, for 783 in total out of 998 that we had originally  (this included some duplicates)
ectos_df<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date)

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex") 

unique(ectos_df$species_jetz )

# ## 2.Data overview -------------------------------------------------------
str( ectos_df)
str(phylo)


# #### 3.Descriptive statistics --------------------------------------------------


# #### 4.Descriptive statistiscs plots -------------------------------------

# ##### 5.Data processing prevalence ----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality ) %>% 
  na.omit() 

#ectos_birds_dff <-ectos_birds_dff  %>% mutate_at("species_jetz", str_replace, " ", "_")
#ectos_birds_dff$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
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

# ##### 5.1.Analyses models prevalence --------------------------------------------------------
###_###_###
  #a) model glmm
ecto_p_glmm <-glmer(ectoparasites_PA~sociality+scale(elevation)+(1|Powder.lvl), #+elevation_midpoint+Powder.lvl
                             data = ectos_birds_dff, 
                             family = "binomial")
###_###_###
  #b) model pglmm
ecto_p_pglmm <-phyr::pglmm(ectoparasites_PA~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl), #+elevation_midpoint+Powder.lvl +(1|Powder.lvl)
                                    data = ectos_birds_dff, # exlude raws with nas this is usefull to check the model assumptions 
                                    family = "binomial",
                                    cov_ranef = list(species_jetz= phylo), #class phylo
                                    #bayes = TRUE,
                                    REML = TRUE, 
                                    verbose = TRUE,
                                    s2.init = .25) # what is this last parameter for
###_###_###
    #c) model pglmm summary
names(ectos_birds_dff)
class(ecto_p_pglmm ) 
summary(ecto_p_pglmm )
print(ecto_p_pglmm )
fixef(ecto_p_pglmm )
predict(ecto_p_pglmm )
rr2::R2(ecto_p_pglmm ) # goodness of fit

    #model pglmm assumptions check DHARMa
    #Overdisperssion
simulationOutput <-DHARMa::simulateResiduals(fittedModel =ecto_p_pglmm, plot = F)
plot(simulationOutput)
plotResiduals(simulationOutput , ecto_p_pglmm$data$sociality)
plotResiduals(simulationOutput , form =ecto_p_pglmm$sociality)

    #Homogenity of variance within groups (Heteroscedasticity) 
testCategorical(simulationOutput, catPred = ectos_birds_dff$sociality)

    # Other test
testUniformity(simulationOutput) #tests if the overall distribution conforms to expectations # 
testOutliers(simulationOutput)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput, catPred = ectos_df_wo_na$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput) ## tests if there are more zeros in the data than expected from the simulations

###_###_###
    #d) model pglmm bayes : hierarchical Bayesian models fitted using integrated nested laplace approximation (INLA)
###_###_###
ecto_p_pglmm_bayes <- phyr::pglmm(ectoparasites_PA~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                     data =ectos_birds_dff, 
                                     cov_ranef = list(species_jetz= phylo), #class phylo can take a phylohgeny or a covariance matrix
                                     family = "binomial", # zeroinflated.binomial in reality I want to use a zeroinflated.binomial but for some reason i can not plot those
                                     bayes = TRUE,
                                     prior = "pc.prior.auto")

summary(ecto_p_pglmm_bayes)
print(ecto_p_pglmm_bayes)
rr2::R2(ecto_p_pglmm_bayes) 
class(ecto_p_pglmm_bayes) 
# plot 
plot_bayes(ecto_p_pglmm_bayes)

# assumptions check model pglmm_ bayes DHARMa ( DHARMa is not working with pgllm Bayes = True)
res_bayes<-DHARMa::simulateResiduals(ecto_prevalence_pglmm_bayes)
class(ecto_prevalence_pglmm_bayes)
ecto_prevalence_pglmm_bayes

## Model Individual  sample with phylognetic + non-phylo effects 
#d) model BRMS bayes
ecto_p_brms_bayes<-brms::brm(ectoparasites_PA~0+sociality+scale(elevation)+
          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
            (1|Powder.lvl)+
            (1|species),
        data=ectos_birds_dff,
        family= bernoulli(), # bernoulli() #zero_inflated_negbinomial()
        data2 = list(phy_cov=phy_cov),
        iter=1000, warmup=500, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
        thin=2,
        control=list(adapt_delta=0.99, max_treedepth=12)) 
saveRDS(ecto_p_brms_bayes, "data/data_analyses/models/model_prevalence_brms_phylo_multiple_obs.RDS")

# Summarize the model
summary (ecto_p_brms_bayes)
levels(ectos_birds_dff$sociality)

fixef() # to get more detailed values for estimates
coef(ecto_p_brms_bayes2) # if you have group-level effects (hierarchical data)

# interpret the model 
#To test whether all regression coefficients are different from zero, we can look at the Credible Intervals that are listed in the summary output or we can visually represent them in density plots.
#If we do so, we clearly see that zero is not included in any of the density plots, meaning that we can be reasonably certain the regression coefficients are different from zero.
#INTERPRETATION:In the model, the parameter for Sociality means the expected difference between non_social(0) and social (1) with all other covariates held constant. we clearly see that zero is included in the density plot for sociality so there is not effect of sociality??
bayes_R2(ecto_p_brms_bayes) 
plot(ecto_p_brms_bayes)
mcmc_plot(ecto_p_brms_bayes) # Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
launch_shinystan()
pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict(ecto_p_brms_bayes)
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=ecto_p_brms_bayes$data$ectoparasites_PA, pp_m[1:200, ])  + 
  coord_cartesian(xlim = c(0, 5))

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
  
# Plots BRMS 
  
color_scheme_set("red")
brmstools::forest(ecto_p_brms_bayes, level = 0.95, show_data = TRUE)

mcmc_plot(ecto_p_brms_bayes ,prob=0.90, prob_outer=0.95) +
  ggtitle("Posterior intervals")+
  theme_minimal(15)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

posterior <- as.array(ecto_p_brms_bayes)
mcmc_areas(posterior  ,prob=0.90, prob_outer=0.95)

names(ecto_p_brms_bayes$formula$resp)

# ####### 6.Data processing abundance lice----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, total_lice) %>% 
  na.omit() %>% 
  filter(total_lice<60)  # removing outliers 

# Finding putliers

#ectos_birds_dff <-ectos_birds_dff  %>% mutate_at("species_jetz", str_replace, " ", "_")
#ectos_birds_dff$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$total_lice<-as.factor(ectos_birds_dff$total_lice)

names(ectos_birds_dff)
is.ultrametric(phylo)

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_dff$species_jetz)) %>% mutate(name=ectos_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylogeny_prevalence, tip$name) 

# ####### 6.Data processing abundance mites----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, total_mites, total_mesostigmatidae, total_no_feathers_mites) %>% 
  na.omit() %>% 
  filter(total_mites<300)  # removing outliers for total mites
#filter(total_mesostigmatidae<100)  # removing outliers for total mites
#filter(total_no_feathers_mites<100)  # removing outliers for total mites

# Outliers
assertr::insist(ectos_birds_dff, within_n_mads(50), total_mites)
assertr::insist(ectos_birds_dff, within_n_mads(2), total_no_feathers_mites)
View(ectos_birds_dff)

#ectos_birds_dff <-ectos_birds_dff  %>% mutate_at("species_jetz", str_replace, " ", "_")
#ectos_birds_dff$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$total_lice<-as.factor(ectos_birds_dff$total_lice)
ectos_birds_dff$total_mites<-as.factor(ectos_birds_dff$total_mites)
ectos_birds_dff$total_mesostigmatidae<-as.factor(ectos_birds_dff$total_mesostigmatidae)
ectos_birds_dff$total_no_feathers_mites<-as.factor(ectos_birds_dff$total_no_feathers_mites)

names(ectos_birds_dff)
is.ultrametric(phylo)

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_dff$species_jetz)) %>% mutate(name=ectos_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylogeny_prevalence, tip$name) 


# ####### 6.1.Analyses models abundance--------------------------------------------------------
# Lice



###_###_###
#d) model pglmm bayes
###_###_###




# Summarize the model
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)

bayes_R2() 
plot()
mcmc_plot()
launch_shinystan()
pp_check(, ndraws = 100)+ xlim(0, 40)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict()
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=m$data$total_lice, pp_m[1:200, ]) + log1 + 
  coord_cartesian(xlim = c(0, 50))

ppc_rootogram(y=m$data$total_lice, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 50), ylim = c(0,30))

ppc_stat(y=m4$data$lice_total, pp_m1, stat = "prop_zero")

#Assumptions check model pglmm_ bayes
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(m, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations


# Mites


###_###_###
#d) model pglmm bayes
###_###_###

# Summarize the model
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)

bayes_R2() 
plot()
mcmc_plot()
launch_shinystan()
pp_check(, ndraws = 100)+ xlim(0, 40)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict()
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=m$data$total_lice, pp_m[1:200, ]) + log1 + 
  coord_cartesian(xlim = c(0, 50))

ppc_rootogram(y=m$data$total_lice, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 50), ylim = c(0,30))

ppc_stat(y=m4$data$lice_total, pp_m1, stat = "prop_zero")

#Assumptions check model pglmm_ bayes
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(m, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations



# ########7.Data processing diversity [In progress]--------------------------------------------------------
ectos_birds_dff_d<-read.csv("data/data_analyses/data_manuscript/", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, total_mites, total_mesostigmatidae, total_no_feathers_mites) %>% 
  na.omit() %>% 
  filter(total_mites<300) =
  
assertr::insist(ectos_birds_dff, within_n_mads(50), total_mites)
assertr::insist(ectos_birds_dff, within_n_mads(2), total_no_feathers_mites)


View(ectos_birds_dff)

#ectos_birds_dff <-ectos_birds_dff  %>% mutate_at("species_jetz", str_replace, " ", "_")
#ectos_birds_dff$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
ectos_birds_dff$total_lice<-as.factor(ectos_birds_dff$total_lice)
ectos_birds_dff$total_mites<-as.factor(ectos_birds_dff$total_mites)
ectos_birds_dff$total_mesostigmatidae<-as.factor(ectos_birds_dff$total_mesostigmatidae)
ectos_birds_dff$total_no_feathers_mites<-as.factor(ectos_birds_dff$total_no_feathers_mites)

names(ectos_birds_dff)
is.ultrametric(phylo)

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_dff$species_jetz)) %>% mutate(name=ectos_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylogeny_prevalence, tip$name)

# ####### 7.1 Analyses models diversity [In progress]--------------------------------------------------------


###_###_###
#d) model pglmm bayes
###_###_###

# Summarize the model
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)

bayes_R2() 
plot()
mcmc_plot()
launch_shinystan()
pp_check(, ndraws = 100)+ xlim(0, 40)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict()
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=m$data$total_lice, pp_m[1:200, ]) + log1 + 
  coord_cartesian(xlim = c(0, 50))

ppc_rootogram(y=m$data$total_lice, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 50), ylim = c(0,30))

ppc_stat(y=m4$data$lice_total, pp_m1, stat = "prop_zero")

#Assumptions check model pglmm_ bayes
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(m, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations


# 0.N Libraries networks --------------------------------------------------------------

# # 1.N Data import-----------------------------------------------------------------

# ####### 8.Data processing networks--------------------------------------------------------
# Lice

# Mites

# ####### 8.1. Data analyses networks--------------------------------------------------------


# ####### 8.1 Analyses models NETWORKS [In progress]--------------------------------------------------------


###_###_###
#d) model pglmm bayes
###_###_###

# Summarize the model
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)

bayes_R2() 
plot()
mcmc_plot()
launch_shinystan()
pp_check(, ndraws = 100)+ xlim(0, 40)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict()
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=m$data$total_lice, pp_m[1:200, ]) + log1 + 
  coord_cartesian(xlim = c(0, 50))

ppc_rootogram(y=m$data$total_lice, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 50), ylim = c(0,30))

ppc_stat(y=m4$data$lice_total, pp_m1, stat = "prop_zero")

#Assumptions check model pglmm_ bayes
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(m, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations


