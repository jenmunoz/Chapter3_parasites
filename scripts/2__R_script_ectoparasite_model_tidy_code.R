#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models cleaned                                                                         ###
### R-code                                                                          ###
### Jenny Munoz      
###
### Last update: March 2023                                                ###
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

# [Presence_absence] Part1_Ectoparasite models_using binary data  1/0


# [Prevalence]Part1 Ectoparasite models_using prevalence data (Corrections after meeting with Jullie et al  using probabilities)_Ectoparasite models_  -----------------------------
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

###_###_###
# [Prevalence_Presence_absence] THE MODELS
###_###_###

#DATA
ectos_df<-read.csv("data/data_analyses/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date)
phylo<- read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex") 

unique(ectos_df$species_jetz )

# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
#ectos_df<-ectos_df %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

# Make sure the tips and the names on the file coincide
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 

a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_df$species_jetz))%>% mutate(name=ectos_df$species_jetz) %>% select(name) %>% arrange(desc(name))

tip<-as.list(setdiff(a,b))
print(tip)
#phylo<-drop.tip (phylogeny_prevalence, tip$name) USE IF NEED TO DROP SOME TIPS

# Important!!!!!Make sure the names  are in the same order in the phylogeny and in the traits 
#( I am not sure if this is relevatnt for the PGLMM cause we have multiple observations)
#rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
#ectos_df<- ectos_df[match(phylo$tip.label,rownames(ectos_df)),]

# Re-strudture the data
# Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?

ectos_df <-ectos_df  %>% mutate_at("species_jetz", str_replace, " ", "_")
ectos_df$elevation_cat<-as.factor(ectos_df$elevation_cat)
ectos_df$foraging_cat<-as.factor(ectos_df$foraging_cat)
ectos_df$species_jetz<-as.factor(ectos_df$species_jetz)
ectos_df$elevation<-as.numeric(ectos_df$elevation)
ectos_df$elevation_midpoint<-as.numeric(ectos_df$elevation_midpoint)
ectos_df$sociality<-as.factor(ectos_df$sociality)

is.ultrametric(phylogeny_prevalence)

#MODEL FITTING AND EVALUALLING MODEL FIT GRAPHICALLY

###_###_###_ Important for model fitting
#predict ( model) # plot predicted vlues and confidence limits on the logit scale, the points are not logit tranformed data, they are "working values to fit the model.....dificult to interpret 
#fitted (model) # willl give us the predicted values transformed to the original scale USE THESE TO PREDICT THE DATA
#visreg(model) # Similarly  plot predicted values and confidence limits on the logit scale, dificult to interpret
#visreg(model, scale="response") # willl give us the predicted values transformed to the original scale 
###_###_###_ ###_###_###_ ###_###_###_ 


###_###_###_ 
# Using PGLMM [ Lets correct for phylogeny ]
###_###_###_ 

#I am not sue about including foraging cat since that some how is included in teh variation per species  ( Removed) 
# we revomed (1|foraging_cat) because it was not significant in individual models 

#lm(proportion_ectoparasites~1, data=ectos_df) # even the average is a linear regression INTERESTING

names( ectos_df)
ecto_prevalence_pglmm <-phyr::pglmm(ectoparasites_PA~sociality+elevation+(1|species_jetz__), #+elevation_midpoint+Powder.lvl
                                      data = ectos_df, 
                                      family = "binomial",
                                      cov_ranef = list(species_jetz= phylo), #class phylo
                                      #bayes = TRUE,
                                      REML = TRUE, 
                                      verbose = TRUE,
                                      s2.init = .25) # what is this last parameter for

names(ectos_df)
summary(ecto_prevalence_pglmm)
print(ecto_prevalence_pglmm)
fixef(ecto_prevalence_pglmm)
predict(ecto_prevalence_pglmm)
rr2::R2(ecto_prevalence_pglmm)
class(ecto_prevalence_pglmm ) 

str(ecto_prevalence_pglmm)


ecto_abundance_pglmm <-phyr::pglmm(total_lice~sociality+elevation+(1|species_jetz__), #+elevation_midpoint +Powder.lvl
                                    data = ectos_df, 
                                    family = "poisson",
                                    cov_ranef = list(species_jetz= phylo), #class phylo
                                    #bayes = TRUE,
                                    REML = TRUE, 
                                    verbose = TRUE,
                                    s2.init = .25) # what is this last parameter for

summary(ecto_abundance_pglmm)



ecto_abundance_pglmm <-phyr::pglmm(total_mites~sociality+elevation+(1|species_jetz__)+Powder.lvl, #+elevation_midpoint
                                   data = ectos_df, 
                                   family = "poisson",
                                   cov_ranef = list(species_jetz= phylo), #class phylo
                                   #bayes = TRUE,
                                   REML = TRUE, 
                                   verbose = TRUE,
                                   s2.init = .25) # what is this last parameter for







png("figures/figures_manuscript/Fig1b.Prevalence_ecto_output_predicted.png", width = 3000, height = 3000, res = 300, units = "px")
plot_data(ecto_prevalence_pglmm.var ="species_jetz", site.var ="sociality",predicted=TRUE)
dev.off()

plot(a)

# Underestanding the summary of the random effects
#The random effect with the largest variance and standard variation is the one with the strongest effect, in our case the phylogenetic effect,
# this implies that the parasites  prevalence is more similar in closely related species 

# Explore if the random effects are important? 
# One way to get an idea is to run a likelihood ratio test on the random effect. 
# This can be achieved by using the pglmm_profile_LRT() function, at least for binomial models ( copied from https://daijiang.github.io/phyr/articles/phyr_example_empirical.html)

phyr::pglmm_profile_LRT(ecto_prevalence_pglmm, re.number = 1) ## here we need to specify qith random effect are we evaluating in the order we put them in the model species is number 1, phylogeny2, and elevation 3

LRTest <- sapply(1:3, FUN = function(x) phyr::pglmm_profile_LRT(ecto_prevalence_pglmm, re.number = x))
colnames(LRTest) <- names(ecto_prevalence_pglmm$ss)
t(LRTest)
