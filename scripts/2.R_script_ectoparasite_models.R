#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Models                                                                           ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
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



# [Presence_absence] Part1_Ectoparasite models_  -----------------------------
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
ectoparasites_df<-read.csv("data/7.ectoparasite_df_presence_absence.csv") # data on presence absence
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




# [Probability_Presence_absence] Part1 B ( Suggestions using probabilities)_Ectoparasite models_  -----------------------------
#All models fitted with pglmm() have class of communityPGLMM. Here is a list of functions that can be used to these models.

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
ectoparasites_df<-read.csv("data/7.ectoparasite_df_presence_absence.csv") # data on presence absence
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


# [Abundance]Part2_Ectoparasite models_  -----------------------------

# lice

#Response variable count 
#Random effects species ( we have multiple measurements per species), account to variation between species other than phylogenetic, maybe redundat with foraging high
#Random effect phylogeny
#Random effect elevation 

library(ape)
library(MCMCglmm)

# Read teh files

lice_df_abundance<-read.csv("data/7.lice_df_abundance.csv")
# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
lice_df_abundance<-lice_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
unique (lice_df_abundance$species_jetz) 
lice_df_abundance<-lice_df_abundance %>% distinct( species_jetz,.keep_all = TRUE)

lice_df_abundance<-read.csv("data/5.lice_df_abundance_manu.csv") # data for manu only


phylo_lice_rooted<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")
# Make sure number in teh data and the phylogenty are consistent

lice_df_abundance <-lice_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
lice_df_abundance$elevation_cat<-as.factor(lice_df_abundance$elevation_cat)
lice_df_abundance$foraging_cat<-as.factor(lice_df_abundance$foraging_cat)
lice_df_abundance$species_jetz<-as.factor(lice_df_abundance$species_jetz)
lice_df_abundance$sociality<-as.factor(lice_df_abundance$sociality)


# Using pglmmm
#I am not sue about including foraging cat since that some how is included in teh variation per species , 
#also sample size seems relevant but not sure how to included since the observations are individual.
# iTHINK THE BEST APPROACH IS TO calculate the mean abundance adn nclude sample size as a random effect

 
l_abun_pglmm <-  phyr::pglmm(total_lice ~ sociality+ (1|elevation_cat)+(1|species_jetz__)+ (1|foraging_cat), 
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

# Using MCMCglmm
#We want to investigate a possible relationship between the abundance and the cofactor : sociality

#The variables nitt and burnin are used to calibrate the MCMCM algorithm:it will iterate for burnin iterations before recording samples (to ensure convergence), and then iterate nitt times. T
#he parameter thin helps us to save memory by saving only every `thin' value and thus, dropping highly auto-correlated values2.
#Note that the use of nodes=`TIPS' or nodes=`ALL' in the inverseA function can have a noticeable impact on auto-correlation: whereas the latter would speed up computation, 
#it can results in higher auto-correlation
# Note the phylo column containing the name of the species in our dataset: it corresponds to the phylogenetic effect we are going to include in our model.
#regarding the random effects, we need to define a set of priors for the variance components of  EACH random effect:

#In we have a continous cofactor We have to add a new random effect taking into account the fact that each species has a ``multiple measurement effect''. 
#Note the new column species, which is identical to the phylo one and will be used for that purpose. 
#First, we will try to fit the model , using the specific mean of the cofactor: # I am not sure i have to do thica casuse i dont have cofactors that are numerical
#lice_df_abundance$spec_mean_cf<-sapply(split(lice_df_abundance$cofactor,data$phylo),mean)[data$phylo] # use this to add the mean random effect

names(lice_df_abundance)
inv.phylo<-inverseA(phylo_lice_rooted,nodes="TIPS",scale=TRUE)
prior2<-list(G=list(G1=list(V=1,nu=0.02,alpha.mu=0, alpha.V=10000),G2=list(V=1,nu=0.02,alpha.mu=0, alpha.V=10000),G3=list(V=1,nu=0.02,alpha.mu=0, alpha.V=10000),G4=list(V=1,nu=0.02,alpha.mu=0, alpha.V=10000)),
             R=list(V=1,nu=0.02))


model_l_a_mcmc<-MCMCglmm(total_lice~sociality,random=~species_jetz+elevation_cat+foraging_cat,
                 family="poisson",ginverse=list(species_jetz=inv.phylo$Ainv),
                 prior=prior2,data=lice_df_abundance,nitt=100000,burnin=1000,thin=500)

# Following are the results for the random effect variances (G-structure, containing the variance of the phylo effect) and the residual variance (R-structure, the residual variance is called units in MCMCglmm). 
#We have information about the posterior mean of the estimate, its 95% credible interval4 and its effective sample size.
#The latter is a measure of the auto-correlation within the parameter sample: it should be close to the MCMC sample size above, or failing that, it should be at least large enough (say more than 1,000). 
#The summary of the fixed effects (intercept and cofactor) are similar, except we also have a ``pMCMC'' value for significance testing if the parameter is different from zero5. # so if it overlaps with zero is not significant? 
#If we are strictly Bayesian, we should not do significance testing because such a concept belongs to the frequentists' paradigm. However, we use ``pMCMC'' as if frequentists' p-values for convenience.
summary(model_l_a_mcmc)
plot(model_l_a_mcmc$Sol)
plot(model_l_a_mcmc$VCV)

lambda <- model_l_a_mcmc$VCV[,'species_jetz']/
  (model_l_a_mcmc$VCV[,'species_jetz']+model_l_a_mcmc$VCV[,'elevation_cat']+model_l_a_mcmc$VCV[,'foraging_cat']+model_l_a_mcmc$VCV[,'units'])

mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda) 
 

#The posterior distribution of the fixed term does not overlap zero. 
#Thus the we can infer that sociality has a statistical effect on lice ectos in this model and is a useful addition to the model

posterior.mode(model_l_a_mcmc$Sol[, "sociality1"])
HPDinterval(model_l_a_mcmc$Sol[, "sociality1"])

#Modeling the mea abundance to take care of the sample size

mean_lice<-lice_df_abundance %>% group_by (species_jetz) %>% 
  summarize(mean_lice=mean(total_lice))

species_atributes<-lice_df_abundance %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean) %>%  filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
species_attributes_distict<-distinct( species_atributes)

mean_lice_abundance<-right_join(species_attributes_distict, mean_lice, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 

###_###_####_###_
#write_csv(mean_lice_abundance,"data/5.lice_df_abundance_means.csv")

#The data
mean_lice_abundance<-read.csv("data/5.lice_df_abundance_means.csv")
mean_lice_abundance<-mean_lice_abundance %>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations

phylo_lice_rooted<-read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")


mean_lice_abundance <-mean_lice_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
mean_lice_abundance$elevation_cat<-as.factor(mean_lice_abundance$elevation_cat)
mean_lice_abundance$foraging_cat<-as.factor(mean_lice_abundance$foraging_cat)
mean_lice_abundance$species_jetz<-as.factor(mean_lice_abundance$species_jetz)
mean_lice_abundance$sociality<-as.factor(mean_lice_abundance$sociality)

# The model 

hist(mean_lice_abundance$mean_lice)
unique(mean_lice_abundance$mean_lice)
names ( mean_lice_abundance)
l_abun_mean_pglmm <-  phyr::pglmm(mean_lice ~ sociality+ (1|elevation_cat)+(1|species_jetz__)+ (1|foraging_cat), 
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

mites_df_abundance<-read.csv("data/7.mites_df_abundance.csv")
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

#phylogeny_for_lice<-read.tree("data/phylo_data/1_host_consensus_tree_lice.tre")
# Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?

mites_df_abundance<-mites_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
mites_df_abundance$elevation_cat<-as.factor(mites_df_abundance$elevation_cat)
mites_df_abundance$foraging_cat<-as.factor(mites_df_abundance$foraging_cat)
mites_df_abundance$species_jetz<-as.factor(mites_df_abundance$species_jetz)
mites_df_abundance$sociality<-as.factor(mites_df_abundance$sociality)

str(mites_df_abundance)

mean(mites_df_abundance$total_no_feathers_mites)
sd(mites_df_abundance$total_no_feathers_mites) # overdispersed variance> mean

# Modeling the data # I would prefer to use a zero inflated model however that is only aviallable in a gassioan approach bt that does no work with my model ( not sure why ye)

names( mites_df_abundance)
m_a_no_f<-  phyr::pglmm(total_no_feathers_mites~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                        data = mites_df_abundance, 
                        family ="poisson", # use when bayes=true "zeroinflated.poisson",
                        cov_ranef = list(species_jetz=phylogeny_for_mites), #class phylo
                        #bayes = TRUE,
                        REML = TRUE, 
                        verbose = TRUE, 
                        s2.init = .25)
summary(m_a_no_f)


m_a_meso<-  phyr::pglmm(total_mesostigmatidae~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
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



mean_m_a_meso<-  phyr::pglmm(mean_mites~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
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
