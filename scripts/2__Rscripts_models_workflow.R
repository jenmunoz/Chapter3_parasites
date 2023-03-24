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

# # 1.Data import-----------------------------------------------------------------
#DATA
# This dataset contains all parasites samples (887) after removing exact duplicated rows, for which we have an assigned elevation, for 783 in total out of 998 that we had originally  (this included some duplicates)
ectos_df<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date)

unique(ectos_df$species_jetz) # this is teh total species that we have samples for 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

unique(ectos_df$species_jetz )

# phylogenetic correlation structure, create a covariance matrix of species
phy_cov<-ape::vcv(phylo, corr=TRUE)


# ## 2.Data overview -------------------------------------------------------
str( ectos_df)
str(phylo)


# #### 3.Descriptive statistics --------------------------------------------------


# #### 4.Descriptive statistiscs plots -------------------------------------

ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, total_lice,total_no_feathers_mites,total_mesostigmatidae ) 

ectos_birds_dff_mean <- ectos_birds_dff %>% select(species_jetz,total_lice,total_no_feathers_mites,total_mesostigmatidae) %>% filter(complete.cases(.)) %>% group_by(species_jetz) %>% mutate(mean_lice = mean(total_lice), mean_nf_mites=mean(total_no_feathers_mites)) %>% droplevels()

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

### double check that the phylo and data match 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_birds_dff$species_jetz)) %>% mutate(name=ectos_birds_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)
tip<-as.list(setdiff(a,b))
print(tip)
# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name)

order<-(as.data.frame(phylo$tip.label))
# Important!!!!!Make sure the names  are in the same order in the phylogeny and in the traits 

rownames(ectos_birds_dff) <- ectos_birds_dff$species_jetz # first make it the row names 
ectos_birds_dff<- ectos_birds_dff[match(phylo$tip.label,rownames(ectos_birds_dff)),]


# releveling
#plot_data <- data %>% select(genus_species, mite_load, ode_mass_g) %>% filter(complete.cases(.)) %>% group_by(genus_species) %>% mutate(avg_mass = mean(ode_mass_g)) %>% droplevels()
#plot_dat$genus_species <- as.factor(plot_dat$genus_species)
#plot_dat <- as.data.frame(plot_dat) %>% mutate(genus_species = fct_relevel(genus_species, plot_tree$tip.label[ordered_tips]))

# it will be nice to include the mean in this figures 
lice_load_plot <- ggplot(ectos_birds_dff, aes(x = total_lice, y =species_jetz)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks=c(0, 1, 5, 10,15, 20, 30,40, 50)) +
  #scale_x_continuous(trans="log1p",breaks=c(0, 1, 2,5, 10, 25, 50, 100, 200, 300)) +
  scale_y_discrete(limits=order$`phylo$tip.label`) +
  ylab("") + xlab("Lice per individual") +
  theme_ridges(center_axis_labels = TRUE) +
  theme() +
geom_jitter(alpha=0.5, col="darkolivegreen4")


lice_load_plot[["data"]][["species_jetz"]]

# The phylogenetic plot with prevalence

#tree_plot <- ggtree(phylo, ladderize=FALSE) + geom_tiplab() + ggplot2::xlim(0, 450)

ColorPalette <- brewer.pal(n = 9, name = "GnBu")
list.names=setNames(ectos_birds_dff$proportion_ectoparasites, ectos_birds_dff$species_jetz)
fmode<-as.factor(setNames(ectos_birds_dff$sociality,ectos_birds_dff$species_jetz))
object = contMap(phylo, list.names, direction = "leftwards", plot=FALSE)
#object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))
object_color<-setMap(object, ColorPalette)


tree_plot<-dotTree(phylo,fmode,colors=setNames(c("red","black"), c("1","0")),ftype="i",fsize=0.5, lwd=4) 
object = contMap(phylo, list.names, direction = "leftwards", plot=FALSE)

object_color<-setMap(object, ColorPalette)

list.names=c(unique(ectos_birds_dff$species_jetz))
scale_y_discrete(limits=list.names)+
  

mites_load_plot <- ggplot(data = ectos_birds_dff, aes(y=species_jetz, x=total_no_feathers_mites)) + 
  geom_jitter(alpha=0.5, col="coral3") + #geom_boxplot(outlier.alpha=0) +
  scale_x_continuous(breaks=c(0, 1, 5, 10, 20, 30,40, 50, 100, 200, 300)) +
  scale_y_discrete(limits=order$`phylo$tip.label`) +
  theme_ridges(center_axis_labels = TRUE) + 
  ylab("") +
  xlab("Mites (non-feather) per individual") +
  theme()
# to eliminate the species names in the axes once we know the order is correct use
theme(
axis.text.y=element_blank(),
axis.title.y=element_blank())

mites_plot[["data"]][["species_jetz"]]

phylo_mite_lice_plot <- lice_load_plot %>% insert_left(tree_plot) %>% insert_right(mites_load_plot) 

phylo_mite_lice_plot <- lice_load_plot %>% insert_left(mites_load_plot) 


ggsave("figures/figures_manuscript/f_Fig0_phylo_mite_lice_plot_descriptive.pdf", plot=phylo_mite_lice_plot, height=10, width=12, units="in")



# Plotting the traits things in the  phylogeny 
#Some examples fr plotting http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html
# Plot phylogenetic tree and  the trait continous trait ( prevalence)
#ContMap #Function plots a tree with a mapped continuous character. 
#The mapping is accomplished by estimating states at internal nodes using ML with fastAnc, and then interpolating the states along each edge using equation [2] of Felsenstein (1985).

#contMap(tree, x, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL,
#lims=NULL, outline=TRUE, sig=3, type="phylogram", direction="rightwards", 
#plot=TRUE, ...)


ectos_pres_abs<-ectos_birds_dff %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation_midpoint=max(elevation_midpoint)) %>% 
  mutate(proportion_ectoparasites=ectoparasites_presence/sample_size)

species_atributes<-ectoparasites_df %>% select( species_jetz, species_clean,species_binomial,BLFamilyLatin, sociality,mass_tidy,elevation_cat, foraging_cat )
species_attributes_distict<-distinct(species_atributes)

ectos_pres_abs_df<-right_join(species_attributes_distict, ectos_pres_abs, by="species_jetz")   %>% arrange(elevation_cat)

ectos_pres_abs_df<-ectos_pres_abs_df %>% distinct(species_jetz,proportion_ectoparasites,.keep_all = TRUE) # Eliminates species taht ocurr at two elevatiosn and keep only one ( in alphabetical order of the stations, e.g for galbula low_elevations is kept and montane is deleted, but the samples size is mantained this is jut to be able to do the phylogenetic analyses properly)

ectos_pres_abs_df<- ectos_pres_abs_df %>% rename(family=BLFamilyLatin)



ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality) 

ectos_birds_dff_PA<-ectos_birds_dff %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation_mean=mean(elevation), sociality=max(sociality))%>% 
  mutate(proportion_ectoparasites=ectoparasites_presence/sample_size)

unique (ectos_birds_dff_PA$species_jetz )
#ectos_birds_dff_mean <- ectos_birds_dff %>% select(species_jetz,total_lice,total_no_feathers_mites,total_mesostigmatidae) %>% filter(complete.cases(.)) %>% group_by(species_jetz) %>% mutate(mean_lice = mean(total_lice), mean_nf_mites=mean(total_no_feathers_mites)) %>% droplevels()
phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include specis form manu and iquitos  so need to rpun the tree in the data processin section

a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_birds_dff$species_jetz)) %>% mutate(name=ectos_birds_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name) 



###_###_###_###_###_###_###_
# Combining both plots
###_###_###_###_###_###_###_

list.names=setNames(ectos_birds_dff_PA$proportion_ectoparasites, ectos_birds_dff_PA$species_jetz)

ColorPalette <- brewer.pal(n = 9, name = "GnBu")

fmode<-as.factor(setNames(ectos_birds_dff_PA$sociality,ectos_birds_dff_PA$species_jetz))
object = contMap(phylo, list.names, direction = "leftwards", plot=FALSE)
#object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))
object_color<-setMap(object, ColorPalette)

png("figures/figures_manuscript/Fig1a.Sociality_and_prevalence_phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
plot(dotTree(phylo,fmode,colors=setNames(c("red","black"),c("1","0")),ftype="i",fsize=0.5, lwd=4),text(x=10,y=-5,"Mixed-species flocks",pos=1))

plot(object_color$tree,colors=object_color$cols,add=TRUE,ftype="off",lwd=5,fsize=0.5,
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)

add.color.bar(10, object_color$cols, title = "", lims = object$lims, digits = 3, prompt=FALSE,x=70,y=-5, lwd=4,fsize=1,subtitle="Ectoparasites Prevalence",pos=4)
dev.off()

### fig 1 
list.names=setNames(ectos_birds_dff_PA$proportion_ectoparasites, ectos_birds_dff_PA$species_jetz)

# Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_birds_dff_PAf) <- ectos_birds_dff_PA$species_jetz # first make it the row names 
ectos_birds_dff_PA<- ectos_birds_dff_PA[match(phylo$tip.label,rownames(ectos_birds_dff_PA)),]

object = contMap(phylo, list.names, direction = "leftwards", plot=FALSE)

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

unique(ectos_df$species_jetz)

# Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_rooted$tip.label,rownames(ectos_df)),]

View(fmode)

png("figures/figures_manuscript/Fig1.sociality_phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
dotTree(phylogeny_rooted,fmode,colors=setNames(c("red","black"),                                              c("1","0")),ftype="i",fsize=0.5, lwd=4)
dev.off()

###_###_###_###_###_###_###_
# Combining both plots
###_###_###_###_###_###_###_

list.names=setNames(ectos_df$proportion_ectoparasites, ectos_df$species_jetz)

ColorPalette <- brewer.pal(n = 9, name = "GnBu")

fmode<-as.factor(setNames(ectos_df$sociality,ectos_df$species_jetz))
object = contMap(phylogeny_rooted, list.names, direction = "leftwards", plot=FALSE)
#object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))
object_color<-setMap(object, ColorPalette)

png("figures/figures_manuscript/Fig1a.Sociality_and_prevalence_phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
plot(dotTree(phylogeny_rooted,fmode,colors=setNames(c("red","black"),c("1","0")),ftype="i",fsize=0.5, lwd=4),text(x=10,y=-5,"Mixed-species flocks",pos=1))

plot(object_color$tree,colors=object_color$cols,add=TRUE,ftype="off",lwd=5,fsize=0.5,
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)

add.color.bar(10, object_color$cols, title = "", lims = object$lims, digits = 3, prompt=FALSE,x=70,y=-5, lwd=4,fsize=1,subtitle="Ectoparasites Prevalence",pos=4)
dev.off()


# ##### 5.Data processing prevalence ectos ----------------------------------------------------
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

phy_cov<-ape::vcv(phylo, corr=TRUE)

# ##### 5.1.Analyses models prevalence ectos --------------------------------------------------------
###_###_###
  #a) model glmm
ectos_p_glmm <-glmer(ectoparasites_PA~sociality+scale(elevation)+(1|Powder.lvl)+(1|species_jetz), #+elevation_midpoint+Powder.lvl
                             data = ectos_birds_dff, 
                             family = "binomial")
###_###_###
  #b) model pglmm
ectos_p_pglmm <-phyr::pglmm(ectoparasites_PA~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl), #+elevation_midpoint+Powder.lvl +(1|Powder.lvl)
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
class(ectos_p_pglmm ) 
summary(ectos_p_pglmm )
print(ectos_p_pglmm )
fixef(ectos_p_pglmm )
predict(ecto_p_pglmm )
rr2::R2(ectos_p_pglmm ) # goodness of fit

    #model pglmm assumptions check DHARMa
    #Overdisperssion
simulationOutput <-DHARMa::simulateResiduals(fittedModel =ectos_p_pglmm, plot = F)
plot(simulationOutput)
plotResiduals(simulationOutput , ectos_p_pglmm$data$sociality)
plotResiduals(simulationOutput , form =ectos_p_pglmm$sociality)

    #Homogenity of variance within groups (Heteroscedasticity) 
testCategorical(simulationOutput, catPred = ectos_birds_dff$sociality)

    # Other test
testUniformity(simulationOutput) #tests if the overall distribution conforms to expectations # 
testOutliers(simulationOutput)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
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
rr2::R2(ecto_p_pglmm_bayes) #0.1940
class(ecto_p_pglmm_bayes) 
# plot 
estimates_plot<-plot_bayes(ecto_p_pglmm_bayes)
png("data/data_analyses/models/model_plots/1.estimates_plot_model_prevalence_pglmm_phylo_multiple_obs_031623.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()

# assumptions check model pglmm_ bayes DHARMa ( DHARMa is not working with pgllm Bayes = True)
res_bayes<-DHARMa::simulateResiduals(ecto_p_pglmm_bayes)
class(ecto_p_pglmm_bayes)

## Model Individual  sample with phylognetic + non-phylo effects 
#d) model BRMS bayes
ecto_p_brms_bayes<-brms::brm(ectoparasites_PA~sociality+scale(elevation)+
          (1|gr(species_jetz, cov = phy_cov))+  #(1|Powder.lvl)
            (1|Powder.lvl)+
            (1|species),
        data=ectos_birds_dff,
        family= bernoulli(), # bernoulli() uses the (link = "logit").#zero_inflated_negbinomial() 
        data2 = list(phy_cov=phy_cov),
        iter=6000, warmup=3000, #First we need the specify how many iteration we want the MCMC to run, We need to specify how many chains we want to run.
        thin=2,
        control=list(adapt_delta=0.99, max_treedepth=12)) 

#saveRDS(ecto_p_brms_bayes, "data/data_analyses/models/1.model_prevalence_brms_phylo_multiple_obs_031623.RDS")
ecto_p_brms_bayes<-readRDS("data/data_analyses/models/1.model_prevalence_brms_phylo_multiple_obs_031623.RDS")
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
   color_scheme_set("teal")

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

png("data/data_analyses/models/model_plots/1.parameters_intervals_plot_model_PREVALENCE_brms_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot_intervals
dev.off()



posterior <- as.array(ecto_p_brms_bayes)
mcmc_areas(posterior  ,prob=0.90, prob_outer=0.95)
names(ecto_p_brms_bayes$formula$resp)

# ##### 5.2 Model predictions plots prevalence ectos----------------------------------



# ####### 6.Data processing abundance lice----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,total_lice,foraging_cat, sociality, total_lice) %>% 
  na.omit() 
  #filter(total_lice<60)  # removing outliers 

str(ectos_birds_dff)
mean(as.numeric(ectos_birds_dff$total_lice), na.rm=T)
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

# ####### 6.1.Analyses models abundance lice--------------------------------------------------------
# Lice

###_###_###
#a) model glmm
###_###_###
 str(ectos_birds_dff)
lice_a_glmm <-glmer(total_lice~sociality+scale(elevation)+(1|Powder.lvl)+(1|species_jetz), #+elevation_midpoint+Powder.lvl
                    data = ectos_birds_dff, 
                    family = "poisson")

# Summary
summary(lice_a_glmm )
rr2::R2(lice_a_glmm)# R2 Predicted 0.3951
fixef(lice_a_glmm)
predict(lice_a_glmm)
class(lice_a_glmm ) 

# Assumptions check
simulationOutput_a_lice<- DHARMa::simulateResiduals(fittedModel =lice_a_glmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_a_lice)
testUniformity(simulationOutput_a_lice) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_a_lice)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_a_lice) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_a_lice) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_a_lice, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_a_lice) ## tests if there are more zeros in the data than expected from the simulations

# zero inflated, overdispersed, with outliers
#b) model pglmm
###_###_###

lice_a_pglmm <-phyr::pglmm(total_lice~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                   data = ectos_birds_dff, 
                                   family = "poisson",
                                   cov_ranef = list(species_jetz= phylo), #class phylo
                                   #bayes = TRUE,
                                   add.obs.re = TRUE,
                                   REML = TRUE, 
                                   verbose = TRUE,
                                   s2.init = .25) # what is this last parameter for

histogram(ectos_birds_dff$total_lice) # id some outliers 
plot_data(lice_a_pglmm)
summary(lice_a_pglmm)
rr2::R2(lice_a_pglmm)
fixef(lice_a_pglmm)
predict(ecto_a_pglmm)

# Assumptions check
simulationOutput_a_lice_2<- DHARMa::simulateResiduals(fittedModel=lice_a_pglmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_a_lice_2)
testUniformity(simulationOutput_a_lice_2) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_a_lice_2)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_a_lice_2) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_a_lice_2) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_a_lice_2, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_a_lice_2) ## tests if there are more zeros in the data than expected from the simulations

# no outliers, but no homogenity of varianse KS and withig groups significant, and zero inflated

#c) model pglmm bayes: hierarchical Bayesian models fitted using integrated nested laplace approximation (INLA)
###_###_###

zip_lice_a_pglmm_bayes <-phyr::pglmm(total_lice~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                         data =ectos_birds_dff, 
                                         family ="zeroinflated.poisson", #POISSON  ="zeroinflated.poisson", #
                                         cov_ranef = list(species_jetz= phylo), #class phylo
                                         bayes = TRUE,
                                         verbose=FALSE,
                                         prior = "inla.default") # consider using    add.obs.re = T
#saveRDS(zip_lice_a_pglmm_bayes, "data/data_analyses/models/2.model_ABUNDANCE_LICE_pglmm_zip_phylo_multiple_obs_17032023.RDS")

#1) Summary of the model 
communityPGLMM.plot.re(x=zip_lice_a_pglmm_bayes ) 
summary(zip_lice_a_pglmm_bayes)   
rr2::R2(zip_lice_a_pglmm_bayes)
class(zip_lice_a_pglmm_bayes)

launch_shinystan(zip_lice_a_pglmm_bayes)

# Plots

estimates_plot<-plot_bayes(zip_lice_a_pglmm_bayes ) # for some reason does not allow me to plot zero inflated poisson 

#png("data/data_analyses/models/model_plots/2.parameters_plot_model_LICE_ABUNDANCE_pglmm_zip_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
#estimates_plot
#dev.off()

# Assumptions check DHARMa does not work with pglm bayes=TRUE so I can not evaluate model fit.
simulationOutput_a_lice_3<- DHARMa::simulateResiduals(fittedModel=zip_lice_a_pglmm_bayes, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_a_lice_3)
testUniformity(simulationOutput_a_lice_3) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_a_lice_3)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_a_lice_3) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_a_lice_3) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_a_lice_3, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_a_lice_3) ## tests if there are more zeros in the data than expected from the simulations


#d) model pglmm bayes : zero_inflated_negbinomial, ind scaled elevation
###_###_###

zinb_lice_a_brms_bayes<-brm(total_lice~sociality+scale(elevation)+
                     (1|gr(species_jetz, cov = phy_cov))+  
                    (1|Powder.lvl) + 
                    (1|species),
                    data=ectos_birds_dff,
                    family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                     data2 = list(phy_cov=phy_cov),
                    iter=4000, warmup=2000,
                    thin=2,
                    control=list(adapt_delta=0.99, max_treedepth=12)) 
zinb_lice_a_brms_bayes<-lice_a_brms_bayes
#saveRDS(lzinb_lice_a_brms_bayes, "data/data_analyses/models/2.model_ABUNDANCE_LICE_brms_zinb_phylo_multiple_obs_17032023.RDS")
zinb_lice_a_brms_bayes<-readRDS("data/data_analyses/models/2.model_ABUNDANCE_LICE_brms_zinb_phylo_multiple_obs_17032023.RDS")


#  exploring zero inflated poisson 
 
zip_lice_a_brms_bayes<-brm(total_lice~sociality+scale(elevation)+
                              (1|gr(species_jetz, cov = phy_cov))+  
                              (1|Powder.lvl) + 
                              (1|species),
                            data=ectos_birds_dff,
                            family=zero_inflated_poisson(),  #zero_inflated_negbinomial()
                            data2 = list(phy_cov=phy_cov),
                            iter=6000, warmup=3000,
                            thin=2,
                            control=list(adapt_delta=0.99, max_treedepth=12)) 

zip_lice_a_brms_bayes<-lice_a_brms_bayes
saveRDS(zip_lice_a_brms_bayes, "data/data_analyses/models/2.model_ABUNDANCE_LICE_brms_zip_phylo_multiple_obs_21032023.RDS")

# some outliers and smaller R2, overall performed porrly compared to zeroinflated negative binomial)_


# Summarize the model
summary(zinb_lice_a_brms_bayes)


fixef(zinb_lice_a_brms_bayes) # to get more detailed values for estimates
coef(zinb_lice_a_brms_bayes) # if you have group-level effects (hierarchical data)

# interpret the model 
#To test whether all regression coefficients are different from zero, we can look at the Credible Intervals that are listed in the summary output or we can visually represent them in density plots.
#If we do so, we clearly see that zero is not included in any of the density plots, meaning that we can be reasonably certain the regression coefficients are different from zero.
#INTERPRETATION:In the model, the parameter for Sociality means the expected difference between non_social(0) and social (1) with all other covariates held constant. we clearly see that zero is included in the density plot for sociality so there is not effect of sociality??
bayes_R2(zinb_lice_a_brms_bayes) #0.32  
bayes_R2(zip_lice_a_brms_bayes) # zip 0.22
plot(zip_lice_a_brms_bayes)
plot(zinb_lice_a_brms_bayes)
mcmc_plot(zinb_lice_a_brms_bayes) # Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
launch_shinystan()
pp_check(lice_a_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(lice_a_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict(zinb_lice_a_brms_bayes)
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=lice_a_brms_bayes$data$total, pp_m[1:200, ])  + 
  coord_cartesian(xlim = c(0, 5))

pp_m<- brms::posterior_predict(zinb_lice_a_brms_bayes)
ppc_rootogram(y=zinb_lice_a_brms_bayes$data$total_lice, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 5), ylim = c(0,30))

ppc_stat(y=lice_a_brms_bayes$data$total_lice, pp_m, stat ="prop_zero")


#Assumptions check model pglmm_ bayes
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(zinb_lice_a_brms_bayes, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations

# Plots BRMS 
color_scheme_set("teal") 
# model fit
png("data/data_analyses/models/model_plots/2.model_fit_ABUNDANCE_LICE_brms_phylo_multiple_obs_032123.png.png",width = 3000, height = 3000, res = 300, units = "px")
pp_m<- brms::posterior_predict(zinb_lice_a_brms_bayes)
ppc_rootogram(y=zinb_lice_a_brms_bayes$data$total_lice, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
dev.off()


# PLOT convergence and parameters distributions
png("data/data_analyses/models/model_plots/2.parameters_distribution_convergence_plot_model_ABUNDANCE_LICE_brms_phylo_multiple_obs_032123.png.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_lice_a_brms_bayes)
dev.off()

# plot posterior distributions 
color_scheme_set("teal")
brmstools::forest(ecto_p_brms_bayes, level = 0.95, show_data = TRUE)

estimates_plot<-mcmc_plot(zinb_lice_a_brms_bayes,prob=0.90, prob_outer=0.95,
                          variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                          type="areas") +
  labs(title="Posterior distributions Lice abundance", subtitle ="Lice abundance with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

estimates_plot_intervals<-mcmc_plot(zinb_lice_a_brms_bayes,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                                    type="intervals") +
  labs(title="Posterior distributions Lice abundance", subtitle ="Lice abundance with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

#png("data/data_analyses/models/model_plots/2.parameters_plot_model_ABUNDANCE_LICE_brms_zinb_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
#estimates_plot
#dev.off()
 
loo(zinb_lice_a_brms_bayes,zip_lice_a_brms_bayes,compare = TRUE)

# ###### 6.2 Model predictions plots abundance lice----------------------------------
#some ideas here: https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html
# this is not workingunsure why
install.packages("modelr")
library(modelr)

ectos_birds_dff %>%
  data_grid(sociality) %>%
  add_epred_draws(zinb_lice_a_brms_bayes, dpar = TRUE) %>%
  ggplot(ectos_birds_dff,aes(x = sociality, y = total_lice, color = sociality)) +
  stat_pointinterval(position = position_dodge(width = .4)) +
  scale_size_continuous(guide = "none") +
  scale_color_manual(values = brewer.pal(6, "Blues")[-c(1,2)])

names(ectos_birds_dff)
grid = ectos_birds_dff %>%
  data_grid(sociality)

means = grid %>%
  add_epred_draws(zinb_lice_a_brms_bayes)

preds = grid %>%
  add_predicted_draws(zinb_lice_a_brms_bayes)

AB %>%
  ggplot(aes(x = response, y = group)) +
  stat_halfeye(aes(x = .epred), scale = 0.6, position = position_nudge(y = 0.175), data = means) +
  stat_interval(aes(x = .prediction), data = preds) +
  geom_point(data = AB) +
  






# ####### 6.Data processing abundance MITES----------------------------------------------------
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality, total_mites, total_mesostigmatidae, total_no_feathers_mites) %>% 
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
#ectos_birds_dff$total_lice<-as.factor(ectos_birds_dff$total_lice)
ectos_birds_dff$total_mites<-as.numeric(ectos_birds_dff$total_mites)
ectos_birds_dff$total_mesostigmatidae<-as.numeric(ectos_birds_dff$total_mesostigmatidae)
ectos_birds_dff$total_no_feathers_mites<-as.numeric(ectos_birds_dff$total_no_feathers_mites)
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

phy_cov<-ape::vcv(phylo, corr=TRUE)


# ####### 6.1.Analyses models abundance MITES--------------------------------------------------------


####
#NON_FEATHER MITES
####

ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality, total_no_feathers_mites) %>% 
  na.omit() %>% 
  filter(total_no_feathers_mites<100)  # removing outliers for total mites

# Removing premnoplex ( Outlier)
ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select(elevation, species_jetz, Powder.lvl,foraging_cat, sociality, total_no_feathers_mites) %>% 
  na.omit() %>% 
  filter(species_jetz!="Premnoplex_brunnescens")  #Removing outliers for total mites


unique (ectos_birds_dff$species_jetz )
#filter(total_mesostigmatidae<100)  # removing outliers for total mites
#filter(total_no_feathers_mites<100)  # removing outliers for total mites

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

# Outliers
assertr::insist(ectos_birds_dff, within_n_mads(50), total_mites)
assertr::insist(ectos_birds_dff, within_n_mads(2), total_no_feathers_mites)

#ectos_birds_dff <-ectos_birds_dff  %>% mutate_at("species_jetz", str_replace, " ", "_")
#ectos_birds_dff$elevation_cat<-as.factor(ectos_birds_dff$elevation_cat)
ectos_birds_dff$foraging_cat<-as.factor(ectos_birds_dff$foraging_cat)
ectos_birds_dff$species_jetz<-as.factor(ectos_birds_dff$species_jetz)
ectos_birds_dff$elevation<-as.numeric(ectos_birds_dff$elevation)
#ectos_birds_dff$elevation_midpoint<-as.numeric(ectos_birds_dff$elevation_midpoint)
ectos_birds_dff$sociality<-as.factor(ectos_birds_dff$sociality)
ectos_birds_dff$Powder.lvl<-as.factor(ectos_birds_dff$Powder.lvl)
#ectos_birds_dff$total_lice<-as.factor(ectos_birds_dff$total_lice)
#ectos_birds_dff$total_mites<-as.numeric(ectos_birds_dff$total_mites)
#ectos_birds_dff$total_mesostigmatidae<-as.numeric(ectos_birds_dff$total_mesostigmatidae)
ectos_birds_dff$total_no_feathers_mites<-as.numeric(ectos_birds_dff$total_no_feathers_mites)
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

phy_cov<-ape::vcv(phylo, corr=TRUE)

#a) model glmm
###_###_###
nf_mites_a_glmm <-glmer(total_no_feathers_mites~sociality+scale(elevation)+(1|Powder.lvl)+(1|species_jetz), #+elevation_midpoint+Powder.lvl
                     data = ectos_birds_dff, 
                     family = "poisson")
# Summary
summary(nf_mites_a_glmm )
rr2::R2(nf_mites_a_glmm)
fixef(nf_mites_a_glmm)
predict(nf_mites_a_glmm)
class(nf_mites_a_glmm) 

# Assumptions check
simulationOutput_nf_mites_a<- DHARMa::simulateResiduals(fittedModel =nf_mites_a_glmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_nf_mites_a)
testUniformity(simulationOutput_nf_mites_a) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_nf_mites_a)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_nf_mites_a) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_nf_mites_a) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_nf_mites_a, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_nf_mites_a) ## tests if there are more zeros in the data than expected from the simulations

# zero inflated?
#b) model pglmm
###_###_###

nf_mites_a_pglmm <-phyr::pglmm(total_no_feathers_mites~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                            data = ectos_birds_dff, 
                            family = "poisson",
                            cov_ranef = list(species_jetz= phylo), #class phylo
                            #bayes = TRUE,
                            add.obs.re = TRUE,
                            REML = TRUE, 
                            verbose = TRUE,
                            s2.init = .25) # what is this last parameter for

histogram(ectos_birds_dff$total_nf_mites) # id some outliers 

summary(nf_mites_a_pglmm)
rr2::R2(nf_mites_a_pglmm)
fixef(nf_mites_a_pglmm)
predict(ecto_a_pglmm)

# Assumptions check
simulationOutput_nf_mites_2<- DHARMa::simulateResiduals(fittedModel=nf_mites_a_pglmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_nf_mites_2)
testUniformity(simulationOutput_nf_mites_2) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_nf_mites_2)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_nf_mites_2) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_nf_mites_2) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_nf_mites_2, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_nf_mites_2) ## tests if there are more zeros in the data than expected from the simulations

# no outliers, but no homogenity of varianse KS and withig groups significant, and zero inflated

#c) model pglmm bayes: hierarchical Bayesian models fitted using integrated nested laplace approximation (INLA)
###_###_###

# zero.inflated poisson did not converge

nf_mites_a_pglmm_bayes <-phyr::pglmm(total_no_feathers_mites~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                  data =ectos_birds_dff, 
                                  family ="zeroinflated.poisson", #POISSON  ="zeroinflated.poisson", #
                                  cov_ranef = list(species_jetz= phylo), #class phylo
                                  bayes = TRUE,
                                  verbose=FALSE,
                                  prior = "inla.default") # consider using    add.obs.re = T


#1) Summary of the model 
communityPGLMM.plot.re(x=nf_mites_a_pglmm_bayes) 
summary(nf_mites_a_pglmm_bayes )   
rr2::R2(nf_mites_a_pglmm_bayes)
class(ecto_abundance_pglmm_bayes)

launch_shinystan(ecto_abundance_pglmm_bayes)

# Plots

estimates_plot<-plot_bayes(nf_mites_a_pglmm_bayes) #  

#png("data/data_analyses/models/model_plots/2m.parameters_plot_model_nf_MITES_ABUNDANCE_zip_pglmm_phylo_multiple_obs_031623.png",width = 3000, height = 3000, res = 300, units = "px")
#estimates_plot
#dev.off()

# Assumptions check DHARMa does not work with pglm bayes=TRUE so I can not evaluate model fit.
simulationOutput_a_nf_mites_3<- DHARMa::simulateResiduals(fittedModel=nf_mites_a_pglmm_bayes, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_a_nf_mites_3)

#d) model pglmm bayes : zero_inflated_negbinomial, ind scaled elevation
###_###_###

zip_nf_mites_a_brms_bayes<-brm(total_no_feathers_mites~sociality+scale(elevation)+
                                  (1|gr(species_jetz, cov = phy_cov))+  
                                  (1|Powder.lvl) + 
                                  (1|species),
                                data=ectos_birds_dff,
                                family=zero_inflated_poisson,  #zero_inflated_negbinomial()
                                data2 = list(phy_cov=phy_cov),
                                iter=6000, warmup=3000,
                                thin=2,
                                control=list(adapt_delta=0.99, max_treedepth=12)) 
#saveRDS(zip_nf_mites_a_brms_bayes,"data/data_analyses/models/2m.model_ABUNDANCE_n_f_MITES_brms_zip_phylo_multiple_obs_wo_premnoplex.RDS")
zip_nf_mites_a_brms_bayes<-readRDS("data/data_analyses/models/2m.model_ABUNDANCE_n_f_MITES_brms_zip_phylo_multiple_obs_wo_premnoplex.RDS")
zip_nf_mites_a_brms_bayes<-readRDS("data/data_analyses/models/2m.model_ABUNDANCE_n_f_MITES_brms_zip_phylo_multiple_obs.RDS")

# removed Premnoplex_brunnescens


zinb_nf_mites_a_brms_bayes<-brm(total_no_feathers_mites~sociality+scale(elevation)+
                               (1|gr(species_jetz, cov = phy_cov))+  
                               (1|Powder.lvl) + 
                               (1|species),
                             data=ectos_birds_dff,
                             family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                             data2 = list(phy_cov=phy_cov),
                             iter=6000, warmup=3000,
                             thin=2,
                             control=list(adapt_delta=0.99, max_treedepth=12)) 

#saveRDS(zinb_nf_mites_a_brms_bayes, "data/data_analyses/models/2m.model_ABUNDANCE_n_f_MITES_brms_zinb_phylo_multiple_obs_wo_premnoplex.RDS")
zinb_nf_mites_a_brms_bayes<-readRDS("data/data_analyses/models/2m.model_ABUNDANCE_n_f_MITES_brms_zinb_phylo_multiple_obs_wo_premnoplex.RDS.RDS")
#zinb_nf_mites_a_brms_bayes<-readRDS("data/data_analyses/models/2m.model_ABUNDANCE_n_f_MITES_brms_zinb_phylo_multiple_obs.RDS")

# removed Premnoplex_brunnescens


# Summarize the model
summary (zinb_nf_mites_a_brms_bayes)
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)

# interpret the model 
#To test whether all regression coefficients are different from zero, we can look at the Credible Intervals that are listed in the summary output or we can visually represent them in density plots.
#If we do so, we clearly see that zero is not included in any of the density plots, meaning that we can be reasonably certain the regression coefficients are different from zero.
#INTERPRETATION:In the model, the parameter for Sociality means the expected difference between non_social(0) and social (1) with all other covariates held constant. we clearly see that zero is included in the density plot for sociality so there is not effect of sociality??
bayes_R2(zinb_nf_mites_a_brms_bayes) 
bayes_R2(zip_nf_mites_a_brms_bayes)
plot(zinb_nf_mites_a_brms_bayes)
mcmc_plot(nf_mites_a_brms_bayes) # Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
launch_shinystan()
pp_check(zinb_nf_mites_a_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
pp_check(zinb_nf_mites_a_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

pp_m<- brms::posterior_predict(zinb_nf_mites_a_brms_bayes)
log1 <- scale_x_continuous(trans="log1p")
ppc_dens_overlay(y=zinb_nf_mites_a_brms_bayes$data$total_no_feathers_mites, pp_m[1:200, ])  + 
  coord_cartesian(xlim = c(0, 5))

ppc_rootogram(y=nf_mites_a_brms_bayes$data$total_no_feathers_mites, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))

ppc_stat(y=nf_mites_a_brms_bayes$data$total_no_feathers_mites, pp_m, stat = "prop_zero")


#Assumptions check model pglmm_ bayes
#remotes::install_github("Pakillo/DHARMa.helpers")

simulate_residuals <- dh_check_brms(zinb_nf_mites_a_brms_bayes, integer = TRUE)
plot( simulate_residuals, form = ectos_birds_dff$sociality)
DHARMa::testDispersion(simulate_residuals)
DHARMa::testZeroInflation(simulate_residuals ) ## tests if there are more zeros in the data than expected from the simulations
testUniformity(simulate_residuals) #tests if the overall distribution conforms to expectations


# Plots BRMS 
color_scheme_set("orange")

# model convergence 
png("data/data_analyses/models/model_plots/2m.model_convergence_ABUNDANCE_nf_MITES_brms_zinb_phylo_multiple_obs_032123_wo_premnoplex.png",width = 3000, height = 3000, res = 300, units = "px")
plot(zinb_nf_mites_a_brms_bayes)
dev.off()


# model fit
png("data/data_analyses/models/model_plots/2m.model_fit_ABUNDANCE_nf_MITES_brms_zinb_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
pp_m<- brms::posterior_predict(zinb_nf_mites_a_brms_bayes)
ppc_rootogram(y=zinb_nf_mites_a_brms_bayes$data$total_no_feathers_mites, pp_m[1:200, ])  +   
  coord_cartesian(xlim = c(0, 100), ylim = c(0,30))
dev.off()

pp_m<- brms::posterior_predict(zinb_nf_mites_a_brms_bayes)

brmstools::forest(zinb_nf_mites_a_brms_bayes, level = 0.95, show_data = TRUE)

estimates_plot<-mcmc_plot(zinb_nf_mites_a_brms_bayes,prob=0.90, prob_outer=0.95,
                          variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                          type="areas") +
  labs(title="Posterior distribution mites (n_f)", subtitle ="with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

estimates_plot_intervals<-mcmc_plot(zinb_nf_mites_a_brms_bayes,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                                    type="intervals") +
  labs(title="Posterior distribution mites(n_f)", subtitle ="with means and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/models/model_plots/2m.parameters_plot_model_ABUNDANCE_nf_MITES_brms_zinb_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()


# Model comparisons 
loo(zip_nf_mites_a_brms_bayes,zinb_nf_mites_a_brms_bayes,compare = TRUE)

zip_nf_mites_a_brms_bayes


# ## ALL MITES ( excluded for now) ------------------------------------------------------------


# Mites
###_###_###
#a) model glmm
###_###_###
mites_a_glmm <-glmer(total_no_feathers_mites~sociality+scale(elevation)+(1|Powder.lvl)+(1|species_jetz), #+elevation_midpoint+Powder.lvl
                     data = ectos_birds_dff, 
                     family ="poisson")
# Summary
summary(mites_a_glmm )
rr2::R2(mites_a_glmm)
fixef(mites_a_glmm)
predict(mites_a_glmm)
class(mites_a_glmm) 

# Assumptions check
simulationOutput_a_mites<- DHARMa::simulateResiduals(fittedModel =mites_a_glmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_a_mites)
testUniformity(simulationOutput_a_mites) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_a_mites)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_a_mites) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_a_mites) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_a_mites, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_a_mites) ## tests if there are more zeros in the data than expected from the simulations

# zero inflated?
#b) model pglmm
###_###_###

mites_a_pglmm <-phyr::pglmm(total_mites~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                            data = ectos_birds_dff, 
                            family = "poisson",
                            cov_ranef = list(species_jetz= phylo), #class phylo
                            #bayes = TRUE,
                            add.obs.re = TRUE,
                            REML = TRUE, 
                            verbose = TRUE,
                            s2.init = .25) # what is this last parameter for

histogram(ectos_birds_dff$total_mites) # id some outliers 

summary(mites_a_pglmm)
rr2::R2(mites_a_pglmm)
fixef(mites_a_pglmm)
predict(ecto_a_pglmm)

# Assumptions check
simulationOutput_a_mites_2<- DHARMa::simulateResiduals(fittedModel=mites_a_pglmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_a_mites_2)
testUniformity(simulationOutput_a_mites_2) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_a_mites_2)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_a_mites_2) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_a_mites_2) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_a_mites_2, catPred = ectos_birds_dff$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_a_mites_2) ## tests if there are more zeros in the data than expected from the simulations

# no outliers, but no homogenity of varianse KS and withig groups significant, and zero inflated

#c) model pglmm bayes: hierarchical Bayesian models fitted using integrated nested laplace approximation (INLA)
###_###_###


zip_mites_a_pglmm_bayes <-phyr::pglmm(total_mites~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                      data =ectos_birds_dff, 
                                      family ="zeroinflated.poisson", #POISSON  ="zeroinflated.poisson", #
                                      cov_ranef = list(species_jetz= phylo), #class phylo
                                      bayes = TRUE,
                                      verbose=FALSE,
                                      prior = "inla.default") # consider using    add.obs.re = T


#1) Summary of the model 
communityPGLMM.plot.re(x=zip_mites_a_pglmm_bayes) 
summary(zip_mites_a_pglmm_bayes)   
rr2::R2(zip_mites_a_pglmm_bayes)
class(ecto_abundance_pglmm_bayes)

launch_shinystan(ecto_abundance_pglmm_bayes)

# Plots

estimates_plot<-plot_bayes(zip_mites_a_pglmm_bayes ) #  

#png("data/data_analyses/models/model_plots/2m.parameters_plot_model_nf_MITES_ABUNDANCE_pglmm_phylo_multiple_obs_032123.png",width = 3000, height = 3000, res = 300, units = "px")
#estimates_plot
#dev.off()

# Assumptions check DHARMa does not work with pglm bayes=TRUE so I can not evaluate model fit.
simulationOutput_a_mites_3<- DHARMa::simulateResiduals(fittedModel=zip_mites_a_pglmm_bayes, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_a_mites_3)

#d) model pglmm bayes : zero_inflated_negbinomial, ind scaled elevation
###_###_###

zinb_mites_a_brms_bayes<-brm(total_mites~sociality+scale(elevation)+
                               (1|gr(species_jetz, cov = phy_cov))+  
                               (1|Powder.lvl) + 
                               (1|species),
                             data=ectos_birds_dff,
                             family=zero_inflated_negbinomial(),  #zero_inflated_negbinomial()
                             data2 = list(phy_cov=phy_cov),
                             iter=6000, warmup=3000,
                             thin=2,
                             control=list(adapt_delta=0.99, max_treedepth=12)) 

#saveRDS(zinb_mites_a_brms_bayes, "data/data_analyses/models/2m.model_ABUNDANCE_MITES_brms_zinb_phylo_multiple_obs.RDS")

zinb_mites_a_brms_bayes<-readRDS("data/data_analyses/models/2m.model_ABUNDANCE_MITES_brms_zinb_phylo_multiple_obs.RDS")

# Summarize the model
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)

# interpret the model 
#To test whether all regression coefficients are different from zero, we can look at the Credible Intervals that are listed in the summary output or we can visually represent them in density plots.
#If we do so, we clearly see that zero is not included in any of the density plots, meaning that we can be reasonably certain the regression coefficients are different from zero.
#INTERPRETATION:In the model, the parameter for Sociality means the expected difference between non_social(0) and social (1) with all other covariates held constant. we clearly see that zero is included in the density plot for sociality so there is not effect of sociality??
bayes_R2(zinb_mites_a_brms_bayes) 
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

#color_scheme_set("red")
brmstools::forest(ecto_p_brms_bayes, level = 0.95, show_data = TRUE)

estimates_plot<-mcmc_plot(lice_a_brms_bayes,prob=0.90, prob_outer=0.95,
                          variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                          type="areas") +
  labs(title="Posterior distributions", subtitle ="with medians and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

estimates_plot_intervals<-mcmc_plot(lice_a_brms_bayes,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    variable = c("b_Intercept", "b_sociality1", "b_scaleelevation","sd_Powder.lvl__Intercept","sd_species__Intercept","sd_species_jetz__Intercept"),
                                    type="intervals") +
  labs(title="Posterior distributions", subtitle ="with means and 95% intervals")+
  theme_minimal(20)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40")+
  xlab("Estimate")

png("data/data_analyses/models/model_plots/2.parameters_plot_model_LICE_ABUNDANCE_brms_zinb_phylo_multiple_obs_031623.png",width = 3000, height = 3000, res = 300, units = "px")
estimates_plot
dev.off()


# ###### 6.2 Model predictions plots abundance MITES----------------------------------

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
# Summarize the model
summary ()
fixef() # to get more detailed values for estimates
coef() # if you have group-level effects (hierarchical data)

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


# #PLOTS FOR BAYESIAN MODELS [FUNCTIONS] ----------------------------------


#### playing with teh functions to create plots for bayesian

# i MODIFIED THE FUNCTION TO BE ABLE TO PLOT IT IN NEGATIVE BINOMIAL


plot_bayes.communityPGLMM <- function(x, n_samp = 1000, sort = TRUE, ...) {
  
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    stop('plot_bayes requires the ggplot2 package but it is unavailable. Use install.packages("ggplot2") to install it.')
  }
  
  if(!x$bayes) {
    stop("plot_bayes only works on communityPGLMM objects fit with bayes = TRUE")
  }
  
  if(!requireNamespace("ggridges", quietly = TRUE)) {
    stop('plot_bayes requires the ggridges package but it is unavailable. Use install.packages("ggridges") to install it.')
  }
  
  r.effect_jen<-c("1|species_jetz", "1|species_jetz__","1|Powder.lvl","1|obs")
  
  re.names <-r.effect_jen
  if (x$family == "zeroinflated.poisson") re.names <- c("residual", re.names)
  random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                         function(x) INLA::inla.rmarginal(1000, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>%
    setNames(re.names) %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "var",
                        values_to = "val") %>%
    dplyr::mutate(effect_type = "Random Effects")
  
  fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(1000, x)) %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "var",
                        values_to = "val") %>%
    dplyr::mutate(effect_type = "Fixed Effects")
  
  samps <- dplyr::bind_rows(random_samps, fixed_samps) %>%
    dplyr::mutate(effect_type = factor(effect_type, 
                                       levels = c("Random Effects", "Fixed Effects")))
  
  ci <- samps %>%
    dplyr::group_by(var, effect_type) %>%
    dplyr::summarise(lower = quantile(val, 0.025),
                     upper = quantile(val, 0.975),
                     mean = mean(val),
                     .groups = "drop_last")
  
  if(sort){
    ci <- dplyr::arrange(ci, mean) %>% dplyr::ungroup() %>% 
      dplyr::mutate(var = factor(as.character(var), levels = as.character(var)))
  }
  
  sig_vars <- ci %>%
    dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                               "CI no overlap with zero",
                               ifelse(sign(lower) == sign(upper),
                                      "CI no overlap with zero",
                                      "CI overlaps zero"))) %>%
    dplyr::select(var, sig)
  
  if(sort){
    samps <- dplyr::mutate(samps, var = factor(var, levels = levels(sig_vars$var)))
  }
  
  samps <- samps %>%
    dplyr::left_join(sig_vars, by = "var") %>%
    dplyr::group_by(var) %>%
    dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
    dplyr::ungroup()
  
  pal <- c("#fc8d62", "#8da0cb")
  
  #samps<-samps %>% filter(var!="1|obs")
  p <- ggplot2::ggplot(samps, ggplot2::aes(val, var, height = ..density..)) +
    ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig), 
                                  stat = "density", adjust = 2, color = "gray70") +
    ggplot2::geom_point(ggplot2::aes(x = mean, y = var), data = ci, inherit.aes = FALSE) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper, y = var), data = ci,
                            inherit.aes = FALSE, height = 0.1) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::scale_alpha_manual(values = c(0.8, 0.2)) +
    ggplot2::scale_fill_manual(values = rev(pal)) +
    ggplot2::facet_wrap(~ effect_type, nrow = 2, scales = "free") +
    ggplot2::ylab("") +
    ggplot2::xlab("Estimate") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(size = 14),
                   strip.text = ggplot2::element_text(size = 16))
  
  p
}

##### Note on interpreting the plots of the posterior distributions with credible intervals
#Remember that the posterior distribution represents our uncertainty (or certainty) in "q"
#after combining the information in the data (the likelihood) with what we knew before collecting data (the prior).
#Bayesian inference we quantify statements like this – that a particular event is “highly likely” – by computing the “posterior probability” of the event,
#which is the probability of the event under the posterior distribution.

#This plots do not tell you the significance of the factors, instead they tell you the effect size of the factors.
# The most common summaries of a posterior distribution are interval estimates and point estimates.
# Point estimates are typically obtained by computing the mean or median (or mode) of the posterior distribution. These are called the “posterior mean” or the “posterior median” (or “posterior mode”).
# Point estimates Essentially this boils down to summarizing the posterior distribution by a single number.
# Interval estimates can be obtained by computing quantiles of the posterior distribution. Bayesian Confidence intervals are often called “Credible Intervals”.
# We can extend this idea to assess the certainty (or confidence) that "q" lies in any interval. For example, from the plot it looks like q
# will very likely lie in the interval [0.2,0.4] because most of the posterior distribution mass lies between these two numbers.
#it is more common to compute Bayesian Confidence Intervals the other way around: specify the level of confidence we want to achieve and find an interval that achieves that level of confidence. 
#This can be done by computing the quantiles of the posterior distribution. For example, the 0.05 and 0.95 quantiles of the posterior would define a 90% Bayesian Confidence Interval.
# Positive estimates are positive effects and  negative values  negative effects, when they overlap with zero you can conclude the effect size does not differ from zero.

# In our models the predictor sociality has two levels 0 sociality and 1 sociality, and so
  #!. The intercept estimate  CORRRESPOND TO  the first  level of sociality ( sociality zero)
  # sociality estimate is the difference between the sociality zero and the sociality one, so if it overlaps with zero means teh difference is close to zero between being social and non social  



