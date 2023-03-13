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
install.packages('brms') # bayesian approach to model phylogenetic data with repides observations
#model assumptions
install.packages("DHARMa")


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
options(repos = c(
  phyr = 'https://daijiang.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
install.packages('phyr')

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
install.packages("TreeTools")
library(TreeTools)
#Libraries ofr plots
library(gridExtra)
library(ggpubr)
library(grid)

# [Presence_absence] Part1_Ectoparasite models_using binary data  1/0 Unsure if this is useful or will give similar results thanthe one below need to rethink this piece

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
# This dataset contains all parasites samples (887) after removing duplicated rows, for which we have an assignated elevation, for 783 in total out of 998 that we had originally  (this included some duplicates)
ectos_df<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date)

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex") 

unique(ectos_df$species_jetz )

# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
#ectos_df<-ectos_df %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

# Make sure the tips and the names on the file coincide and formating of name is consitent
phylo$edge.length  
phylo$tip.label
is.binary(phylo)

# Make sure this two are the same numbers 

a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_df$species_jetz)) %>% mutate(name=ectos_df$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)

# Drop some tips USE IF NEED TO DROP SOME TIPS when using teh full phylogeny
#phylo<-drop.tip (phylogeny_prevalence, tip$name) 

# Important!!!!!Make sure the names  are in the same order in the phylogeny and in the traits 
#( I am not sure if this is relevatnt for the PGLMM cause we have multiple observations)
#rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
#ectos_df<- ectos_df[match(phylo$tip.label,rownames(ectos_df)),]

# Re-strudture the data
# Make sure variables are in the right format, random effects should be factors
#We need to align the variable name and the structure to the names in the column tip.label used for the phylogeny?

ectos_df <-ectos_df  %>% mutate_at("species_jetz", str_replace, " ", "_")
ectos_df$elevation_cat<-as.factor(ectos_df$elevation_cat)
ectos_df$foraging_cat<-as.factor(ectos_df$foraging_cat)
ectos_df$species_jetz<-as.factor(ectos_df$species_jetz)
ectos_df$elevation<-as.numeric(ectos_df$elevation)
ectos_df$elevation_midpoint<-as.numeric(ectos_df$elevation_midpoint)
ectos_df$sociality<-as.factor(ectos_df$sociality)
ectos_df$Powder.lvl<-as.factor(ectos_df$Powder.lvl)

names(ectos_df)

is.ultrametric(phylogeny_prevalence)

#MODEL FITTING AND EVALUALLING MODEL FIT GRAPHICALLY

###_###_###_ Important for model fitting
#predict ( model) # plot predicted vlues and confidence limits on the logit scale, the points are not logit tranformed data, they are "working values to fit the model.....dificult to interpret 
#fitted (model) # willl give us the predicted values transformed to the original scale USE THESE TO PREDICT THE DATA
#visreg(model) # Similarly  plot predicted values and confidence limits on the logit scale, dificult to interpret
#visreg(model, scale="response") # willl give us the predicted values transformed to the original scale 
###_###_###_ ###_###_###_ ###_###_###_ 


###
#PREVALENCE MODELS 
###

###_###_###_ 
# Using PGLMM [ Lets correct for phylogeny ]
###_###_###_ 

#I am not sue about including foraging cat since that some how is included in teh variation per species  ( Removed) 
# we revomed (1|foraging_cat) because it was not significant in individual models 
#notes 
#As with the linear mixed model, it is a very good idea to standardize the predictor (independent) variables to have mean 0 and variance 1. 
#This will make the function more robust and improve the interpretation of the regression coefficients.

#lm(proportion_ectoparasites~1, data=ectos_df) # even the average is a linear regression INTERESTING

# remove ind with NAs for any of the variables of interest 
#   ectos_df_wo_na<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
#   rename(elevation=elevation_extrapolated_date) %>% select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, sociality ) %>% 
#   na.omit()

names( ectos_df)
ecto_prevalence_glmm <-glmer(ectoparasites_PA~sociality+elevation+(1|Powder.lvl), #+elevation_midpoint+Powder.lvl
                             data = ectos_df, 
                             family = "binomial")
ecto_prevalence_glmm <-glmer(ectoparasites_PA~sociality+scale(elevation)+(1|Powder.lvl), #+elevation_midpoint+Powder.lvl
                             data = ectos_df, 
                             family = "binomial")

ecto_prevalence_pglmm <-phyr::pglmm(ectoparasites_PA~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl), #+elevation_midpoint+Powder.lvl
                                      data = ectos_df_wo_na, # exlude raws with nas this is usefull to check the model assumptions 
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
a<-predict(ecto_prevalence_pglmm)
rr2::R2(ecto_prevalence_pglmm) # goodness of fit
class(ecto_prevalence_pglmm ) 

# i am still not sure how to extract and plot confidence intervals for pGLMM that are not bayes


 
ecto_prevalence_pglmm_bayes <- pglmm(ectoparasites_PA~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                         data = ectos_df_wo_na, 
                         cov_ranef = list(species_jetz= phylo), #class phylo
                         family = "binomial", # zeroinflated.binomial in reality I want to use a zeroinflated.binomial but for some reason i can not plot those
                         bayes = TRUE,
                         prior = "pc.prior.auto")

names(ectos_df)
summary(ecto_prevalence_pglmm_bayes)
print(ecto_prevalence_pglmm)
rr2::R2(ecto_prevalence_pglmm_bayes) 
class(ecto_prevalence_pglmm ) 

# plot 
plot_bayes(ecto_prevalence_pglmm_bayes )

#pdf(file="outputs_models/model_prevalence_bayes.pdf")
#plot_bayes(ecto_prevalence_pglmm_bayes )
#dev.off()

#What we are looking for is that the posterior distribution mode is well away from zero, and that it looks relatively symmetrical. 
#If it were skewed and crammed up against the left side of the plot, near zero, we would conclude that the effect is weak (remembering that variance components cannot be less than or equal zero,
#so there will always be some positive probability mass). The most obvious effects (well away from zero) are again the phylogenetic species random effect. 


# Model assumption checks 
#   Read this carefully https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#motivation
#   check is whether assumptions of the model are met by the data. 
#   The typical way to do this is by plotting and/or analyzing the model residuals. 
#   In non-Gaussian models such as this one, this can be less straightforward
#   However, phyr output supports the DHARMa package, which can generated a generalized type of residual known as randomized quantile residuals (or sometimes Dunn-Smyth residuals).
#   These can be calculated and inspected for nearly any error distribution.We can produce standard diagnostic plots for our pglmm model.

install.packages("DHARMa")
library(DHARMa)
citation("DHARMa")
###
# for PREVALENCE PGLMM ( no bayes)
###
# RESIDUAL PLOTS to test for independence and normal distribution of eroors ()
#OVERDISPERSION
# 1. The most important assumption of a linear regression model is that the ERRORS are independent and normally distributed.
# overdispersion is the presence of greater variability in a data set than would be expected based on a given statistical model
#signals of overdispersion QQ: s-shaped QQ plot, distribution test (KS) signficant, QQ: Dispersion test is significant, QQ: Outlier test significant, Res ~ predicted: Quantile fits are spread out too far
#(residuals) should be completely random & unpredictable i.e stochastic. Hence, we want our residuals to follow a normal distribution. 
#And that is exactly what we look for in a residual plot.
#A residual is a measure of how far away a point is vertically from the regression line
#A typical residual plot has the residual values on the Y-axis and the independent variable on the x-axis
###_###_###
# The 'DHARMa' package uses a simulation-based approach to create readily interpretable scaled (quantile) residuals for fitted (generalized) linear mixed models

plot(res)
simulationOutput <- DHARMa::simulateResiduals(fittedModel =ecto_prevalence_pglmm, plot = F)
plot(simulationOutput)

#a qq-plot to detect overall deviations from the expected distribution, by default with added tests for correct distribution (KS test), dispersion and outliers.
#Note that outliers in DHARMa are values that are by default defined as values outside the simulation envelope, not in terms of a particular quantile. Thus, which values will appear as outliers will depend on the number of simulations.
#If you want outliers in terms of a particuar quantile, you can use the outliers() function.

plotQQunif(simulationOutput) # left plot in plot

#plotResiduals (right panel) produces a plot of the residuals against the predicted value (or alternatively, other variable). Simulation outliers (data points that are outside the range of simulated values) are highlighted as red stars. 
#These points should be carefully interpreted, because we actually don’t know “how much” these values deviate from the model expectation. 
plotResiduals(simulationOutput) # right plot in plot

#By default, plotResiduals plots against predicted values. 
#However, you can also use it to plot residuals against a specific other predictors (highly recommend).

#plotResiduals(resids , ecto_prevalence_pglmm)
plotResiduals(simulationOutput , ecto_prevalence_pglmm$data$sociality)
plotResiduals(simulationOutput , form =ecto_prevalence_pglmm$sociality)

plotResiduals(simulationOutput , ecto_prevalence_pglmm$data$elevation) 
plotResiduals(simulationOutput , form =ecto_prevalence_pglmm$elevation)

#
# 2.HOMOGENITY OF VARIANCE  (Heteroscedasticity) to test for independence and normal distribution of eroors
#Heteroskedastic refers to a condition in which the variance of the residual term, or error term, in a regression model varies widely
#
# If the predictor is a factor, or if there is just a small number of observations on the x axis, 
# plotResiduals will plot a box plot with additional tests instead of a scatter plot.
#  under H0 (perfect model), we would expect those boxes to range homogenously from 0.25-0.75. 
# To see whether there are deviations from this expecation, the plot calculates a test for uniformity per box, and a test for homgeneity of variances between boxes. 
#A positive test will be highlighted in red
 # See ?plotResiduas for details, .

plotResiduals(res , ecto_prevalence_pglmm$data$sociality)
testCategorical(simulationOutput, catPred = ectos_df_wo_na$sociality)

#See ?plotResiduas for details under H0 (perfect model), we would expect those boxes to range homogenously from 0.25-0.75. To see whether there are deviations from this expecation, the plot calculates a test for uniformity per box, and a test for homgeneity of variances between boxes. A positive test will be highlighted in red.


# explore teh variables to see that they are same dimensions
length(res$observedResponse)
length(res$fittedPredictedResponse)
length(res$scaledResiduals)
res$simulatedResponse

# OTHER TEST 
testUniformity(simulationOutput) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput, catPred = ectos_df_wo_na$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput) ## tests if there are more zeros in the data than expected from the simulations
testGeneric(simulationOutput) # test if a generic summary statistics (user-defined) deviates from model expectations
testTemporalAutocorrelation(simulationOutput) # tests for temporal autocorrelation in the residuals
testSpatialAutocorrelation(simulationOutput) # tests for spatial autocorrelation in the residuals. Can also be used with a generic distance function, for example to test for phylogenetic signal in the residuals



###
# for PREVALENCE PGLMM ( bayes) 
# Seems like it is not working nsure why since tehre is an example online that suggest it is implemented ni the package 
###

res_bayes<-DHARMa::simulateResiduals(ecto_prevalence_pglmm_bayes)

class(ecto_prevalence_pglmm_bayes)
ecto_prevalence_pglmm_bayes

###
#ABUNDANCE MODELS 
###

# remove ind with NAs for any of the variables of interest 
 ectos_df_abundance_wo_na<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
 rename(elevation=elevation_extrapolated_date) %>%
   select(elevation, species_jetz, Powder.lvl,total_lice, sociality ) %>% 
   na.omit() %>% 
   filter(total_lice<60)  

  ectos_df_abundance_wo_na <- ectos_df_abundance_wo_na  %>% mutate_at("species_jetz", str_replace, " ", "_")
  #ectos_df_abundance_wo_na$elevation_cat<-as.factor( ectos_df_abundance_wo_na$elevation_cat)
  #ectos_df_abundance_wo_na$foraging_cat<-as.factor( ectos_df_abundance_wo_na$foraging_cat)
  ectos_df_abundance_wo_na$species_jetz<-as.factor( ectos_df_abundance_wo_na$species_jetz)
  ectos_df_abundance_wo_na$elevation<-as.numeric( ectos_df_abundance_wo_na$elevation)
  #ectos_df_abundance_wo_na$elevation_midpoint<-as.numeric( ectos_df_abundance_wo_na$elevation_midpoint)
  ectos_df_abundance_wo_na$sociality<-as.factor( ectos_df_abundance_wo_na$sociality)
  ectos_df_abundance_wo_na$Powder.lvl<-as.factor( ectos_df_abundance_wo_na$Powder.lvl)
 

ecto_abundance_pglmm <-phyr::pglmm(total_lice~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                    data = ectos_df_abundance_wo_na, 
                                    family = "poisson",
                                    cov_ranef = list(species_jetz= phylo), #class phylo
                                    #bayes = TRUE,
                                    add.obs.re = TRUE,
                                    REML = TRUE, 
                                    verbose = TRUE,
                                    s2.init = .25) # what is this last parameter for

histogram(ectos_df_abundance_wo_na$total_lice) # id some outliers 

summary(ecto_abundance_pglmm)
rr2::R2(ecto_abundance_pglmm )
fixef(ecto_abundance_pglmm)
predict(ecto_abundance_pglmm)
class(ecto_prevalence_pglmm ) 

# Checking model assumptions
#DHARMa (DHARMa stands for “Diagnostics for HierArchical Regression Models” )
# READ THIS https://stats.stackexchange.com/questions/490680/significant-dispersion-test
#OVERDISPERSION
#signals of overdispersion QQ: s-shaped QQ plot, distribution test (KS) signficant, QQ: Dispersion test is significant, QQ: Outlier test significant, Res ~ predicted: Quantile fits are spread out too far

# R: Our data fitted is not overdisppersed according to the evaluation, only KS significant and when looking at the test [testDispersion()]the value is close to 1, so we wont worry about this 
simulationOutput_abundance<- DHARMa::simulateResiduals(fittedModel =ecto_abundance_pglmm, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T
plot(simulationOutput_abundance)

plotResiduals(simulationOutput_abundance, ecto_abundance_pglmm$ectos_df_abundance_wo_na$sociality)
plotResiduals(simulationOutput_abundance, form =ecto_abundance_pglmm$sociality)
#plotResiduals(simulationOutput_abundance, form =ecto_abundance_pglmm$sociality$"1") # plot residuals for individual group


plotResiduals(simulationOutput , simulationOutput_abundance$data$elevation) 
plotResiduals(simulationOutput , form =simulationOutput_abundance$elevation)

# 2.HOMOGENITY OF VARIANCE  (Heteroscedasticity) to test for independence and normal distribution of eroors
# R: Our data fitted does not have Homogenity of variance within the groups, UNSURE HOW TO SOLVE THIS 
#the test alerts you to the fact that there are some within-group residual distributions that are not uniform,
#i.e. if you plot your residuals for a specific group (which is the one highlighted in red), they don't look uniform -> deviate significantly from your model assumptions.
#What your plots tell you is that there is a misfit that is not random.
# some ideas https://github.com/florianhartig/DHARMa/issues/291

plotResiduals(simulationOutput_abundance , simulationOutput_abundance$data$sociality)
plotResiduals(simulationOutput_abundance , simulationOutput_abundance$data$sociality[["1"]) # individual groups

testCategorical(simulationOutput_abundance, catPred = ectos_df_abundance_wo_na$sociality)


View(simulationOutput_abundance)
# explore the variables to see that they are same dimensions
length(simulationOutput_abundance$observedResponse)
length(simulationOutput_abundance$fittedPredictedResponse)
length(simulationOutput_abundance$scaledResiduals)
res$simulatedResponse
# OTHER TEST 
testUniformity(simulationOutput_abundance) #tests if the overall distribution conforms to expectations
testOutliers(simulationOutput_abundance)#  tests if there are more simulation outliers than expected
testDispersion(simulationOutput_abundance) # tests if the simulated dispersion is equal to the observed dispersion
testQuantiles(simulationOutput_abundance) #fits a quantile regression or residuals against a predictor (default predicted value), and tests of this conforms to the expected quantile
testCategorical(simulationOutput_abundance, catPred = ectos_df_abundance_wo_na$sociality)# tests residuals against a categorical predictor
testZeroInflation(simulationOutput_abundance) ## tests if there are more zeros in the data than expected from the simulations
testGeneric(simulationOutput_abundance) # test if a generic summary statistics (user-defined) deviates from model expectations
testTemporalAutocorrelation(simulationOutput_abundance) # tests for temporal autocorrelation in the residuals
testSpatialAutocorrelation(simulationOutput_abundance) # tests for spatial autocorrelation in the residuals. Can also be used with a generic distance function, for example to test for phylogenetic signal in the residuals


# BAYESIAN _ABUNDANCE
# readhttps://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMaForBayesians.html
#prior	
#Which type of default prior should be used by pglmm? Only used if bayes = TRUE.
#There are currently four options: "inla.default", which uses the default INLA priors; 
#"pc.prior.auto", which uses a complexity penalizing prior (as described in Simpson et al. (2017)) designed to automatically choose good parameters (only available for gaussian and binomial responses); 
#"pc.prior", which allows the user to set custom parameters on the "pc.prior" prior, using the prior_alpha and prior_mu parameters (Run INLA::inla.doc("pc.prec") for details on these parameters); 
#and "uninformative", which sets a very uninformative prior (nearly uniform) by using a very flat exponential distribution. The last option is generally not recommended but may in some cases give estimates closer to the maximum likelihood estimates. 
#"pc.prior.auto" is only implemented for family = "gaussian" and family = "binomial" currently.


ecto_abundance_pglmm_bayes <-phyr::pglmm(total_lice~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                   data = ectos_df_abundance_wo_na, 
                                   family ="zeroinflated.poisson", #POISSON
                                   cov_ranef = list(species_jetz= phylo), #class phylo
                                   bayes = TRUE,
                                   verbose=FALSE,
                                   prior = "inla.default",
                                   add.obs.re = T)

summary(ecto_abundance_pglmm_bayes)   

rr2::R2(ecto_abundance_pglmm_bayes)

 #pdf(file="outputs_models/model_abundance_bayes.pdf")
plot_bayes(ecto_abundance_pglmm_bayes ) # for some reason does not allow me to plot zero inflated poisson 
#dev.off()

plot_bayes.communityPGLMM (ecto_abundance_pglmm_bayes) # IMPORTANT plese note that I modiefied the function ( see the end of the script)cause it was giving me problems with the th random effect ( which i think is the obsevations but somehow it is not being extracted in the model )

DHARMa::checkModel

View(ecto_abundance_pglmm_bayes)
fittedModel=ecto_abundance_pglmm_bayes

simulationOutput_abundance_bayes<-simulateResiduals1(fittedModel, plot= TRUE, re.form = NULL ) #quantreg=T

DHARMa<-plot(simulationOutput_abundance_bayes)

plot(simulationOutput_abundance_bayes)
simulationOutput_abundance_bayes<- DHARMa::simulateResiduals1(fittedModel =ecto_abundance_pglmm_bayes, plot = F,integerResponse = T, re.form = NULL ) #quantreg=T

simulate(ecto_abundance_pglmm_bayes)


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

# Underestanding the summary of the model 
#random effects 
#The random effect with the largest variance and standard variation is the one with the strongest effect, in our case the phylogenetic effect,
# this implies that the parasites  prevalence is more similar in closely related species 

# Explore if the random effects are important? 
# One way to get an idea is to run a likelihood ratio test on the random effect. 
# This can be achieved by using the pglmm_profile_LRT() function, at least for binomial models ( copied from https://daijiang.github.io/phyr/articles/phyr_example_empirical.html)

phyr::pglmm_profile_LRT(ecto_prevalence_pglmm, re.number = 1) ## here we need to specify qith random effect are we evaluating in the order we put them in the model species is number 1, phylogeny2, and elevation 3

LRTest <- sapply(1:3, FUN = function(x) phyr::pglmm_profile_LRT(ecto_prevalence_pglmm, re.number = x))
colnames(LRTest) <- names(ecto_prevalence_pglmm$ss)
t(LRTest)

### other code

plot_data(
  ecto_prevalence_pglmm_bayes,
  sp.var = "species_jetz",
  site.var = "sociality",
  show.sp.names = FALSE,
  show.site.names = FALSE,
  digits = max(3, getOption("digits") - 3),
  predicted = FALSE,
)

###_###_###_###_###_###_###_###_####_####_###_###_###_###_###_###_###_###_###_####_####_###_###_###_###_###_###_###_###_###_####_####_###
###### Using BRMS package   
# see example and documentation here  http://paul-buerkner.github.io/brms/articles/brms_phylogenetics.html
###_###_###_###_###_###_###_###_####_####_######_###_###_###_###_###_###_###_####_####_######_###_###_###_###_###_###_###_####_####_###

## Model 3: Individual mass with phylognetic + non-phylo effects

mean_lice<-ectos_df %>% group_by(species_jetz) %>% 
  summarize(mean_lice=mean(total_lice))

df_ectos_include_mean<-full_join( ectos_df, mean_lice, by="species_jetz")

A <- ape::vcv.phylo(phylo)

  m<- brm ( 
    mean_lice ~ sociality+elevation+ 
      (1 | gr(phylo, cov=A)) + (1|Powder.lvl) + (1|species_jetz), 
       data=df_ectos_include_mean, 
       data2=list(A=A),
       family=zero_inflated_poisson(), 
    prior = c(
      prior(normal(0,10), "b"),
      prior(normal(0,50), "Intercept"),
      prior(student_t(3,0,20), "sd"),
      prior(student_t(3,0,20), "sigma")
    ),
    sample_prior = TRUE, chains = 2, cores = 2, 
    iter = 4000, warmup = 1000
  )
  
  saveRDS(m3, "../results/model_full_physpec.rds")
  # 17 minutes sampling on Dell XPS 13 9343
  

# model5 Removing Phylogenetic covariance structure


# Full Model - no phylo structure

  m5 <- brm(total_lice ~ sociality+ elevation + 
              (1 | species_jetz) + (1|Powder.lvl), 
            data=ectos_df, family=zero_inflated_poisson(), 
            iter=6000, thin=2,
            control=list(adapt_delta=0.90))
  
  mcmc_plot(m5, prob=0.80, prob_outer=0.95, point_est="mean")
  # 2.25 minutes sampling on Dell XPS 13 9343...


#### Checking model assumptions
  ### 
  #OTHER OPTION FOR MODELS NOT SUPPORTED IN DRAMa
  #PGLMM models are currently not supported we need to Using DHARMa for residual checks of unsupported models
  # See details here https://aosmith.rbind.io/2017/12/21/using-dharma-for-residual-checks-of-unsupported-models/#disqus_thread
  
  # i)U se the function simulate if available for your model to create new data
  # ii) Use createDHARMa () 
  #  iii) Need to include three type pf data simulated data, observed data, adn model predictiosn 
  # note if useing count data use integerResponse= TRUE
  
  # need to remove na from dataset s that the dimensions of responde simulated and observed match
  ectos_df_wo_na<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
    rename(elevation=elevation_extrapolated_date) %>% select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, sociality) %>% 
    na.omit()
  
  newData <- simulate(ecto_prevalence_pglmm, nsim=1000)
  
  class(newData)
  
  
  sim_resp_ecto_prevalence_pglmm=createDHARMa(simulatedResponse=newData,
                                              observedResponse = ectos_df_wo_na$ectoparasites_PA,
                                              fittedPredictedResponse = predict(ecto_prevalence_pglmm),
                                              integerResponse = FALSE)
  
  
  dim(simulatedResponse)
  as.matrix( ectos_df_wo_na$ectoparasites_PA)
  dim(observedResponse)
  
  plot(sim_resp_ecto_prevalence_pglmm)
  

    #### playing with teh functions to create plots for bayesian
  
  
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
    
  
  
  ##### Playing with teh Dharma functions 
  
  
  
  simulateResiduals1<- function (fittedModel, n = 250, refit = F, integerResponse = NULL, 
            plot = F, seed = 123, method = c("PIT", "traditional"), rotation = NULL, 
            ...) 
  {
    if (n < 2) 
      stop("error in DHARMa::simulateResiduals: n > 1 is required to calculate scaled residuals")
    checkModel(fittedModel)
    match.arg(method)
    randomState <- getRandomState(seed)
    on.exit({
      randomState$restoreCurrent()
    })
    ptm <- proc.time()
    out = list()
    family = family(fittedModel)
    out$fittedModel = fittedModel
    out$modelClass = class(fittedModel)[1]
    out$additionalParameters = list(...)
    out$nObs = nobs(fittedModel)
    out$nSim = n
    out$refit = refit
    out$observedResponse = getObservedResponse(fittedModel)
    if (out$nObs < length(out$observedResponse)) 
      stop("DHARMA::simulateResiduals: nobs(model) < nrow(model.frame). A possible reason is that you have observation with zero prior weights (ir binomial with n=0) in your data. Calculating residuals in this case wouldn't be sensible. Please remove zero-weight observations from your data and re-fit your model! If you believe that this is not the reason, please file an issue under https://github.com/florianhartig/DHARMa/issues")
    if (is.null(integerResponse)) {
      if (family$family %in% c("zeroinflated.poisson","binomial", "poisson", "quasibinomial", 
                               "quasipoisson", "Negative Binom", "nbinom2", "nbinom1", 
                               "genpois", "compois", "truncated_poisson", "truncated_nbinom2", 
                               "truncated_nbinom1", "betabinomial", "Poisson", "Tpoisson", 
                               "COMPoisson", "negbin", "Tnegbin") | grepl("Negative Binomial", 
                                                                          family$family)) 
        integerResponse = TRUE
      else integerResponse = FALSE
    }
    out$integerResponse = integerResponse
    out$problems = list()
    out$fittedPredictedResponse = getFitted(fittedModel)
    out$fittedFixedEffects = getFixedEffects(fittedModel)
    out$fittedResiduals = getResiduals(fittedModel)
    if (refit == FALSE) {
      out$method = method
      out$simulatedResponse = getSimulations(fittedModel, nsim = n, 
                                             type = "normal", ...)
      checkSimulations(out$simulatedResponse, out$nObs, out$nSim)
      out$scaledResiduals = getQuantile(simulations = out$simulatedResponse, 
                                        observed = out$observedResponse, integerResponse = integerResponse, 
                                        method = method, rotation = rotation)
    }
    else {
      out$method = "traditional"
      out$refittedPredictedResponse <- matrix(nrow = out$nObs, 
                                              ncol = n)
      out$refittedFixedEffects <- matrix(nrow = length(out$fittedFixedEffects), 
                                         ncol = n)
      out$refittedResiduals = matrix(nrow = out$nObs, ncol = n)
      out$refittedPearsonResiduals = matrix(nrow = out$nObs, 
                                            ncol = n)
      out$simulatedResponse = getSimulations(fittedModel, nsim = n, 
                                             type = "refit", ...)
      for (i in 1:n) {
        simObserved = out$simulatedResponse[[i]]
        try({
          refittedModel = getRefit(fittedModel, simObserved)
          out$refittedPredictedResponse[, i] = getFitted(refittedModel)
          out$refittedFixedEffects[, i] = getFixedEffects(refittedModel)
          out$refittedResiduals[, i] = getResiduals(refittedModel)
          out$refittedPearsonResiduals[, i] = residuals(refittedModel, 
                                                        type = "pearson")
        }, silent = TRUE)
      }
      if (anyNA(out$refittedResiduals)) 
        warning("DHARMa::simulateResiduals warning: on refit = TRUE, at least one of the refitted models produced an error. Inspect the refitted model values. Results may not be reliable.")
      dup = sum(duplicated(out$refittedFixedEffects, MARGIN = 2))
      if (dup > 0) {
        if (dup < n/3) {
          warning(paste("There were", dup, "of", n, "duplicate parameter estimates in the refitted models. This may hint towards a problem with optimizer convergence in the fitted models. Results may not be reliable. The suggested action is to not use the refitting procedure, and diagnose with tools available for the normal (not refitted) simulated residuals. If you absolutely require the refitting procedure, try changing tolerance / iterations in the optimizer settings."))
        }
        else {
          warning(paste("There were", dup, "of", n, "duplicate parameter estimates in the refitted models. This may hint towards a problem with optimizer convergence in the fitted models. Results are likely not reliable. The suggested action is to not use the refitting procedure, and diagnose with tools available for the normal (not refitted) simulated residuals. If you absolutely require the refitting procedure, try changing tolerance / iterations in the optimizer settings."))
          out$problems[[length(out$problems) + 1]] = "error in refit"
        }
      }
      out$scaledResiduals = getQuantile(simulations = out$refittedResiduals, 
                                        observed = out$fittedResiduals, integerResponse = integerResponse, 
                                        method = "traditional", rotation = rotation)
    }
    out$time = proc.time() - ptm
    out$randomState = randomState
    class(out) = "DHARMa"
    if (plot == TRUE) 
      plot(out)
    return(out)
  }
  <bytecode: 0x7fc40e5d7538>
    <environment: namespace:DHARMa>
    
    
options(repos = c(
      phyr = 'https://daijiang.r-universe.dev',
      CRAN = 'https://cloud.r-project.org'))
  install.packages('phyr')
    
    
    