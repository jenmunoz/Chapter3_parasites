#######################################################################################
### Chapter 3-parasites and flocking species Flocks in a community context                                         ###
### Part 1 Abundance                                                                             ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
### Last update: Oct 26 2022                                                     ###
################################################################################

# Loading packages --------------------------------------------------------
# libraries for easier manipulation of data
install.packages("tidyr") 
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


#Libraries for data
library(tidyverse)
library(tidyr)
library(plyr)
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

# Part 1 Data --------------------------------------------------------------------
ectoparasites_general<-read.csv("data/7.ectoparasite_raw_pres_abs_07222022.csv")
lice_only<-read.csv("data/7.ectoparsite_raw_lice_abundance_07222022.csv")
mites_only<-read.csv("data/7.ectoparsite_raw_mites_abundance_07222022.csv")
ticks_only<-read.csv("data/7.ectoparsite_raw_ticks_abundance_07222022.csv")

sociality<-

bird_taxonomy 

bird_traits

samples_elevation

phylogenetic_order<-read.csv("additional_info/jetz_information_species_manu.csv",header=TRUE, strip.white=TRUE, na.strings = c("NA","")) 
taxonomy_2021<-read.csv("additional_info/taxonomy_revision_2021.csv",header=TRUE, strip.white=TRUE, na.strings = c("NA",""))




# Part 2  data exploration  AND Visualisation and--------------------------------------------------------------------

#General ectoparasite data
str(ectoparasites_general)
ectoparasites<-ectoparasites_general %>% unite (col = "species_binomial", c('Host.genus','Host.species'), sep='_')
str(ectoparasites)
janitor::get_dupes(ectoparasites) # make sure there are not duplicatesduplicates
host_list<-distinct(ectoparasites, species_binomial , .keep_all = FALSE) # remove duplicated
host_list_general<-  host_list %>%  arrange(desc(species_binomial))

# general lice data
str(lice)
lice<-lice_only%>% unite (col = "species_binomial", c('Host.Genus','Host.Species'), sep='_')
janitor::get_dupes(lice) # make sure there are not duplicatesduplicates
host_list_lice<-distinct(lice, species_binomial , .keep_all = FALSE)
host_list_lice<-  host_list_lice %>%  arrange(desc(species_binomial))

# general mites data
str(mites_only)
mites<-mites_only%>% unite (col = "species_binomial", c('Host.genus','Host.species'), sep='_')
janitor::get_dupes(mites) # make sure there are not duplicatesduplicates
host_list_mites<-distinct(mites, species_binomial , .keep_all = FALSE)
host_list_mites<-host_list_mites %>%  arrange(desc(species_binomial))

# general ticks data
str(ticks_only)
ticks<-ticks_only%>% unite (col = "species_binomial", c('Host.Genus','Host.Species'), sep='_')
janitor::get_dupes(ticks) # make sure there are not duplicatesduplicates
host_list_ticks<-distinct(ticks, species_binomial , .keep_all = FALSE)
host_list_ticks<-host_list_ticks %>%  arrange(desc(species_binomial))


# #Taxonomy Matching ----------------------------------------------------
# Step 1  Taxonomy Remove typos and extra spaces -----------------------
# Remove typos that we found 

# But also there are some species that are not in the general host list, in which we need to reconciliate the taxonomy

dplyr::setdiff( host_list_lice,host_list_general)
dplyr::setdiff( host_list_mites,host_list_general)
dplyr::setdiff( host_list_ticks,host_list_general)

# on the left is how species are in the general list on the righ ow they are in the other list we need to reconciliate this taxonomy
#Phylloscartes_ophthalmicus= Pogonotriccus_ophthalmicus
#Hafferia_fortis= Percnostola_fortis
#Dixiphia_pipra = Pseudopipra pipra
#Chlorospingus_flavopectus = Chlorospingus_opthalmicus
#Cyanocompsa_cyanoides= Cyanoloxia_cyanoides

# Some typos in all the list  general, mites, lice

names(mites)

#mites<-mites %>% mutate(species_clean=species_binomial) #duplicate the species column
names(ectoparasites)
ectoparasites<-ectoparasites %>% mutate(species_clean=str_replace_all(species_binomial,c("Ceratopipra _chloromeros"="Ceratopipra_chloromeros",
                                                                                         "Chlorospingus_opthalmicus"="Chlorospingus_ophthalmicus",
                                                                                         "Galbula_cynescens"="Galbula_cyanescens",
                                                                                         "Mionectes _olivaceus"="Mionectes_olivaceus",
                                                                                         "Rhyncgocyclus_fulvipectus"="Rhynchocyclus_fulvipectus",
                                                                                         "Sclerurus_caudatus"="Sclerurus_caudacutus",
                                                                                         "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                                         "Syndactyla _ucayalae"="Syndactyla_ucayalae",
                                                                                         "XIphorhynchus_elegans"="Xiphorhynchus_elegans",
                                                                                         "Pogonotriccus_ophthalmicus"="Phylloscartes_ophthalmicus",
                                                                                         "Percnostola_fortis"="Hafferia_fortis",
                                                                                         "Pseudopipra_pipra"="Dixiphia_pipra",
                                                                                         "Chlorospingus_opthalmicus"="Chlorospingus_flavopectus",
                                                                                         "Cyanocompsa_cyanoides"="Cyanoloxia_cyanoides")))


                
mites<-mites %>% mutate(species_clean=str_replace_all(species_binomial,c("Ceratopipra _chloromeros"="Ceratopipra_chloromeros",
                                                                  "Chlorospingus_opthalmicus"="Chlorospingus_ophthalmicus",
                                                                  "Galbula_cynescens"="Galbula_cyanescens",
                                                                  "Mionectes _olivaceus"="Mionectes_olivaceus",
                                                                  "Rhyncgocyclus_fulvipectus"="Rhynchocyclus_fulvipectus",
                                                                  "Sclerurus_caudatus"="Sclerurus_caudacutus",
                                                                  "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                  "Syndactyla _ucayalae"="Syndactyla_ucayalae",
                                                                  "XIphorhynchus_elegans"="Xiphorhynchus_elegans",
                                                                  "Pogonotriccus_ophthalmicus"="Phylloscartes_ophthalmicus",
                                                                  "Percnostola_fortis"="Hafferia_fortis",
                                                                  "Pseudopipra_pipra"="Dixiphia_pipra",
                                                                  "Chlorospingus_opthalmicus"="Chlorospingus_flavopectus",
                                                                  "Cyanocompsa_cyanoides"="Cyanoloxia_cyanoides")))

names(lice)
lice<-lice %>% mutate(species_clean=str_replace_all(species_binomial,c("Ceratopipra _chloromeros"="Ceratopipra_chloromeros",
                                                                         "Chlorospingus_opthalmicus"="Chlorospingus_ophthalmicus",
                                                                         "Galbula_cynescens"="Galbula_cyanescens",
                                                                         "Mionectes _olivaceus"="Mionectes_olivaceus",
                                                                         "Rhyncgocyclus_fulvipectus"="Rhynchocyclus_fulvipectus",
                                                                         "Sclerurus_caudatus"="Sclerurus_caudacutus",
                                                                         "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                         "Syndactyla _ucayalae"="Syndactyla_ucayalae",
                                                                         "XIphorhynchus_elegans"="Xiphorhynchus_elegans",
                                                                       "Pogonotriccus_ophthalmicus"="Phylloscartes_ophthalmicus",
                                                                       "Percnostola_fortis"="Hafferia_fortis",
                                                                       "Pseudopipra_pipra"="Dixiphia_pipra",
                                                                       "Chlorospingus_opthalmicus"="Chlorospingus_flavopectus",
                                                                       "Cyanocompsa_cyanoides"="Cyanoloxia_cyanoides")))

names(ticks)
ticks<-ticks %>% mutate(species_clean=str_replace_all(species_binomial,c("Ceratopipra _chloromeros"="Ceratopipra_chloromeros",
                                                                       "Chlorospingus_opthalmicus"="Chlorospingus_ophthalmicus",
                                                                       "Galbula_cynescens"="Galbula_cyanescens",
                                                                       "Mionectes _olivaceus"="Mionectes_olivaceus",
                                                                       "Rhyncgocyclus_fulvipectus"="Rhynchocyclus_fulvipectus",
                                                                       "Sclerurus_caudatus"="Sclerurus_caudacutus",
                                                                       "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                       "Syndactyla _ucayalae"="Syndactyla_ucayalae",
                                                                       "XIphorhynchus_elegans"="Xiphorhynchus_elegans",
                                                                       "Pogonotriccus_ophthalmicus"="Phylloscartes_ophthalmicus",
                                                                       "Percnostola_fortis"="Hafferia_fortis",
                                                                       "Pseudopipra_pipra"="Dixiphia_pipra",
                                                                       "Chlorospingus_opthalmicus"="Chlorospingus_flavopectus",
                                                                       "Cyanocompsa_cyanoides"="Cyanoloxia_cyanoides")))


# Step 2  Reconciliate  reconciliate taxonomy with the current taxonomy of Manu database for traits and sociality -----------------------





We match the taxonomy with jetz taxonomy 
We update the taxonomy to the SACC 2022


 
# Part 3 Modeling abundance~Mice ------------------------------------------------------


# Part 3 Modeling abundance~Lice-------------------------------------------------------------------------
#generalized linear mixed model (GLMMM) :  non-normal data; #is an extension to the generalized linear model (GLM) 
#in which the linear predictor contains random effects in addition to the usual fixed effects
#use lmer() in the lme4 and lmerTest packages or lme() in the nlme package to analyze models containing random effects. T
#hese packages model the variance structure of random effects explicitly.

#pglmm: Phylogenetic Generalized Linear Mixed Model for Community...

# Question when including elevation, we could include elevation for the individual sample when available, but if doing it at the speies level do we include the mean elevation?
#mixed effect binomial models using glmer

# Abundance is counts so we can use  a poisson
#poisson error; fixed effect sociality=categories of variable of direct interest; random effect=foraging type # options calculating the mean abundance per species # or uisng species as a random effect
# species and elevation site has to be factors
#sociality 1, 0
#elevation as a site (as factor) 4 levels bamboo, lowlands, montane, high andes

pglmm(abundance~sociality + (1|species)+ (1|elevationsite)+ (1|foragingstratum)+ (sociality|sp__)) ,  family=poisson #sp

pglmm(abundance~sociality_degree + (1|species)+ (1|elevationsite)+ (1|foragingstratum)+ (sociality|sp__)) ,  family=poisson #sp # 



mod <- phyr::pglmm(ecto_abundance ~ sociality +
                     (1 | sp__) +  #overall phylogenetic effect using sp__, which also automatically includes a nonphylogenetic i.i.d. effect of species. 
                     (1 | siteelevation) + # random effect elevations by site as a factor
                     (sociality | sp__) + #We’ve also included a disturbance-by-phylogenetic species effect ((disturbance | sp__)), which estimates the degree to which social vs. non social  has a phylogenetic signal.Like the main sp__ effect, the sociality-by-sp__ effect also includes an nonphylogenetic species-by-disturbance interaction
                     data = oldfield$data, 
                   cov_ranef = list(sp = oldfield$phy), #we specified the phylogeny in the cov_ranef argument, giving it the name sp which matches sp__ but without the underscores
                   family = "binomial")


# prevalence is a 1 or 0 so we can use binomial # but elevation can not be continuos to be entered as a random effect logit 

mod <- phyr::pglmm(ecto_presence ~ sociality +
                     (1 | sp__) +  #overall phylogenetic effect using sp__, which also automatically includes a nonphylogenetic i.i.d. effect of species. 
                     (1 | siteelevation) + # random effect elevations by site as a factor
                     (sociality | sp__) + #We’ve also included a disturbance-by-phylogenetic species effect ((disturbance | sp__)), which estimates the degree to which social vs. non social  has a phylogenetic signal.Like the main sp__ effect, the sociality-by-sp__ effect also includes an nonphylogenetic species-by-disturbance interaction
                     data = oldfield$data, 
                   cov_ranef = list(sp = oldfield$phy), #we specified the phylogeny in the cov_ranef argument, giving it the name sp which matches sp__ but without the underscores
                   family = "binomial")





pglmm(pres_abs~sociality + (1|species)+ (1|elevationsite)+ (1|foragingstratum)+ (sociality|sp__)) ,  family=binomial /negative binomial? #sp

prevalence~sociality+elevation+(1|species)+(1|foraging)+ (sociality|sp__) #binomial error logit function, species as a random effect, foraging as a random effect # how sociality influences the presence/absence of parasites across the gradient
prevalence~sociality+(1|species)+(1|foraging)+ (1|elevationsite)+(sociality|sp__) #binomial error logit function, species as a random effect, foraging as a random effect # how sociality influences the presence/absence of parasites across the gradient



 glmer(cbind(round(seeds.orig)-round(seeds.intact),round(seeds.orig)) ~ latc +elev.km +seed.sp + (1|date) + (1|siteID),
                      family=binomial,
                      data=seedHav[seedHav$cage.treat=='CT',]) #no warnings

# prevalence is a 1 or 0 so we can use binomial
prevalence~sociality+elevation+ phylogeny  (1|species)+(1|foraging) #binomial error logit function, species as a random effect, foraging as a random effect # how sociality influences the presence/absence of parasites across the gradient

#Including phylogenies

# Model 2 (Eq. 2) example
z <- pglmm(freq ~ sp + X + (1|site) + (X|sp__), data = dat, family = "binomial",
           cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1)

z <- pglmm(freq ~  sociality + elevation+ (1|species) + (X|sp__), data = dat, family = "binomial",
           cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1)

z <- pglmm(freq ~ sp + sociality + (1|site) + (sociality|sp__), data = dat, family = "binomial",
           cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1)


# Part 3 Modeling abundance~Ticks-------------------------------------------------------------------------


# Example

install.packages("phyr")
install.packages("ape")
install.packages("dplyr")



library(phyr)
library(ape)
library(dplyr)

data("oldfield")
View (oldfield)

plot(oldfield$phy)

#With these data we are interested in asking whether there is phylogenetic structure in the distribution of these species, 
# as well as whether disturbance has any overall effects.

mod <- phyr::pglmm(pres ~ disturbance +
                     (1 | sp__) +  #overall phylogenetic effect using sp__, which also automatically includes a nonphylogenetic i.i.d. effect of species. 
                     (1 | site) + 
                     (disturbance | sp__) + #We’ve also included a disturbance-by-phylogenetic species effect ((disturbance | sp__)), which estimates the degree to which occurrence in disturbed vs. undisturbed habitat has a phylogenetic signal.Like the main sp__ effect, the disturbance-by-sp__ effect also includes an nonphylogenetic species-by-disturbance interaction
                     (1 | sp__@site), #This is a “nested” phylogenetic effect. This means that we are fitting an effect that had covariance proportional to the species phylogeny, but independently for each site (or “nested” in sites)
                   data = oldfield$data, 
                   cov_ranef = list(sp = oldfield$phy), #we specified the phylogeny in the cov_ranef argument, giving it the name sp which matches sp__ but without the underscores
                   family = "binomial")

summary(mod)

#The results of the random effects imply that the strongest effect is an overall nonphylogenetic species effect, 
#followed closely by a disturbance-by-species effect. This implies that species vary strongly in their occurrence
#in disturbed or undisturbed sites, that there is a clear community difference between these treatments. 
#The next strongest effect is the nested phylogenetic effect. But how can we know if this effect is strong enough to take seriously?
#Well one way to get an idea is to run a likelihood ratio test on the random effect. This can be achieved by using the pglmm_profile_LRT() function, at least for binomial models.

#The other result from this model is that there is a strong fixed effect of disturbance.
#In the context of a binomial multivariate model such as pglmm, this means there is an overall increase in the probability of occurrence in disturbed sites. 
#In other words, disturbed sites have a higher species richness at the site level (noting that expected alpha species richness of a site can be expressed as Gamma richness * E(prob_site(occurrence))).
