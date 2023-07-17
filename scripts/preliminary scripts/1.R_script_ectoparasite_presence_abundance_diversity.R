#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Taxonomy cleaning Abundance                                                                             ###
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
library(phytools)
library(metafor)
library (phangorn) # to reconstruct a maximum clade credibility tree

library(tidyverse)
library(skimr)

#Libraries ofr plots
library(gridExtra)
library(ggpubr)
library(grid)

# Part 1 Data --------------------------------------------------------------------
ectoparasites_general<-read.csv("data/data_raw/7.ectoparasite_raw_pres_abs_07222022.csv")
lice_only<-read.csv("data/data_raw/7.ectoparsite_raw_lice_abundance_07222022.csv")
unique(lice_only$total_lice)
mites_only<-read.csv("data/data_raw/7.ectoparsite_raw_mites_abundance_07222022.csv")
ticks_only<-read.csv("data/data_raw/7.ectoparsite_raw_ticks_abundance_07222022.csv")

flocks<-read.csv ("data/0.flocks_manu_complete_18052022.csv")
bird_traits_manu<-read.csv("data/4.df_traits_manu_birds.csv") # this is t he final file with traits of Manu  this files includes sociality data binomial 1_0 AND including diet and foraging from PCoA
jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv")
manu_detections<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022.csv") # with elevations

###_###_###_###_###_###_###_###_###_###_###_###_###_###_
#traits data selected 
df_traits_selected<-bird_traits_manu %>% 
  select(sociality,species,mass_tidy,ForStrat_ground,ForStrat_understory,ForStrat_midhigh,ForStrat_canopy,ForStrat_aerial)
###_###_###_###_###_###_###_###_###_###_###_###_###_###_


# Other data --------------------------------------------------------------


###_###_###_###_###_###_###_###_###_###_###_###_###_###_
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

#Including the number of detections
detection_flocks<-flocks %>% group_by(species_clean) %>% 
  summarize(detections_flocks=n()) %>% 
  arrange(desc(detections_flocks)) 
  
#Checking for species with very low ocurrence in flcoks 

View(flocks)
###_###_###_###_###_###_###_###_###_###_###_###_###_###_


phylogenetic_order<-read.csv("additional_info/jetz_information_species_manu.csv",header=TRUE, strip.white=TRUE, na.strings = c("NA","")) 
taxonomy_2021<-read.csv("additional_info/taxonomy_revision_2021.csv",header=TRUE, strip.white=TRUE, na.strings = c("NA",""))

####  elevation midpoint per species

manu_detections<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022.csv") # with elevations

elevation_midpoint<-manu_detections %>% group_by(species_taxonomy_SACC_2021) %>% 
  summarize(elevation_midpoint=median(elevation))

#write_csv(elevation_midpoint,"data/1.elevation_midpoint_manu_species.csv")

anti_join(ectoparasite_list_clean, df_traits_list, by=c("species_clean"="species")) # speceis that are in the ectoparasite list that do not have a matcj in b 



# Part 2  data restructuring and exploration--------------------------------------------------------------------

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

# Error_typo in ecto list name= corrected typo or manu list names


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
                                                                                         "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                                         "Chlorothraupis_turdina"="Schiffornis_turdina", # typo
                                                                                         "Pogonotriccus_opthalmicus"="Phylloscartes_ophthalmicus",
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
                                                                  "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                  "Chlorothraupis_turdina"="Schiffornis_turdina", # typo
                                                                  "Pogonotriccus_opthalmicus"="Phylloscartes_ophthalmicus",
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
                                                                       "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                       "Chlorothraupis_turdina"="Schiffornis_turdina", # typo
                                                                       "Pogonotriccus_opthalmicus"="Phylloscartes_ophthalmicus",
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
                                                                       "XIphorhynchus_elegans"="Xiphorhynchus_elegans", #typo
                                                                       "Sporathraupis _cyanocephala"="Sporathraupis_cyanocephala",
                                                                       "Chlorothraupis_turdina"="Schiffornis_turdina", # typo
                                                                       "Pogonotriccus_opthalmicus"="Phylloscartes_ophthalmicus",
                                                                       "Pogonotriccus_ophthalmicus"="Phylloscartes_ophthalmicus", #taxo matching
                                                                       "Percnostola_fortis"="Hafferia_fortis",#taxo matching
                                                                       "Pseudopipra_pipra"="Dixiphia_pipra",#taxo matching
                                                                       "Chlorospingus_opthalmicus"="Chlorospingus_flavopectus",#taxo matching
                                                                       "Cyanocompsa_cyanoides"="Cyanoloxia_cyanoides")))



# Make sure all listare consistent now 

dplyr::setdiff( ticks$species_clean,ectoparasites$species_clean)
dplyr::setdiff( lice$species_clean,lice$species_clean)
dplyr::setdiff( mites$species_clean,mites$species_clean)

#Remove the underscore

ectoparasites<-ectoparasites %>% mutate_at("species_clean", str_replace, "_", " ")
lice<-lice %>% mutate_at("species_clean", str_replace, "_", " ")
mites<-mites %>% mutate_at("species_clean", str_replace, "_", " ")
ticks<-ticks %>% mutate_at("species_clean", str_replace, "_", " ")

# Step 2  Reconciliate taxonomy with the current taxonomy of Manu database for traits and sociality imported from chapter 1  -----------------------

bird_traits_manu<-read.csv("data/4.df_traits_manu_birds.csv") # this is t he final file with traits of Manu including diet and foraging from PCoA 

###_###_###_###_###
#Warning make sure you use the column species and not "species checked" from the traits data base. Species is the 2021 taxonomy
###_###_###_###_###

df_traits_selected<-bird_traits_manu %>% 
  select(sociality,species,mass_tidy,ForStrat_ground,ForStrat_understory,ForStrat_midhigh,ForStrat_canopy,ForStrat_aerial)

# Validate that the data is align with taxonomy
#Make sure taxonomy is the same that in the Manu data base
df_traits_list<-df_traits_selected %>% distinct(species)
ectoparasite_list_clean<-ectoparasites %>% distinct(species_clean)
#check if when joining them we get all the 184 species ofr which we have ectots

#dplyr::intersect( df_traits_list$species,ectoparasite_list_clean$species_clean) # that is not the case we are recovering only 150 
# there are 33 species in which the taxonomy differs from Manu 2021

# we do na antijoin to recover the lsit of species that need taxonomy matching 

anti_join(ectoparasite_list_clean, df_traits_list, by=c("species_clean"="species")) # speceis that are in the ectoparasite list that do not have a matcj in b 

# some of the spcies recovered are form differenet communities and are of no interest because we do not have sociality information 
# but some of the species recovered are of importance
# Ectos list name= manu list names

#"Veniliornis affinis"="Dryobates affinis"
#"Frederickena unduligera"="Frederickena unduliger"
#"Hylophilus ochraceiceps"="Tunchiornis ochraceiceps"
#"Phaetornis superciliosus"="Phaethornis superciliosus"
#"Chlorospingus ophthalmicus"="Chlorospingus flavopectus"
#"Percnostola goeldii"="Akletos goeldii"

# other species 
#Thamnomanes caesius NEW
#Oneillornis lunulatus NEW
#Pithys albifrons NEW
#Gymnopithys leucaspis NEW
#Phlegopsis erythroptera NEW
#Tachyphonus surinamus NEW

ectoparasites_manu<-ectoparasites %>% mutate(species_clean=str_replace_all(species_clean,c("Veniliornis affinis"="Dryobates affinis",
                                                                                              "Frederickena unduligera"="Frederickena unduliger",
                                                                                              "Hylophilus ochraceiceps"="Tunchiornis ochraceiceps",
                                                                                              "Phaetornis superciliosus"="Phaethornis superciliosus",
                                                                                              "Chlorospingus ophthalmicus"="Chlorospingus flavopectus",
                                                                                              "Percnostola goeldii"="Akletos goeldii")))


mites_manu<-mites %>% mutate(species_clean=str_replace_all(species_clean,c("Veniliornis affinis"="Dryobates affinis",
                                                                              "Frederickena unduligera"="Frederickena unduliger",
                                                                              "Hylophilus ochraceiceps"="Tunchiornis ochraceiceps",
                                                                              "Phaetornis superciliosus"="Phaethornis superciliosus",
                                                                              "Chlorospingus ophthalmicus"="Chlorospingus flavopectus",
                                                                              "Percnostola goeldii"="Akletos goeldii")))

lice_manu<-lice %>% mutate(species_clean=str_replace_all(species_clean,c("Veniliornis affinis"="Dryobates affinis",
                                                                            "Frederickena unduligera"="Frederickena unduliger",
                                                                            "Hylophilus ochraceiceps"="Tunchiornis ochraceiceps",
                                                                            "Phaetornis superciliosus"="Phaethornis superciliosus",
                                                                            "Chlorospingus ophthalmicus"="Chlorospingus flavopectus",
                                                                            "Percnostola goeldii"="Akletos goeldii")))

ticks_manu<-ticks %>% mutate(species_clean=str_replace_all(species_clean,c("Veniliornis affinis"="Dryobates affinis",
                                                                              "Frederickena unduligera"="Frederickena unduliger",
                                                                              "Hylophilus ochraceiceps"="Tunchiornis ochraceiceps",
                                                                              "Phaetornis superciliosus"="Phaethornis superciliosus",
                                                                              "Chlorospingus ophthalmicus"="Chlorospingus flavopectus",
                                                                              "Percnostola goeldii"="Akletos goeldii")))

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# Double check that the species thata re not shared are species thata re not from manu
lice_list_clean<-lice_manu%>% distinct(species_clean)
anti_join(lice_list_clean, df_traits_list, by=c("species_clean"="species")) # speceis that are in the ectoparasite list that do not have a matcj in b 
mites_list_clean<-mites_manu%>% distinct(species_clean)
anti_join(mites_list_clean, df_traits_list, by=c("species_clean"="species")) # speceis that are in the ectoparasite list that do not have a matcj in b 
ticks_list_clean<-ticks_manu%>% distinct(species_clean)
anti_join(ticks_list_clean, df_traits_list, by=c("species_clean"="species")) # speceis that are in the ectoparasite list that do not have a matcj in b 

# Step 3  Merge ectoparasites data with traits data  ------------------------------------------------------------------------
# make sure the numbers of rows is euqal to the number of rows in eah parasite file 

ectoparasites_manu_traits<-inner_join(ectoparasites_manu,bird_traits_manu, by = c("species_clean"="species"))
lice_manu_traits<-inner_join(lice_manu,bird_traits_manu, by = c("species_clean"="species"))
mites_manu_traits<-inner_join(mites_manu,bird_traits_manu, by = c("species_clean"="species"))
ticks_manu_traits<-inner_join(ticks_manu,bird_traits_manu, by = c("species_clean"="species"))

# Step 4  Match the taxonomy with Jetz taxonomy  ------------------------------------------------------------------------
# We match the taxonomy with Jetz taxonomy that we aleady curated form manu data( see list of species that were manually matched)
jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv")
# Species that were manually matched 

#species_jetz	species_taxonomy_SACC_2021
#Aglaiocercus kingi	Aglaiocercus kingii
#Amazona mercenaria	Amazona mercenarius
#Aramides cajanea	Aramides cajaneus
#Clypicterus oseryi	Cacicus oseryi
#Strix albitarsis	Ciccaba albitarsis
#Dendroplex picus	Dendroplex picus
#Campylorhamphus pucherani	Drymotoxeres pucheranii
#Epinecrophylla erythrura	Epinecrophylla erythrura
#Frederickena unduligera	Frederickena unduliger
#Myrmotherula hauxwelli	Isleria hauxwelli
#Micrastur buckleyi	Micraster buckleyi
#Micrastur semitorquatus	Micraster semitorquatus
#Phaeothlypis fulvicauda	Myiothlypis fulvicauda
#Ochthoeca frontalis	Ochthoeca frontalis
#Ochthoeca pulchella	Ochthoeca pulchella
#Notiochelidon flavipes	Orochelidon flavipes
#Thryothorus genibarbis	Pheugopedius genibarbis
#Phylloscartes ophthalmicus	Phylloscartes ophthalmicus
#Phylloscartes orbitalis	Phylloscartes orbitalis
#Phylloscartes poecilotis	Phylloscartes poecilotis
#Picumnus subtilis	Picumnus subtilus
#Premnornis guttuligera	Premnornis guttuliger
#Pteroglossus beauharnaesii	Pteroglossus beauharnaisii
#Pygochelidon cyanoleuca	Pygochelidon cyanoleuca

#lets make sure things align create the list

ectoparasites_manu_traits_list<-ectoparasites_manu_traits %>% distinct(species_clean)
lice_manu_traits_list<-lice_manu_traits%>% distinct(species_clean)
mites_manu_traits_list<-mites_manu_traits%>% distinct(species_clean)
ticks_manu_traits_list<-ticks_manu_traits%>% distinct(species_clean)

jetz_taxonomy_manu_only_list<-jetz_taxonomy_manu_only%>% distinct(species_taxonomy_SACC_2021)

# Make sure the list  of traits "ectoparasites_manu_traits_list"  and jetz alings and non species is out of it 
anti_join(ectoparasites_manu_traits_list,jetz_taxonomy_manu_only_list, by=c("species_clean"="species_taxonomy_SACC_2021")) # species that are in the ectoparasite trait list that do not have a match on jetzz
anti_join(lice_manu_traits_list,jetz_taxonomy_manu_only_list, by=c("species_clean"="species_taxonomy_SACC_2021")) # species that are in the ectoparasite trait list that do not have a match on jetzz
anti_join(mites_manu_traits_list,jetz_taxonomy_manu_only_list, by=c("species_clean"="species_taxonomy_SACC_2021")) # species that are in the ectoparasite trait list that do not have a match on jetzz
anti_join(ticks_manu_traits_list,jetz_taxonomy_manu_only_list, by=c("species_clean"="species_taxonomy_SACC_2021")) # species that are in the ectoparasite trait list that do not have a match on jetzz

# Good databases are matching correctly, there are bot names in the ectos data sets that are not present on manu jetz

# Step 5  Generate the databases for analyses  ------------------------------------------------------------------------
#Select the pieces that we need from the ectoparasite traits and the jettz dataset

#ectoparasites_manu_traits
#View(jetz_taxonomy_manu_only)
#names(jetz_taxonomy_manu_only)
#names(ectoparasites_manu_traits)

#SELEC WHAT WE WANT FRM EACH DATA SET

jetz_tax_manu_selected<-jetz_taxonomy_manu_only%>% 
  select(species_taxonomy_SACC_2021,
         species_jetz,
         WJSpecID,
         TipLabel,
         PatchClade,
         Hackett_FineClades,
         Hackett_CoarseClades,
         English,
         FileName,
         Taxo,
         BLFamilyLatin,
         BLFamilyEnglish,
         FamSequID,
         IOCOrder,
         PassNonPass,
         OscSubOsc)

ectoparasites_manu_traits_selected<-ectoparasites_manu_traits %>% select(-species_checked,
                                                                         -Diet.Inv,
                                                                         -Diet.Vend,
                                                                         -Diet.Vect,
                                                                         -Diet.Vfish,
                                                                         -Diet.Vunk,
                                                                         -Diet.Scav,
                                                                         -Diet.Fruit,
                                                                         -Diet.Nect,
                                                                         -Diet.Seed,
                                                                         -Diet.PlantO,
                                                                         -ForStrat_watbelowsurf,
                                                                         -ForStrat_wataroundsurf,
                                                                         -ForStrat_aerial,
                                                                         -PelagicSpecialist,
                                                                         -A1.x,
                                                                         -A2.x,
                                                                         -A3.x,
                                                                         -A4.x,
                                                                         -A1.y,
                                                                         -A2.y,
                                                                         -A3.y,
                                                                         -A4.y)


lice_manu_traits_selected<-lice_manu_traits %>% select(-species_checked,
                                                                         -Diet.Inv,
                                                                         -Diet.Vend,
                                                                         -Diet.Vect,
                                                                         -Diet.Vfish,
                                                                         -Diet.Vunk,
                                                                         -Diet.Scav,
                                                                         -Diet.Fruit,
                                                                         -Diet.Nect,
                                                                         -Diet.Seed,
                                                                         -Diet.PlantO,
                                                                         -ForStrat_watbelowsurf,
                                                                         -ForStrat_wataroundsurf,
                                                                         -ForStrat_aerial,
                                                                         -PelagicSpecialist,
                                                                         -A1.x,
                                                                         -A2.x,
                                                                         -A3.x,
                                                                         -A4.x,
                                                                         -A1.y,
                                                                         -A2.y,
                                                                         -A3.y,
                                                                         -A4.y)

mites_manu_traits_selected<-mites_manu_traits %>% select(-species_checked,
                                                       -Diet.Inv,
                                                       -Diet.Vend,
                                                       -Diet.Vect,
                                                       -Diet.Vfish,
                                                       -Diet.Vunk,
                                                       -Diet.Scav,
                                                       -Diet.Fruit,
                                                       -Diet.Nect,
                                                       -Diet.Seed,
                                                       -Diet.PlantO,
                                                       -ForStrat_watbelowsurf,
                                                       -ForStrat_wataroundsurf,
                                                       -ForStrat_aerial,
                                                       -PelagicSpecialist,
                                                       -A1.x,
                                                       -A2.x,
                                                       -A3.x,
                                                       -A4.x,
                                                       -A1.y,
                                                       -A2.y,
                                                       -A3.y,
                                                       -A4.y)

ticks_manu_traits_selected<-ticks_manu_traits %>% select(-species_checked,
                                                         -Diet.Inv,
                                                         -Diet.Vend,
                                                         -Diet.Vect,
                                                         -Diet.Vfish,
                                                         -Diet.Vunk,
                                                         -Diet.Scav,
                                                         -Diet.Fruit,
                                                         -Diet.Nect,
                                                         -Diet.Seed,
                                                         -Diet.PlantO,
                                                         -ForStrat_watbelowsurf,
                                                         -ForStrat_wataroundsurf,
                                                         -ForStrat_aerial,
                                                         -PelagicSpecialist,
                                                         -A1.x,
                                                         -A2.x,
                                                         -A3.x,
                                                         -A4.x,
                                                         -A1.y,
                                                         -A2.y,
                                                         -A3.y,
                                                         -A4.y)


# Generate the data sets 
# Warning use the species clean and not the species binomial, and the species_taxonomy_Sac_2021 NO species jetz

ectoparasites_df_jetz<-inner_join(ectoparasites_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
lice_df_jetz<-inner_join(lice_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
mites_df_jetz<-inner_join(mites_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
ticks_df_jetz<-inner_join(ticks_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))

#We match the taxonomy with jetz taxonomy 

#We update the taxonomy to the SACC 2022


# Part 3  Structuring files as required for the models------------------------------------------------------------------

# create factors for the elevation and foraging strata
ectoparasites_df_jetz<-inner_join(ectoparasites_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
lice_df_jetz<-inner_join(lice_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
mites_df_jetz<-inner_join(mites_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
ticks_df_jetz<-inner_join(ticks_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))

# #Part 3.a  Presence absence data------------------------------------------------------------------
#ticks<-ticks %>% mutate_at("species_clean", str_replace, "_", " ")
###_###_###_###
#  presence absence analyses~ General ectoparasites
###_###_###_###
names(ectoparasites_df_jetz)

ectoparasites_df<-ectoparasites_df_jetz %>% mutate( foraging_cat = case_when(
                                                (ForStrat_ground==50&ForStrat_understory==50)~"ground-understory",
                                                (ForStrat_midhigh==50&ForStrat_understory==50)~"undersory-midhigh",
                                                (ForStrat_midhigh==50&ForStrat_canopy==50)~"midhigh-canopy",
                                                (ForStrat_midhigh>=30&ForStrat_understory>=30&ForStrat_ground>=30)~"ground-midhigh",
                                                ForStrat_ground>=50 ~ "ground", 
                                                ForStrat_understory>=50~"understory",
                                                ForStrat_midhigh>=50~"midhigh",
                                                ForStrat_canopy>=50~"canopy",
                                                TRUE~  "other"))%>% mutate(elevation_site=Full_Label) %>% 
  mutate_at("elevation_site", str_replace_all,(c(  "W"="high_montane_manu__",
                                                  "SP"="montane_manu__", "TL"="montane_manu__", 
                                                  "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                  "P"="lowland_manu__",
                                                  "M"="lowland_iquitos__",
                                                  "V"="lowland_iquitos__" ))) %>% 
              separate(elevation_site, c("elevation_cat","sample"),"__")

#View(ectoparasites_df)

ectoparasites_df$elevation_cat<-replace(x = ectoparasites_df$elevation_cat,
                   list =  !ectoparasites_df$elevation_cat %in% c('lowland_iquitos', 'lowland_manu', 'high_montane_manu','montane_manu'), 
                   values =  'other_iquitos')

unique(ectoparasites_df$elevation_cat)

# i have not figure it out this part so will change those number to Iquits manually later on


# We need the presence absence of ectosand sociality to be O and 1 
ectoparasites_df <-ectoparasites_df %>% mutate_at("sociality", str_replace, "yes", "1")
ectoparasites_df <-ectoparasites_df %>% mutate_at("sociality", str_replace, "no", "0")

ectoparasites_df <-ectoparasites_df %>% mutate_at("Lice", str_replace, "Yes", "1")
ectoparasites_df <-ectoparasites_df %>% mutate_at("Lice", str_replace, "No", "0")

ectoparasites_df <-ectoparasites_df %>% mutate_at("Mites", str_replace, "No", "0")
ectoparasites_df <-ectoparasites_df %>% mutate_at("Mites", str_replace, "no", "0")
ectoparasites_df <-ectoparasites_df %>% mutate_at("Mites", str_replace, "Yes", "1")
ectoparasites_df <-ectoparasites_df %>% mutate_at("Mites", str_replace, "Yes?", "1")

ectoparasites_df <-ectoparasites_df %>% mutate_at("Ticks", str_replace, "Yes", "1")
ectoparasites_df <-ectoparasites_df %>% mutate_at("Ticks", str_replace, "No", "0")

#unique(ectoparasites_df$Ticks)

# We need the random factors  to be factors (elevation and faraging stratum,and species )

ectoparasites_df$elevation_cat<-as.factor(ectoparasites_df$elevation_cat)
ectoparasites_df$foraging_cat<-as.factor(ectoparasites_df$foraging_cat)
ectoparasites_df$species_jetz<-as.factor(ectoparasites_df$species_jetz)

write_csv(ectoparasites_df, "data/7.ectoparasite_df_presence_absence.csv")
ectoparasites_df<-read.csv("data/7.ectoparasite_df_presence_absence.csv")
unique(ectoparasites_df$foraging_cat)




# #Part 3.b  Abundance data------------------------------------------------------------------
length(ectoparasites_df)
unique( ectoparasites_df$Full_Label)
###_###_###_###
# Abundance analyses_Lice
###_###_###_###

lice_df<-lice_df_jetz %>% mutate( foraging_cat = case_when(
  (ForStrat_ground==50&ForStrat_understory==50)~"ground-understory",
  (ForStrat_midhigh==50&ForStrat_understory==50)~"undersory-midhigh",
  (ForStrat_midhigh==50&ForStrat_canopy==50)~"midhigh-canopy",
  (ForStrat_midhigh>=30&ForStrat_understory>=30&ForStrat_ground>=30)~"ground-midhigh",
  ForStrat_ground>=50 ~ "ground", 
  ForStrat_understory>=50~"understory",
  ForStrat_midhigh>=50~"midhigh",
  ForStrat_canopy>=50~"canopy",
  TRUE~  "other"))%>% mutate(elevation_site=Full_Label) %>% 
  mutate_at("elevation_site", str_replace_all,(c(  "W"="high_montane_manu__",
                                                   "SP"="montane_manu__", "TL"="montane_manu__", 
                                                   "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                   "P"="lowland_manu__",
                                                   "M"="lowland_iquitos__",
                                                   "V"="lowland_iquitos__" ))) %>% 
  separate(elevation_site, c("elevation_cat","sample"),"__")

lice_df$elevation_cat<-replace(x = lice_df$elevation_cat,
                                        list =  !lice_df$elevation_cat %in% c('lowland_iquitos', 'lowland_manu', 'high_montane_manu','montane_manu'), 
                                        values =  'other_iquitos')

unique(lice_df$elevation_cat)
# We need the presence sociality to be O and 1 
lice_df <-lice_df %>% mutate_at("sociality", str_replace, "yes", "1")
lice_df <-lice_df %>% mutate_at("sociality", str_replace, "no", "0")

# We need the random factors  to be factors (elevation and faraging stratum,and species )

lice_df$elevation_cat<-as.factor(lice_df$elevation_cat)
lice_df$foraging_cat<-as.factor(lice_df$foraging_cat)
lice_df$species_jetz<-as.factor(lice_df$species_jetz)
lice_df$total_lice<-as.numeric(lice_df$total_lice)

unique (lice_df$total_lice)

str(lice_df)
###_###_###_###
write_csv(lice_df, "data/7.lice_df_abundance.csv")
lice_df<-read.csv( "data/7.lice_df_abundance.csv")

###_###_###_###


###_###_###_###
# Abundance analyses_mites Sample.Full..
###_###_###_###
names(mites_df_jetz)

mites_df<-mites_df_jetz %>% mutate( foraging_cat = case_when(
  (ForStrat_ground==50&ForStrat_understory==50)~"ground-understory",
  (ForStrat_midhigh==50&ForStrat_understory==50)~"undersory-midhigh",
  (ForStrat_midhigh==50&ForStrat_canopy==50)~"midhigh-canopy",
  (ForStrat_midhigh>=30&ForStrat_understory>=30&ForStrat_ground>=30)~"ground-midhigh",
  ForStrat_ground>=50 ~ "ground", 
  ForStrat_understory>=50~"understory",
  ForStrat_midhigh>=50~"midhigh",
  ForStrat_canopy>=50~"canopy",
  TRUE~  "other"))%>% mutate(elevation_site=Sample.Full..) %>% 
  mutate_at("elevation_site", str_replace_all,(c(  "W"="high_montane_manu__",
                                                   "SP"="montane_manu__", "TL"="montane_manu__", 
                                                   "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                   "P"="lowland_manu__",
                                                   "M"="lowland_iquitos__",
                                                   "V"="lowland_iquitos__" ))) %>% 
  separate(elevation_site, c("elevation_cat","sample"),"__")

mites_df$elevation_cat<-replace(x = mites_df$elevation_cat,
                               list =  !mites_df$elevation_cat %in% c('lowland_iquitos', 'lowland_manu', 'high_montane_manu','montane_manu'), 
                               values =  'other_iquitos')
unique(mites_df$elevation_cat)


# We need the presence sociality to be O and 1 
mites_df <-mites_df %>% mutate_at("sociality", str_replace, "yes", "1")
mites_df <-mites_df %>% mutate_at("sociality", str_replace, "no", "0")

# We need the random factors  to be factors (elevation and faraging stratum,and species )

mites_df$elevation_cat<-as.factor(mites_df$elevation_cat)
mites_df$foraging_cat<-as.factor(mites_df$foraging_cat)
mites_df$species_jetz<-as.factor(mites_df$species_jetz)
mites_df$total_mites<-as.numeric(mites_df$total_mites)

unique (mites_df$total_mites)

###_###_###_###
write_csv(mites_df, "data/7.mites_df_abundance.csv")
mites_df<-read.csv( "data/7.mites_df_abundance.csv")
str(mites_df)
###_###_###_###


###_###_###_###
# Abundance analyses_ticks
###_###_###_###

names(ticks_df_jetz)
ticks_df<-ticks_df_jetz %>% mutate( foraging_cat = case_when(
  (ForStrat_ground==50&ForStrat_understory==50)~"ground-understory",
  (ForStrat_midhigh==50&ForStrat_understory==50)~"undersory-midhigh",
  (ForStrat_midhigh==50&ForStrat_canopy==50)~"midhigh-canopy",
  (ForStrat_midhigh>=30&ForStrat_understory>=30&ForStrat_ground>=30)~"ground-midhigh",
  ForStrat_ground>=50 ~ "ground", 
  ForStrat_understory>=50~"understory",
  ForStrat_midhigh>=50~"midhigh",
  ForStrat_canopy>=50~"canopy",
  TRUE~  "other"))%>% mutate(elevation_site=Sample.Full..) %>% 
  mutate_at("elevation_site", str_replace_all,(c(  "W"="high_montane_manu__",
                                                   "SP"="montane_manu__", "TL"="montane_manu__", 
                                                   "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                   "P"="lowland_manu__",
                                                   "M"="lowland_iquitos__",
                                                   "V"="lowland_iquitos__" ))) %>% 
  separate(elevation_site, c("elevation_cat","sample"),"__")

ticks_df$elevation_cat<-replace(x = ticks_df$elevation_cat,
                                list =  !ticks_df$elevation_cat %in% c('lowland_iquitos', 'lowland_manu', 'high_montane_manu','montane_manu'), 
                                values =  'other_iquitos')

unique(ticks_df$elevation_cat)

# We need the presence sociality to be O and 1 
ticks_df <-ticks_df %>% mutate_at("sociality", str_replace, "yes", "1")
ticks_df <-ticks_df %>% mutate_at("sociality", str_replace, "no", "0")
# We need the random factors  to be factors (elevation and faraging stratum,and species )
ticks_df$elevation_cat<-as.factor(ticks_df$elevation_cat)
ticks_df$foraging_cat<-as.factor(ticks_df$foraging_cat)
ticks_df$species_jetz<-as.factor(ticks_df$species_jetz)
ticks_df$total_lice<-as.numeric(ticks_df$Total)

unique (ticks_df$Total)
###_###_###_###
write_csv(ticks_df, "data/7.ticks_df_abundance.csv")
ticks_df<-read.csv( "data/7.ticks_df_abundance.csv")

###_###_###_###

# #Part 3.c  Diversity Lice data------------------------------------------------------------------
### Lice 
#At the genus level
str(lice_df)

lice_df_wide<-lice_df %>% group_by(Host.Family,BLFamilyLatin,species_clean,species_jetz,Lice.Genus, Lice.Species, Lice.Genus2, Lice.Species2, Lice.Genus3, Lice.Species3,
                                   sociality,foraging_cat, elevation_cat,TipLabel, Full_Label ) %>% 
  summarise()

keycol <- "column"
valuecol <- "ectos_genus"
gathercols <- c("Lice.Genus", "Lice.Genus2", "Lice.Genus3")

lice_df_diversity_genus<-gather(lice_df_wide, keycol, valuecol, gathercols)%>% arrange(desc(species_jetz)) 

###_###_###_###
#write_csv(lice_df_diversity_genus,"data/7.lice_df_diversity_genus.csv")
###_###_###_###


# At the species level 


lice_df<-lice_df %>% unite (col = "ectos_binomial", c('Lice.Genus','Lice.Species'), sep='_')
lice_df<-lice_df %>% unite (col = "ectos_binomial2", c('Lice.Genus2','Lice.Species2'), sep='_')
lice_df<-lice_df %>% unite (col = "ectos_binomial3", c('Lice.Genus3','Lice.Species3'), sep='_')


lice_df_wide2<-lice_df %>% group_by(Host.Family,BLFamilyLatin,species_clean,species_jetz,ectos_binomial,ectos_binomial2, ectos_binomial3,
                                   foraging_cat, elevation_cat,TipLabel, Full_Label, sociality) %>% 
  summarise()

keycol <- "column"
valuecol <- "ectos_species"
gathercols <- c("ectos_binomial","ectos_binomial2", "ectos_binomial3")


lice_df_diversity_species<-gather(lice_df_wide2, keycol, valuecol, gathercols) %>% arrange(desc(species_jetz))

  ###_###_###_###
  write_csv(lice_df_diversity_species,"data/7.lice_df_diversity_species.csv")
  ###_###_###_###               
                
# # Part 3.c Diversity Mites data------------------------------------------------------------------
  
# Mites
  
  #At the group level
 # view(mites_df)
  
 names(mites_df)
  mites_df_wide<-mites_df%>% group_by(Host.Family,BLFamilyLatin,species_clean,species_jetz,Mite.Group, Mite.Genus, Mite.Group2, Mite.Genus2,                                    sociality,foraging_cat, elevation_cat,TipLabel, Sample.Full.. ) %>% dplyr::summarise()
  
  keycol <- "column"
  valuecol <- "mites_group"
  gathercols <- c("Mite.Group", "Mite.Group2")
  
  mites_df_diversity_group<-gather(mites_df_wide, keycol, valuecol, gathercols)%>% arrange(desc(species_jetz)) 
  
  #View(mites_df_diversity_group)
  ###_###_###_###
  write_csv(mites_df_diversity_group,"data/7.mites_df_diversity_group.csv")
  ###_###_###_###
  
  # At the genus level 
  
 # view(mites_df)
  names(mites_df)
  
  mites_df_wide2<-mites_df%>% group_by(Host.Family,BLFamilyLatin,species_clean,species_jetz,Mite.Group, Mite.Genus, Mite.Group2, Mite.Genus2,
                                      sociality,foraging_cat, elevation_cat,TipLabel, Sample.Full.. ) %>% dplyr::summarise()
  
  keycol <- "column"
  valuecol <- "mites_genus"
  gathercols <- c("Mite.Genus", "Mite.Genus2")
  
  mites_df_diversity_genus<-gather(mites_df_wide, keycol, valuecol, gathercols)%>% arrange(desc(species_jetz)) 
  
  #View(mites_df_diversity_genus)
  ###_###_###_###
  write_csv(mites_df_diversity_genus,"data/7.mites_df_diversity_genus.csv")
  ###_###_###_###
  # Diversity level at in ticks(we dont have enogth information)
  

# Part 4 Phylogenetic analyses ---------------------------------------------------

  # Master taxonomy form jetz
  
  taxonomy_jetz<-read.csv( "data/PhyloMasterTax_jetz.csv")
  ectoparasite_general_jetz<-read.csv("data/7.ectoparasite_df_presence_absence.csv") # data of general ectoparasites with jetz taxonomy
  
  names (taxonomy_jetz)
 # View(taxonomy_jetz)
  names (ectoparasite_general_jetz)
  #intersect(taxonomy_jetz, ectoparasite_general_jetz, by=c("Scientific"="species_jetz"))
  
  # Make sure there are not differences in the lsit of spcies with the master taxonomy from jetz
  anti_join(ectoparasite_general_jetz,taxonomy_jetz, by=c("species_jetz"="Scientific")) # speceis that are in the ectoparasite list that do not have a matcj in b 
  
  #Upload the list of species in bird.org
  
  # Step 1 Extractig 1000 trees from  birdtree.org-----------------------------------------------------------------------
  #Import phylogeny from birdtree.org
  #warning we will get the trees only for species for which we have some parasite samples( general parasite list)
  #Go to birdtree.org and get a Hackettbackbone tree of the species of interest.
  #filter species from the list then they will send you a output.nex element that you will need to use here
  
  # Read the tree
  # there are two options for reading the three
  #ape::read.nexus(filename, multiPhylo=TRUE)
  host_species_tree <- read.nexus("data/phylo_data/tree_pruner/output_bird_parasites.nex") #This is the tree with all species
  #host_species_tree <- read.nexus("data/phylo_data/tree_pruner/output_bird_parasites_genetic.nex")  #This is the tree with species with genetic info # I ma nost sure whe one tree has a 1 tip?
 # when using the one ony with genetic information we have to reduce the number of species with ectoparasites samples to match the ones for which genetci info is avialable # so i will use th eother
  is.rooted.multiPhylo(host_species_tree) # the trees are rooted
  is.ultrametric.multiPhylo(host_species_tree) # when using genetic data is not ultrametric
  print (host_species_tree, details=TRUE)
  #host_species_tree[9789]<-NULL # for some reason this genetic tree has less tips (weird I will just remove it)
  random_host_tree<-sample(host_species_tree,size=1)[[1]] # select one tree
  random_host_tree1<-sample( host_species_tree,size=1000)# select a thounsand trees tahta re rooted
  class(random_host_tree)
  is.rooted.multiPhylo(random_host_tree) # the trees are rooted
  is.ultrametric.multiPhylo(random_host_tree)
  edge.widthMap()
  
  example <- read.nexus("data/phylo_data/tree_pruner/output_example_genetic.nex") 
  is.ultrametric.multiPhylo(example)
  
  class( host_species_tree )# Must be multiPhylo
  
  #
  ##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__
  # Opt1 Generate a consensus tree -----------------------------------------------
  
  # Opt 1  Generate a consensus tree  with phytotools -------------------------------
  # no sure what is the difference between phytools and ape
  ###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
  # Methods to compute consensus edge lengths (branch lengths for a consensus topology()
  
  #Method 3: library phytotools Compute the non-negative least squares edge lengths on the consensus tree 
  #using the mean patristic distance matrix. (Function argument method="least.squares".) 
  #the phangorn function nnls.tree is used and the option rooted will be set to is.rooted(tree) for the consensus tree.
  #If the input trees are rooted & ultrametric, this can be used to produce a consensus tree that is also ultrametric.
  example_tree<-phytools::consensus.edges(example,method="least.squares")# generates a consensus tree with branch lenghts
  example_ape<- ape::consensus(example, p = 1, check.labels = TRUE, rooted = TRUE) # othier method for extracting the consensus tree . p=0.5 to get the majority-rule concesnsus tree
  is.rooted.phylo(example_ape) # Is this one rotted because we have the genetic data only???
  
  # consensus edges function computes for a tree under some criterion
  #If the input trees are rooted & ultrametric, this can be used to produce a consensus tree that is also ultrametric.
  host_consensus_tree<-phytools::consensus.edges(host_species_tree,method="least.squares")# generates a consensus tree with branch lenghts
  is.rooted.phylo(host_consensus_tree) # why the consensus tree is not rooted?
  is.ultrametric.phylo(host_consensus_tree)# but it is ultrametric
  # generatinga consensus tree with ape
  host_consensus_tree_ape<- ape::consensus(random_host_tree1, p =1, check.labels = TRUE, rooted = TRUE) # othier method for extracting the consensus tree . p=0.5 to get the majority-rule concesnsus tree
  # wHY Is the consensus tree not rooted
  is.rooted.phylo(host_consensus_tree_ape) # I am not sure why the consesus tree is not rooted i n the original wer rooted, quesion for julie!!!
  
  #Branch lenght 
  host_species_tree_rooted<-ape::compute.brlen(  host_consensus_tree_ape,method = "Grafen", power = 1)
  class(  branches)
  is.rooted.phylo( host_species_tree_rooted)
  
  #plots
  plotTree(host_consensus_tree,fsize=0.01, type="fan")
  plot(host_consensus_tree,type="fan",no.margin = TRUE, cex=0.5)
  host_consensus_tree
  
  # save phylogeny
  saveRDS(host_consensus_tree, file='data/phylo_data/1_host_consensus_tree_Manuspecies.rds')
  str(birdtree)
  
  ## Save tree
  write.tree(host_consensus_tree, file="data/phylo_data/1_host_consensus_tree_Manuspecies.tre")
  write.nexus(host_consensus_tree, file="data/phylo_data/1_host_consensus_tree_Manuspecies.nex")
  write.tree(host_consensus_tree, file="data/phylo_data/1_host_consensus_tree_Manuspecies.txt") 
  
  write.tree(random_host_tree,file="data/phylo_data/1_host_consensus_tree_Manuspecies_onerooted.tre") # this is rooted and hs brach lenghts 
  write.tree(random_host_tree,file="data/phylo_data/1_host_consensus_tree_Manuspecies_onerooted.nex") # this is rooted and hs brach lenghts 
  
  # NOTE IMPORTANT: for some reason my consensus tree is not rooted, and does not have brancj leght so for practicity I extracted one tree that has brach lenghts nad is rooted to run the analysis until I figure this out   

  
# Part 5 [ANALISES] #####Modeling###### ectoparasites presence absence -------------------------
  # To underestand teh models and out puts better see this paper: https://www.iecolab.org/wp-content/uploads/2020/10/phyr_2020.pdf
  # also see this blog with an example https://daijiang.github.io/phyr/articles/phyr_example_empirical.html
  #All models fitted with pglmm() have class of communityPGLMM. Here is a list of functions that can be used to these models.
  
 # pglmm_matrix_structure(): produce the whole covariance matrix
  #pglmm_plot_re(): plot images of random term matrix, see vignettes plot-re
  #pglmm_predicted_values() or fitted(): extract fitted values
  #pglmm_profile_LRT(): to test significance of random terms; only works with binomial models
  #plot_data(): plot data used (and optionally predicted) by fitted model
  #plot_bayes(): plot posterior distributions of random and fixed effects
 # summary(): summary of model fit
  #print(): summary of model fit
  #residuals(): residuals values
  #fixef(): estimates of fixed effects
  #ranef(): estimates of random terms (variance and standard deviation)
  
  library(phyr)
  
# prevalence is a 1 or 0 so we can use binomial # but elevation can not be continuos to be entered as a random effect logit 
 
  ###_###_###_###_##
  # The data
  ###_###_###_###_##
ectoparasites_df<-read.csv("data/7.ectoparasite_df_presence_absence.csv") # data on presence absence
names( ectoparasites_df)
unique(ectoparasites_df$Mites)
phylogeny<- read.nexus("data/phylo_data/1_host_consensus_tree_Manuspecies.nex") 

phylogeny_rooted<- read.tree("data/phylo_data/1_host_consensus_tree_Manuspecies_onerooted.tre")   # include all species not only genetic is roooted nad ultrametric

View(phylogeny)
class(phylogeny)
  
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
  
  # WARNING STILL NEED TO INCLUDE SAMPLE SIZE IN THE ANALYSES AS A RANDOM EFFECT but maybe not here because it is ndividual samples
  # this will be required in the diversity analyses
  ectoparasites_df<- ectoparasites_df %>% mutate(ectoparasites_PA=Lice+Mites+Ticks)
  ectoparasites_df$ectoparasites_PA[ectoparasites_df$ectoparasites_PA>=1]<-1   # convert the numerical values that we have without core to lowland iquitos
  unique(ectoparasites_df$ectoparasites_PA)
  
  names(ectoparasites_df)
  
 
  
  
  
  ecto_PA <-  phyr::pglmm(ectoparasites_PA ~ sociality+ (1|foraging_cat)+ (1|elevation_cat)+(1|species_jetz__), 
                                    data = ectoparasites_df, 
                                    family = "binomial",
                                    cov_ranef = list(species_jetz= phylogeny_rooted), #class phylo
                                    #bayes = TRUE,
                                    REML = TRUE, 
                                    verbose = TRUE,
                                    s2.init = .25) # what is this last parameter for
  
  binaryPGL
  
  summary(ecto_PA)
 predict(ecto_PA)
 
 
 
 ecto_PA_model_glmm<-  lme4::glmer(ectoparasites_PA ~ sociality+ (1|foraging_cat)+ (1|elevation_cat)+(1|species_jetz), 
                               data = ectoparasites_df, 
                               family = "binomial",
                               #bayes = TRUE,
                               verbose = TRUE)
 
 
 install.packages("MCMCglmm")
 library(MCMCglmm)
 
 # GOODNESS OF FIT  
 # Read this paper to underestand better https://academic.oup.com/sysbio/article/68/2/234/5098616?login=false
 # An additional complication is created by model with random effects. Given that random effects are very flexible model components 
 #(for example, nothing stops you from fitting a random effect for each observation in your dataset), a straight-up calculation of variance explains isn’t meaningful. T
 That said, methods that can produce a useful R2 metric in the complex situation have been developed. 
 #The package rr2 is able to calculate several flavors of R2, and supports phyr’s pglmm model object. Let’s try it!
 #
 install.packages("rr2")
 library(rr2)
 rr2::R2(ecto_PA)
 rr2::R2(ecto_PA_model_glmm)
 rr2::R2(MCMC)
 
 
  
  names( lice_diversity)
  
  
  
  binaryPGLMM(Y ~ X1, phy=phy, data=sim.dat)
  
  #bayes = TRUE
  # Also exploring with a glmm ( not phylogenetically corrected)
  ecto_PA_glmm<-lme4::glmer (ectoparasites_PA~sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz), 
                     data = ectoparasites_df, 
                     family ="binomial")
                     #bayes = TRUE)
  summary(ecto_PA_glmm)
  

  #The results of the random effects imply that the strongest effect ( the variables withthe higher variance) is an overall nonphylogenetic and phylogenetic species effect,
  #This implies that species vary strongly in their presence absence of parasites.
  
  #The other result from this model is that there is not a strong fixed effect of  sociality ( as indicated by teh significance * p-value? on the model summary)
  #In the context of a binomial multivariate model such as pglmm, this means there is NOT an overall increase in the probability of occurrence of parasite in social MSF birds. 

  # lice only
  
  l <-  phyr::pglmm(Lice ~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                    data = ectoparasites_df, 
                    family = "binomial",
                    cov_ranef = list(species_jetz= phylogeny), #class phylo
                    #bayes = TRUE,
                    REML = TRUE, 
                    verbose = TRUE, 
                    s2.init = .25)
  summary( l )
  
  rr2::R2(l)
  
  
  # mites only
  
  
m <-  phyr::pglmm(Mites ~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                    data = ectoparasites_df, 
                    family = "binomial",
                    cov_ranef = list(species_jetz= phylogeny), #class phylo
                    #bayes = TRUE,
                    REML = TRUE, 
                    verbose = TRUE, 
                    s2.init =.25)
  summary( m)
  rr2::R2(m)
  
  
  
  # ticks only

  t <-  phyr::pglmm(Ticks ~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                    data = ectoparasites_df, 
                    family = "binomial",
                    cov_ranef = list(species_jetz= phylogeny), #class phylo
                    #bayes = TRUE,
                    REML = TRUE, 
                    verbose = TRUE,
                    s2.init = .25)                   # an array of initial estimates of s2 for each random effect that scales the variance. # If s2.init is not provided for family="binomial", these are set to 0.25.

  summary( t)
  
  
# Part 6 [ANALISES] Modeling ~Lice Abundance-------------------------------------------------------------------------
  #generalized linear mixed model (GLMMM) :  non-normal data; #is an extension to the generalized linear model (GLM) 
  #in which the linear predictor contains random effects in addition to the usual fixed effects
  #use lmer() in the lme4 and lmerTest packages or lme() in the nlme package to analyze models containing random effects. T
  #hese packages model the variance structure of random effects explicitly.
   
  #pglmm: Phylogenetic Generalized Linear Mixed Model for Community...
  
  # Abundance is counts so we can use  a poisson but it is zero infladed 
  #poisson error; fixed effect sociality=categories of variable of direct interest; random effect=foraging type # options calculating the mean abundance per species # or uisng species as a random effect
  # species and elevation site has to be factors
  #sociality 1, 0
  #elevation as a site (as factor) several levels bamboo, lowlands, montane, high andes
  
  # Lets check for zero infalted data 
  
  100*sum(lice_df_abundance$total_lice== 0)/nrow(lice_df_abundance)
  
  # 42 % of our data is zeros( truth zeros)? i guess yes cause we collected the sample for teh individual

  #### Abundance
  # Modeling the individual abundances
  
  lice_df_abundance<-read.csv("data/7.lice_df_abundance.csv")
  length( lice_df_abundance)
  names(lice_df_abundance)
  phylogeny_for_lice<- read.nexus("data/phylo_data/1_host_consensus_tree_lice.nex")
  
  #phylogeny_for_lice<-read.tree("data/phylo_data/1_host_consensus_tree_lice.tre")

  # Filter only manu species 
  lice_df_abundance <-lice_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
  unique(lice_df_abundance$elevation_cat)
 
  # Make sure variables are in teh right format, random effects should be factors
  #We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?
  lice_df_abundance <-lice_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
  lice_df_abundance$elevation_cat<-as.factor(lice_df_abundance$elevation_cat)
  lice_df_abundance$foraging_cat<-as.factor(lice_df_abundance$foraging_cat)
  lice_df_abundance$species_jetz<-as.factor(lice_df_abundance$species_jetz)
  lice_df_abundance$sociality<-as.factor(lice_df_abundance$sociality)
  
  mean(lice_df_abundance$total_lice)
  sd(lice_df_abundance$total_lice)
  
str(lice_df_abundance)  

lice_df_abundance1<-lice_df_abundance %>% filter(total_lice!=0) %>% view()

unique(lice_df_abundance1$total_lice)


  # Modeling the data # I would prefer to use a zero inflated model however that is only aviallable in a gassioan approach bt that does no work with my model ( not sure why ye)
  

  l_a<-  phyr::pglmm(total_lice~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                    data = lice_df_abundance1, 
                    family ="poisson" , # use when bayes=true "zeroinflated.poisson",
                    cov_ranef = list(species_jetz=phylogeny_for_lice), #class phylo
                    #bayes = TRUE,
                    REML = TRUE, 
                    verbose = TRUE, 
                    s2.init = .25)
  summary( l_a )
  
  #pglmm_compare()  #It simultaneously estimates the strength of phylogenetic signal in the residuals and gives an approximate conditional likelihood ratio test for the hypothesis that there is no signal.
  
coef(l_a)
  # Modeling the mean abundances 
  
mean_lice<-lice_df_abundance %>% group_by (species_jetz) %>% 
  summarize(mean_lice=mean(total_lice))

species_atributes<-lice_df_abundance %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean)
species_attributes_distict<-distinct( species_atributes)

mean_lice_abundance<-right_join(species_attributes_distict, mean_lice, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 

###_###_####_###_
write_csv(mean_lice_abundance,"data/7.lice_df_abundance_means.csv")
###_###_####_###_
l_a_mean<-  phyr::pglmm(mean_lice~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                   data = mean_lice_abundance, 
                   family ="poisson", # use when bayes=true "zeroinflated.poisson",
                   cov_ranef = list(species_jetz=phylogeny_for_lice), #class phylo
                   #bayes = TRUE,
                   REML = TRUE, 
                   verbose = TRUE, 
                   s2.init = .25)

summary( l_a_mean)
### Trying with a glmm intend without the phylogenetic correction
  l_a_glmm<-  lme4::glmer(total_lice~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz), 
                     data = lice_df_abundance, 
                     family=poisson (link = "log"), # use when bayes=true "zeroinflated.poisson",
                     #REML = TRUE, 
                     verbose = TRUE)
  
  Anova(l_a_glmm, type=3)
  
  summary(l_a_glmm)
  
  anova(l_a, l_a_glmm) # not applicabel to PGLM 
  
#### Modeling infestation   
### Creating an infestation variable
  
lice_df_abundance<-lice_df_abundance %>% mutate( infestation_lice = case_when(
total_lice>20 ~ "1", 
TRUE~  "0"))


#write_csv(lice_df_abundance,"data/7.lice_df_abundance.csv")
  
  
names(lice_df_abundance)
  lice_df_abundance <-lice_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
  lice_df_abundance$elevation_cat<-as.factor(lice_df_abundance$elevation_cat)
  lice_df_abundance$foraging_cat<-as.factor(lice_df_abundance$foraging_cat)
  lice_df_abundance$species_jetz<-as.factor(lice_df_abundance$species_jetz)
  lice_df_abundance$sociality<-as.factor(lice_df_abundance$sociality)
  lice_df_abundance$infestation_lice<-as.numeric(lice_df_abundance$infestation_lice)
  
  
  unique(lice_df_abundance$infestation_lice)
  
  inf_l_model <-  phyr::pglmm(infestation_lice ~ sociality+ (1|foraging_cat)+(1|elevation_cat)+(1|species_jetz__), 
                          data = lice_df_abundance, 
                          family = "binomial",
                          cov_ranef = list(species_jetz= phylogeny_for_lice), #class phylo
                          #bayes = TRUE,
                          REML = TRUE, 
                          verbose = TRUE,
                          s2.init = .25) # what is this last parameter for
  
  summary(inf_l_model)
  
  class(inf_l_model)
  print(inf_l_model)
  pglmm_profile_LRT(inf_l_model)
  ranef(inf_l_model)
  
  
  # modeling without the phylogenetic component
  
  inf_l_model_glmm<-  lme4::glmer(infestation_lice~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz), 
                          data = lice_df_abundance, 
                          family=binomial, # use when bayes=true "zeroinflated.poisson",
                          #REML = TRUE, 
                          verbose = TRUE)
  summary(inf_l_model_glmm)
  

# Part 7 [ANALISES] Modeling abundance~Mite ------------------------------------------------------
  
 mites_df_abundance<-read.csv("data/data_analyses/7.dff_mites_abundance.csv",na.strings=c(""))
  names(mites_df_abundance)
  phylogeny_abundance_mites<- read.nexus("data/phylo_data/1_host_consensus_tree_mites.nex")
  phylogeny_abundance_mites<- read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_mites_abundance.nex")
  # Abundance is counts so we can use  a poisson but it is zero infladed (no opticon in PGLMM to take care of this)
  #poisson error; fixed effect sociality=categories of variable of direct interest; random effect=foraging type
  # species and elevation site has to be factors
  #sociality 1, 0
  #elevation as a site (as factor) several levels bamboo, lowlands, montane, high andes
  
  # Filter only manu species 
  
  mites_df_abundance <-mites_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")
  unique(mites_df_abundance$elevation_cat)
  
  # Lets check for zero infalted data 
  
  100*sum(mites_df_abundance$total_mites== 0)/nrow(mites_df_abundance)
  
  # 21 % of our data is zeros( truth zeros)? i guess yes cause we collected the sample for teh individual
  
  #### Abundance
  # Modeling the individual abundances
  
  #phylogeny_for_lice<-read.tree("data/phylo_data/1_host_consensus_tree_lice.tre")
  # Make sure variables are in teh right format, random effects should be factors
  #We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?
  
 mites_df_abundance <-mites_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
 mites_df_abundance$elevation_cat<-as.factor(mites_df_abundance$elevation_cat)
 mites_df_abundance$foraging_cat<-as.factor(mites_df_abundance$foraging_cat)
 mites_df_abundance$species_jetz<-as.factor(mites_df_abundance$species_jetz)
 mites_df_abundance$sociality<-as.factor(mites_df_abundance$sociality)
  
  mean(mites_df_abundance$total_mites)
  sd(mites_df_abundance$total_mites) # overdispersed variance> mean
  
  # Modeling the data # I would prefer to use a zero inflated model however that is only aviallable in a gassioan approach bt that does no work with my model ( not sure why ye)
  
  names( mites_df_abundance)
  m_a<-  phyr::pglmm(total_mites~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                     data = mites_df_abundance, 
                     family ="poisson", # use when bayes=true "zeroinflated.poisson",
                     cov_ranef = list(species_jetz=phylogeny_for_mites), #class phylo
                     #bayes = TRUE,
                     REML = TRUE, 
                     verbose = TRUE, 
                     s2.init = .25)
  summary(m_a)  
  
  # Modeling the mean abundances 
  
   mean_mites<-mites_df_abundance %>% group_by (species_jetz) %>% 
    summarize(mean_mites=mean(total_mites))
  
  species_atributes<-mites_df_abundance %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean)
  species_attributes_distict<-distinct( species_atributes)
  
  mean_mites_abundance<-right_join(species_attributes_distict, mean_mites, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 
  
  ###_###_####_###_
  write_csv(mean_mites_abundance,"data/7.mites_df_abundance_means.csv",header=TRUE,na.strings=c(""))
  ###_###_####_###_
  

  m_a_mean<-  phyr::pglmm(mean_mites~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__), 
                          data = mean_mites_abundance, 
                          family ="poisson", # use when bayes=true "zeroinflated.poisson",
                          cov_ranef = list(species_jetz=phylogeny_for_mites), #class phylo
                          #bayes = TRUE,
                          REML = TRUE, 
                          verbose = TRUE, 
                          s2.init = .25)
  
  summary( m_a_mean)
  
  ###_###_###_###_###_###_###_###
  # Modeling no feather mites only (cause feathers mites can be misleding)
  ###_###_###_###_###_###_###_###
  
  
  mites_df_abundance<-read.csv("data/7.mites_df_abundance.csv")
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
  
  #### mean non feather mites 
  
  mean_mites_mesostigmatidae<-mites_df_abundance %>% group_by (species_jetz) %>% 
    summarize(mean_mites=mean(total_mesostigmatidae))
  
  species_atributes<-mites_df_abundance %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean)
  species_attributes_distict<-distinct( species_atributes)
  
  mean_mites_abundance_mesostigmatidae<-right_join(species_attributes_distict, mean_mites_mesostigmatidae, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 
  
  
  
# Part 8 [ANALISES] Modeling abundance~Ticks-------------------------------------------------------------------------


#
# Note for underestanding the results

#p values for the fixed effects are given by a Wald test and for the random effects by profile likelihood, 
#al- though we recommend bootstrap-based tests when computation- ally feasible.

phyr::pglmm(Lice ~ sociality + (1 | Scientific__) +   (1 | elevation_cat) + (1|foraging_cat), 
            data=ectoparasites_df,
            family = "binomial",
            cov_ranef=phylogeny)

phyr::pglmm(Lice ~ sociality + (1 | species_jetz__) +   (1 | elevation_cat) + (sociality | sp__), data=ectoparasites_df,family = "binomial",
            cov_ranef=phylogeny)

ectoparasites_df$elevation_cat<-as.factor(ectoparasites_df$elevation_cat)
ectoparasites_df$foraging_cat<-as.factor(ectoparasites_df$foraging_cat)
ectoparasites_df$species_jetz<-as.factor(ectoparasites_df$species_jetz)

write_csv(ectoparasites_df, "data/7.ectoparasite_df_presence_absence.csv")
presence_absence

mod <- phyr::pglmm(Lice ~ sociality +
                     (1 | species_jetz__) +  #overall phylogenetic effect using sp__, which also automatically includes a nonphylogenetic i.i.d. effect of species. 
                     (1 | elevation_cat) + # random effect elevations by site as a factor
                     (sociality | sp__) + #We’ve also included a disturbance-by-phylogenetic species effect ((disturbance | sp__)), which estimates the degree to which social vs. non social  has a phylogenetic signal.Like the main sp__ effect, the sociality-by-sp__ effect also includes an nonphylogenetic species-by-disturbance interaction
                     data = presence_absence, 
                   cov_ranef =phylogeny, #we specified the phylogeny in the cov_ranef argument, giving it the name sp which matches sp__ but without the underscores
                   family = "binomial")


glmer(cbind(round(seeds.orig)-round(seeds.intact),round(seeds.orig)) ~ latc +elev.km +seed.sp + (1|date) + (1|siteID),
      family=binomial,
      data=seedHav[seedHav$cage.treat=='CT',]) #no warnings




# Part 9 [ANALISES] Modeling diversity~Lice --------------------------------------------------------------
#WARNING IN THIS ANALYSES WE ARE NOT INTERESTED IN THE ZEROS (BECAUSE HAVING A DIVERSITY OF ZERO DOES NOT HAVE ECOLOGICL MEANING FOR US, SO WE CAN DROP THE ZEROS OR NAS)

# at the genus level 
lice_df_diversity<-read.csv("data/7.lice_df_diversity_genus.csv",header=TRUE,na.strings=c(""))
names(lice_df_diversity)
phylogeny_for_lice_diversity<- read.nexus("data/phylo_data/1_host_tree_lice_diversity.nex")

# Filtering only manu data 

lice_df_diversity <-lice_df_diversity %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

## Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?
lice_df_diversity <-lice_df_diversity  %>% mutate_at("species_jetz", str_replace, " ", "_")
lice_df_diversity$elevation_cat<-as.factor(lice_df_diversity$elevation_cat)
lice_df_diversity$foraging_cat<-as.factor(lice_df_diversity$foraging_cat)
lice_df_diversity$species_jetz<-as.factor(lice_df_diversity$species_jetz)
lice_df_diversity$sociality<-as.factor(lice_df_diversity$sociality)

# REMOVE EMPTY ROWS (for some samples we do not have Lice genus2 and tree but when summarizing the table the field were created)
lice_df_diversity<-lice_df_diversity %>% drop_na(valuecol)
str(lice_df_diversity)

#richness at the genus level by host species
lice_richness_genus<-lice_df_diversity %>% group_by(species_jetz) %>% 
  summarise(richness=n_distinct(valuecol))

species_atributes<-lice_df_diversity %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean)
species_attributes_distict<-distinct( species_atributes)

lice_df_richness_genus<-right_join(species_attributes_distict, lice_richness_genus, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 

#lice_df_richness<-lice_df_diversity %>% group_by(species_jetz) %>% 
  #summarise(richness=n_distinct(valuecol), across()) 
#
#richness per sample

lice_df_diversity %>% group_by(species_jetz, Full_Label) %>% 
  summarise(richness=n_distinct(valuecol)) %>% 
  View()

#lice_df_richness_summary<-lice_df_richness %>% distinct(species_jetz,richness, elevation_cat, foraging_cat, sociality )

mean(lice_df_diversity$total_lice)
sd(lice_df_diversity$total_lice)

# We are worried that diversity increases with sample size, lets check is that is the case
# summarize the number of samples per species  
n_samples_lice<-lice_df_diversity %>% group_by(species_jetz) %>% 
  summarise(n_samples_lice=n_distinct(Full_Label))

is.na(lice_df_diversity$valuecol)

####Correlation sample size and diversity #WARNING

lice_diversity_samplesize<-full_join(lice_df_richness_genus, n_samples_lice, by="species_jetz")

plot(richness~n_samples_lice, data=lice_diversity_samplesize)

model<-lm(richness~n_samples_lice, data=lice_diversity_samplesize)
abline(model)

# let's include the sample size as a random effect? # is this a correct approach, or how do we correct for sample size 

lice_df_richness_genera<-right_join(lice_df_richness_genus, n_samples_lice, by="species_jetz")


#Modeling diversity without phylogeny

l_r_model_glmm<-  lme4::glmer(richness~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz) +(1|n_samples_lice), 
                                data = lice_df_richness_genera, 
                                family=poisson(), # use when bayes=true "zeroinflated.poisson",
                                #REML = TRUE, 
                                verbose = TRUE)
summary(l_r_model_glmm)
list_host_lice_diversity<-as.data.frame(unique(lice_df_richness_summary$species_jetz))


# Model
names( lice_abundance)
l_r<-  phyr::pglmm(richness~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__)+(1|n_samples_lice), 
                   data = lice_df_richness_genera, 
                   family ="poisson", # use when bayes=true "zeroinflated.poisson",
                   cov_ranef = list(species_jetz=phylogeny_for_lice_diversity), #class phylo
                   #bayes = TRUE,
                   REML = TRUE, 
                   verbose = TRUE, 
                   s2.init = .25)
summary( l_r)

### Exploratory plots

ggplot(lice_df_richness_genera, aes(x = sociality, y=richness, fill= sociality, group=sociality)) +
  geom_boxplot() +
  labs(title="a) Lice diversity (total genera)",y="Total genera", x = "sociality")+
  theme_classic(20)
###_###_###_###_
# at the species level 
###_###_###_###_

lice_df_diversity_species<-read.csv("data/7.lice_df_diversity_species.csv",header=TRUE,na.strings=c(""))
names(lice_df_diversity_species)
phylogeny_for_lice<- read.nexus("data/phylo_data/1_host_consensus_tree_lice.nex")

# Filtering only manu data 
lice_df_diversity_sp<-lice_df_diversity_species %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

names(lice_df_diversity_sp)
# REMOVE EMPTY ROWS (for some samples we do not have Lice genus2 and tree but when summarizing the table the field were created)
lice_df_diversity_sp<-lice_df_diversity_sp %>% drop_na(valuecol)

#richness at the genus level by host species
lice_df_diversity_sp<-lice_df_diversity_sp %>% group_by(species_jetz) %>% 
  summarise(richness_sp=n_distinct(valuecol))

species_atributes<-lice_df_diversity_species %>% select(elevation_cat, sociality, foraging_cat, species_jetz, species_clean)
species_attributes_distict<-distinct( species_atributes)

lice_df_richness_sp<-right_join(species_attributes_distict, lice_df_diversity_sp, by="species_jetz")  # speceis that are in the ectoparasite list that do not have a matcj in b 

# We are worried that diversity increases with sample size, lets check is that is the case
# summarize the number of samples per species  
n_samples_lice_sp<-lice_df_diversity_species %>% group_by(species_jetz) %>% 
  summarise(n_samples_lice=n_distinct(Full_Label))

is.na(lice_df_diversity$valuecol)

####Correlation sample size and diversity #WARNING

lice_diversity_sp_samplesize<-inner_join(lice_df_richness_sp, n_samples_lice_sp, by="species_jetz")

plot(richness_sp~n_samples_lice, data=lice_diversity_sp_samplesize)

model<-lm(richness_sp~n_samples_lice, data=lice_diversity_sp_samplesize)
abline(model)

str(lice_diversity_sp_samplesize)

#The model

l_d_sp_model_glmm<-  lme4::glmer(richness_sp~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz) +(1|n_samples_lice), 
                              data = lice_diversity_sp_samplesize, 
                              family=poisson(), # use when bayes=true "zeroinflated.poisson",
                              #REML = TRUE, 
                              verbose = TRUE)
summary(l_d_sp_model_glmm)



names( lice_abundance)
l_r_sp<-  phyr::pglmm(richness_sp~ sociality + (1|elevation_cat) + (1|foraging_cat)+(1|species_jetz__)+(1|n_samples_lice), 
                   data = lice_diversity_sp_samplesize, 
                   family ="poisson", # use when bayes=true "zeroinflated.poisson",
                   cov_ranef = list(species_jetz=phylogeny_for_lice_diversity), #class phylo
                   #bayes = TRUE,
                   REML = TRUE, 
                   verbose = TRUE, 
                   s2.init = .25)
summary( l_a )

###
ggplot(lice_diversity_sp_samplesize, aes(x = sociality, y=richness_sp, group=sociality, color=sociality)) +
  geom_boxplot() +
  scat
  ylim(0,10)+
  theme_classic(20)

View(lice_diversity_sp_samplesize)





# Part 10 PGLS ------------------------------------------------------------

fit <- gls(response~explanatory, correlation=corPagel(1, phy=my.tree), data=dat)
  res <- residuals(fit, type="normalized")
