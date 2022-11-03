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
unique(lice_only$total_lice)
mites_only<-read.csv("data/7.ectoparsite_raw_mites_abundance_07222022.csv")
ticks_only<-read.csv("data/7.ectoparsite_raw_ticks_abundance_07222022.csv")

flocks<-read.csv ("data/0.flocks_manu_complete_18052022.csv")
bird_traits_manu<-read.csv("data/4.df_traits_manu_birds.csv") # this is t he final file with traits of Manu  this files includes sociality data binomial 1_0 AND including diet and foraging from PCoA
jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv")
manu_detections<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022.csv") # with elevations

###_###_###_###_###_###_###_###_###_###_###_###_###_###_
#traits data selected 
df_traits_selected<-bird_traits_manu %>% 
  select(sociality,species,mass_tidy,ForStrat_ground,ForStrat_understory,ForStrat_midhigh,ForStrat_canopy,ForStrat_aerial)
###_###_###_###_###_###_###_###_###_###_###_###_###_###_


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
###_###_###_###_###_###_###_###_###_###_###_###_###_###_



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

# Step 5 Generate the databases for analyses  ------------------------------------------------------------------------
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


# Part 3 Structuring files as required for the models------------------------------------------------------------------

# create factors for the elevation and foraging strata
ectoparasites_df_jetz<-inner_join(ectoparasites_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
lice_df_jetz<-inner_join(lice_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
mites_df_jetz<-inner_join(mites_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))
ticks_df_jetz<-inner_join(ticks_manu_traits_selected,jetz_tax_manu_selected, by = c("species_clean"="species_taxonomy_SACC_2021"))

# Part 3.a Presence absence------------------------------------------------------------------
#ticks<-ticks %>% mutate_at("species_clean", str_replace, "_", " ")
###_###_###_###
#  presence absence analyses~ General ectoparasites
###_###_###_###

ectoparasites_df<-ectoparasites_df_jetz %>% mutate( foraging_cat = case_when(ForStrat_ground>50 ~ "ground", 
                                                ForStrat_understory>50~"understory",
                                                ForStrat_midhigh>50~"midhigh",
                                                ForStrat_canopy>50~"canopy",
                                                (ForStrat_ground=50)&(ForStrat_understory=50)~"ground-understory",
                                                (ForStrat_midhigh=50)&(ForStrat_understory=50)~"understory-midhigh",
                                                (ForStrat_midhigh=50)&(ForStrat_canopy=50)~"midhigh-canopy",
                                                (ForStrat_midhigh>=30)&(ForStrat_understory>=30)&(ForStrat_ground>=30)~"ground-midhigh",
                                                TRUE~  "other"))%>% 
  mutate(elevation_site=Full_Label ) %>% 
  mutate_at("elevation_site", str_replace_all,(c( "P"="lowland_manu__", "W"="highmontane_manu__",
                                                  "SP"="montane_manu__", "TL"="montane_manu__", 
                                                  "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                  "M"="lowland_iquitos__",
                                                  "V"="lowland_iquitos__" ))) %>% 
              separate(elevation_site, c("elevation_cat","sample"),"__")

ectoparasites_df$elevation_cat[ectoparasites_df$elevation_cat>=0]<-"lowland_iquitos"   # convert the numerical values that we have without core to lowland iquitos

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

# Part 3.b Abundance------------------------------------------------------------------

###_###_###_###
# Abundance analyses_Lice
###_###_###_###

lice_df<-lice_df_jetz %>% mutate( foraging_cat = case_when(ForStrat_ground>50 ~ "ground", 
                                                                             ForStrat_understory>50~"understory",
                                                                             ForStrat_midhigh>50~"midhigh",
                                                                             ForStrat_canopy>50~"canopy",
                                                                             (ForStrat_ground=50)&(ForStrat_understory=50)~"ground-understory",
                                                                             (ForStrat_midhigh=50)&(ForStrat_understory=50)~"understory-midhigh",
                                                                             (ForStrat_midhigh=50)&(ForStrat_canopy=50)~"midhigh-canopy",
                                                                             (ForStrat_midhigh>=30)&(ForStrat_understory>=30)&(ForStrat_ground>=30)~"ground-midhigh",
                                                                             TRUE~  "other"))%>% 
  mutate(elevation_site=Full_Label ) %>% 
  mutate_at("elevation_site", str_replace_all,(c( "P"="lowland_manu__", "W"="highmontane_manu__",
                                                  "SP"="montane_manu__", "TL"="montane_manu__", 
                                                  "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                  "M"="lowland_iquitos__",
                                                  "V"="lowland_iquitos__" ))) %>% 
  separate(elevation_site, c("elevation_cat","sample"),"__")

lice_df$elevation_cat[lice_df$elevation_cat>=0]<-"lowland_iquitos"   # convert the numerical values that we have without core to lowland iquitos

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
###_###_###_###


###_###_###_###
# Abundance analyses_mites
###_###_###_###
names(mites_df_jetz)
mites_df<-mites_df_jetz %>% mutate( foraging_cat = case_when(ForStrat_ground>50 ~ "ground", 
                                                           ForStrat_understory>50~"understory",
                                                           ForStrat_midhigh>50~"midhigh",
                                                           ForStrat_canopy>50~"canopy",
                                                           (ForStrat_ground=50)&(ForStrat_understory=50)~"ground-understory",
                                                           (ForStrat_midhigh=50)&(ForStrat_understory=50)~"understory-midhigh",
                                                           (ForStrat_midhigh=50)&(ForStrat_canopy=50)~"midhigh-canopy",
                                                           (ForStrat_midhigh>=30)&(ForStrat_understory>=30)&(ForStrat_ground>=30)~"ground-midhigh",
                                                           TRUE~  "other"))%>% 
  mutate(elevation_site=Sample.Full..) %>% 
  mutate_at("elevation_site", str_replace_all,(c( "P"="lowland_manu__", "W"="highmontane_manu__",
                                                  "SP"="montane_manu__", "TL"="montane_manu__", 
                                                  "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                  "M"="lowland_iquitos__",
                                                  "V"="lowland_iquitos__" ))) %>% 
  separate(elevation_site, c("elevation_cat","sample"),"__")

mites_df$elevation_cat[mites_df$elevation_cat>=0]<-"lowland_iquitos"   # convert the numerical values that we have without core to lowland iquitos

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
ticks_df<-ticks_df_jetz %>% mutate( foraging_cat = case_when(ForStrat_ground>50 ~ "ground", 
                                                             ForStrat_understory>50~"understory",
                                                             ForStrat_midhigh>50~"midhigh",
                                                             ForStrat_canopy>50~"canopy",
                                                             (ForStrat_ground=50)&(ForStrat_understory=50)~"ground-understory",
                                                             (ForStrat_midhigh=50)&(ForStrat_understory=50)~"understory-midhigh",
                                                             (ForStrat_midhigh=50)&(ForStrat_canopy=50)~"midhigh-canopy",
                                                             (ForStrat_midhigh>=30)&(ForStrat_understory>=30)&(ForStrat_ground>=30)~"ground-midhigh",
                                                             TRUE~  "other"))%>% 
  mutate(elevation_site=Sample.Full..) %>% 
  mutate_at("elevation_site", str_replace_all,(c( "P"="lowland_manu__", "W"="highmontane_manu__",
                                                  "SP"="montane_manu__", "TL"="montane_manu__", 
                                                  "TU"="high_montane_manu__", "VC"="lowland_manu__",
                                                  "M"="lowland_iquitos__",
                                                  "V"="lowland_iquitos__" ))) %>% 
  separate(elevation_site, c("elevation_cat","sample"),"__")

ticks_df$elevation_cat[ticks_df$elevation_cat>=0]<-"lowland_iquitos"   # convert the numerical values that we have without core to lowland iquitos
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

# Part 3.c Diversity Lice------------------------------------------------------------------
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
write_csv(lice_df_diversity_genus,"data/7.lice_df_diversity_genus.csv")
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
                
  # Part 3.c Diversity Mites------------------------------------------------------------------
  
# Mites
  
  #At the group level
  view(mites_df)
  
 names(mites_df)
  mites_df_wide<-mites_df%>% group_by(Host.Family,BLFamilyLatin,species_clean,species_jetz,Mite.Group, Mite.Genus, Mite.Group2, Mite.Genus2,                                    sociality,foraging_cat, elevation_cat,TipLabel, Sample.Full.. ) %>% dplyr::summarise()
  
  keycol <- "column"
  valuecol <- "mites_group"
  gathercols <- c("Mite.Group", "Mite.Group2")
  
  mites_df_diversity_group<-gather(mites_df_wide, keycol, valuecol, gathercols)%>% arrange(desc(species_jetz)) 
  
  View(mites_df_diversity_group)
  ###_###_###_###
  write_csv(mites_df_diversity_group,"data/7.mites_df_diversity_group.csv")
  ###_###_###_###
  
  # At the genus level 
  
  view(mites_df)
  names(mites_df)
  
  mites_df_wide2<-mites_df%>% group_by(Host.Family,BLFamilyLatin,species_clean,species_jetz,Mite.Group, Mite.Genus, Mite.Group2, Mite.Genus2,
                                      sociality,foraging_cat, elevation_cat,TipLabel, Sample.Full.. ) %>% dplyr::summarise()
  
  keycol <- "column"
  valuecol <- "mites_genus"
  gathercols <- c("Mite.Genus", "Mite.Genus2")
  
  mites_df_diversity_genus<-gather(mites_df_wide, keycol, valuecol, gathercols)%>% arrange(desc(species_jetz)) 
  
  View(mites_df_diversity_genus)
  ###_###_###_###
  write_csv(mites_df_diversity_genus,"data/7.mites_df_diversity_genus.csv")
  ###_###_###_###
  
# Diversit level at the ticks leves( we dont have enight information)

                
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
