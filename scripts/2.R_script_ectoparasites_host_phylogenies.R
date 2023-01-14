#######################################################################################
### Chapter 3-parasites                                         ###
### Part 3 Phylogent trees estimations                                                                         ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
### Last update: Oct 26 2022                                                     ###
################################################################################

# Loading packages --------------------------------------------------------
#rm(list=ls()) # function to delete all objects
#1#Loading packages
install.packages("ape")
#install.packages("here")
install.packages("phytools")
install.packages("tidyverse")
install.packages("metafor")
install.packages("phangorn") # to reconstruct a maximum clade credibility tree

library(ape)
#library(here)
library(phytools)
library(metafor)
library (phangorn) # to reconstruct a maximum clade credibility tree

library(tidyverse)
library(skimr)


# Getting the data ready  -------------------------------------------------
# Master taxonomy form jetz
taxonomy_jetz<-read.csv( "data/PhyloMasterTax_jetz.csv")
ectoparasite_general_jetz<-read.csv("data/7.ectoparasite_df_presence_absence.csv") # data of general ectoparasites with jetz taxonomy
names (taxonomy_jetz)
names (ectoparasite_genaral_jetz)
intersect(taxonomy_jetz, ectoparasite_genaral_jetz, by=c(Scientific,species_jetz))

# Make sure there are not differences in the lsit of spcies with the master taxonomy from jetz
anti_join(ectoparasite_general_jetz,taxonomy_jetz, by=c("species_jetz"="Scientific")) # speceis that are in the ectoparasite list that do not have a matcj in b 
#Upload the list of species in bird.org


# [Presence absence all parasites]  -----------------------------------------------------------------------
#Step 1 Extractig 10000 trees from  birdtree.org
#Import phylogeny from birdtree.org
#warning we will get the trees only for species for which we have some parasite samples( general parasite list)
#Go to birdtree.org and get a Hackettbackbone tree of the species of interest.
#filter species from the list then they will send you a output.nex element that you will need to use here
# You can use species only with genetic material or teh overall list, we will be using the overall list for now
# becausewhen using the one ony with genetic information we have to reduce the number of species with ectoparasites samples to match the ones for which genetci info is avialable # so i will use th eother


# Read the tree
# there are two options for reading the trees
#ape::read.nexus(filename, multiPhylo=TRUE) #opt 1
host_species_tree <- read.nexus("data/phylo_data/tree_pruner/output_bird_parasites.nex") #This is the tree with all species  of host including Iquitos and Manu
host_species_tree_manu<- read.nexus("data/phylo_data/tree_pruner/1.output_ectos_pres_abs_100_trees_manu_only.nex") # this is for all species with samples from Manu

# Check the properties of the trees
is.rooted.multiPhylo(host_species_tree_manu) # the trees are rooted
is.ultrametric.multiPhylo(host_species_tree_manu) # when using genetic data  only is not ultrametric not sure why
print (host_species_tree_manu, details=TRUE)

#host_species_tree[9789]<-NULL # for some reason this genetic tree has less tips (weird I will just remove it)
random_host_tree_manu<-sample(host_species_tree_manu,size=1)[[1]] # select one tree
#random_host_tree1<-sample( host_species_tree,size=1000)# select a thounsand trees tahta re rooted
class(random_host_tree_manu)
is.rooted.multiPhylo(random_host_tree) # the trees are rooted
is.ultrametric.multiPhylo(random_host_tree)
class( host_species_tree )# Must be multiPhylo

write.tree(random_host_tree_manu, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_ectos_pres_abs.tre")
write.nexus(random_host_tree_manu, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_ectos_pres_abs.nex")


#  [Abundance] ------------------------------------------------------------
####_####_
#For lice only
###_###

#Double check the names are consistent 
lice_df_abundance<-read.csv("data/7.lice_df_abundance.csv") # make sure ou filtered the iquitos data  if desired

lice_df_abundance<-lice_df_abundance %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

unique(lice_df_abundance$species_jetz)
taxonomy_jetz<-read.csv( "data/PhyloMasterTax_jetz.csv")
# Make sure there are not differences in the lsit of spcies with the master taxonomy from jetz
anti_join(lice_df_abundance,taxonomy_jetz, by=c("species_jetz"="Scientific")) # speceis that are in the ectoparasite list that do not have a matcj in b 

# Use the list to extract species from bird.org
# Read the tree
#host_species_tree_lice <- ape::read.nexus("data/phylo_data/tree_pruner/output_bird_lice_abundance.nex") # this is for all species with samples including iquitos and Manu
host_species_tree_lice_manu<-ape::read.nexus("data/phylo_data/tree_pruner/1. output_lice_pres_abs_100_trees_manu_only.nex") # this is amanu species only 
class(host_species_tree_lice_manu)# Must be multiPhylo

# Check the properties of the trees
is.rooted.multiPhylo(host_species_tree_lice_manu) # the trees are rooted
is.ultrametric.multiPhylo(host_species_tree_lice) # when using genetic data  only is not ultrametric not sure why
print (host_species_tree_lice_manu, details=TRUE)

# Lets extract ne tree ( this is just before I figure it out how to extract a concensus rooted tree if that is neccesary)
random_host_tree_lice_manu<-sample(host_species_tree_lice_manu,size=1)[[1]] # select one tree
is.rooted(random_host_tree_lice_manu)
write.tree(random_host_tree_lice_manu, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.tre")
write.nexus(random_host_tree_lice_manu, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")


#  [Diversity] ------------------------------------------------------------
####_####_
#For lice only
###_###

#Double check the names are consistent 
lice_richness_manu_sp<-read.csv("data/5.lice_richness_sp_df_manu.csv") # make sure ou filtered the iquitos data  if desired
lice_richness_manu_sp<-lice_richness_manu_sp %>% distinct( species_jetz,.keep_all = TRUE)
unique(lice_richness_manu_sp)
taxonomy_jetz<-read.csv( "data/PhyloMasterTax_jetz.csv")
# Make sure there are not differences in the lsit of spcies with the master taxonomy from jetz
anti_join(lice_richness_manu_sp,taxonomy_jetz, by=c("species_jetz"="Scientific")) # speceis that are in the ectoparasite list that do not have a matcj in b 

# Use the list to extract species from bird.org
# Read the tree
#host_species_tree_lice <- ape::read.nexus("data/phylo_data/tree_pruner/output_bird_lice_abundance.nex") # this is for all species with samples including iquitos and Manu
host_species_tree_lice_manu<-ape::read.nexus("data/phylo_data/tree_pruner/1. output_lice_pres_abs_100_trees_manu_only.nex") # this is manu species only 
class(host_species_tree_lice_manu)# Must be multiPhylo

# Check the properties of the trees
is.rooted.multiPhylo(host_species_tree_lice_manu) # the trees are rooted
is.ultrametric.multiPhylo(host_species_tree_lice_manu) # when using genetic data  only is not ultrametric not sure why
print (host_species_tree_lice_manu, details=TRUE)

# Lets extract one tree ( this is just before I figure it out how to extract a concensus rooted tree if that is neccesary)
random_host_tree_lice_manu<-sample(host_species_tree_lice_manu,size=1)[[1]] # select one tree
is.rooted(random_host_tree_lice_manu)
write.tree(random_host_tree_lice_manu, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.tre")
write.nexus(random_host_tree_lice_manu, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun.nex")


# [Networks] ----------------------------------------------------------------

# For analyzing Abundance for  individuals
# This tree include only social species because we only have degree for social species
# Read the tree
#host_species_tree_lice <- ape::read.nexus("data/phylo_data/tree_pruner/output_bird_lice_abundance.nex") # this is for all species with samples including iquitos and Manu
host_species_tree_lice_manu_abundance_networks<-ape::read.nexus("data/phylo_data/tree_pruner/1.output_ectos_abundance_100_trees_manu_only_social_networks.nex") # this is manu species only 
class(host_species_tree_lice_manu_abundance_networks)# Must be multiPhylo

# Check the properties of the trees
is.rooted.multiPhylo(host_species_tree_lice_manu_abundance_networks) # the trees are rooted
is.ultrametric.multiPhylo(host_species_tree_lice_manu_abundance_networks) # when using genetic data  only is not ultrametric not sure why
print (host_species_tree_lice_manu_abundance_networks, details=TRUE)

# Lets extract one tree ( this is just before I figure it out how to extract a concensus rooted tree if that is neccesary)
random_host_tree_lice_manu_abun_networks<-sample(host_species_tree_lice_manu_abundance_networks,size=1)[[1]] # select one tree
is.rooted(random_host_tree_lice_manu_abun_networks)
write.tree(random_host_tree_lice_manu_abun_networks, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abun_manu_networks.tre")
write.nexus(random_host_tree_lice_manu_abun_networks, file="data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abunmanu_networks.nex")

