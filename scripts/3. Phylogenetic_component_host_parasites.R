#######################################################################################
### Chapter 3-parasites                                         ###
### Part 3 Phylogenetic component                                                                           ###
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



##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__

# Master taxonomy form jetz

taxonomy_jetz<-read.csv( "data/PhyloMasterTax_jetz.csv")
ectoparasite_general_jetz<-read.csv("data/7.ectoparasite_df_presence_absence.csv") # data of general ectoparasites with jetz taxonomy

names (taxonomy_jetz)
names (ectoparasite_genaral_jetz)
intersect(taxonomy_jetz, ectoparasite_genaral_jetz, by=c(Scientific,species_jetz))

# Make sure there are not differences in the lsit of spcies with the master taxonomy from jetz
anti_join(ectoparasite_general_jetz,taxonomy_jetz, by=c("species_jetz"="Scientific")) # speceis that are in the ectoparasite list that do not have a matcj in b 

#Upload the list of species in bird.org


# Part1_Phylogeny for all bird with Ectoparasites samples (presence_absence) ------------------------------------
# Step 1 Extractig 1000 trees from  birdtree.org
#Import phylogeny from birdtree.org
#warning we will get the trees only for species for which we have some parasite samples( general parasite list)
#Go to birdtree.org and get a Hackettbackbone tree of the species of interest.
#filter species from the list then they will send you a output.nex element that you will need to use here

# Read the tree
# there are two options for reading the three
ape::read.nexus(filename, multiPhylo=TRUE)
host_species_tree <- read.nexus("data/phylo_data/tree_pruner/output_bird_parasites.nex") # this is for all species with samples 
class(multi_species_tree)# Must be multiPhylo
__
# Opt 1  Generate a consensus tree  with phytotools -------------------------------

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# Methods to compute consensus edge lengths (branch lengths for a consensus topology

#Method 3: library phytotools Compute the non-negative least squares edge lengths on the consensus tree 
#using the mean patristic distance matrix. (Function argument method="least.squares".) 
#If the input trees are rooted & ultrametric, this can be used to produce a consensus tree that is also ultrametric.

# consensus edges function computes for a tree under some criterion
host_consensus_tree<-consensus.edges(host_species_tree,method="least.squares")# generates a consensus tree with branch lenghts

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

#### Extra  to calculate phylogeneticdistances

# Calculate phylogenetic distance of species pairs 
ph_distance<-cophenetic(host_consensus_tree)
ph_distance[1:10,1:10]
View(ph_distance)

# Converting matrix into long data frame 
ph_distance<-as.data.frame(ph_distance)
data<-ph_distance  
data<-tibble::rownames_to_column(data,"species") # Apply rownames_to_column
df_phylo_distance<-as.data.frame.table(`row.names<-`(as.matrix(data)[,-1],data$species))
View(df_phylo_distance)
#Rename and filter diagnals
df_phylo_distance<-df_phylo_distance%>% dplyr::rename(species1=Var1, species2=Var2, phylo_distance=Freq) %>% 
  filter(species1!=species2) 

# lets unite the columns
df_phylo_distance<-df_phylo_distance %>% unite(species_pair,1:2, sep ="&", remove=FALSE)
df_phylo_distance$phylo_distance<-as.numeric(df_phylo_distance$phylo_distance)

write.csv(df_phylo_distance, "phylo_analyses/outputs/df_phylo_distance.csv")

# Part2_Phylogeny for all bird with lice samples (Abundance) ------------------------------------

####_####_
#For lice only
###_###

#Double check the names are consistent 
lice_df_abundance<-read.csv("data/7.lice_df_abundance.csv")
unique(lice_df_abundance$species_jetz)

taxonomy_jetz<-read.csv( "data/PhyloMasterTax_jetz.csv")
# Make sure there are not differences in the lsit of spcies with the master taxonomy from jetz
anti_join(lice_df_abundance,taxonomy_jetz, by=c("species_jetz"="Scientific")) # speceis that are in the ectoparasite list that do not have a matcj in b 

# U se the list to extract species from bird.org
# Read the tree
# there are two options for reading the three
host_species_tree_lice <- ape::read.nexus("data/phylo_data/tree_pruner/output_bird_lice_abundance.nex") # this is for all species with samples 
class(host_species_tree_lice)# Must be multiPhylo

#Opt 1  Generate a consensus tree  with phytotools 

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# Methods to compute consensus edge lengths (branch lengths for a consensus topology

#Method 3: library phytotools Compute the non-negative least squares edge lengths on the consensus tree 
#using the mean patristic distance matrix. (Function argument method="least.squares".) 
#If the input trees are rooted & ultrametric, this can be used to produce a consensus tree that is also ultrametric.

# consensus edges function computes for a tree under some criterion
host_species_tree_lice<-consensus.edges(host_species_tree_lice,method="least.squares")# generates a consensus tree with branch lenghts

#plots
plotTree(host_species_tree_lice,fsize=0.01, type="fan")
plot(host_species_tree_lice,type="fan",no.margin = TRUE, cex=0.5)
host_species_tree_lice

# save phylogeny
saveRDS(host_species_tree_lice, file='data/phylo_data/1_host_consensus_tree_lice.rds')
str(birdtree)

## Save tree
write.tree(host_species_tree_lice,  file="data/phylo_data/1_host_consensus_tree_lice.tre")
write.nexus(host_species_tree_lice,  file="data/phylo_data/1_host_consensus_tree_lice.nex")
write.tree(host_species_tree_lice,  file="data/phylo_data/1_host_consensus_tree_lice.txt") 

# Part2_Phylogeny for all bird with Mites samples (Abundance) ------------------------------------

####_####_
#For Mites only
###_###

#Double check the names are consistent 
mites_df_abundance<-read.csv("data/7.mites_df_abundance.csv")
unique(mites_df_abundance$species_jetz)

taxonomy_jetz<-read.csv( "data/PhyloMasterTax_jetz.csv")


# Make sure there are not differences in the lsit of spcies with the master taxonomy from jetz
anti_join(mites_df_abundance,taxonomy_jetz, by=c("species_jetz"="Scientific")) # speceis that are in the ectoparasite list that do not have a matcj in b 

# Use the list to extract species from bird.org
# Read the tree
# there are two options for reading the three
host_species_tree_mites<- ape::read.nexus("data/phylo_data/tree_pruner/output_bird_mites_abundance.nex") # this is for all species with samples 
class(host_species_tree_mites)# Must be multiPhylo

#Opt 1  Generate a consensus tree  with phytotools 

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# Methods to compute consensus edge lengths (branch lengths for a consensus topology

#Method 3: library phytotools Compute the non-negative least squares edge lengths on the consensus tree 
#using the mean patristic distance matrix. (Function argument method="least.squares".) 
#If the input trees are rooted & ultrametric, this can be used to produce a consensus tree that is also ultrametric.

# consensus edges function computes for a tree under some criterion
host_species_tree_mites<-consensus.edges(host_species_tree_mites,method="least.squares")# generates a consensus tree with branch lenghts

View(host_species_tree_mites)
#plots
plotTree(host_species_tree_mites,fsize=0.01, type="fan")
plot(host_species_tree_mites,type="fan",no.margin = TRUE, cex=0.5)
host_species_tree_mites

# save phylogeny
saveRDS(host_species_tree_mites, file='data/phylo_data/1_host_consensus_tree_mites.rds')
str(birdtree)

## Save tree
write.tree(host_species_tree_mites,  file="data/phylo_data/1_host_consensus_tree_mites.tre")
write.nexus(host_species_tree_mites,  file="data/phylo_data/1_host_consensus_tree_mites.nex")
write.tree(host_species_tree_mites,  file="data/phylo_data/1_host_consensus_tree_mites.txt") 



##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__
#Other methods  to extract phylogeny and branch lenghts --------------------------------------------------------

#Method 2: Compute the mean edge length for each edge in the consensus tree setting the length
#for each tree in which the edge is absent to zero. (Default setting. Function arguments method="mean.edge" and if.absent="zero".)
#manu.consensus<-consensus.edges(multi_species_tree)
#plotTree(manu.consensus,fsize=0.4)

#Method 3: Compute the mean edge length, but ignore trees in which the edge is absent.
#(Function arguments method="mean.edge" and if.absent="ignore".)
#manu.consensus<-consensus.edges(multi_species_tree,if.absent="ignore")

# We can use the library ape, function "concensus" or the library phangorn, function "maxCladeCred"
# majority-rule consensus tree (p = 0.5)
consensus_tree <- consensus(multi_species_tree, p = 0.5, check.labels = TRUE)
consensus_tree 
# consensus_tree is unrooted, so compute branch lengths 
birdtree <- compute.brlen(consensus_tree)
#Plot the tree
plot(birdtree, no.margin = TRUE, direction = 'upwards', cex = 4, srt = -20)
###___######___######___######___###
#need to loook at this part more carefully
# check that tips match occurrence data
birdlist2 <- read.csv(file = here('Data/1_bird_species_list.csv'))
birdlist2$Jetz_taxonomy <- gsub(" ", "_", birdlist2$Jetz_taxonomy)

setdiff(birdlist2$Jetz_taxonomy, birdtree$tip.label) # none
all(birdtree$tip.label %in% birdlist2$Jetz_taxonomy) # TRUE
###___######___######___######___######___######___###

# export at width = 14000 
plot(birdtree, no.margin = TRUE, direction = 'upwards', cex = 4, srt = -20)
plot(birdtree, type="fan",no.margin = TRUE, cex = 0,srt = -20)


#Export (write) trees
#To save a tree to a text file, use ape::write.tree(tree, file='filename.txt') for Newick format (widely supported by most phylogenetic software), 
#or ape::write.nexus(tree, file='filename.nex') for Nexus format.
#also check
tidytree()

#Other methods_Option 3 Reconstructing a maximun clade credibility tree -----------------------------------------------------
#library (phangorn) # to reconstruct a maximum clade credibility tree
#maxCladeCred() # FUNCTION to reconstruct a maximum credibility tree

#consensus topologies can be created from downloaded samples using scripts supplied with the R-packages ‘ape’ or ‘RPhylip,’ 
#though these do not supply edge lengths. For that, one could use MrBayes. For “All species” bird trees, 
#consensus representations can be misleading (for both edge lengths and topologies), and we advise that analyses using 
#full trees are done across a reasonable number (»100) of draws from the distributions we supply.


# Diversity_analyses_consensus_trees --------------------------------------

#lice consensus tree for diversity analyses


host_species_tree_lice_diversity<- ape::read.nexus("data/phylo_data/tree_pruner/output_lice_diversity.nex") # this is for all species with samples 
class(host_species_tree_lice_diversity)# Must be multiPhylo

#Opt 1  Generate a consensus tree  with phytotools 


# consensus edges function computes for a tree under some criterion
host_tree_lice_diversity<-consensus.edges(host_species_tree_lice_diversity,method="least.squares")# generates a consensus tree with branch lenghts

#plots
plot(host_tree_lice_diversity,type="fan",no.margin = TRUE, cex=0.5)

# save phylogeny
saveRDS(host_tree_lice_diversity, file='data/phylo_data/1_host_tree_lice_diversity.rds')
str(birdtree)

## Save tree
write.tree(host_tree_lice_diversity,  file="data/phylo_data/1_host_tree_lice_diversity.tre")
write.nexus(host_tree_lice_diversity,  file="data/phylo_data/1_host_tree_lice_diversity.nex")

