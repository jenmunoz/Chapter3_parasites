##############################################################################-# 
## Analyses chapter 2 PhD_ Life history traits 
##Phylogenetic tree of Manu species
## R-code 
## Jenny Munoz
#### last update: Apr 21 2021
###############################################################################-#
# Assumptions":

#Notes: we will be using the  phylogenetic tree downloaded from https://birdtree.org

# Setup--------------------------------------------------------
rm(list=ls()) # function to delete all objects
#1#Loading packages
install.packages("ape")
install.packages("here")
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

# ------Data--------------------------------------------------------

#species_manu<-read_csv("1.Manu_bird_species_detections_tidy_tax&elev_method.csv") 
#taxonomy_2021<-read.csv("taxonomy_revision_2021.csv")  # list of changes 
#jetz_manu<-read_csv("jetz_information_species_manu.csv")

getwd()
list_manu<-read_csv("1.Manu_bird_species_list_tax_changes_tracking_tidy.csv") # list of species with 2021 taxonomy and revisited taxonomy
jetz_list<-read_csv("phylo_analyses/PhyloMasterTax_Jetz.csv")
jetz_manu<-read_csv("phylo_analyses/list_manu_jetz_tax_missmatches_corrected.csv")

##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__
# Data cleaning_Combinig the files ------------------------------------------------------
#Matching the files
str(jetz_list)
str(list_manu)

# These are the specie that are common in the two files 
intersect( jetz_list$species,list_manu$species) #658 species are common between the two files, in total we have 689 in the Manu list so we are missing 32 species 
n_unique(list_manu$species)

# These are the species that are different between the twofiles
setdiff(list_manu$species,jetz_list$species) 

# Let's create a list of those species and try to resolve the taxonomy
missmatches<-anti_join(list_manu,jetz_list,by="species")
missmatches<-missmatches%>%mutate(missmatches= TRUE)

list_manu_jetz_tax<-full_join(list_manu,jetz_list, by="species") 
list_manu_jetz_tax<-full_join(list_manu_jetz_tax,missmatches, by="species") 

# Export_file to resolve some of  the taxonomy manually
#write_csv(list_manu_jetz_tax, "list_manu_jetz_tax_w_missmatches_1.csv")

# Read file with the taxonomy align with jetz
jetz_manu<-read_csv("list_manu_jetz_tax_missmatches_corrected.csv")

n_unique(jetz_manu$species)
n_unique(jetz_manu$species_jetz)

#I want to make sure we have the names are consistent with jetz (double check)
setdiff(jetz_manu$species_jetz,jetz_list$species) 

# "Pyrrhyra roseifrons" is a typo but also is a species that does not have a record in jetz


# Extract the list of species with unknown taxonomy 
unk_taxonomy<-anti_join(list_sp_flocks_SAmerica_taxonomy, taxonomy_higher_level, by="species")
n_unique(unk_taxonomy$species)
 
#Other calculations
#---------------#
list_sp_taxo_changes<-taxonomy_2021 %>% 
  filter (notes_taxonomy=="REV_TAX") %>%
  select(species, species_taxonomy_SACC_2021, notes_taxonomy)

list_manu_w_tax_changes<-full_join(list_manu, list_sp_taxo_changes, by="species_taxonomy_SACC_2021")
write_csv(list_manu_w_tax_changes, "1.Manu_bird_species_list_tax_changes_tracking_tidy.csv")

# Check that non matching files do not exist
anti_join( list_sp_taxo_changes, list_manu,by="species_taxonomy_SACC_2021")

#---------------#

#
##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__
# Extractig 1000 trees -----------------------------------------------------
#Import phylogeny from birdtree.org
#Go to birdtree.org and get a Hackettbackbone tree of the species of interest,
#filter species from the list then they will send you a output.nex element that you will need to use here

ape::read.nexus(filename, multiPhylo=TRUE)
multi_species_tree <- read.nexus('phylo_analyses/tree_pruner/output.nex') 
class(multi_species_tree)# Must be multiPhylo
#multitree <- read.nexus(here('Data/1_jetz_birdtree.nex')) # From many, but I never found a data file qith taht name so I asumed is the one named output.nex

#
##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__
# Opt1 Generate a consensus tree -----------------------------------------------

# Opt 1  Generate a consensus tree  with phytotools -------------------------------

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
# Methods to compute consensus edge lengths (branch lengths for a consensus topology

#Method 3: library phytotools Compute the non-negative least squares edge lengths on the consensus tree 
#using the mean patristic distance matrix. (Function argument method="least.squares".) 
#If the input trees are rooted & ultrametric, this can be used to produce a consensus tree that is also ultrametric.

# consensus edges function computes for a tree under some criterion
manu_consensus_tree<-consensus.edges(multi_species_tree,method="least.squares")# generates a consensus tree with branch lenghts

#plots
plotTree(manu_consensus_tree,fsize=0.01, type="fan")
plot(manu_consensus_tree,type="fan",no.margin = TRUE, cex=0.01)
manu.consensus

# Calculate phylogenetic distance of species pairs 
ph_distance<-cophenetic(manu.consensus)
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

# save phylogeny
saveRDS(manu_consensus_tree, file='phylo_analyses/outputs/1_birdtreeManuspecies.rds')
str(birdtree)

## Save tree
write.tree(manu_consensus_tree, file="phylo_analyses/outputs/1_birdtreeManuspecies.tre")
write.nexus(manu_consensus_tree, file="phylo_analyses/outputs/MyNexusTreefilespecies.nex")
write.tree(manu_consensus_tree, file='phylo_analyses/outputs/birdtreeManuspecies.txt') 

##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__##__
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

# save phylogeny
saveRDS(birdtree, file='1_birdtreeManu.rds')
str(birdtree)

# Save tree
write.tree(birdtree, file="1_birdtreeManu.tre")
write.nexus(MyTree, file="MyNexusTreefile.nex")
write.tree(birdtree, file='birdtreeManu.txt') 
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





# create list of taxa for birdtree.org
compiled_tbl <- readRDS(file = here('RDS_files/2_Table1_compiled_HOF_parameters.rds'))

birdlist <- unique(compiled_tbl$species)

#Extra notes
#consensus.edges---Compute consensus edges for a tree under some criterion
#contMap Map continuous trait evolution on the tree
data<-data(anole.data)
class(anole.data)
plotTree(anoletree)

