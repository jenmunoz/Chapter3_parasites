####################################################################
####################################################################
##  Network analyses for parasite host across te gradient
## 
## Explorartory analyses for Manu Networks
## R-code Adapted from SNA workshop (Farine 2020)
#### last updated: Nov 10 2022
####################################################################

# 1. initial set up ----------------------------------------------------------

## clear the decks
# rm(list = ls())

## where am I working?
here::here()

## load required packages
install.packages("asnipe")
install.packages("assortnet")
install.packages("sna")
install.packages("igraph")
install.packages("vegan")
install.packages("network")
install.packages("ndtv") #helpful to create dynamic plots, and dynamic networks
install.packages("bipartite") # To compare environmental gradients networks see paper  Jochen Fründ  2021
#install.packages("betalinkr")
install.packages("tidyverse")  
install.packages("purrr")
install.packages("skimr")
install.packages("igraph")

library(igraph)
library(asnipe)
library(assortnet)
library(sna)
library(igraph)
library(vegan)
library(network)
library(ndtv) #helpful to create dynamic plots, and dynamic networks
library(bipartite)
#library(betalinkr)

library(tidyverse)   
library(purrr)
library(skimr)

## set a plotting theme
theme_set(theme_few())


# 2. Read the data --------------------------------------------------------
#My data

flocks_manu<-read.csv ("data/0.flocks_manu_complete_18052022.csv", header=T)

jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv")


# Checking the data strcuture
flocks_manu<-flocks_manu %>% filter(database_decision=="include") 
#Look for duplicates there are 755 observations, 4 repited entrees
{flocks_manu} %>%
  dplyr::group_by( flock_id, species_taxonomy_SACC_2021) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 
#flocks %>% group_by( flock_id, species_taxonomy_SACC_2021) %>%filter(n() > 1)
#eliminate duplicates 
flocks_manu<-flocks_manu %>% distinct( flock_id, species_taxonomy_SACC_2021, .keep_all = TRUE)
unique(flocks_manu$flock_id)
# create list of flocking species
flocks_list<-flocks_manu %>% 
  distinct( species_taxonomy_SACC_2021, keep_all=FALSE) %>% 
  add_column(flocking="yes") %>% 
  select(-keep_all)

n_unique(flocks_manu$flock_id) ##  ## how many unique flocks?



# 3. Networks along the gradient (nont separated data by communities) ------------------------------------------
#(i) mapping species interactions, 
#(ii) quantifying network structure, and 
#(iii) quantifying specific interaction patterns 
#(iv) exammining the contribution of community cahnegs and species interactions


####Create matrices

#create a column with values 
flocks_manu_net<-flocks_manu%>% mutate(value=1)
flocks_manu_net<-flocks_manu_net%>% group_by(species_taxonomy_SACC_2021, flock_id, elevation) %>% 
  summarise(presence=sum(value))
#Change to presence absence
flocks_manu_net$presence[flocks_manu_net$presence>0] <- 1 #(if needed to use presence absence)

detection_matrix_flocks_manu<-flocks_manu_net %>%select(flock_id,presence,species_taxonomy_SACC_2021) %>% 
pivot_wider(names_from=species_taxonomy_SACC_2021, values_from=presence)

as.data.frame(detection_matrix_flocks_manu)# save as data frame

#replaces na with zeros
detection_matrix_flocks_manu[is.na(detection_matrix_flocks_manu)] <- 0
detection_matrix_flocks_manu<-arrange(detection_matrix_flocks_manu, flock_id) # arrange by flock number
View(detection_matrix_flocks_manu)

# create the matrix
#write.csv(detection_matrix_flocks_manu, "8.detection_matrix_all_flocks_manu_gradient.csv")


# 2.C. Create the networks -----------------------------------------------------
###
View(network)
matrix<-read.csv("8.detection_matrix_all_flocks_manu_gradient.csv",header=TRUE) # Full network manu
matrix<-arrange(matrix, flock_id) # arrange by flcok number


attributes<-read.csv<-read.csv("data/8.attributes_all_flocks_manu.csv", header=T,stringsAsFactors=FALSE) ####Remember that this MUST be in teh same order that the matrix of  detections our case by Flock site
attributes<-arrange(attributes, flock_id) # arrange by flcok number

##-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_
#Warning!!! #Delete the first column with the flock id  and check that the first row and the first column actually contain data rather than column names or row names
###-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_

#Using  "asnipe" we create the network of species associations[matrix of degrees, association strengh or co-ocurrence between species pairs] the fuction get_network() is using the SRI (Simple ration Index that is already accounting for opportunity). See help for other options

# Structuring the network 
matrix<-(matrix[,-1]) # run two times!!! #Delete the first column with the flock id  and check that the first row and the first column actually contain data rather than column names or row names
rownames(matrix) <-matrix[,1] # make the firs column the raw names
matrix<-(matrix[,-1]) # !!! #Delete the first column with the flock id  and check that the first row and the first column actually contain data rather than column names or row names
View(matrix)

##-##-####-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_
# Warning make sure attributes and matrix for the network ane in the same order!!!!!
##-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_##-##-##-##_##_##_##_##_##_

# Convert into a matrix
matrix<-as.matrix(matrix)
class(matrix)

# # Create the network 
#Using  "asnipe" we create the network of species associations[matrix of degrees, association strengh or co-ocurrence between species pairs] 
#the fuction get_network() is using the SRI (Simple ration Index that is already accounting for opportunity). See help for other options

network_manu<-get_network(association_data=matrix,data_format = "GBI",association_index = "SRI",) #lall manu
net_manu<-graph.adjacency(network_manu, mode="undirected", weighted=TRUE, diag=FALSE)
class(network_manu)
##_###_
#measure weighted degree [or node’s strength , that is, the sum of the weights of all its links (Barrat, Barthelemy, et al., 2004 #Calculate the degree weighthed [ the sum of the weight of the edges for each node()] using igraph
##_###_
deg_weighted <- graph.strength(net_manu)
deg_weighted<-as.data.frame(deg_weighted)
#######
##Calculating degree#
#######
network_binary<-network_manu
network_binary[network_binary> 0] <- 1
deg_binary <- rowSums(network_binary)
deg_binary<-as.data.frame(deg_binary)

write.csv(deg_binary, "data/8.network_outputdegree_network_all_sp_manu.csv")
write.csv(deg_weighted, "data/8.network_degree_weighted_network_all_sp_manu.csv")

# Adding a column for  jetz taxonomy

degree_manu<-read.csv("data/8.network_outputdegree_network_all_sp_manu.csv") %>% 
  rename("species_taxonomy_SACC_2021"="X")


degree_w_manu<-read.csv("data/8.network_degree_weighted_network_all_sp_manu.csv") %>% 
  rename("species_taxonomy_SACC_2021"="X")


jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv")%>% 
  mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_") %>% 
  select(species_taxonomy_SACC_2021,species,species_jetz,TipLabel) %>% 
  distinct(species_taxonomy_SACC_2021,.keep_all = TRUE )



degree_manu_jetz<-left_join (degree_manu, jetz_taxonomy_manu_only, by="species_taxonomy_SACC_2021")
degree_w_manu_jetz<-left_join (degree_w_manu, jetz_taxonomy_manu_only, by="species_taxonomy_SACC_2021")

write.csv(degree_manu_jetz, "data/8.network_outputdegree_network_all_sp_manu_jetz_tax.csv")
write.csv(degree_w_manu_jetz, "data/8.network_degree_weighted_network_all_sp_manu_jetz_tax.csv")

                                

View(deg_binary)
# Other ways of calcualting metrics

# Degree igraph 
#network has to be converted to a graph variable (my_network)
detach("package:sna")
library(igraph)

deg_binary_sna<- degree(network_manu, gmode="graph",ignore.eval=TRUE)  #Number of conexion in each node

# Binary degree 
deg_binary_igraph<-degree(network_manu) 

#no sure if this graoh is infomrtive?
hist(degree_distribution(net_manu, cumulative = FALSE, mode=c("total")))

hist(degree_distribution())
#Weighted degree or strenght!!!
deg_igraph <- graph.strength(net_manu) %>% view()

