#######################################################################################
### Chapter 3-parasites and flocking species                                       ###
### Part 1 Models cleaned                                                          ###
### R-code                                                                         ###
### Jenny Munoz      
###
### Last update: March 2023                                                        ###
################################################################################ ### ###


#networks
install.packages("asnipe")
install.packages("assortnet")
install.packages("sna")
install.packages("igraph")
install.packages("vegan")
install.packages("network")
install.packages("ndtv") #helpful to create dynamic plots, and dynamic networks
install.packages("bipartite") # To compare environmental gradients networks see paper  Jochen Fründ  2021
#install.packages("betalinkr")

#networks
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

flocks_networks<-read.csv("data/9.flocks_for_networks_all_flocks.csv")
pairwise_range_limits<-read.csv("data/9.pairwise_range_limits_thamnophilidae.csv")

# Creating a function to extract the networks
#Warning make sure that the names in the flocks_networks and parwiseran limits/parameters have the same structure (e.g space or dot between genus and species name)
#Also make sure that there are not extra columns in the file, if there are extracolumns then male sure you exclude that column e.g   network1<-(network0[,-1])

# this works pretty well NOW
association.db<-NULL

for (i in 1:nrow(pairwise_range_limits)) {
  low<-pairwise_range_limits$low_limit[i]
  high<-pairwise_range_limits$high_limit[i]
  flock_list<-flocks_networks[flocks_networks$elevation>(low)& flocks_networks$elevation<(high),]
  network0<-(flock_list[,-1])
  network1<-(network0[,-1])
  rownames(network1) <- network1[,1] # make the first column the raw names(flock id shoudl be the row names)
  network2<-(network1[,-1])
  network3<-as.matrix(network2)
  detected<-apply(network3, 2, sum) # filter species that have more than 2 detections
  network3<-network3[,detected>3]
  if(!is.null(nrow(network3))){ if(nrow(network3)>20){
    my_network<-get_network(association_data=network3,data_format = "GBI",association_index = "SRI")
    species<-unlist(strsplit(pairwise_range_limits$species_pair[i],"&"))
    row_vector<-dimnames(my_network)[[1]]
    col_vector<-dimnames(my_network)[[2]]
    association<-my_network[row_vector==species[1], col_vector==species[2]] 
    if (length(association)==0) {association<-NA}# when in the range there are ot enoug bandadas to compute
    association_pair<-data.frame(species1=species[1], species2=species[2],association=association,n=nrow(flock_list))
    association.db<-rbind(association.db, association_pair)}}
}

View(association.db)

str(network3)
class(network3)
View(network3)
association.db<-association.db %>% mutate (species_pair= paste0(species1,"&",species2))
#Filter duplicated pairs
association.db<-association.db %>% 
  mutate(species_pair = paste(pmin(species1, species2), pmax(species1, species2), sep = "&")) %>%
  distinct(species_pair,.keep_all = TRUE) 


#### Rebuilding For the ranges analyses 

# Getting the data ready 

manu_detections<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022.csv")# with elevations
jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv") %>% 
  distinct(species_taxonomy_SACC_2021, .keep_all = TRUE)
list_manu_jetz<-read.csv("data/0.list_manu_species_jetz_taxonomy.csv") %>% select( species_jetz, species_from_detections,species_taxonomy_SACC_2021) %>% 
  distinct(species_taxonomy_SACC_2021, .keep_all = TRUE)

#list_manu<-manu_detections %>% distinct(species_clean,species_taxonomy_SACC_2021)
#list_manu_jetz<-full_join(list_manu,jetz_taxonomy_manu_only, by="species_taxonomy_SACC_2021") %>% rename(species_from_detections=species_clean)
#write.csv(list_manu_jetz, "data/0.list_manu_species_jetz_taxonomy.csv")

manu_detections_jetz<-left_join(manu_detections, list_manu_jetz, by=("species_taxonomy_SACC_2021"), multiple='all')
#write.csv(manu_detections_jetz, "data/1.Manu_bird_species_detections_tidy_taxonomy_18052022_jetz_included.csv")

###_###_###_###_
# The data 
###_###_###_###_

manu_detections_jetz<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022_jetz_included.csv")
samples<-read.csv("data/data_analyses/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv")

# the ranges observed
observed_species_range_limits<-manu_detections_jetz %>% group_by(species_jetz) %>% 
  summarise(low_limit=min(elevation), high_limit=max(elevation)) %>% 
  mutate_at("species_jetz", str_replace, " ", "_")

#View(species_range_limits)
flocks_networks<-read.csv("data/9_flocks_for_networks_all_flocks.csv")
#flocks_networks_file<-read.csv("data/9_flocks_for_networks_all_flocks_reformatted.csv")


  




association.db<-NULL
for(i in 2:nrow(samples)){
  elevation<-samples$elevation[i]
  species_list<-species_range_limits[species_range_limits$low_limit<(elevation)&species_range_limits$high_limit>(elevation),]
  flock_list<-flocks_networks %>% select(c( elevation,species_list$species))
  rownames(flock_list) <- flock_list[,1] # make the first column the raw names(flock id shoudl be the row names)
  flock_list_1<-(flock_list[,-1])
  #my_network10<-get_network(association_data= network_matrix,data_format = "GBI",association_index = "SRI")
}


association.db<-NULL
for(i in 1:nrow(samples)){
  elevation<-samples$elevation[i]
  species_list<-species_range_limits[observed_species_range_limits$low_limit<(elevation)&observed_species_range_limits$high_limit>(elevation),]
  flock_list<-flocks_networks_file %>% filter(species %in% species_list$species)
  network_matrix<- t(flock_list)
  names(network_matrix)<-network_matrix[1,]
  #my_network10<-get_network(association_data= network_matrix,data_format = "GBI",association_index = "SRI")
}





flocks_networks %>% select(c(species_list$species))


step 1: Copy 1st row to header:
  
  names(dat) <- dat[1,]
step 2: Delete 1st row :
  
  dat <- dat[-1,]




network3
View( network_matrix)

class(network_matrix)
View(network_matrix)

flock_list3<-flocks_networks_file[flocks_networks_file$species==(species_list$species)]

flocks_networks_file %>% filter(species %in% species_list$species)

flock_list<-flocks_networks_file[flocks_networks_file$species%in%(species_list$species)]

flock_list<-match_df(flocks_networks_file,species_list)

network3


flock_list<-flocks_networks_file[a]

match_d

flock_list<-flocks_networks_file["flocks_networks_file$species"=="species_list$species"]

View(flock_list)
net<-igraph::graph.adjacency(  my_network5, mode="undirected", weighted=TRUE, diag=FALSE)
plot(net)

class(network_matrix)
View(flock_list)

flocks_networks_file$species


a<-species_list$species


flock_list<-flocks_networks_file[species_list$species,]

View(flock_list)


for (i in 1:nrow(samples)) {
  low<-species_range_limits$low_limit[i]
  high<-species_range_limits$high_limit[i]
  flock_list<-flocks_networks[flocks_networks$elevation>(low)& flocks_networks$elevation<(high),]
  network0<-(flock_list[,-1])
  network1<-(network0[,-1])
  rownames(network1) <- network1[,1] # make the first column the raw names(flock id shoudl be the row names)
  network2<-(network1[,-1])
  network3<-as.matrix(network2)
  detected<-apply(network3, 2, sum) # filter species that have more than 2 detections
  network3<-network3[,detected>3]
  if(!is.null(nrow(network3))){ if(nrow(network3)>20){
    my_network<-get_network(association_data=network3,data_format = "GBI",association_index = "SRI")
    species<-unlist(strsplit(pairwise_range_limits$species_pair[i],"&"))
    row_vector<-dimnames(my_network)[[1]]
    col_vector<-dimnames(my_network)[[2]]
    association<-my_network[row_vector==species[1], col_vector==species[2]] 
    if (length(association)==0) {association<-NA}# when in the range there are ot enoug bandadas to compute
    association_pair<-data.frame(species1=species[1], species2=species[2],association=association,n=nrow(flock_list))
    association.db<-rbind(association.db, association_pair)}}
}

for (i in 1:nrow(species_range_limits)) {
elevation<-species_range_limits$low_limit[i]
  high<-species_range_limits$high_limit[i]
flock_list<-flocks_networks[flocks_networks$elevation>(low)& flocks_networks$elevation<(high),]




