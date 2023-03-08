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

flocks_networks<-read.csv("data/9_flocks_for_networks_all_flocks.csv")
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
#write.csv(list_manu_jetz, "data/1.list_manu_species_jetz_taxonomy.csv")

manu_detections_jetz<-left_join(manu_detections, list_manu_jetz, by=("species_taxonomy_SACC_2021"), multiple='all')
#write.csv(manu_detections_jetz, "data/1.Manu_bird_species_detections_tidy_taxonomy_18052022_jetz_taxo_included.csv")


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


###_###_###_###_
# The data 
###_###_###_###_

# Detections
manu_detections_jetz<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022_jetz_taxo_included.csv")
# ectoparasite samples for social species only 
samples<-read.csv("data/data_analyses/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv") %>% 
  rename(elevation=elevation_extrapolated_date) %>% 
  filter(sociality=="1") %>% 
  filter(elevation!="NA") # remove ectoparasite samples for which we dont have elevation
unique(samples$elevation)


# the ranges observed
observed_species_range_limits<-manu_detections_jetz %>% group_by(species_jetz) %>% 
  summarise(low_limit=min(elevation), high_limit=max(elevation)) %>% 
  mutate_at("species_jetz", str_replace, " ", "_")

#write.csv(observed_species_range_limits,"data/data_analyses/1.observed_species_range_limits.csv")

# the flocks lists
flocks_list<-read.csv("data/data_analyses/1.manu_bird_species_list_flocks.csv")

list_manu_jetz<-read.csv("data/1.list_manu_species_jetz_taxonomy.csv") %>% select( species_jetz, species_from_detections,species_taxonomy_SACC_2021) %>% 
  distinct(species_taxonomy_SACC_2021, .keep_all = TRUE) %>% 
  mutate_at("species_jetz", str_replace, " ", "_") %>% 
  mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_")

flock_list_jetz<-inner_join(flocks_list,list_manu_jetz, by="species_taxonomy_SACC_2021" )


# the ranges observed flock members
observed_species_range_limits_flocks<-inner_join(observed_species_range_limits,flock_list_jetz, by=("species_jetz"))
#write.csv(observed_species_range_limits_flocks,"data/data_analyses/1.observed_species_range_limits_flocks.csv")

#View(species_range_limits)
flocks_networks<-read.csv("data/9_flocks_for_networks_all_flocks.csv")
#flocks_networks_file<-read.csv("data/9_flocks_for_networks_all_flocks_reformatted.csv")

###_###_###_###_###_###_###_###_###_
# Steps to extract degree for each indivuald sample
###_###_###_###_###_###_###_###_###_

#i)Extract the elevation of a sample
#ii) Extract the species name and sample level
#iii) Extract the flocking species list of species co-ocurring at that elevation
#iv) Filter in the individual flocks, the flpcking species only
# v) Buld the network 
# vi) Extract the parameters for each individual species at that elevation 

# ALternatively 
# Repit i) to iii), then build one network and subset by the species list 
#shoudl give the same result 

network_metrics_ectos_samples.db<-NULL
for(i in 1:nrow(samples)){
  elevation<-samples$elevation[i]
  species<-samples$species_taxonomy_SACC_2021[i]
  sample_full_label<-samples$Full_Label[i]
  species_list<-observed_species_range_limits_flocks[observed_species_range_limits_flocks$low_limit<(elevation)&observed_species_range_limits_flocks$high_limit>(elevation),]
  flock_list<-flocks_networks %>% select(c(flock_id,species_list$species_taxonomy_SACC_2021)) # 
  rownames(flock_list)<- flock_list[,1] # make the first column the raw names(flock id shoudl be the row names)
  flocks_formatted<-(flock_list[,-1]) # delete first row
  network_flock_elevation<-asnipe::get_network(association_data=flocks_formatted,data_format = "GBI",association_index = "SRI")
  net_flocks_elevation<-igraph::graph.adjacency(network_flock_elevation, mode="undirected", weighted=TRUE, diag=FALSE)
  degree_centrality<-as.data.frame(igraph::degree(net_flocks_elevation))
  degree_centrality$species_network<-row.names(degree_centrality )
  degree_centrality_species<-degree_centrality %>% filter(species_network==species) %>% rename(degree=`igraph::degree(net_flocks_elevation)`)
  #w_degree<-rowSums(network_flock_elevation) # manually but unsure his is correct
  deg_weighted<-as.data.frame(igraph::graph.strength(net_flocks_elevation))
  deg_weighted$species_network<-row.names(deg_weighted)
  deg_weighted_species<-deg_weighted %>% filter(species_network==species)%>% rename(w_degree=`igraph::graph.strength(net_flocks_elevation)`)
  network_metrics<-as.data.frame(bind_cols(sample_full_label,species, degree_centrality_species$degree, deg_weighted_species$w_degree, elevation))
  network_metrics_ectos_samples.db<-rbind(network_metrics_ectos_samples.db,network_metrics) # this is an odd wa of writing it but it works
}

network_metrics_ectoparasite_samples<-as.data.frame(network_metrics_ectos_samples.db) %>% 
  rename(sample=...1, species=...2,degree=...3, elevation=...5,w_degree=...4)
#network_metrics<-data.frame(sample_full_label,species ,degree_centrality_species$degree, deg_weighted_species$w_degree, elevation)

write.csv( network_metrics_ectoparasite_samples, "data/data_analyses/ 7_dff_network_metrics_ectoparasite_samples_FILE.csv")

plot(network_metrics_ectoparasite_samples$elevation,data$degree)
plot(network_metrics_ectoparasite_samples$elevation,data$w_degree)

ggplot( network_metrics_ectoparasite_samples, aes(x=elevation, y=degree, by=species,color=species))+
  geom_point(alpha=0.5)+
  stat_smooth(method=glm,by=species)+
  guides(colour=FALSE)
  


# the other approach 


network_metrics_ectos_samples.db<-NULL
for(i in 1:nrow(samples)){
  elevation<-samples$elevation[i]
  species<-samples$species_taxonomy_SACC_2021[i]
  sample_full_label<-samples$Full_Label[i]
  species_list<-observed_species_range_limits_flocks[observed_species_range_limits_flocks$low_limit<(elevation)&observed_species_range_limits_flocks$high_limit>(elevation),]
  flock_list<-flocks_networks %>% select(c(-elevation,-X)) # 
  rownames(flock_list)<- flock_list[,1] # make the first column the raw names(flock id shoudl be the row names)
  flocks_formatted<-(flock_list[,-1]) # delete first row
  network_flock_elevation<-asnipe::get_network(association_data=flocks_formatted,data_format = "GBI",association_index = "SRI")
  net_flocks_elevation<-igraph::graph.adjacency(network_flock_elevation, mode="undirected", weighted=TRUE, diag=FALSE)
  subset<-subgraph(net_flocks_elevation,species_list$species_taxonomy_SACC_2021 )
  degree_centrality<-as.data.frame(igraph::degree(subset))
  degree_centrality$species_network<-row.names(degree_centrality )
  degree_centrality_species<-degree_centrality %>% filter(species_network==species) %>% rename(degree=`igraph::degree(subset)`)
  #w_degree<-rowSums(network_flock_elevation) # manually but unsure his is correct
  deg_weighted<-as.data.frame(igraph::graph.strength(subset))
  deg_weighted$species_network<-row.names(deg_weighted)
  deg_weighted_species<-deg_weighted %>% filter(species_network==species)%>% rename(w_degree=`igraph::graph.strength(subset)`)
  network_metrics<-as.data.frame(bind_cols(sample_full_label,species, degree_centrality_species$degree, deg_weighted_species$w_degree, elevation))
  network_metrics_ectos_samples.db<-rbind(network_metrics_ectos_samples.db,network_metrics) # this is an odd wa of writing it but it works
}

network_metrics_ectoparasite_samples<-as.data.frame(network_metrics_ectos_samples.db) %>% 
  rename(sample=...1, species=...2,degree=...3, elevation=...5,w_degree=...4)

write.csv( network_metrics_ectoparasite_samples, "data/data_analyses/ 7_dff_network_metrics_ectoparasite_samples_FILE_method2.csv")


##### just exploring the change in degre by species

manu_detections_jetz<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022_jetz_taxo_included.csv") %>% 
  filter(method=="flock") %>% 
  distinct(species_taxonomy_SACC_2021,elevation)


network_metrics_gradient.db<-NULL
for(i in 1:nrow(manu_detections_jetz)){
  elevation<-manu_detections_jetz$elevation[i]
  species<-manu_detections_jetz$species_taxonomy_SACC_2021[i]
  method<-manu_detections_jetz$methodl[i]
  species_list<-observed_species_range_limits_flocks[observed_species_range_limits_flocks$low_limit<(elevation)&observed_species_range_limits_flocks$high_limit>(elevation),]
  flock_list<-flocks_networks %>% select(c(-elevation,-X)) # 
  rownames(flock_list)<- flock_list[,1] # make the first column the raw names(flock id shoudl be the row names)
  flocks_formatted<-(flock_list[,-1]) # delete first row
  network_flock_elevation<-asnipe::get_network(association_data=flocks_formatted,data_format = "GBI",association_index = "SRI")
  net_flocks_elevation<-igraph::graph.adjacency(network_flock_elevation, mode="undirected", weighted=TRUE, diag=FALSE)
  subset<-subgraph(net_flocks_elevation,species_list$species_taxonomy_SACC_2021 )
  degree_centrality<-as.data.frame(igraph::degree(subset))
  degree_centrality$species_network<-row.names(degree_centrality )
  degree_centrality_species<-degree_centrality %>% filter(species_network==species) %>% rename(degree=`igraph::degree(subset)`)
  #w_degree<-rowSums(network_flock_elevation) # manually but unsure his is correct
  deg_weighted<-as.data.frame(igraph::graph.strength(subset))
  deg_weighted$species_network<-row.names(deg_weighted)
  deg_weighted_species<-deg_weighted %>% filter(species_network==species)%>% rename(w_degree=`igraph::graph.strength(subset)`)
  network_metrics<-as.data.frame(bind_cols(sample_full_label,species, degree_centrality_species$degree, deg_weighted_species$w_degree, elevation))
  network_metrics_gradient.db<-rbind(network_metrics_ectos_samples.db,network_metrics) # this is an odd wa of writing it but it works
}






create a file with changes every 50 m of elevation and see how degree changes

network_metrics<-as.data.frame(bind_cols(sample_full_label,species, degree_centrality_species$degree, deg_weighted_species$w_degree, elevation))

association_pair<-data.frame(species1=species[1], species2=species[2],association=association,n=nrow(flock_list))

association.db<-rbind(association.db, association_pair)

network_metrics_ectos_samples.db

network_metrics_ectos_samples

View(network_metrics_ectos_samples.db)
species

degree_centrality$HeaderName<- row.names(degree_centrality )



class(deg_weighted)

View(degree_centrality)

unique(deg_weighted$name)
##_###_
#measure weighted degree [or node’s strength , that is, the sum of the weights of all its links (Barrat, Barthelemy, et al., 2004 #Calculate the degree weighthed [ the sum of the weight of the edges for each node()] using igraph
##_###_
deg_weighted <- graph.strength(net_manu)
deg_weighted<-as.data.frame(deg_weighted)

net_manu<-igraph::graph.adjacency(network_manu, mode="undirected", weighted=TRUE, diag=FALSE)

#######
##Calculating degree#
#######
network_binary<-network_manu
network_binary[network_binary> 0] <- 1
deg_binary <- rowSums(network_binary)
deg_binary<-as.data.frame(deg_binary)



  #my_network10<-get_network(association_data= network_matrix,data_format = "GBI",association_index = "SRI")
}

species_taxonomy_SACC_2021

View(flock_list1)


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




