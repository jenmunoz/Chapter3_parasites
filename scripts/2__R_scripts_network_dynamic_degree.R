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

association.db<-association.db %>% mutate (species_pair= paste0(species1,"&",species2))
#Filter duplicated pairs
association.db<-association.db %>% 
  mutate(species_pair = paste(pmin(species1, species2), pmax(species1, species2), sep = "&")) %>%
  distinct(species_pair,.keep_all = TRUE) 
