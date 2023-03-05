####################################################################
####################################################################
##  Network analyses for parasite host across the gradient
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

flocks_manu<-read.csv ("data/0.flocks_manu_complete_18052022.csv", header=T) %>% 
  mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_") 
  

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


jetz_taxonomy_manu_only<-read.csv("data/4.list_manu_jetz_tax_missmatches_corrected.csv")%>% 
  mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_") %>% 
  select(species_taxonomy_SACC_2021,species,species_jetz,TipLabel) %>% 
  distinct(species_taxonomy_SACC_2021,.keep_all = TRUE )


flocks_manu_jetz<-left_join (flocks_manu, jetz_taxonomy_manu_only, by="species_taxonomy_SACC_2021")

flocks_list_jetz<-flocks_manu_jetz %>% 
  distinct( species_taxonomy_SACC_2021, keep_all=FALSE) %>% 
  select(-keep_all)


# 3. Networks along the gradient (nont separated data by communities) ------------------------------------------
#(i) mapping species interactions, 
#(ii) quantifying network structure, and 
#(iii) quantifying specific interaction patterns 
#(iv) exammining the contribution of community cahnegs and species interactions


####Create matrices

#create a column with values 
flocks_manu_net<-flocks_manu_jetz%>% mutate(value=1)
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
write.csv(detection_matrix_flocks_manu, "data/8.detection_matrix_all_flocks_manu_gradient_jetz.csv")


# 2.C. Create the networks -----------------------------------------------------
###
View(network)
#matrix<-read.csv("data/8.detection_matrix_all_flocks_manu_gradient.csv",header=TRUE) # Full network manu
matrix<-read.csv("data/8.detection_matrix_all_flocks_manu_gradient_jetz.csv",header=TRUE) # Full network manu with jetz taxonomy 

matrix<-arrange(matrix, flock_id) # arrange by flcok number

View(matrix)

attributes<-read.csv("data/8.attributes_all_flocks_manu.csv", header=T,stringsAsFactors=FALSE) ####Remember that this MUST be in teh same order that the matrix of  detections our case by Flock site
attributes<-arrange(attributes, flock_id) # arrange by flcok number

atributes_species<-read.csv("data/7.dff_lice_abundance_means.csv", header=T,stringsAsFactors=FALSE)

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

network_manu<-asnipe::get_network(association_data=matrix,data_format = "GBI",association_index = "SRI",) #lall manu
net_manu<-igraph::graph.adjacency(network_manu, mode="undirected", weighted=TRUE, diag=FALSE)
plot(net_manu)
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

#write.csv(deg_binary, "data/8.network_outputdegree_network_all_sp_manu.csv")
#write.csv(deg_weighted, "data/8.network_degree_weighted_network_all_sp_manu.csv")

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

#write.csv(degree_manu_jetz, "data/8.network_outputdegree_network_all_sp_manu_jetz_tax.csv")
#write.csv(degree_w_manu_jetz, "data/8.network_degree_weighted_network_all_sp_manu_jetz_tax.csv")

                                

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

# Plotting ectoparasites in thenexworks 
network_manu<-get_network(association_data=matrix,data_format = "GBI",association_index = "SRI",) #lall manu
net_manu<-graph.adjacency(network_manu, mode="undirected", weighted=TRUE, diag=FALSE)


# add attributesto the igrah object 

flocks_list_jetz<-flocks_manu_jetz %>% 
  distinct( species_jetz, keep_all=FALSE) %>% 
  select(-keep_all) %>% 
  mutate_at("species_jetz", str_replace, " ", "_") 
  


atributes_species<-read.csv("data/data_analyses/7.dff_lice_abundance_means.csv", header=T,stringsAsFactors=FALSE) %>% select(species_jetz, mean_lice)
flock_species_ectos<-left_join(flocks_list_jetz,atributes_species, by="species_jetz")
#flock_species_ectos[is.na(flock_species_ectos)] <- 0

#Set attributes in the igraph object V is for vertex 

unique(flock_species_ectos$mean_lice )

# if names matches you can add attributes to the igraph object
V(net_manu)$lice<-flock_species_ectos$mean_lice # adding lice abundance o plot it later
#V(net_manu)$color<-col[V(net_manu)$lice]

V(net_manu)$degree<-deg_binary # adding degree
V(net_manu)$strength<-deg_weighted # adding weoghted degree

# degree comes from here
network_binary<-network_manu
network_binary[network_binary> 0] <- 1
deg_binary <- rowSums(network_binary)
class(deg_binary)

# the plot
png("figures/figures_manuscript/Fig5.Network_degree_lice_abundance_NET.png", width = 2500, height = 3100, res = 300, units = "px")
plot (net_manu,layout=layout_nicely(net_manu),vertex.size=V(net_manu)$degree*0.09,vertex.label.cex=0.2, edge.width=0.2,  edge.lty=c("solid"), edge.color="#00798c", vertex.color=ColorPalette[as.numeric(V(net_manu)$lice)])
dev.off()

png("figures/figures_manuscript/Fig5.Network_strength_lice_abundance_NET.png", width = 2500, height = 3100, res = 300, units = "px")
plot (net_manu,layout=layout_nicely(net_manu),vertex.size=V(net_manu)$strength,vertex.label.cex=0.2, edge.width=0.2,  edge.lty=c("solid"), edge.color="#00798c", vertex.color=ColorPalette[as.numeric(V(net_manu)$lice)])
dev.off()


# THE COLOR PALLETS

palf <- colorRampPalette(c("gray80", "dark red"))
col.1 <- adjustcolor("orange red", alpha=0.9)
col.2 <- adjustcolor("orange", alpha=0.9)
node.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE) 
node.col <- node.pal (10)

ColorPalette <- brewer.pal(n = 20, name = "Reds")

ColorPalette <- brewer.pal(n = 20, name = "Greens")


#slategray3","dodgerblue","darkolivegreen3","goldenrod1"))
##object_color<-setMap(fmode,ColorPalette)

#png("figures/figures_manuscript/Fig1a.Sociality_and_prevalence_phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
#plot(dotTree(phylogeny_rooted,fmode,colors=setNames(c("red","black"),c("1","0")),ftype="i",fsize=0.5, lwd=4),text(x=10,y=-5,"Mixed-species flocks",pos=1))



# # The models for networks -----------------------------------------------

# data prep
ectos_dfff<-read.csv("data/data_analyses/7.dff_ectos_pres_abs.csv")# data on prevalence # should have same number of rows than teh phylogeny 
phylogeny_prevalence<- read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex") 
mean_lice_abundance<-read.csv("data/data_analyses/7.dff_lice_abundance_means.csv") %>% select(species_jetz, mean_lice)
mean_mites_abundance<-read.csv("data/data_analyses/7.dff_mites_abundance_means_non_feathers.csv") %>% select(species_jetz, mean_mites)
ectos_dff<-left_join(ectos_dfff, mean_lice_abundance, by="species_jetz")
ectos_dff$mean_lice[is.na(ectos_dff$mean_lice)] <- 0
ectos_dff<-ectos_dff %>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations
write.csv(ectos_dff, "data/data_analyses/1_dff_ectoparasites_prevalence_abundance_w_zeros.csv")
#

# DATA
network_parameters<-read.csv("data/data_analyses/1.dff_degree_n_strength_network_jetz_manu.csv", header=T,stringsAsFactors=FALSE) %>% select(species_jetz, deg_binary, deg_weighted,species_taxonomy_SACC_2021)
lice_means<-read.csv("data/data_analyses/7.dff_lice_abundance_means.csv")
pasite_parameters<-read.csv("data/data_analyses/1_dff_ectoparasites_prevalence_abundance_w_zeros.csv")
network_parasite_parameters<-inner_join(network_parameters,pasite_parameters, by="species_jetz")

unique(network_parasite_parameters$species_jetz)

phylo=phylo_lice_rooted_abund_networks_manu<- read.nexus("data/phylo_data/1_host_tree_Manuspecies_onetree_rooted_lice_abunmanu_networks.nex")
str(phylo_lice_rooted_abund_networks_manu)
phylo_lice_rooted_abund_networks_manu$tip.label
# Make sure number in teh data and the phylogenty are consistent

# Reformating 
#lice_df_abundance <-lice_df_abundance  %>% mutate_at("species_jetz", str_replace, " ", "_")
network_parasite_parameters$elevation_cat<-as.factor(network_parasite_parameters$elevation_cat)

network_parasite_parameters$species_jetz<-as.factor(network_parasite_parameters$species_jetz)
network_parasite_parameters$sociality<-as.factor(network_parasite_parameters$sociality)

network_parasite_parameters$deg_binary <-as.numeric(network_parasite_parameters$deg_binary )
network_parasite_parameters$deg_weighted  <-as.numeric(network_parasite_parameters$deg_weighted)


# Aling phylogeny and atributes

a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(network_parasite_parameters$species_jetz))%>% mutate(name=network_parasite_parameters$species_jetz) %>% select(name) %>% arrange(desc(name))

tip<-as.list(setdiff(a,b))
print(tip)
phylo<-drop.tip (phylo, tip$name)

network_parasite_parameters %>% filter(species_jetz==c(a$species_jetz))

network_parasite_parameters_filtered<-inner_join(network_parasite_parameters,a, by=c("species_jetz"="name"))



# THE MODEL

names(sociality_continous_abundance_lice)

l_abun_degree<-  phyr::pglmm(mean_lice~deg_binary +(1|elevation_cat)+(1|species_jetz__),
                                 data = network_parasite_parameters_filtered,
                                 family = "gaussian", 
                                 cov_ranef = list(species_jetz=phylo), #class phylo
                                 #bayes = TRUE,
                                 REML= FALSE,  # NOT SURE WHEN TO USE ML
                                 verbose = TRUE,
                                 s2.init = .25) # what is this last parameter for

l_abun_degree_w<-  phyr::pglmm(mean_lice~deg_weighted +(1|elevation_cat)+(1|species_jetz__),
                             data = network_parasite_parameters_filtered,
                             family = "gaussian", 
                             cov_ranef = list(species_jetz=phylo), #class phylo
                             #bayes = TRUE,
                             REML= FALSE,  # NOT SURE WHEN TO USE ML
                             verbose = TRUE,
                             s2.init = .25) # what is this last parameter for

rr2::R2(l_abun_degree_w)

ggplot(network_parasite_parameters_filtered, aes(x=deg_binary, y= mean_lice)) +
  geom_point()+
  geom_smooth(method = glm)+
  geom_jitter(height = 0.01)+
  scale_y_continuous("mean_lice", limits = c(0,20)) +
  labs(title="b) Degree per species ")+
  theme_classic(20)


ggplot(network_parasite_parameters_filtered, aes(x=deg_weighted, y= mean_lice)) +
  geom_point()+
  geom_smooth(method = gam)+
  geom_jitter(height = 0.01)+
  scale_y_continuous("mean_lice", limits = c(0,20)) +
  labs(title="b) Degree_W_per species ")+
  theme_classic(20)


predictions<-predict(l_abun_degree_w, type="response")
data<- cbind(network_parasite_parameters_filtered, predictions)
 
 data %>% 
  ggplot(aes(x=deg_weighted, y=mean_lice)) +
  geom_smooth(aes(x=deg_weighted, y=Y_hat))+
  scale_y_continuous("mean_lice", limits = c(0,10)) 
   
                 
 data %>% 
   ggplot(aes(x=deg_binary, y=mean_lice)) +
   geom_smooth(aes(x=deg_binary, y=Y_hat))+
   scale_y_continuous("mean_lice", limits = c(0,10)) 
 
  

  
  ggplot(aes(x2, fit)) +
  geom_smooth_ci(f1
  

predict_gam(model_2, values = list(f1 = c(0.5, 1, 1.5))) %>%
  ggplot(aes(x2, fit)) +
  geom_smooth_ci(f1)

# Important!!!!!Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_prevalence$tip.label,rownames(ectos_df)),]


