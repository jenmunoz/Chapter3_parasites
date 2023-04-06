#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Figures                                                                       ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
### Last update: Nov 2022                                                     ###
################################################################################
# #### 4.Descriptive statistiscs plots -------------------------------------

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

ectos_birds_dff<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, total_lice,total_no_feathers_mites,total_mesostigmatidae,date ) %>% 
  filter(species_jetz!="Premnoplex_brunnescens")

# just keeping the prevalences 
ectos_birds_dff_PA<-read.csv("data/data_analyses/data_manuscript/7.dff_all_ectos_prevalence_abundance_individual_elevation_FILE.csv", na.strings =c("","NA")) %>% 
  rename(elevation=elevation_extrapolated_date) %>%
  select( species_jetz, ectoparasites_PA, sociality ) %>% 
  filter(species_jetz!="Premnoplex_brunnescens")

unique (ectos_birds_dff_PA$species_jetz)

str(ectos_birds_dff_PA)
ectos_birds_dff_PA$ectoparasites_PA<-as.numeric(ectos_birds_dff_PA$ectoparasites_PA)
ectos_birds_dff_PA$species_jetz<-as.factor(ectos_birds_dff_PA$species_jetz)

# creating summaries per species
ectos_birds_dff_PA_species<-ectos_birds_dff_PA %>% group_by(species_jetz) %>% 
  summarise(total_presences=sum(ectoparasites_PA), sample_size=n(), sociality=max(sociality)) %>% 
  mutate(ectos_prevalence=total_presences/sample_size) %>% 
  filter(species_jetz!="Premnoplex_brunnescens") %>% 
  #filter(sample_size>3)


#ectos_birds_dff_mean <- ectos_birds_dff %>% select(species_jetz,total_lice,total_no_feathers_mites,total_mesostigmatidae) %>% filter(complete.cases(.)) %>% group_by(species_jetz) %>% mutate(mean_lice = mean(total_lice), mean_nf_mites=mean(total_no_feathers_mites)) 

ectos_birds_dff_PA_species
### double check that the phylo and data match 
a<-(as.data.frame(phylo$tip.label))%>% mutate(name=phylo$tip.label) %>% select(name) %>% arrange(desc(name))
b<-(as.data.frame(ectos_birds_dff$species_jetz)) %>% mutate(name=ectos_birds_dff$species_jetz) %>% select(name) %>% arrange(desc(name)) %>% distinct(name)

tip<-as.list(setdiff(a,b))
print(tip)
# Drop some tips USE IF NEED TO DROP SOME TIPS when using the full phylogeny
phylo<-drop.tip (phylo, tip$name)

order<-(as.data.frame(phylo$tip.label)) # list of specie sin phylogenetic order
# Important!!!!!Make sure the names  are in the same order in the phylogeny and in the traits 
#rownames(ectos_birds_dff) <- ectos_birds_dff$species_jetz # first make it the row names 
#ectos_birds_dff<- ectos_birds_dff[match(phylo$tip.label,rownames(ectos_birds_dff)),]


# releveling
#plot_data <- data %>% select(genus_species, mite_load, ode_mass_g) %>% filter(complete.cases(.)) %>% group_by(genus_species) %>% mutate(avg_mass = mean(ode_mass_g)) %>% droplevels()
#plot_dat$genus_species <- as.factor(plot_dat$genus_species)
#plot_dat <- as.data.frame(plot_dat) %>% mutate(genus_species = fct_relevel(genus_species, plot_tree$tip.label[ordered_tips]))

# it will be nice to include the mean in this figures 
lice_load_plot <- ggplot(ectos_birds_dff, aes(x = total_lice, y =species_jetz)) +
  coord_cartesian(clip = "off") +
  geom_jitter(alpha=0.4, col="darkcyan")+
  scale_x_continuous(breaks=c(0, 1, 5, 10,15, 20, 30,40, 50)) +
  geom_point(data=ectos_birds_dff_mean, shape=124, size=1.5,aes(y=species_jetz, x=mean_nf_mites))+
  #scale_x_continuous(trans="log1p",breaks=c(0, 1, 2,5, 10, 25, 50, 100, 200, 300)) +
  scale_y_discrete(limits=order$`phylo$tip.label`) +
  ylab("") + xlab("Lice per individual") +
  theme_ridges(center_axis_labels = TRUE) +
  theme_classic(10)+
  theme(panel.grid.major.y = element_line( size=.05, color="grey"))

lice_load_plot[["data"]][["species_jetz"]]

mites_load_plot <- ggplot(data = ectos_birds_dff, aes(y=species_jetz, x=total_no_feathers_mites)) + 
  geom_jitter(alpha=0.4, col="olivedrab") + #geom_boxplot(outlier.alpha=0) +
  geom_point(data=ectos_birds_dff_mean, shape=124, size=1.5,aes(y=species_jetz, x=mean_nf_mites))+
  scale_x_continuous(breaks=c(0, 1, 5, 10, 20, 30,40, 50, 100, 200, 300)) +
  scale_y_discrete(limits=order$`phylo$tip.label`) +
  theme_ridges(center_axis_labels = TRUE) + 
  ylab("") +
  xlab(" Mites per individual (Non_feather mites)") +
  theme_classic(10)+
  theme( axis.text.y=element_blank(),
         axis.title.y=element_blank(),
         panel.grid.major.y = element_line( size=.05, color="grey"))


#  remove this part to make sure teh species are in the same orderto eliminate the species names in the axes once we know the order is correct use
theme(
  axis.text.y=element_blank(),
  axis.title.y=element_blank())


# The phylogenetic plot with prevalence 

#tree_plot <- ggtree(phylo, ladderize=FALSE) + geom_tiplab() + ggplot2::xlim(0, 450)

ColorPalette <- brewer.pal(n = 9, name = "YlGnBu")
list.names=setNames(ectos_birds_dff_PA_species$ectos_prevalence, ectos_birds_dff_PA_species$species_jetz)
fmode<-as.factor(setNames(ectos_birds_dff_PA_species$sociality,ectos_birds_dff_PA_species$species_jetz))
tree_plot_ectos<- contMap(phylo, list.names, direction = "leftwards", plot=FALSE)
#object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))
object_color<-setMap(tree_plot_ectos, ColorPalette)
tree_plot_sociality<-dotTree(phylo,fmode,colors=setNames(c("yellow","black"), c("1","0")),ftype="i",fsize=0.5, lwd=4) 
plot(tree_plot_sociality)
# 
png("figures/figures_manuscript/Fig1.PhyloTree_prevalence>4.png", width = 2500, height = 3100, res = 300, units = "px")
plot(dotTree(phylo,fmode,colors=setNames(c("coral3","black"),c("1","0")),ftype="i",fsize=0.5, lwd=4),text(x=10,y=-5,"Mixed-species flocks",pos=1))
plot(object_color$tree,colors=object_color$cols,add=TRUE,ftype="off",lwd=5,fsize=0.5,
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
add.color.bar(10, object_color$cols, title = "", lims = object$lims, digits = 3, prompt=FALSE,x=70,y=-5, lwd=4,fsize=1,subtitle="Ectoparasites Prevalence",pos=4)
dev.off()

unique( ectos_birds_dff_PA_species$species_jetz)

phylo_mite_lice_plot <- lice_load_plot %>% insert_right(mites_load_plot) 


#ggsave("figures/figures_manuscript/Fig1b__mite_lice_plot.png", plot=phylo_mite_lice_plot, height=10, width=12, units="in")




