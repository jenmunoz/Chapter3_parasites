#######################################################################################
### Chapter 3-parasites and flocking species                                     ###
### Part 1 Figures                                                                       ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
### Last update: Nov 2022                                                     ###
################################################################################
#  # Figure 1####  -------------------------------------
# Phylogeny and ectoparasites plots for figure 1
# ### #### #### #### #### #### #### ### 

ectos_birds_df<-read.csv("data/data_manuscript/3_dff_all_ectos_prevalence_abundance_diversity_individual_elevation_mass_FILE_TIDY.csv", na.strings =c("","NA")) %>% 
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality, year_seasonality, mass_tidy_species, total_lice, total_no_feathers_mites ) %>% 
  na.omit() 

names(ectos_birds_df)

#%>%  filter(species_jetz!="Premnoplex_brunnescens")


str(ectos_birds_df)

ectos_pres_abs<-ectos_birds_df %>% group_by(species_jetz ) %>% 
  summarise(ectoparasites_presence=(sum(ectoparasites_PA)), sample_size=(n()), elevation=mean(elevation), sociality=max(sociality), mass=max(mass_tidy_species), foraging_cat=first(foraging_cat))%>% 
  mutate(ectos_prevalence=ectoparasites_presence/sample_size) %>% 
  na.omit() 

ectos_birds_dff<-ectos_pres_abs %>% distinct(species_jetz,ectos_prevalence,.keep_all = TRUE)%>% 
  filter(sample_size>3) 

unique(ectos_birds_dff$ectos_prevalence)
#View(ectos_birds_dff)

### the phylogeny 

phylo<-read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex")  # This include speceis form manu and iquitos  so need to rpun the tree in the data processin section

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



# The phylogenetic plot with prevalence 

# ### ### ###
### Figure 1a
# ### ### ###
#tree_plot <- ggtree(phylo, ladderize=FALSE) + geom_tiplab() + ggplot2::xlim(0, 450)

ColorPalette <- brewer.pal(n = 9, name = "YlGnBu")
list.names=setNames(ectos_birds_dff$ectos_prevalence, ectos_birds_dff$species_jetz)
fmode<-as.factor(setNames(ectos_birds_dff$sociality,ectos_birds_dff$species_jetz))
tree_plot_ectos<- contMap(phylo, list.names, direction = "leftwards", plot=FALSE)
#object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))
object_color<-setMap(tree_plot_ectos, ColorPalette)
tree_plot_sociality<-dotTree(phylo,fmode,colors=setNames(c("yellow","black"), c("1","0")),ftype="i",fsize=0.5, lwd=4) 
plot(tree_plot_sociality)
# 
#png("figures/figures_pdf_manuscript/Fig1.PhyloTree_prevalence.pdf", width = 2500, height = 3100, res = 300, units = "px")

pdf(file="figures/figures_pdf_manuscript/Figure1a.PhyloTree_ectos_prevalence.pdf", width =10, height =10)
plot(dotTree(phylo,fmode,colors=setNames(c("#FF6633","#CCCCCC"),c("1","0")),ftype="i",fsize=0.8, lwd=4),text(x=10,y=-3,"Mixed-species flocks",pos=1,lwd=8,fsize=2))
plot(object_color$tree,colors=object_color$cols,add=TRUE,ftype="off",lwd=5,fsize=0.8,
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
add.color.bar(15, object_color$cols, title = "Ectoparasite Prevalence", lims = object$lims, digits = 3, prompt=FALSE,x=70,y=-3, lwd=10,fsize=1, subtitle="Ectoparasite Prevalence")
dev.off()

unique( ectos_birds_dff_PA_species$species_jetz)

#ggsave("figures/figures_manuscript/Fig1b__mite_lice_plot.png", plot=phylo_mite_lice_plot, height=10, width=12, units="in")


### individual plots 

# releveling
#plot_data <- data %>% select(genus_species, mite_load, ode_mass_g) %>% filter(complete.cases(.)) %>% group_by(genus_species) %>% mutate(avg_mass = mean(ode_mass_g)) %>% droplevels()
#plot_dat$genus_species <- as.factor(plot_dat$genus_species)
#plot_dat <- as.data.frame(plot_dat) %>% mutate(genus_species = fct_relevel(genus_species, plot_tree$tip.label[ordered_tips]))

# it will be nice to include the mean in this figures 
 
unique(ectos_birds_df$total_lice)
log1p(62)

names(ectos_birds_df)

lice_load_plot <- ggplot(ectos_birds_df, aes(x =log1p(total_lice), y =species_jetz)) +
  coord_cartesian(clip = "off") +
  geom_jitter(width = 0.02,alpha=0.6,size=3, col="darkcyan")+
  #scale_x_continuous(trans="log1p",breaks=c(0, 1, 5, 10,15, 20,50))+ # log1p allows the calcualtion of log of values >0
  #scale_x_continuous(trans="log1p",breaks=c(0, 1, 2,5, 10, 25, 50, 100, 200, 300)) +
  scale_y_discrete(limits=order$`phylo$tip.label`) +
  ylab("") + xlab("log1p (Lice per individual host)") +
  theme_ridges(center_axis_labels = TRUE) +
  theme_classic(10)+
  theme(panel.grid.major.y = element_line( size=.01, color="grey100"))

lice_load_plot[["data"]][["species_jetz"]]

mites_load_plot<- ggplot(data = ectos_birds_df, aes(x=log1p(total_no_feathers_mites), y=species_jetz)) + 
 coord_cartesian(clip = "off") +
 geom_jitter(width = 0.02, alpha=0.6,size=3, col="seagreen") + 
  #geom_boxplot(outlier.alpha=0) +
  #geom_point(data=ectos_birds_dff_mean, shape=124, size=1.5,aes(y=species_jetz, x=mean_nf_mites))+
  #scale_x_continuous(breaks=c(0, 1, 5, 10, 20, 30,40, 50, 100, 200, 300)) +
  scale_y_discrete(limits=order$`phylo$tip.label`) +
  theme_ridges(center_axis_labels = TRUE) + 
  ylab("") +
  xlab("log1p (Mites per individual host)") +
  theme_classic(10)+
  theme( axis.text.y=element_blank(),
         axis.title.y=element_blank(),
         panel.grid.major.y = element_line( size=.01, color="grey100"))


dff_lice_diversity<-read.csv("data/data_manuscript/3_dff_ectos_diversity_species_prevalence_abundance_diversity_elevation_mass_FILE_TIDY.csv",na.strings =c("","NA"))%>% filter(total_sample_size>9)
dff_lice_diversity$cumulative_richness[is.na(dff_lice_diversity$cumulative_richness)] = 0

richness_load_plot<- ggplot(data = dff_lice_diversity, aes(x=(cumulative_richness), y=species_jetz)) + 
  coord_cartesian(clip = "off") +
  geom_jitter(width = 0.01, alpha=0.6, size=3, col="orange3") + 
  #geom_boxplot(outlier.alpha=0) +
  #geom_point(data=ectos_birds_dff_mean, shape=124, size=1.5,aes(y=species_jetz, x=mean_nf_mites))+
  #scale_x_continuous(breaks=c(0, 1, 5, 10, 20, 30,40, 50, 100, 200, 300)) +
  scale_y_discrete(limits=order$`phylo$tip.label`) +
  theme_ridges(center_axis_labels = TRUE) + 
  ylab("") +
  xlab("Lice genera per host species") +
  theme_classic(10)+
  theme( axis.text.y=element_blank(),
         axis.title.y=element_blank(),
         panel.grid.major.y = element_line( size=.01, color="grey100"))
#  remove this part to make sure teh species are in the same orderto eliminate the species names in the axes once we know the order is correct use
theme(
  axis.text.y=element_blank(),
  axis.title.y=element_blank())
# the plot integrating the three plots

# ### ### ###
### Figure 1b
# ### ### ###
phylo_mite_lice_richness_plot <- mites_load_plot %>% insert_left(lice_load_plot) %>% insert_right(richness_load_plot)

ggsave("figures/figures_pdf_manuscript/Figure1b.phylo_mite_lice_richness_plot.pdf", plot=phylo_mite_lice_richness_plot , height=10, width=10, units="in")
phylo_mite_lice_richness_plot <- mites_load_plot %>% insert_left(lice_load_plot) %>% insert_right(richness_load_plot)


# # Figure 2 #### --------------------------------------------------------------
###_###_###_###_###

# Figure 2 a PREVALENCE ####### 1.3 ***Selected*** model prevalence ectos INFECTION (included mass) ----------------------------

# model 
selected_ecto_infection_brms_bayes_no_int<-readRDS("results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED.RDS")

# plots
color_scheme_set("blue")

#ECreible intervals
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot_intervals<-mcmc_plot(selected_ecto_infection_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS INFECTION ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3.Prevalence_sociality_credible)intervals.pdf", width =10, height =10)
estimates_plot_intervals
dev.off()

bayes_R2(selected_ecto_infection_brms_bayes_no_int)

# half student only allows positive values

# # Figure 3  Suplementary#### --------------------------------------------------------------
###_###_###_###_###
# Figure 3 a PREVALENCE ####### 1.3 ***Selected*** model prevalence ectos INFECTION (included mass) ----------------------------

# plots
color_scheme_set("blue")
# model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3.Prevalence_sociality_convergence.pdf", width =10, height =10)
plot(selected_ecto_infection_brms_bayes_no_int)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3.Prevalence_sociality_fit.pdf", width =10, height =10)
pp_check(selected_ecto_infection_brms_bayes_no_int, type = "dens_overlay", ndraws = 100) 
dev.off()

# model estimates density

estimates_plot<-mcmc_plot(selected_ecto_infection_brms_bayes_no_int,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS INFECTION ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3.Prevalence_sociality_estimates.pdf", width =10, height =10)
estimates_plot
dev.off()






conditional_effects(selected_ecto_infection_brms_bayes_no_int)
marginal_effects(selected_ecto_infection_brms_bayes_no_int)
plot( conditional_effects(selected_ecto_infection_brms_bayes_no_int), 
      points = TRUE, 
      point_args = list(width = .05, shape = 1))
bayes_R2() # R2 0.1529









