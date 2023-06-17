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
  select(elevation, species_jetz, Powder.lvl,ectoparasites_PA, foraging_cat,sociality,year_seasonality, mass_tidy_species, total_lice, total_no_feathers_mites ) %>% 
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


# #  Figure 2 #### --------------------------------------------------------------
###_###_###_###_###

# # Figure 2 a  INFECTION ***Selected***  ----------------------------

# model 
selected_ecto_infection_brms_bayes_no_int<-readRDS("results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED_antfollowers_included.RDS")

# plots
color_scheme_set("brightblue")

#Credible intervals
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot_intervals<-mcmc_plot(selected_ecto_infection_brms_bayes_no_int,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS INFECTION ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2a1.Infection_sociality_credible_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# Conditional effects for varaibles tahta are important predictors
#prevalence_seasonality<-conditional_effects(selected_ecto_infection_brms_bayes_no_int, "year_seasonality",points=TRUE, rug=TRUE)
conditional<-conditional_effects(selected_ecto_infection_brms_bayes_no_int)

prevalence_seasonality<-plot(conditional, plot = FALSE)[["year_seasonality"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Day of the year")+
  ylab("Ectoparasite infection")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2a4.Infection_seasonality.pdf", plot=prevalence_seasonality , height=10, width=10, units="in")

#plot(conditional_effects(selected_ecto_infection_brms_bayes_no_int), ask = FALSE)

bayes_R2(selected_ecto_infection_brms_bayes_no_int) # R2 0.1529

# # # Figure 2 a Infection Networks ---------------------------------------
# model degree
selected_ecto_p_brms_bayes_no_int_degree_prior<-readRDS("results/selected_models/4_1_M2PND_model_INFECTION_b_brms_phylo_multiple_obs_no_interactions_degree_prior_SELECTED_ind_mass_scaled.RDS")


estimates_plot_intervals<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION ~DEGREE ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2a2.Infection_sociality_degree_credible_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# Conditional effects for varaibles that are important predictors
conditional<-conditional_effects(selected_ecto_p_brms_bayes_no_int_degree_prior)

prevalence_seasonality_degree<-plot(conditional, plot = FALSE)[["year_seasonality"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Day of the year")+
  ylab("Ectoparasite infection")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2a4.Infection_seasonality_degree.pdf", prevalence_seasonality_degree , height=10, width=10, units="in")

prevalence_degree<-plot(conditional, plot = FALSE)[["degree"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Sociality degree")+
  ylab("Ectoparasite infection")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2a4.Infection_degree.pdf", prevalence_degree , height=10, width=10, units="in")

#plot(conditional_effects(selected_ecto_infection_brms_bayes_no_int), ask = FALSE)

bayes_R2()
# Model strength
selected_ecto_p_brms_bayes_no_int_strength_prior<-readRDS("results/selected_models/4_1_2_M2PNS_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior_SELECTED_ind_mass_scaled.RDS")
bayes_R2(selected_ecto_p_brms_bayes_no_int_strength_prior)

#PLOTS
color_scheme_set("brightblue")

estimates_plot_intervals<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION ~STRENGTH ")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2a3.Infection_sociality_strength_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# Conditional effects for varaibles that are important predictors
conditional<-conditional_effects(selected_ecto_p_brms_bayes_no_int_strength_prior)

prevalence_seasonality_strength<-plot(conditional, plot = FALSE)[["year_seasonality"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Day of the year")+
  ylab("Ectoparasite infection")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2a4.Infection_seasonality_strength.pdf", prevalence_seasonality_strength , height=10, width=10, units="in")

prevalence_stregth<-plot(conditional, plot = FALSE)[["degree"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Sociality strength")+
  ylab("Ectoparasite infection")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2a4.Infection_strenght.pdf", prevalence_stregth , height=10, width=10, units="in")



# # # Figure 2 b Lice abundance -------------------------------------------
# the model

selected_zinb_a_lice_brms_bayes_no_int_priors<-readRDS("results/selected_models/3_M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_ind_mass_scaled_SELECTED_antbirds_included.RDS")
###_###_###_##
#PLOTS
###_###_###_##

color_scheme_set("teal")

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ZINB LICE ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2b1.Lice_abundance_sociality_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# Conditional effects for varaibles tahta are important predictors
#prevalence_seasonality<-conditional_effects(selected_ecto_infection_brms_bayes_no_int, "year_seasonality",points=TRUE, rug=TRUE)
conditional<-conditional_effects(selected_zinb_a_lice_brms_bayes_no_int_priors)

lice_body_mass<-plot(conditional, plot = FALSE)[["mass_tidy_species"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Host body mass")+
  ylab("Lice abundance")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2b4.lice_and_body_mass.pdf", plot=lice_body_mass , height=10, width=10, units="in")

#plot(conditional_effects(selected_ecto_infection_brms_bayes_no_int), ask = FALSE)

bayes_R2() # R2 0.1529

# # # Figure 2 b Lice abundance Networks-------------------------------------------

####
# Model degree
###
selected_zinb_a_lice_brms_bayes_no_int_degree_prior<-readRDS ("results/selected_models/4_3_M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior_ind_mass_scaled.RDS")
###The best model is the zero inflated negatve binomialr 

# plots 
color_scheme_set("teal") 

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE DEGREE ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2b2.Lice_abundance_sociality_degree_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()
####
# Model strength
###

selected_zinb_a_lice_brms_bayes_no_int_strength_prior<-readRDS("results/selected_models/4_3_2_M2LNS.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_STRENGTH_prior_ind_mass_scaled.RDS")
###_###_###_###_###_###_###_###_###_###_
###The best model is the zero inflated negatve binomialr 

# plots 
color_scheme_set("teal") 

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE STRENGTH ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2b3.Lice_abundance_sociality_strength_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()
# #  Figure 2 c Mites (no-feather) abundance ---------------------------
# model 
selected_zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_SELECTED_antfollowers_included.RDS")

# plots 
color_scheme_set("green") 

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2c1.Mites_abundance_sociality_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# # Conditional effects for varaibles tahta are important predictors
#prevalence_seasonality<-conditional_effects(selected_ecto_infection_brms_bayes_no_int, "year_seasonality",points=TRUE, rug=TRUE)
conditional<-conditional_effects(selected_zinb_a_nf_mites_brms_bayes_no_int_prior)

mites_elevation<-plot(conditional, plot = FALSE)[["elevation"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Elevation")+
  ylab("Mites abundance")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2c4.Mice_and_elevation.pdf", plot=mites_elevation , height=10, width=10, units="in")

#plot(conditional_effects(selected_ecto_infection_brms_bayes_no_int), ask = FALSE)

# # # Figure 2 c Mites (no-feather) abundance Networks---------------------------

# model degree
selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-readRDS( "results/selected_models/4_3_3_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior_individual_mass_scaled.RDS")
###_###_###_##
#PLOTS
###_###_###_##

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~DEGREE")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2c2.Mites_abundance_sociality_degree_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# # Conditional effects for varaibles tahta are important predictors
#prevalence_seasonality<-conditional_effects(selected_ecto_infection_brms_bayes_no_int, "year_seasonality",points=TRUE, rug=TRUE)
conditional<-conditional_effects(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior)

mites_elevation_degree<-plot(conditional, plot = FALSE)[["elevation"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Elevation")+
  ylab("Mites abundance")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2c4.Mice_and_elevation_degreemodel.pdf", plot=mites_elevation , height=10, width=10, units="in")


###_###_###_###_###_###_###_###_###_###_
#Model strength 
###_###_###_###_###_###_###_###_###_###_

selected_zinb_a_nf_mites_brms_bayes_no_int_strength_prior<-readRDS( "results/selected_models/4_3_4_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_strength_prior.RDS")

#PLOTS
###_###_###_##
color_scheme_set("green")

estimates_plot_intervals<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~STRENGTH")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2c3.Mites_abundance_sociality_degree_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()


# # Figure 2 d  PREVALENCE ***Selected***--------------------------------------------------------------
# model
ecto_p_brms_bayes_no_int_species_priors_zobi<-readRDS("results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_zobi_antbirds_included.RDS")

# Plots
color_scheme_set("orange") 

# Credible intervals
estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_zobi,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2d1.Prevalence_sociality_credible_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# Conditional effects [ No predictor was considered important ]
conditional<-conditional_effects(ecto_p_brms_bayes_no_int_species_priors_zobi)


# # #Figure 2 d Prevalence Networks  ----------------------------------------
# model degree
ecto_p_brms_bayes_no_int_species_priors_degree_zobi<-readRDS("results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_DEGREE_zobi.RDS")

# Plots
color_scheme_set("orange") 

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_degree_zobi,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~DEGREE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2d2.Prevalence_socialitydegree_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

# Conditional effects [ No predictor was considered important ]
conditional<-conditional_effects(ecto_p_brms_bayes_no_int_species_priors_zobi)

prevalence_conditional<-plot(conditional, plot = FALSE)[[""]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Day of the year")+
  ylab("Ectoparasite infection")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2d5.Prevalence_conditional.pdf", plot=prevalence_seasonality , height=10, width=10, units="in")

###_###_###_###
# model strength 
ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi<-readRDS("results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_STRENGHT_zobi.RDS")

# Plots
color_scheme_set("orange") 

estimates_plot_intervals<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~STRENGHT")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2d3.Prevalence_socialitystrength_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()



# # # Figure 2e Lice richness --------------------------------------------

selected_poisson_lice_diversity_sociality_no_int_priors<-readRDS( "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_NO_truncated.RDS")

estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_sociality_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY ")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2e1.LiceRichness_sociality_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()




# # # Figure 2e Lice richness Networks--------------------------------------------

# I am using the non truncated model because truccated keepp giving me many divergent transitions
# model degree
selected_poisson_lice_diversity_degree_no_int_priors<-readRDS("results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_degree.RDS")

color_scheme_set("red") 

estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_degree_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY~DEGREE ")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2e2.LiceRichness_sociality_degree_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

bayes_R2(selected_poisson_lice_diversity_degree_no_int_priors)

# model strength
###
#Strength 
###
#the model

selected_poisson_lice_diversity_w_degree_no_int_priors<-readRDS( "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_w_degree.RDS")

color_scheme_set("red") 

# intervals
estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_w_degree_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY ~W_degree")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2e3.LiceRichness_socialitySTRENGTH_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()




# # Figure 3  Supplementary mat #### --------------------------------------------------------------
###_###_###_###_###
# Figure 3 a PREVALENCE/Infection ####### 1.3 ***Selected*** model prevalence ectos INFECTION (included mass) ----------------------------

# plots
color_scheme_set("brightblue")
# model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3a.Infection_sociality_convergence.pdf", width =20, height =10)
plot(selected_ecto_infection_brms_bayes_no_int)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3a.Infection_sociality_fit.pdf", width =10, height =10)
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

pdf(file="figures/figures_pdf_manuscript/FigureS3a.Infection_sociality_estimates.pdf", width =10, height =10)
estimates_plot
dev.off()


#  #Figure 3 a Infection Networks -------------------------------------------

#PLOTS
color_scheme_set("brightblue")

# model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3a2.Infection_sociality_degree_convergence.pdf", width =15, height =10)
plot(selected_ecto_p_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3a2.Infection_sociality_degree_fit.pdf", width =15, height =10)
pp_check(selected_ecto_p_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION~DEGREE")+
  theme_minimal(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3a2.Infection_sociality_degree_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()

###_###_###_
# MODEL STRENGHT
###_###_###_

color_scheme_set("brightblue")

# model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3a3.Infection_sociality_strength_convergence.pdf", width =15, height =10)
plot(selected_ecto_p_brms_bayes_no_int_strength_prior)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3a3.Infection_sociality_strength_fit.pdf", width =10, height =10)
pp_check(selected_ecto_p_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 100)+ xlim(0,5)
dev.off()


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model
estimates_plot<-mcmc_plot(selected_ecto_p_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="INFECTION~STRENGTH")+
  theme_minimal(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3a3.Infection_sociality_strength_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()


# Figure 3 b Lice abundance  ------------------------------------------
selected_zinb_a_lice_brms_bayes_no_int_priors<-readRDS ("results/selected_models/3_M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_ind_mass_scaled_SELECTED_antbirds_included.RDS")

color_scheme_set("teal")

# model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3b1.Lice_abundance_sociality_convergence.pdf", width =15, height =10)
plot(selected_zinb_a_lice_brms_bayes_no_int_priors)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3b1.Lice_abundance_sociality_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_lice_brms_bayes_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

# estimates
estimates_plot<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ZINB LICE ABUNDANCE ")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3b1.Lice_abundance_sociality_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()

# # Figure 3 b Lice abundance Networks -----------------------------------------------

# DEGREE 

# plots 
color_scheme_set("teal") 
# poisson 
#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3b2.Lice_abundance_sociality_degree_convergence.pdf", width =15, height =10)
plot(selected_zinb_a_lice_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3b2.Lice_abundance_sociality_degree_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_lice_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE DEGREE ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3b2.Lice_abundance_sociality_degree_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()

#conditional effects

conditional<-conditional_effects(selected_zinb_a_lice_brms_bayes_no_int_degree_prior)

lice_body_mass_degree<-plot(conditional, plot = FALSE)[["mass_tidy_species"]] +
  scale_color_grey() +
  scale_fill_grey() +
  xlab("Host body mass")+
  ylab("Lice abundance")+
  theme_classic(30)
ggsave("figures/figures_pdf_manuscript/Figure2b4.lice_and_body_mass_degreemodel.pdf", lice_body_mass_degree , height=10, width=10, units="in")


# STRENGTH 

#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3b3.Lice_abundance_sociality_strength_convergence.pdf", width =15, height =10)
plot(selected_zinb_a_lice_brms_bayes_no_int_strength_prior)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3b3.Lice_abundance_sociality_strength_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_lice_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 50)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_zinb_a_lice_brms_bayes_no_int_strength_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle =" ZINB LICE ABUNDANCE STRENGTH ")+
  xlim(-5,5)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3b3.Lice_abundance_sociality_strength_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()




# # # Figure 3 c Mites abundance ------------------------------------------

color_scheme_set("green") 

#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3c1.Mites_abundance_sociality_convergence.pdf", width =15, height =10)
plot(selected_zinb_a_nf_mites_brms_bayes_no_int_prior)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3c1.Mites_abundance_sociality_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES
estimates_plot<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals ", subtitle ="ZINB MITES ABUNDANCE  ")+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3c1.Mites_abundance_sociality_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()

# # # Figure 3 Mites abundance networks -----------------------------------

# degree
selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior<-readRDS( "results/selected_models/4_3_3_M2MND_model_MITES_ABUNDANCE_nb_brms_phylo_multiple_obs_no_interactions_degree_prior_individual_mass_scaled.RDS")

color_scheme_set("green")

# model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3c2.Mites_abundance_sociality_degree_convergence.pdf", width =15, height =10)
plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3c2.Mites_abundance_sociality_degree_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES
estimates_plot<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~DEGREE")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3c2.Mites_abundance_sociality_degree_intervals.pdf", width =15, height =10)
estimates_plot
dev.off()

# strength
# model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3c3.Mites_abundance_sociality_degree_convergence.pdf", width =15, height =10)
plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3c3.Mites_abundance_sociality_degree_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES

estimates_plot<-mcmc_plot(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title=" Posterior distributions with medians and 95% intervals", subtitle ="ZINB MITES ABUNDANCE ~STRENGTH")+
  theme_classic(30)+
  xlim(-5,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey20")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3c3.Mites_abundance_sociality_degree_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()



# Figure 3 d PREVALENCE ---------------------------------------------------
#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3d.Prevalence_sociality_convergence.pdf", width =10, height =10)
plot(ecto_p_brms_bayes_no_int_species_priors_zobi)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3d.Prevalence_sociality_fit.pdf", width =10, height =10)
pp_check(ecto_p_brms_bayes_no_int_species_priors_zobi, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_zobi,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3d.Prevalence_sociality_estimates.pdf", width =10, height =10)
estimates_plot
dev.off()

# Figure 3 d Prevalence networks -------------------------------------------

# degree

#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3d2.Prevalence_socialitydegree_convergence.pdf", width =15, height =10)
plot(ecto_p_brms_bayes_no_int_species_priors_degree_zobi)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3d2.Prevalence_socialitydegree_fit.pdf", width =15, height =10)
pp_check(ecto_p_brms_bayes_no_int_species_priors_degree_zobi, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()


#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model


estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_degree_zobi,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~DEGREE")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3d2.Prevalence_socialitydegree_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()

###
# Strength
###
#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3d2.Prevalence_socialitystrength_convergence.pdf", width =15, height =10)
plot(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3d2.Prevalence_socialitystrength_fit.pdf", width =15, height =10)
pp_check(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi, type = "dens_overlay", ndraws = 100)+ xlim(0, 5)
dev.off()

#pp_check(ecto_p_brms_bayes, ndraws = 100)+ xlim(0, 5)  #  test for the model fit to the data .need to modify the scale of this plot posterior predictive checks, 100 random draws or distributions created by the model 
#pp_check(ecto_p_brms_bayes, type="bars", ndraws = 100)+ xlim(0, 20) 

#MODEL ESTIMATES
# Dots represent means of posterior distribution along with 95% CrIs, as estimated by the bmod5 model

estimates_plot<-mcmc_plot(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals", subtitle ="ECTOS PREVALENCE~STRENGHT")+
  theme_classic(30)+
  xlim(-2,5)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3d2.Prevalence_socialitystrength_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()















# # # FIgure 3 e Lice richness  -------------------------------------------
color_scheme_set("red") 

#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3e1.LiceRichness_sociality_convergence.pdf", width =15, height =10)
plot(selected_poisson_lice_diversity_sociality_no_int_priors)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigurS3e1.LiceRichness_sociality_fit.pdf", width =15, height =10)
pp_check(selected_poisson_lice_diversity_sociality_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()


#ESTIMATES
estimates_plot<-mcmc_plot(selected_poisson_lice_diversity_sociality_no_int_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="LICE DIVERSITY ")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3e1.LiceRichness_sociality_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()

# # # Figure 3 e Lice Richness NetworKs --------------------------------------------

#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3e2.LiceRichness_socialitydegree_convergence.pdf", width =15, height =10)
plot(selected_poisson_lice_diversity_degree_no_int_priors_trunc)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3e2.LiceRichness_socialitydegree_fit.pdf", width =15, height =10)
pp_check(selected_poisson_lice_diversity_degree_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_poisson_lice_diversity_degree_no_int_priors_trunc,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="LICE DIVERSITY ~DEGREE")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3e2.LiceRichness_socialitydegree_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()

###
#Strength 
###
#the model

selected_poisson_lice_diversity_w_degree_no_int_priors<-readRDS( "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_w_degree.RDS")

color_scheme_set("red") 

# intervals
estimates_plot_intervals<-mcmc_plot(selected_poisson_lice_diversity_w_degree_no_int_priors,prob=0.90, prob_outer=0.95,point_est = "mean",
                                    type="intervals") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="  LICE DIVERSITY ~W_degree")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/Figure2e3.LiceRichness_socialitySTRENGTH_intervals.pdf", width =15, height =10)
estimates_plot_intervals
dev.off()

#model convergence 
pdf(file="figures/figures_pdf_manuscript/FigureS3e3.LiceRichness_socialityStrength_convergence.pdf", width =15, height =10)
plot(selected_poisson_lice_diversity_w_degree_no_int_priors)
dev.off()

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3e3.LiceRichness_socialityStrength_fit.pdf", width =15, height =10)
pp_check(selected_poisson_lice_diversity_w_degree_no_int_priors, type = "dens_overlay", ndraws = 100)+ xlim(0, 20)
dev.off()

#ESTIMATES

estimates_plot<-mcmc_plot(selected_poisson_lice_diversity_w_degree_no_int_priors,prob=0.90, prob_outer=0.95,
                          type="areas") +
  labs(title="Posterior distributions with medians and 95% intervals ", subtitle ="LICE DIVERSITY ~W DEGREE")+
  xlim(-4,4)+
  theme_classic(30)+
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10")+
  xlab("Estimate")

pdf(file="figures/figures_pdf_manuscript/FigureS3e3.LiceRichness_socialityStrength_estimates.pdf", width =15, height =10)
estimates_plot
dev.off()


