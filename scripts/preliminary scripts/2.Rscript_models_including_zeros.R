
###_###_###
# ALL THE MODELS
###_###_###

#DATA
ectos_df<-read.csv("data/data_analyses/7.dff_ectos_pres_abs.csv")# data on prevalence FOR THE MODEL
phylogeny_prevalence<- read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex") 

###_####_####
ectos_dfff<-read.csv("data/data_analyses/7.dff_ectos_pres_abs.csv")# data on prevalence # should have same number of rows than teh phylogeny 
phylogeny_prevalence<- read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex") 
mean_lice_abundance<-read.csv("data/data_analyses/7.dff_lice_abundance_means.csv") %>% select(species_jetz, mean_lice)
mean_mites_abundance<-read.csv("data/data_analyses/7.dff_mites_abundance_means_non_feathers.csv") %>% select(species_jetz, mean_mites)
ectos_dff<-left_join(ectos_dfff, mean_lice_abundance, by="species_jetz")
ectos_dff$mean_lice[is.na(ectos_dff$mean_lice)] <- 0
ectos_df<-left_join(ectos_dff, mean_mites_abundance, by="species_jetz")
ectos_df$mean_mites[is.na(ectos_df$mean_mites)] <- 0

###_####_####

###_
# Keep data from Manu only ( Since I am not sure about iquitos metodology of parasite extraction)
ectos_df<-ectos_df %>% filter(elevation_cat!="lowland_iquitos",elevation_cat!="other_iquitos")

# Make sure the tips and the names on the file coincide
phylogeny_prevalence$edge.length  
phylogeny_prevalence$tip.label
is.binary(phylogeny_prevalence)

tips<- as.data.frame(phylogeny_prevalence$tip.label)
names<-as.data.frame(ectos_df$species_jetz)
anti_join(tips,names, by=c("phylogeny_prevalence$tip.label"="ectos_df$species_jetz"))

ectos_df<-ectos_df %>% distinct( species_jetz,.keep_all = TRUE) # remove the species that are duplicated because they ocur at two elevations


# Important!!!!!Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_prevalence$tip.label,rownames(ectos_df)),]

# Re-strudture the data
# Make sure variables are in teh right format, random effects should be factors
#We need to aling the variable name and the structure to the names in the column tip.label used for the phylogeny?

ectos_df <-ectos_df  %>% mutate_at("species_jetz", str_replace, " ", "_")

ectos_df$elevation_cat<-as.factor(ectos_df$elevation_cat)
ectos_df$foraging_cat<-as.factor(ectos_df$foraging_cat)
ectos_df$species_jetz<-as.factor(ectos_df$species_jetz)
ectos_df$sociality<-as.numeric(ectos_df$sociality)
ectos_df$elevation_midpoint<-as.numeric(ectos_df$elevation_midpoint)
ectos_df$sociality<-as.factor(ectos_df$sociality)

#as.data.frame(unique(ectos_df$species_jetz)) %>% View()

str(ectos_df) # should have same number of rows than teh phylogeny 

is.ultrametric(phylogeny_prevalence)

#MODEL FITTING AND EVALUALLING MODEL FIT GRAPHICALLY

###_###_###_ Important for model fitting
#predict ( model) # plot predicted vlues and confidence limits on the logit scale, the points are not logit tranformed data, they are "working values to fit the model.....dificult to interpret 
#visreg(model) # Similarly  plot predicted values and confidence limits on the logit scale, dificult to interpret
#fitted (model) # willl give us the predicted values transformed to the original scale 
#visreg(model, scale="response") # willl give us the predicted values transformed to the original scale 
###_###_###_ ###_###_###_ ###_###_###_ 


###_###_###_ 
# Using PGLMM [ Lets correct for phylogeny ]
###_###_###_ 

#I am not sue about including foraging cat since that some how is included in teh variation per species  ( Removed) 
#elevation_cat # only hs three categories so i am not sure I can used as a ranmod effect
# we revomed (1|foraging_cat) because it was not significant in individual models 

unique(ectos_df$species_jetz)
mean(ectos_df$proportion_ectoparasites) # mean prevalnece

#lm(proportion_ectoparasites~1, data=ectos_df) # even the average is a linear regression INTERESTING

str( ectos_df)
ectos_lice_abundance_mean <-  phyr::pglmm(mean_lice~sociality+sample_size+(1|species_jetz__)+(1|elevation_cat), #+elevation_midpoint
                                      data = ectos_df, 
                                      family = "gaussian",
                                      cov_ranef = list(species_jetz= phylogeny_prevalence), #class phylo
                                      #bayes = TRUE,
                                      REML = TRUE, 
                                      verbose = TRUE,
                                      s2.init = .25) # what is this last parameter for



summary(ectos_lice_abundance_mean)
predict(ecto_prevalence_pglmm, type="response")



#png("figures/figures_manuscript/Fig1b.Prevalence_ecto_output_predicted.png", width = 3000, height = 3000, res = 300, units = "px")
#plot_data(ectos_lice_abundance_mean,sp.var ="species_jetz", site.var ="sociality",predicted=TRUE)
#dev.off()



# Underestanding the summary of the random effects
#The random effect with the largest variance and standard variation is the one with the strongest effect, in our case the phylogenetic effect,
# this implies that the parasites  prevalence is more similar in closely related species 

# Explore if the random effects are important? 
# One way to get an idea is to run a likelihood ratio test on the random effect. 
# This can be achieved by using the pglmm_profile_LRT() function, at least for binomial models ( copied from https://daijiang.github.io/phyr/articles/phyr_example_empirical.html)

phyr::pglmm_profile_LRT(ectos_lice_abundance_mean, re.number = 1) ## here we need to specify qith random effect are we evaluating in the order we put them in the model species is number 1, phylogeny2, and elevation 3

LRTest <- sapply(1:3, FUN = function(x) phyr::pglmm_profile_LRT(ectos_lice_abundance_mean, re.number = x))
colnames(LRTest) <- names(ectos_lice_abundance_mean$ss)
t(LRTest)



###_###_###
# PLOTS PREVALENCE [Presence Absence]
###_###_###

# PLotting the model predictions 
#### How well does PGLMM it predict the data
newdata <- data.frame(elevation_cat = ectos_df$elevation_cat,
                      species_jetz = ectos_df$species_jetz,
                      sociality = ectos_df$sociality,
                      sample_size = ectos_df$sample_size)

predictions<- predict(ectos_lice_abundance_mean,newdata = newdata, type = "response" ) ##newdata = newdata

# or use fitedd instead
#predictions<- fitted(ecto_prevalence_pglmm,newdata = newdata) ##newdata = newdata

ectos_df_predicted <- cbind(ectos_df, predictions)

str(ectos_df_predicted)

# lets calculae the mean of the predited values on the untransformed scale
predictions_summary<- ectos_df_predicted %>% 
  group_by(sociality) %>%      
  dplyr::summarise(mean_lice_predicted= mean(Y_hat), sd=sd(Y_hat), n = n()) %>% 
  mutate(se= sd/(sqrt(n)))

colnames(predictions_summary) <- c("sociality", "mean_lice_predicted", "sd", "n", "se")

head(predictions_summary)

# make a plot of model predictions (that also shows data)
png("figures/figures_manuscript/Fig2_lice_abundance_means_pglmm_w_zeros.png", width = 3000, height = 3000, res = 300, units = "px")
ggplot(data = ectos_df_predicted, aes(x = sociality, y = mean_lice))+
  # geom_point(data = ectos_df, aes(x=sociality, y = proportion_ectoparasites),color="grey",size=2)+
  geom_jitter(data = ectos_df_predicted, aes(x=sociality, y = mean_lice),color="grey",size=3,width = 0.07)+
  geom_segment(data = predictions_summary, aes(x = sociality, y = mean_lice_predicted, xend = sociality, yend =mean_lice_predicted+sd, color="red"),show_guide = FALSE)+
  geom_segment(data = predictions_summary, aes(x = sociality, y = mean_lice_predicted, xend = sociality, yend =mean_lice_predicted-sd, color="red"),show_guide = FALSE)+
  #geom_jitter(data = ectos_df_predicted, aes(x=sociality, y = Y_hat), color="red", size=4,shape=19,width = 0.07)+
  geom_point(data = predictions_summary, aes(x=sociality, y = mean_lice_predicted), color="red", size=4,shape=19)+
  scale_y_continuous("Ectoparasites prevalence", limits = c(0,50)) +
  scale_x_discrete("Sociality")+
  geom_hline(yintercept = 3.3, linetype = "dashed")+
  theme_classic(40)
dev.off()

mean(ectos_df$mean_lice)


###_###_###_###_###_###_###_
# Combining both plots
###_###_###_###_###_###_###_


phylo=phylogeny_prevalence<- read.nexus("data/phylo_data/consensus/1_consensus_birdtreeManu_ectos_prevalence.nex") 

##### Finding incongruences
tips<- as.data.frame(phylo$tip.label)
names<-as.data.frame(ectos_df$species_jetz)
unique(mites_df_abundance$species_jetz)

a<-tips %>% arrange(desc(phylo$tip.label)) %>% mutate(name=phylo$tip.label) %>% select(name)
b<-names %>%arrange(desc(ectos_df$species_jetz))%>% mutate(name=ectos_df$species_jetz)%>% select(name)

tip<-as.list(setdiff(a,b))
print(tip)
phylo<-drop.tip (phylo, tip$name)



# Make sure the names  are in the same order in the phylogeny and in the traits
rownames(ectos_df) <- ectos_df$species_jetz # first make it the row names 
ectos_df<- ectos_df[match(phylogeny_mites$tip.label,rownames(ectos_df)),]

list.names=setNames(ectos_df$mean_lice, ectos_df$species_jetz)


###_###_###_###_###_###_###_
# Combining both plots
###_###_###_###_###_###_###_

#ColorPalette <- brewer.pal(n = 4, name = "YlGnBu")
ColorPalette <- brewer.pal(n = 8, name = "Paired")

fmode<-as.factor(setNames(ectos_df$sociality,ectos_df$species_jetz))
object = contMap(phylo, list.names, direction = "leftwards", plot=FALSE)
#object_color<-setMap(object, c("snow3","darkslategray3","dodgerblue","darkolivegreen3","goldenrod1"))
object_color<-setMap(object, ColorPalette)

png("figures/figures_manuscript/Fig2a.Sociality_and_lice_abundance_phylotree_mean_w_zeros.png", width = 2500, height = 3100, res = 300, units = "px")
plot(dotTree(phylo,fmode,colors=setNames(c("red","black"),c("1","0")),ftype="i",fsize=0.5, lwd=4),text(x=10,y=-5,"Mixed-species flocks",pos=1))

plot(object_color$tree,colors=object_color$cols,add=TRUE,ftype="off",lwd=5,fsize=0.5,
     xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
     ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)

add.color.bar(9, object_color$cols, title = "", lims = object$lims, digits = 3, prompt=FALSE,x=70,y=-5, lwd=4,fsize=1,subtitle="mean lice abundance",pos=4)
dev.off()


