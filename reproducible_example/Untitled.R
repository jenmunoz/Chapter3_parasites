#DATA


data_df<-read.csv("7.df_data_example.csv", na.strings =c("","NA")) %>% na.omit() 
phylo<-read.nexus("1_consensus_birdtree.nex") 
unique(data_df$species_jetz )


#STRUCTURE
ectos_df$species_jetz<-as.factor(ectos_df$species_jetz)
ectos_df$elevation<-as.numeric(ectos_df$elevation)
ectos_df$sociality<-as.factor(ectos_df$sociality)
ectos_df$Powder.lvl<-as.factor(ectos_df$Powder.lvl)

# MODEL 
pglmm_bayes <-phyr::pglmm(total_parasites~sociality+scale(elevation)+(1|species_jetz__)+(1|Powder.lvl),
                                         data = data_df, 
                                         family ="zeroinflated.poisson", #POISSON
                                         cov_ranef = list(species_jetz= phylo), #class phylo
                                         bayes = TRUE,
                                         verbose=FALSE,
                                         prior = "inla.default") # consider using  add.obs.re = T

#ASSUMPTIONS CHECK ( )
simulationOutput<-DHARMa::simulateResiduals(fittedModel= pglmm_bayes , plot= TRUE, re.form = NULL ) #quantreg=T


#PLOT 
summary(pglmm_bayes)
phyr::plot_bayes(pglmm_bayes)


