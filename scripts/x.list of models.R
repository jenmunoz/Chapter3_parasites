#models

#Prevalence

ecto_p_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions.RDS")
ecto_p_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/M1P_model_prevalence_b_brms_phylo_multiple_obs_no_interactions_priors_default_intercept.RDS")

ecto_p_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/M2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE.RDS")
ecto_p_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/2PND.model_prevalence_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")


#lICE
zinb_a_lice_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions.RDS") # lice abundance default priors             !!!! warrning 20 pareto k are high!!!! what to do?
zinb_a_lice_brms_bayes_all_interactions_priors<-readRDS("data/data_analyses/model_selection/M1L.model_prevalence_zinb_brms_LICE_ABUNDANCE_phylo_multiple_obs_all_interactions_priors.RDS")# Lice abundance weakly informative priors !!!! warrning 20 pareto k are high!!!! what to do?

zinb_a_lice_brms_bayes_no_int_degree<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_zinb_brms_ABUNDANCE_LICE_phylo_multiple_obs_no_interactions_DEGREE.RDS")
zinb_a_lice_brms_bayes_no_int_degree_prior<-readRDS("data/data_analyses/model_selection/M2LND.model_lICE_ABUNDANCE_b_brms_phylo_multiple_obs_no_interactions_DEGREE_prior.RDS")


#MITES
zinb_a_nf_mites_brms_bayes_no_int<-readRDS("data/data_analyses/model_selection/M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions.RDS")# Mites abundance default priors 
zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("data/data_analyses/model_selection/1NFM.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior.RDS")# Mites abundance weakly informative priors 


 # Mites abundance weakly by degree informative priors 

