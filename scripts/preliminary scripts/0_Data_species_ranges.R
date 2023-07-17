# Extracting data ranges for Santi 

# Detections
manu_detections_jetz<-read.csv("data/1.Manu_bird_species_detections_tidy_taxonomy_18052022_jetz_taxo_included.csv") %>% filter(database_decision=="include")

unique(manu_detections_jetz$database_decision)

list_manu_jetz<-read.csv("data/1.list_manu_species_jetz_taxonomy.csv") %>% select( species_jetz, species_from_detections,species_taxonomy_SACC_2021) %>% 
  distinct(species_taxonomy_SACC_2021, .keep_all = TRUE) %>% 
  mutate_at("species_jetz", str_replace, " ", "_") %>% 
  mutate_at("species_taxonomy_SACC_2021", str_replace, " ", "_") %>% 
  mutate_at("species_from_detections", str_replace, " ", "_") 


# the ranges observed
observed_species_range_limits<-manu_detections_jetz %>% group_by(species_jetz) %>% 
  summarise(low_limit=min(elevation), high_limit=max(elevation),n=n()) %>% 
  mutate_at("species_jetz", str_replace, " ", "_") 

# Add more taxonomy fields

observed_species_range_limits_taxo<-inner_join(observed_species_range_limits,list_manu_jetz, by="species_jetz")
write.csv(observed_species_range_limits_taxo,"data/data_analyses/1.observed_species_range_limits_Manu_birds.csv")
