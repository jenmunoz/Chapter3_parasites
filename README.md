# Chapter3_parasites
This is the code for the chapter 3 of my PhD


## Sociality ##
Sociality status is based on flcok observations from Manu only, some particular species have changes in teh status

Lophotriccus pileatus::This species has been obserevd in flocks but those flcoks were incomplete and tehrefore not included in the analyses, I changed 0 to 1 in sociality manually

##File names and description ##

Raw data files:

[1]The raw files of bird doversity form Manu diversity 

0.flocks_manu_complete_18052022.csv
0.netting_manu_complete_18052022.csv

[2] The tidy files of manu diversity
1.Manu_bird_species_detections_tidy_taxonomy_18052022.csv
1.manu_bird_species_list_flocks.csv
1.Manu_bird_species_list_tidy_18052022.csv

[3] The functional traits of Manu birds
4.df_traits_manu_birds.csv

[4]  Files to correct taxonomy 
# We match the taxonomy with Jetz taxonomy that we aleady curated form manu data( see list of species that were manually matched)

0.species_manual_taxonomy_match.csv
4.list_manu_jetz_tax_missmatches_corrected.csv


These are the raw files of the ectoparasies without taxonomy corrected or any calculation
EctoparsiteRawData_All_07.22.2022.xlsx
7.ectoparasite_raw_pres_abs_07222022.csv
7.ectoparsite_raw_lice_abundance_07222022.csv
7.ectoparsite_raw_mites_abundance_07222022
7.ectoparsite_raw_ticks_abundance_07222022.csv

The above files were used to create the datafiles (df)  using script 1.
5.ectos_pres_abs_df copy.csv # this includes data from manu only
5.ectos_pres_abs_df.csv # this files includes data from manu and iquitos

5.lice_df_abundance_manu.csv # lice manu only 
5.lice_df_abundance_means    # lice abundance manu onlu
5.lice_richness_sp_df_manu.csv  # lice richness manu only 

These files are the ones used for the analyses it  includes the parasite information with traits information, jetz taxonomy  and taxonomy corrected to 2022

7.ectoparasite_df_presence_absence.csv 
7.lice_df_abundance_means.csv
7.lice_df_abundance.csv
7.lice_df_diversity_genus.csv
7.lice_df_diversity_species.csv
7.mites_df_abundance_means.csv
7.mites_df_abundance.csv
7.mites_df_diversity_genus.csv
7.mites_df_diversity_group.csv
7.ticks_df_abundance.csv

# 




