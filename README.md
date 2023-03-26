# Chapter3_parasites
This is the code for the chapter 3 of my PhD


## Sociality ##
Sociality status is based on flcok observations from Manu only, some particular species have changes in thh status

Lophotriccus pileatus::This species has been obserevd in flocks but those floks were incomplete and therefore not included in the analyses, I changed 0 to 1 in sociality manually

##File names and description ##

Raw data files:

[1]The raw files of bird diversity form Manu diversity 

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
5.lice_df_abundance_means    # lice abundance manu only
5.lice_richness_sp_df_manu.csv  # lice richness manu only 

These files are the ones used for the analyses it  includes the parasite information with traits information, jetz taxonomy  and taxonomy corrected to 2022: 

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

# From those files we generated the datasets for data_analyses, selecting only relevant collumns with the following characeterustics: Script 2 

Jetz taxonomy
only Manu data, 
excluded iquitos, 
included sociality_binomial
included and selected some bird traits
included elevation of the sample
included elevation midpoint ofthe species 


# files 
7.dff_all_ectos_prevalence_abundance_individual.csv # includes  Abundance and prevalence 

###
We the included the individual elevation for the samples using teh metadata from samples 
#0.ectoparasite_samples_elevation_metadata.csv

And for some ( few smples nputed data manually)

#P121006	Rupicola_peruvianus 1362 The correct name of this is SP121006 AND THE ELEVATION 1362
#"P120889" is Leptopogon superciliaris CORRRECT NAME SP120889  AND THE ELEVATION IS 1395
#P120085 Formicarius rufipectus   correct name  SP120085 elev 1362
#P115912 Leptopogon amaurocephalus elevation 381 


### Created collumn with extrapolated date 
This column take the NA for elevation and use the information about date and station to see which other species were captured the same data in the same station and used that elevation 


included degree_general_network (at the species level)
included refined degree at the individual level





Notes on analyses with PGLMM phyr

Goodness of Fit ( package rr2)
Another question about our model is simply how well it fits the data. A metric appealing in its simplicity is the classic R2 metric, which purports to estimate the proportion of the variance in the response explained by the model predictors. Once again, though, when using non-Gaussian errors this can become a bit tricky. An additional complication is created by model with random effects. Given that random effects are very flexible model components (for example, nothing stops you from fitting a random effect for each observation in your dataset), a straight-up calculation of variance explains isn’t meaningful. That said, methods that can produce a useful R2 metric in the complex situation have been developed. The package rr2 is able to calculate several flavors of R2, and supports phyr’s pglmm model object. Let’s try it!

rr2::R2(mod)

rr2::R2(mod_bayes)


# notes on analyses with brms 

#priors 

For the zero-inflated Poisson model in brms, the default prior distributions are as follows:

Fixed effects:
The default prior for fixed effects is a Student's t-distribution with 3 degrees of freedom, a mean of 0, and a scale of 10. This prior is relatively vague and is intended to be a default that works well in many situations.
Random effects:
The default prior for the intercepts of the random effects is a Student's t-distribution with 3 degrees of freedom, a mean of 0, and a scale of 10.
The default prior for the standard deviation of the random effects is a half Student's t-distribution with 3 degrees of freedom, a mean of 0, and a scale of 1.
The default prior for the correlation matrix of the random effects is a LKJ prior with a shape parameter of 1, which is a weakly informative prior that encourages shrinkage towards zero.
These default priors are intended to be relatively uninformative and to allow the data to drive the inference. However, it is important to remember that these priors may not always be appropriate for your specific analysis and that it is often beneficial to use informative priors based on prior knowledge or previous research.


The default prior distributions for zero-inflated negative binomial in brms are:

Fixed effects:
The default prior for fixed effects is a Student's t-distribution with 3 degrees of freedom, a mean of 0, and a scale of 10.
Random effects:
The default prior for the intercepts of the random effects is a Student's t-distribution with 3 degrees of freedom, a mean of 0, and a scale of 10.
The default prior for the standard deviation of the random effects is a half Student's t-distribution with 3 degrees of freedom, a mean of 0, and a scale of 1.
The default prior for the correlation matrix of the random effects is a LKJ prior with a shape parameter of 1, which is a weakly informative prior that encourages shrinkage towards zero.
In addition to these default priors, brms allows you to specify other prior distributions for the model parameters, including more informative priors based on prior knowledge or previous research. It is often a good idea to check the appropriateness of the default priors for your specific analysis and data, and to adjust them if necessary.



example of non informative priors for my model

non_informative_prior <- prior(normal(0, 10), class = "b") + 
                        prior(student_t(3, 0, 10), class = "Intercept") + 
                        prior(normal(0, 10), class = "sd") + 
                        prior(lkj(1), class = "cor")

# Example in how to interpret the model outputs of brms 

# This is an output 

#Family: zero_inflated_negbinomial 
  Links: mu = log; shape = identity; zi = identity 
Formula: total_lice ~ scale(degree, center = TRUE) + scale(elevation) + (1 | gr(species_jetz, cov = phy_cov)) + (1 | Powder.lvl) + (1 | species) 
   Data: dff_ectos_network_individual_metrics (Number of observations: 535) 
  Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 2;
         total post-warmup draws = 6000

Group-Level Effects: 
~Powder.lvl (Number of levels: 6) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     0.92      0.75     0.08     2.93 1.00     3625     4045

~species (Number of levels: 82) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     0.71      0.38     0.05     1.45 1.01     1140     2547

~species_jetz (Number of levels: 82) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     2.33      0.70     0.91     3.68 1.01     1534     2467

Population-Level Effects: 
                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  -1.21      1.23    -3.80     1.04 1.00     4778     5082
scaledegreecenterEQTRUE     0.17      0.19    -0.21     0.55 1.00     5543     5583
scaleelevation             -0.35      0.25    -0.85     0.14 1.00     5459     5581

Family Specific Parameters: 
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.40      0.08     0.29     0.58 1.00     5031     5040
zi        0.08      0.06     0.00     0.23 1.00     5116     5320

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).

HERE SOME INTERPRETATION 

This is the summary output of a zero-inflated negative binomial model fitted using the brms package in R. Here are some interpretations of the different sections:

Formula: This section shows the formula used for the model, including the response variable (total_lice) and the predictor variables (scale(degree, center = TRUE), scale(elevation), and the grouping variables (1 | gr(species_jetz, cov = phy_cov)), (1 | Powder.lvl), and (1 | species)).
Data: This section shows the number of observations in the data set (Number of observations: 535).
Draws: This section shows the number of chains (4), the number of iterations (iter = 6000), the number of warmup iterations (warmup = 3000), and the thinning rate (thin = 2).

Group-Level Effects: This section shows the estimated standard deviation (sd) of the random intercept for each grouping variable. For example, the estimated standard deviation of the random intercept for Powder.lvl is 0.92. The Bulk_ESS and Tail_ESS values are measures of the effective sample size for each parameter.
Population-Level Effects: This section shows the estimated coefficients (Estimate) and standard errors (Est.Error) for each predictor variable. For example, the estimated coefficient for scale(degree, center = TRUE) is 0.17, indicating that, on average, the expected value of total_lice increases by 0.17 units for every one unit increase in degree after centering. The Bulk_ESS and Tail_ESS values are measures of the effective sample size for each parameter.
Family Specific Parameters: This section shows the estimated parameters for the distribution of the response variable (total_lice). For a zero-inflated negative binomial model, there are two parameters: shape and zi. The shape parameter is the dispersion parameter and the zi parameter is the probability of excess zeros. The Bulk_ESS and Tail_ESS values are measures of the effective sample size for each parameter.

# Other example 

zero_inflated_negbinomial 
  Links: mu = log; shape = identity; zi = identity 
Formula: total_lice ~ sociality + scale(elevation) + (1 | gr(species_jetz, cov = phy_cov)) + (1 | Powder.lvl) + (1 | species) 
   Data: ectos_birds_dff (Number of observations: 782) 
  Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 2;
         total post-warmup draws = 4000

Group-Level Effects: 
~Powder.lvl (Number of levels: 7) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     0.85      0.70     0.07     2.70 1.00     2075     2662

~species (Number of levels: 131) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     0.88      0.33     0.17     1.48 1.01      658     1109

~species_jetz (Number of levels: 131) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     2.00      0.63     0.76     3.21 1.01      675     1289

Population-Level Effects: 
               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         -1.12      1.01    -3.22     0.73 1.00     3161     3190
sociality1         0.16      0.45    -0.72     1.05 1.00     3427     3418
scaleelevation    -0.20      0.21    -0.59     0.21 1.00     3005     3450

Family Specific Parameters: 
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.39      0.06     0.29     0.53 1.00     3167     3373
zi        0.07      0.06     0.00     0.20 1.00     3595     3397

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).

# Interpretation 
The population-level effects estimates represent the average effect of each predictor on the response variable, after accounting for the other predictors in the model.

In this model, the population-level effects estimates are:

Intercept: -1.12
sociality1: 0.16
scaleelevation: -0.20
Intercept: The intercept represents the expected value of the response variable when all predictors are equal to zero. In this model, the intercept estimate is -1.12, which means that the expected value of the response variable (total_lice) is -1.12 when sociality and elevation are both zero.

sociality1: The sociality predictor is a categorical variable that takes on two values (0 and 1). The estimate for sociality1 is 0.16, which means that, on average, the expected value of the response variable increases by 0.16 when sociality is equal to 1 (compared to when sociality is equal to 0), while holding elevation constant.

scaleelevation: The scaleelevation predictor is a continuous variable. The estimate for scaleelevation is -0.20, which means that, on average, the expected value of the response variable decreases by 0.20 for every one-unit increase in elevation, while holding sociality constant.

It's important to note that the interpretation of these estimates assumes that all other variables in the model are held constant. The effect of a predictor may change depending on the values of the other predictors in the model.

### Additionally 
Yes, of course! In this model, the population-level effect estimate for the "sociality1" predictor variable is 0.16. This means that, on average, social species tend to have 0.16 more total lice than non-social species.

Since sociality is coded as a binary variable (0 or 1), this effect is in comparison to the reference group, which is non-social species (coded as 0). Therefore, this effect estimate of 0.16 indicates that, on average, social species have 0.16 more total lice than non-social species after accounting for other variables in the model.




