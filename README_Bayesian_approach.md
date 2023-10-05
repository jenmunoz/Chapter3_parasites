# Chapter3_parasites
This is a cummary in my own words ofthe workflow for bayesian in the way I have aapplied to this manuscript 
It also includes refrecences to papers and links to blogs to underestand better bayesian.

# Underestanding parameters in brms models https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup


Bayesian workflow [in progress]
Visualizing the bayesian workflow: https://www.monicaalexander.com/posts/2020-28-02-bayes_viz/
paper: https://academic.oup.com/jrsssa/article/182/2/389/7070184
a more complex but complete workflow https://arxiv.org/pdf/2011.01808.pdf

if run into rpoblems check this out: https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

#1 Explore your data with some simple graphs

#2 Prior predictive checks

#2.1 Weakly informative priors

#3 Run the models

#4 Model posterior predictive checks 
  The idea of posterior predictive checks is to compare our observed data to replicated data from the model. 
  If our model is a good fit, we should be able to use it to generate a dataset that resembles the observed data.

# Use Cross validation for model selections (loo)
The last piece of this workflow is comparing models using leave-one-out cross validation. We are interested in estimating the out-of-sample predictive accuracy at each point  when all we have to fit the model is data that without point . We want to estimate the leave-one-out (LOO) posterior predictive densities and a summary of these across all points, which is called the LOO expected log pointwise predictive density (
elpd. The bigger the numbers, the better we are at predicting the left out point 
By comparing the (elpd_looO)  across models, we can choose the model that has the higher value. By looking at values for the individual points, we can see which observations are particularly hard to predict.

# See detais for K-pareto and other diagnosis below

**Model selection workflow :

1)Analyze the posterior, do they make sense
2Cross validation checking (use loo)
2.1) if random effect present, high pareto K could indicate  misspecification but also flexibility of the model. importance sampling in PSIS-LOO can fail for “random effect” model. 
2.2) Option 1, refi the model using reloo, but takes very long time
2.3) Use k-fold and re-fit the model 10 times, each time leaving out 10 % of teh observations and compare the k-models with loo_compare()
Use posterior predictive checking pp_check  to see if the model can predict the proportion of zeros well



------
**Prior selection workflow ( my learning process)

I am reading this now very useful: http://svmiller.com/blog/2021/02/thinking-about-your-priors-bayesian-analysis/
aslo this : https://github.com/paul-buerkner/brms/issues/131
and this :https://discourse.mc-stan.org/t/default-student-t-priors-in-brms/17197/7


FIRST: Prior selection matter!!!

i) At first I started running the models with the defaul priors, because I thoug that will be reasnable without known more about priors.
Also the defauls look ok, for the fixed efefct predictors the priors were non-informative priors (flat) with I though it was an k decision.
However when I started reading more, it seem that on informative priors are not as hrmfull as I though. Sevral papers and blogs suggest that this is not a good ideam and  Most of the priors that people think are uninformative turn out not to be, and can be harmfull when saple size is not huge and when effec sizes are intrnsically small.

Default priors in brms change with teh model so it tries t o addapt to your model ???
This ia a blog about it https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations, another https://statmodeling.stat.columbia.edu/2018/04/03/justify-my-love/ and some papers https://projecteuclid.org/journals/bayesian-analysis/volume-1/issue-3/Prior-distributions-for-variance-parameters-in-hierarchical-models-comment-on/10.1214/06-BA117A.full

ii) My conclusion after reading all this and without much knowledge of teh system that i ma working on, It seems taht using weakly informative priors is a better option since this 

Weakly informative prior should contain enough information to regularize: the idea is that the prior rules out unreasonable parameter values but is not so strong as to rule out values that might make sense

Here are some notes in weakly informative priors:We characterize a prior distribution as weakly informative if it is proper but is set upso that the information it does provide is intentionally weaker than whatever actual prior knowledge is available. In general any problem has some natural constraints that would allow aweakly-informative model. dor regression models on the logarithmic orlogit scale, with predictors that are binary or scaled to have standard deviation 1, wecan be sure for most applications that e®ect sizes will be less than 10, or certainly lessthan 100.

In our case the data are scaled and centered ( mean zero,a dn sd of 1), therefore we are  we are using a t-student for the predictors with 3 degrees of freedom, a mean of zero and 


ALTHOUGH I am now happy to use the weakly informative priors, interestingly I went to see max farell paper ad recreate their models and for some reason they have flat priors for the fixed effects !!!! ( I am assumnig that was ust an error uploadding the model but it made me freak out). ANyways my main concern now is how to decide the prior for teh intercept cause i have no idea what to do there. 
ALSO FOND NA INTERESTING PAPAER THAT LOOS AT BRMS TOO NAD JUST USE the default parameters (Barrow et al 2019.Deeply+conserved+susceptibility+in+a+multi-host,+m)




#5 levels of priors
#Flat prior (not usually recommended);
#Super-vague but proper prior: normal(0, 1e6) (not usually recommended);
#Weakly informative prior, very weak: normal(0, 10);
#Generic weakly informative prior: normal(0, 1);
#Specific informative prior: normal(0.4, 0.2) or whatever. Sometimes this can be expressed as a scaling followed by a generic prior: theta = 0.4 + 0.2*z; z ~ normal(0, 1);


Notes on priors: ( Gelman 2006)
1)non informative priors =Improper priors, uniform the uniform(0;A) model yields a limiting proper posterior densityas A ! 1, as long as the number of groups J is at least 3. Thus, for a ¯nite butsu±ciently large A, inferences are not sensitive to the choice of A.
Noninformative prior distributions are intended to allow Bayesian inference for param-eters about which not much is known beyond the data included in the analysis at hand.
 Uniform prior distributions are possible (e.g. by setting stan_glm's prior argument to NULL) but, unless the data is very strong, they are not recommended and are not non-informative, giving the same probability mass to implausible values as plausible ones.



2) weakly informative priors:We characterize a prior distribution as weakly informative if it is proper but is set upso that the information it does provide is intentionally weaker than whatever actualprior knowledge is available. Iin general any problem has some natural constraints that would allow aweakly-informative model. dor regression models on the logarithmic orlogit scale, with predictors that are binary or scaled to have standard deviation 1, wecan be sure for most applications that e®ect sizes will be less than 10, or certainly lessthan 100.

Weakly informative prior should contain enough information to regularize: the idea is that the prior rules out unreasonable parameter values but is not so strong as to rule out values that might make sense

Super-weak priors to better diagnose major problemsFor example, normal(0,100) priors added for literally every parameter in the model. The assumption is that everything's on unit scale so these priors will have no effect--unless the model is blowing up from nonidentifiablity. The super-weak prior allows you to see problems without the model actually blowing up. We could even have a setting in Stan that automatically adds this information for every parameter.

We view any noninformative or weakly-informative prior distribution as inherentlyprovisional|after the model has been ¯t, one should look at the posterior distributionand see if it makes sense. If the posterior distribution does not make sense, this impliesthat additional prior knowledge is available that ha

Some links to prior selection 

# Prior selection https://paul-buerkner.github.io/brms/reference/set_prior.html
# very informative blog https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
#other interesting read with  many reference that i explored https://statmodeling.stat.columbia.edu/2018/04/03/justify-my-love/

## other deils on priors that I find less informative 
Recommendations: for standard deviation parameters ( Reffects?) 
a)In fitting hierarchical models, we recommend starting with a noninformative uniformprior density on standard deviation parameters, We expect this will generally workwell unless the number of groups J is low (below 5, say).
b)For a noninformative but proper prior distribution, we recommend approximating theuniform density on ¾® by a uniform on a wide range (for example, U(0; 100) in the SATcoaching example) or a half-normal centered at 0 with standard deviation set to a highvalue such a 100. The latter approach is particularly easy to program as a N(0; 1002)prior distribution 

c)When more prior information is desired, for instance to restrict ¾® away from verylarge values, we recommend working within the half-t family of prior distributions,which are more °exible and have better behavior near 0, compared to the inverse-gamma family. A reasonable starting point is the half-Cauchy family, with scale set toa value that is high but not o® the scale;




**Model selection workflow :

1)Analyze the posterior, do they make sense
2Cross validation checking (use loo)
2.1) if random effect present, high pareto K could indicate  misspecification but also flexibility of the model. importance sampling in PSIS-LOO can fail for “random effect” model. 
2.2) Option 1, refi the model using reloo, but takes very long time
2.3) Use k-fold and re-fit the model 10 times, each time leaving out 10 % of teh observations and compare the k-models with loo_compare()
Use posterior predictive checking pp_check  to see if the model can predict the proportion of zeros well

---
**Model selection notes:


# Cross validation for model selection https://mc-stan.org/loo/articles/online-only/faq.html#modelselection

#Cross-validation is a family of techniques that try to estimate how well a model would predict previously unseen data by using fits of the model to a subset of the data to predict the rest of the data.
#Cross-validation can be used to:
#Asses the predictive performance of a single model
#Asses model misspecification or calibration of the predictive distribution of a single model
#Compare multiple models
#Select a single model from multiple candidates
#Combine the predictions of multiple models
#Even if the goal of the model is not to make predictions, a model which makes bad or badly calibrated predictions is less likely to provide useful insights to a phenomenon studied.

#Using cross-validation for a single model
#Two basic cases why to use cross-validation for one model are:
# We want to know how good predictions the model can make for future or otherwise unseen observations.
 # We want to know if the model describes the observed data well, but we are not going make any predictions for the future.

#How to use cross-validation for model selection?
#First avoid model selection by using the model which includes all predictors and includes all uncertain things. Then optimal thing is to integrate over all the uncertainties. 
#When including many components to a model, it is useful to think more carefully about the prior. For example, if there are many predictors, it is useful to use priors that
#a) state that only some of the effects are big, or b) many effects are big and correlating (it is not possible to have a large number of big independent effects Tosh et al. (2021)).
#If there is explicit utility or loss for observing future predictor values (e.g. medical tests) use decision theory.
#If there is implicit cost for bigger models (e.g. bigger model more difficult to explain or costs of feature measurements are unknown), 
#choose a smaller model which similar predictive performance as the biggest model. If there are only a small number of models, overfitting due to selection process is small. 
#If there are a large number of models, as for example often in variable selection, then the overfitting due to the selection process can be a problem (Piironen and Vehtari, 2017) and more elaborate approaches, such as projection predictive variable selection is recommended.
#If there is application specific utility or loss function, use that to assess practically relevant difference in predictive performance of two models.
#If there is no application specific utility or loss function, use log score, ie elpd. If elpd difference (elpd_diff in loo package) is less than 4, the difference is small (Sivula, Magnusson and Vehtari, 2020)).
#If elpd difference (elpd_diff in loo package) is larger than 4, then compare that difference to standard error of elpd_diff (provided e.g. by loo package) (Sivula, Magnusson and Vehtari, 2020). See also Section How to interpret in Standard error (SE) of elpd difference (elpd_diff)?.
#If there is a large number of models compared, there is possibility of overfitting in model selection.

#Can cross-validation be used for hierarchical / multilevel models?
#The short answer is “Yes”. Hierarchical model is useful, for example, if there are several subjects and for each subject several trials. As discussed in When is cross-validation valid?, it is useful to think of the prediction task or generalizability over different exchangeable entities. We can use different types of cross-validation to choose the focus. This means that also different forms of cross-validation are valid for hierarchical models

# workflow for model selection after reading some of the demo studies
#A case study Roaches cross-validation demo with “random effects” models Cross validation demo https://avehtari.github.io/modelselection/roaches.html

#Model selection note: 
# 1) Analyze the posterior, do they make sense
#2) Cross validation checking (use loo)
#2.1) if random effect present, high pareto K could indicate  misspecification but also flexibility of the model. importance sampling in PSIS-LOO can fail for “random effect” model. 
#2.2) Option 1, refi the model using reloo, but takes very long time
#2.3) Use k-fold and re-fit the model 10 times, each time leaving out 10 % of teh observations and compare the k-models with loo_compare()
#Use posterior predictive checking pp_check  to see if the model can predict the proportion of zeros well


#Steps to check model misspecification 
# i) Calculate loo 
loo1<-loo(ecto_p_brms_bayes_no_int,save_psis = TRUE)
#ii) Plot pareto k diagnostics
plot(loo1)
#iii) leave-one-out cross-validation marginal posterior predictive checks Gabry et al (2018) in combination with bayesplot()
# For a good model, the distribution of LOO-PIT values should be uniform
#LOO-PIT values for our model (thick curve) is compared to many independently generated samples (each the same size as our dataset) 
yrep <- posterior_predict(ecto_p_brms_bayes_no_int)

ppc_loo_pit_overlay(
  y = ectos_birds_dff$ectoparasites_PA,
  yrep = yrep,
  lw = weights(loo1$psis_object)
)

#What to do if I have many high Pareto 𝑘?
#If 𝑘̂<0.5 k  then the corresponding component of elpd_loo is estimated with high accuracy. If 0.5<𝑘̂<0.7the accuracy is lower, but still OK.
#If 𝑘̂>0.7, then importance sampling is not able to provide useful estimate for that component/observation. Pareto-𝑘̂is also useful as a measure of influence of an observation. 
#Highly influential observations have high 𝑘̂values. Very high 𝑘̂values often indicate model misspecification, outliers or mistakes in data processing. 
#If p_loo <𝑝 and the number of parameters 𝑝 is relatively large compared to the number of observations (e.g., 𝑝>𝑁/5), it is likely that the model is so flexible or the population prior so weak that it’s difficult to predict the left out observation (even for the true model). This happens, for example, in the simulated 8 schools (Vehtari, Gelman and Gabry, 2017), random effect models with a few observations per random effect, and Gaussian processes and spatial models with short correlation lengths.

# Solution 1 

#One option is to re-fit the model without these problematic observations, 
#and directly calculate the loo statistic directly for them.

#Exclude those with high values 
loo2<-loo(MODEL, moment_match = TRUE)

#loo2 <- loo(model, k_threshold=0.7)

# i THINK THIS IS MORE OF OUR CASES
#If p_loo < p and the number of parameters p is relatively large compared to the number of observations (e.g., p>N/5),
#it is likely that the model is so flexible or the population prior so weak that it’s difficult to predict the left out observation (even for the true model). This happens, for example, in the simulated 8 schools (in VGG2017), random effect models with a few observations per random effect, and Gaussian processes and spatial models with short correlation lengths.

# Using ELPD for model selection 
# ELPD is good for model comparison as it measures the goodness of the whole predictive distribution
#As quick rule: If elpd difference (elpd_diff in loo package) is less than 4, the difference is small (Sivula, Magnusson and Vehtari, 2020). 
#If elpd difference (elpd_diff in loo package) is larger than 4, then compare that difference to standard error of elpd_diff (provided e.g. by loo package) (Sivula, Magnusson and Vehtari, 2020).
#The value for deciding what is small or large can be based on connection to Pseudo-BMA+-weights (Yao et al., 2018). See also How to interpret in Standard error (SE) of elpd difference (elpd_diff)?. https://mc-stan.org/loo/articles/online-only/faq.html#se_diff
#SE tends to be underestimated especially if the number of observations is small or the models are badly misspecified.
# When the difference (elpd_diff) is larger than 4, the number of observations is larger than 100 and the model is not badly misspecified then normal approximation and SE are quite reliable description of the uncertainty in the difference. 
#Differences smaller than 4 are small and then the models have very similar predictive performance and it doesn’t matter if the normal approximation fails or SE is underestimated (Sivula, Magnusson and Vehtari, 2020).
----



Model selection links to relevan information and papers

# Checking model overfitting and missbehavior and model comparison  ----------------------------------------------

# they promote the use of  WAIC (Widely applicable information criterion), and LOO (leave-one-out cross-validation). These use calculations of the log-likelihood across the entire posterior.
#For predictive models, a common method for picking models is to test the predictive accuracy of a model on new data, or some held-out portion of  the data you already have (but not the portion you used to build the model).
#If your model does well predicting a random datapoint excluded from your model when you build it, it is likely not overfit. However, re-running your model  many times, dropping a different data point every time would take a lot of time making it unpractical for complex models.
#The "loo" package in R employs a algorithm to approximate what the performance of your model would be if you performed full leave-one-out cross-validation.
#If you are interested in how the approximation is calculated,  and how to determine if you can reasonably make this approximation  for you model, you can read the original paper by Aki Vehtari which  explains WAIC and LOO: https://arxiv.org/abs/1507.04544 and check out the example on CRAN: https://cran.r-project.org/web/packages/loo/index.html


#As the proportion of zeros is quite high in the data, it is worthwhile to test also a zero-inflated negative-binomial model, 
#which is a mixture of two models - logistic regression to model the proportion of extra zero counts - negative-binomial model
# notes from Max farell https://github.com/maxfarrell/qcbs_stan_workshop/blob/master/QCBS_stan.Rmd
# some step by step in brms and bayesian https://ourcodingclub.github.io/tutorials/brms/
#Notes check this blog for discussion in STAN BRMS situations: https://discourse.mc-stan.org
# check this blog for discussion in model selection https://discourse.mc-stan.org/t/model-selection-in-brms/30492
# about loo package interpretation https://mc-stan.org/loo/reference/loo-glossary.html
# Using the loo package https://mc-stan.org/loo/articles/loo2-example.html
# some examples here : https://avehtari.github.io/modelselection/roaches.html
# Infereing from P_loo https://discourse.mc-stan.org/t/a-quick-note-what-i-infer-from-p-loo-and-pareto-k-values/3446
# Avoidng model refit with moment match when pareko is high https://mc-stan.org/loo/articles/loo2-moment-matching.html
#how to use the loo package to carry out Pareto smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO) for purposes of model checking and model comparison. https://mc-stan.org/loo/articles/loo2-example.html

# Notes for model selection with BRMS, [Work in progress]
#loo gives us warnings about the Pareto diagnostics, which indicate that for some observations the leave-one-out posteriors are different enough from the full posterior that importance-sampling is not able to correct the difference. We can see more details by printing the loo object.
# definition [elpd_loo] The ELPD is the theoretical expected log pointwise predictive density for a new dataset 
# definition [elpd_loo SE] This standard error is a coarse description of our uncertainty about the predictive performance for unknown future data. 
# In well behaving cases p_loo < N and p_loo < p, where p is the total number of parameters in the model and n ( number of observations)
#elpd_loo differences of more than 4 magnitude are consider lareg enough to make a difference in teh model
# elpd_loo large compared to se_diff indicate changes in the model selection too
# notes more model overfiiting in BRMS 
# Look at the number of observations (226 is shown in the first line), 
# The effective number of parameters, p_loo ()  p_loo < N (number of observations) and p_loo < p (parameters of the model)
# and the number of high Pareto k’s ( All Pareto k’s are good, so there are no highly influential individual observations, This indicates that loo computation is reliable)

#To improve the accuracy of the loo() result above, we could perform leave-one-out cross-validation by explicitly leaving out single observations 
#and refitting the model using MCMC repeatedly. However, the Pareto 𝑘diagnostics indicate that there are 19 observations which are problematic. 
#This would require 19 model refits which may require a lot of computation time.
loo(zip_a_nf_mites_brms_bayes_no_int,reloo = TRUE) # refits the model leaving the observation with hgh pareto out , it is an slow alternative 
#Instead of refitting with MCMC, we can perform a faster moment matching correction to the importance sampling for the problematic observations. 
loo_2<-loo(zip_a_nf_mites_brms_bayes_no_int, moment_match = TRUE)
plot(loo_2)








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










