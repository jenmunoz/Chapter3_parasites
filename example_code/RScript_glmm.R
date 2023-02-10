setwd("~/Dropbox/aaDocuments/Research/3_Writing_up/YUS_malaria")

library(dplyr)
library(nlme)
library(lme4)
library(phyr)
library(ape)
library(ggplot2)
theme_set(theme_classic())

mydata <- read.csv("example_code/full_png.csv")

str(mydata)

# fit a generalized linear mixed model (good suggestion by reviewer)
z <- glmer(status ~ elevation + (1|lineage), family = binomial(link = "logit"), data = mydata)
           
summary(z)


newdata <- data.frame(elevation = mydata$elevation,
  lineage = mydata$lineage)

newdata$Predicted_Response <- predict(z, newdata = newdata, type = "response")


# model  predicts overall patterns pretty well

summary_predictions <- newdata %>%   # Specify data frame
  group_by(elevation) %>%      # Specify group indicator
  dplyr::summarise(prevalence <- mean(Predicted_Response))

colnames(summary_predictions) <- c("elevation", "prevalence")

ggplot(summary_predictions, aes(x = elevation, y = prevalence))+
  geom_point()


# fit a glmm that incorporates phylogenetic information
# use phyr package; function is pglmm

# read in phylogeny
phy <- read.nexus("example_code/mcc.nex")

phy$tip.label

splist <- unique(mydata$lineage)
splist %in%phy$tip.label
# make a couple names match

mydata$lineage <- as.character (mydata$lineage)
mydata$lineage[mydata$lineage == "Pachycephala_rufinucha"] <- "Aleadryas_rufinucha"
mydata$lineage[mydata$lineage == "Amblyornis_germanus"] <- "Amblyornis_macgregoriae"
mydata$lineage[mydata$lineage == "Peneothello_bimaculatus"] <- "Peneothello_bimaculata"
mydata$lineage[mydata$lineage == "Pitohui_dicrous"] <- "Pitohui_dichrous"
mydata$lineage[mydata$lineage == "Ptilorrhoa_castanonotus"] <- "Ptilorrhoa_castanonota"
mydata$lineage[mydata$lineage == "Oedistoma_iliolophus"] <- "Toxorhamphus_iliolophus"
mydata$lineage[mydata$lineage == "Tregallasia_leucops"] <- "Tregellasia_leucops"

splist <- unique(mydata$lineage)
splist %in%phy$tip.label
# they match now

# see this website for information on using pglmm
# https://daijiang.github.io/phyr/articles/phyr_example_empirical.html

mod <- phyr::pglmm (status ~ elevation + (1 | lineage__), family = "binomial", cov_ranef = list(lineage = phy), data = mydata)
# note! this takes ~ 10 minutes to run
summary(mod)

# how to know if random effects are strong enough to take seriously? 
# Well one way to get an idea is to run a likelihood ratio test on the random effect. 
# This can be achieved by using the pglmm_profile_LRT() function, at least for binomial models.
# (this code from link above)

LRTest <- sapply(1:6, FUN = function(x) phyr::pglmm_profile_LRT(mod, re.number = x))
colnames(LRTest) <- names(mod$ss)
t(LRTest)

newdata <- data.frame(elevation = mydata$elevation,
                      lineage = mydata$lineage)

predictions<- predict(mod, newdata = newdata, type = "response")

mydata2 <- cbind(mydata, predictions)

summary <- mydata2 %>%   # Specify data frame
  group_by(elevation) %>%      # Specify group indicator
  dplyr::summarise(prevalence <- mean(status), sd <- sd(status), n <- n())
colnames(summary) <- c("elevation", "prevalence", "sd", "n")

summary$se <- summary$sd/sqrt(summary$n)
head(summary)

# make a plot of model predictions (that also shows data)
ggplot(data = mydata2, aes(x = elevation, y = Y_hat))+
  geom_smooth()+
  geom_point(data = summary, aes(x = elevation, y = prevalence))+
  geom_segment(data = summary, aes(x = elevation, y = prevalence, xend = elevation, yend = prevalence + se))+
  geom_segment(data = summary, aes(x = elevation, y = prevalence, xend = elevation, yend = prevalence - se))+
  scale_y_continuous("infection prevalence", limits = c(0,1)) +
  geom_hline(yintercept = 0.4535, linetype = "dashed")+ # show average prevalence for whole dataset?
  scale_x_continuous("elevation (m)")
# figure out how to add site means from raw data to this
ggsave("Figure1.pdf", height = 4, width = 5)

# can also make a plot that gives trendlines for different species
species20 <- newdata %>%   # Specify data frame
  group_by(lineage) %>%      # Specify group indicator
  dplyr::summarise(number <- n())
colnames(species20) <- c("species", "n")
species20 <- species20[species20$n >19,]

# show trends with elevation for all 35 species with 20 or more samples in dataset
ggplot(data = mydata2[mydata2$lineage%in%species20$species,], aes(x = elevation, y = Y_hat, color = lineage))+
  # geom_point() +
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_viridis_d(option = "viridis")+
  scale_y_continuous("infection prevalence", limits = c(0,1)) +
  scale_x_continuous("elevation (m)")+
  theme(legend.position="none")
ggsave("Figure2.pdf", height = 4, width = 5)




trial <- mydata %>%   # Specify data frame
  group_by(lineage) %>%      # Specify group indicator
  dplyr::summarise(infected <- sum(status), n <- n())
colnames(trial) <- c("lineage", "infected", "n")

trial$prevalence <- trial$infected / trial$n

# make a figure showing jewel babblers and berrypeckers

jbs <- c("Ptilorrhoa_caerulescens", "Ptilorrhoa_castanonota", "Ptilorrhoa_leucosticta")
mydata_jb <- mydata[mydata$lineage %in% jbs,]

ggplot(data = mydata_jb, aes(x = elevation, y = status, color = lineage))+
  geom_jitter(height = 0.05, width = 10)+
  scale_y_continuous("infection", limits = c(-0.1,1.1)) +
  scale_x_continuous("elevation (m)", limits = c(0, 3000))+
  scale_color_viridis_d(begin = 0.95, end= 0.1)+
  theme(legend.position="none")
ggsave("jewelbabblers.pdf", height = 4, width = 5)


bps <- c("Melanocharis_longicauda", "Melanocharis_nigra", "Melanocharis_versteri")
mydata_bps <- mydata[mydata$lineage %in% bps,]
mydata_bps$lineage <- as.factor(mydata_bps$lineage)
mydata_bps$lineage <- relevel(mydata_bps$lineage, "Melanocharis_nigra")

ggplot(data = mydata_bps, aes(x = elevation, y = status, color = lineage))+
  geom_jitter(height = 0.05, width = 10)+
  scale_y_continuous("infection", limits = c(-0.1,1.1)) +
  scale_x_continuous("elevation (m)", limits = c(0, 3000))+
  scale_color_viridis_d(begin = 0.85, end= 0)+
  theme(legend.position="none")
ggsave("berrypeckers.pdf", height = 4, width = 5)



wh <- c("Pachycephala_hyperythra", "Pachycephala_soror", "Pachycephala_schlegelii")
mydata_wh <- mydata[mydata$lineage %in% wh,]
mydata_wh$lineage <- as.factor(mydata_wh$lineage)
mydata_wh$lineage <- factor(mydata_wh$lineage, c("Pachycephala_hyperythra", "Pachycephala_soror", "Pachycephala_schlegelii"))

ggplot(data = mydata_wh, aes(x = elevation, y = status, color = lineage))+
  geom_jitter(height = 0.05, width = 10)+
  scale_y_continuous("infection", limits = c(-0.1,1.1)) +
  scale_x_continuous("elevation (m)", limits = c(0, 3000))+
  scale_color_viridis_d(begin = 0.85, end= 0)+
  theme(legend.position="none")
ggsave("whistlers.pdf", height = 4, width = 5)

############make the phylocomp map
tennex = read.nexus("example_code/mcc.nex")
tenprev = read.csv("example_code/10plus.csv")
tenmid = read.csv("example_code/10mid.csv")
tenlist = setNames(tenprev$prevalence, tenprev$Jetz_name)
tenlist.v2 = setNames (tenmid$midpoint, tenmid$..Jetz_name)

object = contMap(tennex, tenlist, direction = "leftwards", plot=FALSE)
n=length(object$cols)
object$cols[1:n]=viridis(n)

tenprev



png("phylotree.png", width = 2500, height = 3100, res = 300, units = "px")
plot(object, mar=c(5.1,1,1,1), ylim=c(1-0.09*(Ntip(object$tree)-1), Ntip(object$tree)), legend=FALSE)
add.color.bar(40, object$cols, title = "Infection Prevalence", lims = object$lims, digits = 3, prompt=FALSE,x=0,y=1-0.08*(Ntip(object$tree)-1), lwd=4,fsize=1,subtitle="")
dev.off()
title(xlab = "Time from root (Ma)")


object = contMap(tennex, tenlist.v2, direction = "leftwards", plot=FALSE)
n=length(object$cols)
object$cols[1:n]=viridis(n)

png("phylotree_midpoint.png", width = 2500, height = 3100, res = 300, units = "px")
plot(object, mar=c(5.1,1,1,1), ylim=c(1-0.09*(Ntip(object$tree)-1), Ntip(object$tree)), legend=FALSE)
add.color.bar(40, object$cols, title = "Midpoint elevation (m)", lims = object$lims, digits = 3, prompt=FALSE,x=0,y=1-0.08*(Ntip(object$tree)-1), lwd=4,fsize=1,subtitle="")
dev.off()

###pagels lambda
names = match(prev.mid$Jetz_name, tennex$tip.label)
z = nlme::gls(prevalence~midpoint, data = prev.mid2, correlation = corPagel(0.5, tennex, fixed=FALSE))
summary(z)
