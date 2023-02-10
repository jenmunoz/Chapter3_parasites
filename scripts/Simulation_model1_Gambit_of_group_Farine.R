

####### 

# Set parameters
N <- c(10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120)
n.obs <- c(5, 10, 15, 20, 25, 30, 35, 40)
params <- expand.grid(N=N,n.obs=n.obs)

# effect present?
effect <- FALSE

# control for group_size as a location?
location <- FALSE

# metric ("DEGREE", "EVCENT", or "BETWEEN")
metric <- "DEGREE"

## repeats
n.reps <- 100 	# repeated networks for each parameter set
n.perm <- 1000	# number of swaps in the pre-network permutation
n.perm2 <- 1000	# number of swaps in the node permutation


# Somewhere to store the results
file_out <- paste("/home/farine/Desktop/Damien/results_effect",effect,"_spatial",location,"_metric",metric,".RData",sep="")


# Run code

library(asnipe)
library(sna)

# empty dataframe for the results
results <- data.frame()

# loop through each set of parameters
for (zzz in 1:nrow(params)) {


## Storage of final measures
sds <- rep(NA, n.reps) # store sd of degree distribution

effect.size <- coef.controlled <- rep(NA, n.reps)
effect.sizes.node.upper <- effect.sizes.node.lower <- rep(NA, n.reps)
coef.trait <- coef.trait.obs.control <- t.trait <- rep(NA, n.reps)	# storage for coefficient values of the model without and with controls for pre-network permutations

coefs.data <- coefs.data.controlled <- coefs.node <- t.data <- matrix(NA, nrow=n.reps, ncol=n.perm) # storage for coefficient values from pre-network permutations and node permutations
effect.sizes.node <- matrix(NA, nrow=n.reps, ncol=n.perm2)
p.data <- p.node <- p.node.control <- p.node.obs.control <- p.t <- rep(NA, n.reps) #P values for: datastream perms, node perms, node perms from double perm, node perms controlling for observations (and location if applicable)
p.model <- p.model.obs.control <- p.model.controlled <- p.model.informed.prenetwork <- rep(NA, n.reps) # P values for model-based tests: uninformed lm, lm controlling for observations, lm controlled using pre-network perm


# loop through to create replicated networks
for (zz in c(1:n.reps)) {

# create a data frame with IDs
ids <- data.frame(ID=1:params$N[zzz])

# Generate some possible groups
groups <- 0.5*params$N[zzz]*params$n.obs[zzz]

# individual trait values
ids$TRAIT <- rnorm(params$N[zzz],0,2)

# make these all positive
ids$TRAIT <- ids$TRAIT-min(ids$TRAIT) + 1


# assign number of observations per individual (sums to N*n.obs)
ids$OBS <- rpois(params$N[zzz],params$n.obs[zzz])
m <- params$N[zzz]*params$n.obs[zzz]
d <- sum(ids$OBS)-m
if (d > 0) {
	rem <- sample(1:params$N[zzz],d,prob=ids$OBS,replace=TRUE)
	rem <- table(rem)
	ids$OBS[as.numeric(names(rem))] <- ids$OBS[as.numeric(names(rem))]-rem
}
if (d < 0) {
	add <- sample(1:params$N[zzz],abs(d),prob=ids$OBS,replace=TRUE)
	add <- table(add)
	ids$OBS[as.numeric(names(add))] <- ids$OBS[as.numeric(names(add))]+add
}
ids$OBS[ids$OBS <= 0] <- 1

# create a vector of group sizes (larger number = more individuals)
group_size <- sample(c(1:10),groups,replace=TRUE)

# Create blank observation matrix
obs <- matrix(0,nrow=groups,ncol=params$N[zzz])


# Allocate individuals to groups, giving smaller trait values a higher probability to be in smaller groups
if (effect) {
	which.ind <- order(ids$TRAIT, decreasing=FALSE) # allocate in increasing order of TRAIT
} else {
	which.ind <- 1:nrow(ids) 	# if random, simply select individuals in random order
}

for (i in which.ind) {

	group_size.tmp <- group_size-rowSums(obs)
	probs <- 1/(1+(group_size.tmp)^2)
	probs[group_size.tmp==0] <- 0
	g <- sample(1:groups,ids$OBS[i],prob=probs)
	obs[g,i] <- 1

}

# remove empty groups
if (any(rowSums(obs) == 0)) {
	group_size <- group_size[-which(rowSums(obs)==0)]
	obs <- obs[-which(rowSums(obs)==0),]
}

# calculate most observed location per ind
max_loc <- function(groups,group_size) {
	grps <- table(group_size[groups==1])
	max <- names(grps)[which(grps==max(grps))]
	if(length(max) > 1){
		max <- sample(max, 1)
	}
	return(max)
}

ids$LOCATION <- apply(obs,2,max_loc,group_size)
ids$OBS <- colSums(obs)


## NOW DO NETWORK ANALYSIS ON THESE DATA

# Calculate network
network <- get_network(obs)

get_network_metric <- function(network, metric) {# Calculate Metric
	if (metric == "DEGREE") {
		METRIC <- rowSums(network)
	} else if (metric == "EVCENT") {
		METRIC <- sna::evcent(network)
	} else if (metric == "BETWEEN") {
		METRIC <- sna::betweenness(1/network, ignore.eval=FALSE)
	} else {
		stop("INCORRECT METRIC SPECIFIED")
	}
	return(METRIC)
}
ids$METRIC <- get_network_metric(network, metric)


### Fix effect relationship (mostly for betweenness which is hard to model)
# do this by sampling individuals in order from largest to smallest network metric
# and select traits proportionally to their number
if (effect & cor(ids$TRAIT,ids$METRIC) < 0.2) {
	metric.order <- order(ids$METRIC,decreasing=FALSE)
	metric.order <- metric.order[order(rnorm(length(metric.order),mean=1:length(metric.order),sd=length(metric.order)^(1/3)))]
	ids <- ids[order(ids$TRAIT,decreasing=FALSE),]
	obs <- obs[,metric.order]
	network <- get_network(obs)
	ids$METRIC <- get_network_metric(network, metric)
	ids$LOCATION <- apply(obs,2,max_loc,group_size)
	ids$OBS <- colSums(obs)
}


# Calculate effects
model <- lm(METRIC~TRAIT,data=ids)
coef.trait[zz] <- coefficients(model)[2]
p.model[zz] <- summary(model)$coefficients[2,4]
t.trait[zz] <- summary(model)$coefficients[2,3]

# Calculate effects with informed controls
if (location == FALSE) {
	model <- lm(METRIC~TRAIT+OBS,data=ids)
	coef.trait.obs.control[zz] <- coefficients(model)[2]
	p.model.obs.control[zz] <- summary(model)$coefficients[2,4]
} else {
	model <- lm(METRIC~TRAIT+OBS+LOCATION,data=ids)
	coef.trait.obs.control[zz] <- coefficients(model)[2]
	p.model.obs.control[zz] <- summary(model)$coefficients[2,4]
}

## Data permutations
# Create random networks
if (location) {
	networks_rand <- network_permutation(obs, association_matrix=network, permutations=n.perm, locations=group_size, within_location=TRUE)
} else {
	networks_rand <- network_permutation(obs, association_matrix=network, permutations=n.perm)
}

# Calculate degrees for each network
deg_rand <- apply(networks_rand,1,function(x) { get_network_metric(x, metric)})

# function to get coefficients for each randomisation
get_res <- function(x) {
	model <- lm(x~TRAIT,data=ids)
	coef <- coefficients(model)[2]
	return(coef)
}
# function to get coefficients for each randomisation
get_res2 <- function(x) {
	model <- lm(x~TRAIT,data=ids)
	coef <- summary(model)$coefficients[2,3]
	return(coef)
}

# function to get coefficients for each randomisation
get_res3 <- function(x) {
if (location == FALSE) {
	model <- lm(x~TRAIT+OBS,data=ids)
} else {
	model <- lm(x~TRAIT+OBS+LOCATION,data=ids)
}
	coef <- coefficients(model)[2]
	return(coef)
}

# coefficients based on original degree
coefs.data[zz,] <- apply(deg_rand,2,function(x) { get_res(x) })
t.data[zz,] <- apply(deg_rand,2,function(x) { get_res2(x) })
coefs.data.controlled[zz,] <- apply(deg_rand,2,function(x) { get_res3(x) })

# get controlled effect size
effect.size[zz] <- coef.trait[zz] - median(coefs.data[zz,])

# controlled degree (after pre-network permutation)
ids$METRIC.cont <- ids$METRIC - apply(deg_rand,1,median)

# Calculate effect controlling using pre-network permutation
model <- lm(METRIC.cont~TRAIT,data=ids)
coef.controlled[zz] <- coefficients(model)[2]
p.model.controlled[zz] <- summary(model)$coefficients[2,4]


## Node permutations (store these in temporary arrays)
coefs.tmp <- rep(NA,n.perm2)	# straight-up node permutation
coefs.tmp.control <- rep(NA,n.perm2)	# node permutation controlling using pre-network permutations
coefs.tmp.obs.control <- rep(NA,n.perm2)	# node permutation controlling for number of observations
for (i in 1:n.perm2) {
	ids$TRAIT.tmp <- sample(ids$TRAIT)
	model.tmp <- lm(METRIC~TRAIT.tmp,data=ids)
	coefs.tmp[i] <- coefficients(model.tmp)[2]
	model.tmp <- lm(METRIC.cont~TRAIT.tmp,data=ids)
	coefs.tmp.control[i] <- coefficients(model.tmp)[2]
	if (location == FALSE) {
		model.tmp <- lm(METRIC~TRAIT.tmp+OBS,data=ids)
		coefs.tmp.obs.control[i] <- coefficients(model.tmp)[2]
	} else {
		model.tmp <- lm(METRIC~TRAIT.tmp+OBS+LOCATION,data=ids)
		coefs.tmp.obs.control[i] <- coefficients(model.tmp)[2]
	}
}

# P values
p.data[zz] <- sum(coef.trait[zz]<=coefs.data[zz,])/n.perm	# P value from the pre-network permutation
p.node[zz] <- sum(coef.trait[zz]<=coefs.tmp)/n.perm2		# P value from the simple node permutation
p.node.control[zz] <- sum(coef.controlled[zz] <= coefs.tmp.control)/n.perm2	# P value from the double permutation test
p.node.obs.control[zz] <- sum(coef.trait.obs.control[zz] <= coefs.tmp.obs.control)/n.perm2	# P value from the node permutation controlling for #obs and location (if applicable)
p.t[zz] <- sum(t.trait[zz] <= t.data[zz,])/n.perm	# P value based on the pre-network permutations using t values
p.model.informed.prenetwork[zz] <- sum(coef.trait.obs.control[zz]<=coefs.data.controlled[zz,])/n.perm	# prenetwork permutation with lm controlling for OBS

# sd
sds[zz] <- sd(ids$DEGREE)

}

r <- data.frame(N=params$N[zzz], n.obs=params$n.obs[zzz], coef.trait=coef.trait, coef.trait.obs.control=coef.trait.obs.control, coef.controlled=coef.controlled, p.model=p.model, p.model.obs.control=p.model.obs.control, p.model.controlled=p.model.controlled, p.node.control=p.node.control, p.node.obs.control=p.node.obs.control, effect.size=effect.size, p.data=p.data, p.node=p.node, p.t, sds=sds, p.data.controlled=p.model.informed.prenetwork)
results <- rbind(results, r) 

}

save(results, params, N, n.obs, file=file_out)


