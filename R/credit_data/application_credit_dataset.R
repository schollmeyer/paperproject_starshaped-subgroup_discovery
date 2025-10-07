# set working directory
setwd("C:/Git/paperproject_starshaped-subgroup_discovery")
# needed libraries
library(foreign)
library(rsubgroup)
library(ggplot2)
# load the data set credit ( https://archive.ics.uci.edu/dataset/144/statlog+german+credit+data)
# Thisis already the curated dataset, compare TODO
# The variable sex (column 9) will be later deleted because it is actually not
# the variable sex but instead a combination of sex and marital status
data(credit.data)
dat <- credit.data
dim(dat)
# [1] 1000   21

#generate formal context of covariates (sex (column 9) and credit status (column 21) are excluded)
# TODO: comment

whole_context <- oofos:::get_auto_conceptual_scaling(dat[,-c(9,21)])
# For the classical subgroup discoery we make a sample split approach because
# optimizing under H0 is computationally very expensive
set.seed(1234567)
indexs <- sample((1:1000),size =500)
training_context <- whole_context[indexs,]
test_context <- whole_context[-indexs,]

# generate objective vector for optimization. This corresponds
# (up to a multiplicative constant) to the
# Piatetsky-Shapiro quality function for the target value credit-state==good
objective <- oofos:::compute_objective(dat[indexs,],"class","good")
whole_objective <- oofos:::compute_objective(dat,"class","good")

table(objective)

# objective
# -0.00628930817610063  0.00293255131964809
# 159                  341
#
# 341 good states and 159 bad states


############################################################################
# Classical Subgroup Discovery on a subset of size n=500
############################################################################


# generate gurobi model for classical subgroup discovery
model <- oofos:::optimize_on_context_extents(context=training_context,objective=objective,binary_variables="all" )
# For later computations it is more effective to treat the constraint matrix
# not as a simple triplet matrix
model$A <- as.matrix(model$A)
# add additional constraints that may help to speed up the optimization
model <- oofos:::add_attr_antiimplications(model)
#compute a preliminary result with timelimit of 10 minits
# result_1 <- gurobi::gurobi(model,list(timelimit=60*10))



# save result
# saveRDS(result_1,"results_credit_data/result_1.RDS")

result_1 <- readRDS("results_credit_data/result_1.RDS")

result_1$objval
# [1] 0.3223593


result_1$runtime
# [1] 601.01

# use semioptimal solution to tighten the search space via optimistic estimates
# (in the style of M. Boley and H. Grosskreutz. Non-redundant subgroup discovery us-
# ing a closure system.)
model_2 <- oofos:::add_sos_constraints(model,result_1$objval)
# add semioptimal solution as a start vector
model_2$start <- round(result_1$x,2)
# result_2 <- gurobi::gurobi(model_2,list(PoolSolutions=20,PoolSearchMode=2))


#saveRDS(result_2,"results_credit_data/result_2.RDS")

result_2 <- readRDS("results_credit_data/result_2.RDS")

#overall runtie in hours
(result_1$runtime  + result_2$runtime)/3600
# [1] 9.11968

result_2$objval

# [1] 0.3572364

#Piatetsky-Shapiro quality value
oofos:::quality(model_2,result_2)$piatetsky_shapiro
#[1] 38.738

# For comparison: subgroup discovery with R package rsubgroup:
task <- CreateSDTask(dat[indexs,], as.target("class", "good"), SDTaskConfig(attributes=colnames(dat),qf="ps",method="sdmap",k=1,maxlen=100,discretize=TRUE,nbins=3,nodefaults=FALSE))
result_rsubgroup <- DiscoverSubgroupsByTask(task)
(result_rsubgroup[[1]])@quality
#[1] 38.328
# value is a little bit smaller than for exhaustive subgroup discovery via MILP




# Subgroup with largest value of the Piatetsky-Shapiro quality function
extent <- round(result_2$x[(1:500)],2)
intent <- round(result_2$x[-(1:500)],2)
colnames(training_context)[which(intent==1)]
# ...
# manually read out extreme attributes:

# checking_status: no checking
# credit_amount in [522,11054]
# duration in [4,60]  # *
# age in [20,65] # min*
# existing_credits in [1,3] # min*
# installment_commitment in [1,4] #*
# residence_since in [1,4] #*
# num_dependents in[1,2] #*



# Sample Splitting analysis of statistical significance
 
pre_extent <- rep(0,1000)
pre_extent[indexs[which(result_2$x>=0.5)]] <- 1
intent <- oofos:::compute_psi(pre_extent,whole_context)
extent <- oofos:::compute_phi(intent,whole_context)

n_1 <- sum((extent*(whole_objective>0))[-indexs] )
N_1 <- sum(extent[-indexs])

n_2 <- sum(((1-extent)*(whole_objective>0))[-indexs] )
N_2 <- sum((1-extent)[-indexs])

X <- rbind(c(n_1, N_1-n_1),c(n_2,N_2-n_2))
fisher.test(X,alternative="greater")


# For comparison: rsubgroup:
# TODO


############################################################################
# Sample splitting test
############################################################################
# TODO
extent <- round(result_2$x[(1:500)],2)
intent <- round(result_2$x[-(1:500)],2)
extent_test <- oofos:::compute_phi(intent,test_context)
test_objective <- oofos:::compute_objective(dat[-indexs,],"class","good")
observed_test_statistic <- sum(test_objective*extent_test)
observed_test_statistic

n_rep <-1000000
h0 <- rep(0, n_rep)
set.seed(1234567)
for(k in (1:n_rep)){
  h0[k] <- sum(sample(test_objective)*extent_test)
}
plot(ecdf(h0))
mean(h0 >= observed_test_statistic)

#0

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#  END OF CURATED PART

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

############################################################################
# Classical Subgroup Discovery on a subset of size n=500
############################################################################
# Laufzeit ca 40 min
#unter H0: Laufzeit ca. > 83864s
# unter H0 mit gleichem Trick: 13040.84 s
# mit bestem SOS constraint: 3274s

result_2$objval
#[1] 0.1326288
set.seed(1234567)
indexs <- sample((1:1000),size=500)
sampled_context <- context[indexs,]

sampled_objective <- oofos:::compute_objective(dat[indexs,],"class","good")

model_sampled <- oofos:::optimize_on_context_extents(context=sampled_context,objective=sampled_objective,binary_variables="all" )
model_2_sampled <- oofos:::add_attr_antiimplications(model_sampled)
model_2_sampled <- oofos:::add_sos_constraints(model_2_sampled,0.278)
result_sampled <- gurobi::gurobi(model_2_sampled,list(timelimit=60*40))


############################################################################
# Starshaped Subgroup Discovery with angle based stylized betweenness
############################################################################

indexs_numerical_variables <- c(2,5,8,11,13,16,18)

absb <- get_absb(as.matrix(dat[,indexs_numerical_variables]))

saveRDS(absb,"results_credit_data/absb.RDS")

## obsb

#obsb <- get_obsb(training_context)
obsb <- readRDS("results_credit_data_final_n_500/obsb.RDS")
ssd_vc_5 <- oofos:::discover_starshaped_subgroups(stylized_betweenness=obsb,objective=objective,complexity_measure=oofos:::compute_width,complexity_control=5,params=list(outputflag=1))
saveRDS(ssd_vc_5,"results_credit_data_ssd_vc_5")
test_vc_5 <- oofos::compute_starshaped_distr_test(ssd_vc_5)



#sample splitting
objective2 <- oofos:::compute_objective(dat[indexs,],"class","good",weights=c(rep(1,250),rep(0,250)))
gc();ssd_vc_70_sample_splitting <- oofos:::discover_starshaped_subgroups(stylized_betweenness=obsb,objective=objective2,complexity_measure=oofos:::compute_width,complexity_control=70,params=list(outputflag=0))


#TODO
mask <- c(rep(0,250),rep(1,250))
star <- ssd_vc_70_sample_splitting$star
pre_extent[indexs[which(result_2$x>=0.5)]] <- 1
intent <- oofos:::compute_psi(pre_extent,whole_context)
extent <- oofos:::compute_phi(intent,whole_context)

n_1 <- sum(mask * star *(objective>0))
N_1 <- sum(mask * star )

n_2 <- sum(mask * (1-star) *(objective>0))
N_2 <- sum(mask *(1- star) )

X <- rbind(c(n_1, N_1-n_1),c(n_2,N_2-n_2))
fisher.test(X,alternative="greater")