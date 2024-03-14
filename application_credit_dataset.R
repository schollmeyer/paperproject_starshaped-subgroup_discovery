# set working directory
setwd("C:/Git/paperproject_starshaped-subgroup_discovery")
library(foreign)
library(rsubgroup)#
dat_old <- read.arff("datasets/dataset_31_credit-g.arff")
data(credit.data)
dat <- credit.data
dim(dat)
# 

#generate formal context


whole_context <- oofos:::get_auto_conceptual_scaling(dat[,-c(9,21)])
#

indexs <- sample((1:1000),size =500)
training_context <- oofos:::get_auto_conceptual_scaling(dat[indexs,-c(9,21)])#whole_context[indexs,]
test_context <- whole_context[-indexs,]
# remove duplicated columns
#context <- t(unique(t(context)))
dim(whole_context)

# generate objective vector for optimization. This corresponds to the
# Piatetsky-Shapiro quality function for the target value credit-state==good
objective <- oofos:::compute_objective(dat[indexs,],"class","good")#class","good")

table(objective)
#
# 700 good states and 300 bad states


############################################################################
# Classical Subgroup Discovery
############################################################################


# generate gurobi model for classical subgroup discovery
model <- oofos:::optimize_on_context_extents(context=training_context,objective=objective,binary_variables="all" )

model$A <- as.matrix(model$A)# optimize model under time constraint of 20 minutes for a first semioptimal solution
result_1 <- gurobi::gurobi(model,list(timelimit=60*10))

result_1$objval
# [1] 0.257619


result_1$runtime
# [1] 2400.381

# save result
saveRDS(result_1,"results_credit_data/result_1.RDS")




# add additional contsraints that are always valid for inetrordinally scaled variables x (implications of the form <= x \rightarrow not >= x+1 etc. 
model_2 <- oofos:::add_attr_antiimplications(model)

# use semioptimal solution to tighten the search space via optimistic estimates (in the style of ***)
model_2 <- oofos:::add_sos_constraints(model_2,result_1$objval)
# add semioptimal solution as a start vector
model_2$start <- round(result_1$x,2)
result_2 <- gurobi::gurobi(model_2)


result_2$objval

# [1] 0.3519048
# no improvement in 20 minutes

# Optimization without timelimit
result_2 <- gurobi::gurobi(model_2)

# save result
# noch nicht geschen
saveRDS(result_2,"results_credit_data/result_2.RDS")

############################################################################
# Classical Subgroup Discovery on a subset of size n=500
############################################################################

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

discovery_absb <- oofos:::discover_starshaped_subgroups(stylized_betweenness=absb,objective=objective,local_vc_dimension=100,params=list(outputflag=1)) 
