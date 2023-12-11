# set working directory
setwd("C:/Git/paperproject_starshaped-subgroup_discovery")
library(foreign)
dat <- read.arff("dataset_31_credit-g.arff")
dim(dat)
# 

#generate formal context
context <- oofos:::get_auto_conceptual_scaling(dat[,-21])

# remove duplicated columns
context <- t(unique(t(context)))
dim(context)

# generate objective vector for optimization. This corresponds to the
# Piatetsky-Shapiro quality function for the target value credit-state==good
objective <- oofos:::compute_objective(dat,"class","good")

table(objective)
#
# 300 good states and 700 bad states


############################################################################
# Classical Subgroup Discovery
############################################################################


# generate gurobi model for classical subgroup discovery
model <- oofos:::optimize_on_context_extents(context=context,objective=objective,binary_variables="all" )

# optimize model under time constraint of 20 minutes for a first semioptimal solution
result_1 <- gurobi::gurobi(model,list(timelimit=60*20))

result_1$objval

result_1$runtime
# save result
saveRDS(result_1,"results_credit_data/result_1.RDS")




# add additional contsraints that are always valid for inetrordinally scaled variables x (implications of the form <= x \rightarrow not >= x+1 etc. 
model_2 <- oofos:::add_attr_antiimplications(model)

# use semioptima solution to tighten tthe search space via optimistic estimates (in the style of ***)
model_2 <- oofos:::add_sos_constraints(model_2,result_1$objval)
# add semioptimal solution as a start vector
model_2$start <- round(result_1$x,2)
result_2 <- gurobi::gurobi(model_2,list(timelimit=60*20))
