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


# generate gurobi model for classical subgroup discovery
model <- oofos:::optimize_on_context_extents(context=context,objective=objective,binary_variables="all" )

# optimize model under time constraint of 20 minutes for a first semioptimal solution
result_1 <- gurobi::gurobi(model,list(timelimit=60*20))

saveRDS(result_1,"result_1.RDS")

model_2 <- oofos:::add_attr_antiimplications(model)
model_2 <- oofos:::add_sos_constraints(model_2,0.35)#result_1$objval)
model_2$start <- round(result_1$x,2)
result_2 <- gurobi::gurobi(model_2)#,list(timelimit=60*5))
