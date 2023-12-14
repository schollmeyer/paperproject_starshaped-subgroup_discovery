# working directory setzen !!!:

# setwd()

# needed packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SC3")

library(SingleCellExperiment)
library(SC3)

## additional functions needed:

source("additional_functions_for_gene_expression_data.R")
source("used_stylized_betweenness_functions.R")


sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(yan),
        logcounts = log2(as.matrix(yan) + 1)
    ), 
    colData = ann
)

# preprocessing counts (gene filter plus scaling)
x <- t(counts(sce))
x <- gene_filter(x)
x <- log2(1+x)
x <- scaling(x)

#
# generate formal context (with interordinal scaling) for some nalysis of x
context <- oofos:::get_auto_conceptual_scaling(x)

dim(context)
# remove duplicated columns
context <- t(unique(t(context)))
dim(context)

set.seed(1234567)
indexs <- sample(seq_len(ncol(context)),size=30000)
vc_model <- oofos::compute_extent_vc_dimension(context[,indexs])

vc=gurobi::gurobi(vc_model,list(timelimit=60*20))


# compute objective
y <- (colData(sce))[,1]
table(y)/length(y)

# y
#     16cell      2cell      4cell      8cell      blast     zygote 
# 0.17777778 0.06666667 0.13333333 0.22222222 0.33333333 0.06666667 

# Analysis of target "blast"
# objective <- oofos:::compute_objective(data.frame(y=y%in%c("blast","8cell")),"y","TRUE")

 objective <- oofos:::compute_objective(data.frame(y=y%in%c("blast")),"y","TRUE")
 table(objective)


 absb <- get_absb(x)
 max_absb <- max(absb)
 local_vc_dimensions <- rep(0,nrow(x))
 for(k in seq_len(nrow(x))){local_vc_dimensions[k] <- oofos:::compute_width(absb[k,,]>=max_absb/2.5)$width}
 
 table(local_vc_dimensions)
 
# local_vc_dimensions
# 89 
# 90 


starshaped_discovery <- oofos::discover_starshaped_subgroups(stylized_betweenness=absb,objective=objective,local_vc_dimension=9)
starshaped_discovery$objval
test <- oofos:::compute_starshaped_distr_test(starshaped_discovery, n_rep=100)
