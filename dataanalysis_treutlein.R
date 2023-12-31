##########################################
# Import data (code from hemberg-lab): https://github.com/hemberg-lab/scRNA.seq.datasets/blob/master/R/treutlein.R
#########################################

## set working directory accordingly !!!

# load additional functions for computing stylized betweenness:
source("used_stylized_betweenness_functions.R")
# load additional functions needed for gene expression data (gene filter and scaling has to be applied):
source("additional_functions_for_gene_expression_data.R")

### DATA
d <- read.table("datasets/nature13173-s4.txt")
d <- t(d)
genes <- d[,1][5:nrow(d)]
# remove genes and bulk samples
d <- d[,2:(ncol(d) - 2)]
exprs_data <- as.data.frame(matrix(as.numeric(d[5:nrow(d),]), ncol = ncol(d)))
rownames(exprs_data) <- genes
colnames(exprs_data) <- d[1,]

### ANNOTATIONS
ann <- data.frame(cell_type1 = d[4,])
rownames(ann) <- d[1,]

# This code does not work anymore
### SINGLECELLEXPERIMENT
#source("../utils/create_sce.R")
#sceset <- create_sce_from_normcounts(exprs_data, ann)
#saveRDS(sceset, file = "treutlein.rds")

############################################

# manual construction of x and y:
 x <- t(exprs_data)
 x <- gene_filter(x)
 x <- log2(1+x)
 dim(x)
 x <- scaling(x)
 
 context <- oofos:::get_auto_conceptual_scaling(x)
 
 # This context has VC dimension 80!!!
 
 y <- ann$cell_type1
 table(y)/length(y)
 
# cell_type1
#     AT1      AT2       BP ciliated    Clara 
#      41       12       13        3       11

objective <- oofos:::compute_objective(data.frame(y=y),"y","AT1") 
obsb <- get_obsb(context)
absb <- get_absb(context)
gbsb <- get_gbsb(x)

saveRDS(obsb,"results_treutlein/obsb.RDS")
saveRDS(absb,"results_treutlein/absb.RDS")
saveRDS(gbsb,"results_treutlein/gbsb.RDS")


# absb:
starshaped_discovery_absb <- oofos::discover_starshaped_subgroups(stylized_betweenness=absb,objective=objective,local_vc_dimension=8)
starshaped_discovery_absb$objval
test_absb <- oofos:::compute_starshaped_distr_test(starshaped_discovery_absb, n_rep=100)
saveRDS(starshaped_discovery_absb,"results_treutlein/starshaped_discovery_absb")
saveRDS(test_absb,"results_treutlein/test_absb")

# gbsb:

starshaped_discovery_gbsb <- oofos::discover_starshaped_subgroups(stylized_betweenness=gbsb,objective=objective,local_vc_dimension=8)
starshaped_discovery_gbsb$objval
test_gbsb <- oofos:::compute_starshaped_distr_test(starshaped_discovery_gbsb, n_rep=100)
saveRDS(starshaped_discovery_gbsb,"results_treutlein/starshaped_discovery_gbsb")
saveRDS(test_gbsb,"results_treutlein/test_gbsb")



## gene expression dataset (Treutlein: Reconstructing lineage hierarchies of the distal lung epithelium using single-cell RNA-seq
## https://www.nature.com/articles/nature13173


# needed packages

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scRNAseq")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

library(SingleCellExperiment)

# read treulein data

dat <- readRDS("datasets/treutlein.rds")

y <- colData(dat)$cell_type1
x <- counts(dat)
x <- t(x)

# apply gene filter like described here:  ????

i <-NULL

for(k in seq_len(ncol(x))){

if(mean(x[,k]> 2) <= 0.06 | mean(x[,k]>0) >=.94){i <- c(i,k)}
}

dim(x[,-i])

x <- x[,-i]
x <- t(unique(t(x)))
dim(x)
X <- get_auto_conceptual_scaling(x)

dim(X)

X <- t(unique(t(X)))

dim(X)

CS <- colSums(X)

i <- which(CS==79)
length(i)
# 80
# VC dimension: 80


# define context
context <- X

# compute objective (target is "AT1")
objective <- compute_objective(as.data.frame(y),"y","AT1")
stylized_betweenness <- ocompute_attribute_counting_betweenness(t(context))

discovery <- discover_starshaped_subgroups(stylized_betweenness = A*B,objective=objective,local_vc_dimension=40)

discovery$objval
discovery_test <- compute_starshaped_distr_test(discovery)
discovery_test$p_value
