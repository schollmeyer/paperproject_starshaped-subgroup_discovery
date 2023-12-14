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
