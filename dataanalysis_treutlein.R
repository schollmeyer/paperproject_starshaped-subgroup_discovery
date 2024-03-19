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
 
 # y
#  AT1      AT2       BP ciliated    Clara 
#  0.5125   0.1500   0.1625   0.0375   0.1375 
 
objective <- oofos:::compute_objective(data.frame(y=y),"y","AT1") 

gbsb <- get_gbsb(x)



saveRDS(gbsb,"results_treutlein/absb.RDS")


# gbsb:
starshaped_discovery_gbsb <- oofos::discover_starshaped_subgroups(stylized_betweenness=gbsb,objective=objective,local_vc_dimension=8)
starshaped_discovery_gbsb$objval
test_gbsb <- oofos:::compute_starshaped_distr_test(starshaped_discovery_gbsb, n_rep=10000)
saveRDS(starshaped_discovery_gbsb,"results_treutlein/starshaped_discovery_gbsb")
saveRDS(test_gbsb,"results_treutlein/test_gbsb")





