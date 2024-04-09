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

absb <- get_absb(context)
saveRDS(absb,"results_treutlein_absb/absb.RDS")

gbsb <- readRDS("results_treutlein_gbsb/gbsb.RDS")

plot(absb,gbsb)




# absb:
discovery <- oofos::discover_starshaped_subgroups(stylized_betweenness=absb,objective=objective,complexity_control=8)
discovery$objval
test_gbsb <- oofos:::compute_starshaped_distr_test(discovery, n_rep=20000)
saveRDS(starshaped_discovery_gbsb,"results_treutlein_absb/starshaped_discovery.RDS")
saveRDS(test_gbsb,"results_treutlein_absb/test.RDS")


##### extensive analysis for gbsb
set.seed(1234567)
vc_dimensions <- seq(1,80,length.out=500)
p_values <- rep(0,length(vc_dimensions))
p_values_param <- rep(0,length(vc_dimensions))
objvalues <- array(0,c(500,500))
objval <- rep(0,500)
for(k in (1:500)){
  discovery <- oofos::discover_starshaped_subgroups(absb,objective=objective,complexity_control = vc_dimensions[k])
  test <- oofos::compute_starshaped_distr_test(discovery,n_rep=500)
  objvalues[k,] <- test$objvalues
  objval[k] <- discovery$objval
  p_values[k] <- test$p_value
  p_values_param[k] <-  (discovery$objval-mean(test$objvalues))/sd(test$objvalues)
  plot(vc_dimensions,(p_values_param))
}

# saveRDS(p_values_param,"results_treutlein_absb/p_values_param.RDS")
# saveRDS(objvalues,"results_treutlein_absb/objvalues.RDS")
# saveRDS(objval,"results_treutlein_absb/objval.RDS")
# saveRDS(vc_dimensions,"results_treutlein_absb/vc_dimensions.RDS")


pplot <- ggplot(data=data.frame(x=vc_dimensions,y=p_values_param), aes(x=x,y=y))


pplot + layer(mapping = NULL,   position = "identity",   stat="identity",  geom = "point") +labs(x = "VC dimension",y="significance in standard deviations")



