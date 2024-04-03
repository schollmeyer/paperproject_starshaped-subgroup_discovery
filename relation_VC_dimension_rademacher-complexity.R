compute_rademacher_complexity <-function (incidence,n_rep=1000) {
  if (is.null(incidence)) {
    return(list(complexity = 0))
  }
  n_rows <- nrow(incidence)
  if (n_rows != ncol(incidence)) {
    print("No square matrix")
    return(NULL)
  }
  if (any(!incidence %in% c(0, 1))) {
    print("Not all entries of the incidence matrix are in{0,1}")
    return(NULL)
  }
  incidence <- oofos:::compute_transitive_hull(incidence)
  diag(incidence) <- 1
  incidence <- oofos:::compute_pseudoreduction(incidence)
  
  if (n_rows == 0) {
    return(list(complexity = 0))
  }
  if (n_rows == 1) {
    return(list(complexity = 1))
  }
  result <- 0
  
  
  model <- oofos:::get_model_from_quasiorder(t(incidence))
  if (is.null(model)) {
    model <- list(A = matrix(0, nrow = 1, ncol = n_rows), 
                  rhs = 1, sense = "<=")
  }
  
  model$lb <- rep(0, n_rows)
  model$ub <- rep(1, n_rows)
  
  model$modelsense <- "max"
  
  
  
  
  for(k in seq_len(n_rep)){
    rademacher_variables <- sample(c(-1,1),size=n_rows,replace=TRUE)
    model$obj <- rademacher_variables#objective
    b <- gurobi::gurobi(model, params = list(outputflag=0))
    result <- result+b$objval
    
  }
  
  return(list(complexity=result/n_rep))
}


cut_incidence <- function(incidence, cut_value, complexity_measure=width,interval = stats::quantile(unique(as.vector(incidence)), 
                                                                                                    c(0, 1))) 
{
  vc <- (complexity_measure(compute_quotient_order(compute_transitive_hull(incidence >= 
                                                                             interval[2]))))[[1]]
  if (vc <= cut_value) {
    return(incidence >= interval[2])
  }
  f <- function(C, incidence) {
    complexity <- (complexity_measure(compute_quotient_order(compute_transitive_hull(incidence >= 
                                                                                       C))))[[1]]
    return(complexity - cut_value)
  }
  ans <- stats::uniroot(f, interval = interval, incidence = incidence)
  return(list(root=ans$root,incidence=compute_transitive_hull(incidence >= ans$root)))
}


n_rep <- 10000
N <- 500
rademacher_complexity <- rep(0,n_rep)
vc_dimension <- rep(0,n_rep)
cutting_value <- rep(0,n_rep)


F_obsb <- ecdf(obsb)
obsb <- F_obsb(obsb)
dim(obsb) <- c(N,N,N)
for(k in (1:n_rep)){
  l <- sample((1:N),size=1)

  
  cut_value <- runif(1)
  
  vc_dimension[k] <- oofos::compute_width(obsb[l,,]>=cut_value)$width
  #vc_dimension[k] <- sample((10:300),size=1)
  #result <- cut_incidence(obsb[k,,],cut_value=vc_dimension[k], complexity_measure=width)
  #rademacher_complexity[k] <- compute_rademacher_complexity(obsb[l,,]>=cut_value)$complexity#result$incidence)$complexity
  cutting_value[k] <- cut_value
  plot(log2(vc_dimension),(cutting_value))#rademacher_complexity)
   
}

cor(rademacher_complexity[(1:k)],vc_dimension[(1:k)]) 

# [1] 0.9993353

summary(lm(rademacher_complexity[(1:k)]~vc_dimension[(1:k)]))
  
# Call:
#   lm(formula = rademacher_complexity[(1:k)] ~ vc_dimension[(1:k)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -6.7397 -0.3556  0.0203  0.2680  3.3644 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.185801   0.049516   3.752 0.000207 ***
#   vc_dimension[(1:k)] 0.546162   0.001093 499.493  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8016 on 332 degrees of freedom
# Multiple R-squared:  0.9987,	Adjusted R-squared:  0.9987 
# F-statistic: 2.495e+05 on 1 and 332 DF,  p-value: < 2.2e-16
