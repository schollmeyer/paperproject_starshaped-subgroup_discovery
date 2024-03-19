compute_rademacher_complexity <-function (incidence,n_rep=10) {
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
    incidence <- oofos:::compute_quotient_order(incidence)
    n_rows <- nrow(incidence)
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
	  b <- gurobi::gurobi(model, params = list(outputflag=0))#params)
	  result <- result+b$objval
	
	}
	
    return(list(complexity=result/n_rep))
}


######


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
