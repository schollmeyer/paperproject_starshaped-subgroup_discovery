compute_rademacher_complexity <-function (incidence,n_rep=1000) {
  # This function computes the (empirical) rademacher complexity of the family
  # of all downssets w.r.t. the binary relation incidence
  # If R subseteq G^3 is a ternary relation on G, then for fixed center point
  # g_i the family of all starshaped sets w.r.t. R and with center point g_i
  # is exactly the family of all downsets w.r.t. the relation
  # R(g_i,\cdot,\cdot) (Under certain conditions on R this relation is a
  # partial order. But note that the function works for arbitrary binary
  # relations)

  # handling of special cases:
  if (is.null(incidence)) {
    return(list(complexity = 0))
  }
  n_rows <- nrow(incidence)

  # exlusion of pathological cases
  if (n_rows != ncol(incidence)) {
    print("No square matrix")
    return(NULL)
  }
  if (any(!incidence %in% c(0, 1))) {
    print("Not all entries of the incidence matrix are in{0,1}")
    return(NULL)
  }

  # instead of the downsets w.r.t. the relation incidence we look ath the upsets
  # w.r.t. the transitive reduction of the relation incidence This family
  # coincides with the family of all downsets w.r.t. the relation incidence
  incidence <- oofos:::compute_transitive_hull(incidence)
  diag(incidence) <- 1
  incidence <- oofos:::compute_pseudoreduction(incidence)

  # handling of special cases:
  if (n_rows == 0) {
    return(list(complexity = 0))
  }
  if (n_rows == 1) {
    return(list(complexity = 1))
  }
  result <- 0

  # build gurobi model for optimization on downsets
  # The function oofos:::get_model_from_quasiorder computes the model for
  # optimization on upsets. Therefore we use t(incidence) instead of incidence
  # to optimize on downsets
  model <- oofos:::get_model_from_quasiorder(t(incidence))

  # handling of special case:
  if (is.null(model)) {
    model <- list(A = matrix(0, nrow = 1, ncol = n_rows),
                  rhs = 1, sense = "<=")
  }
  # add the bounds for the decision variables
  model$lb <- rep(0, n_rows)
  model$ub <- rep(1, n_rows)

  # maximization instead of minimization
  model$modelsense <- "max"



  # generate Rademacher random variables and optimize
  # n_rep times to estimate the empirical Rademacher complexity
  for(k in seq_len(n_rep)){
    rademacher_variables <- sample(c(-1,1),size=n_rows,replace=TRUE)
    # The model is always the same. Only the objective changes for different
    # realizations of the Rademacher random variables
    model$obj <- rademacher_variables#objective
    # optimization with gurobi
    b <- gurobi::gurobi(model, params = list(outputflag=0))
    result <- result+b$objval

  }
  # return mean sup value
  return(list(complexity=result/n_rep))
}
