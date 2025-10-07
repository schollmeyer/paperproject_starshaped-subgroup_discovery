# needed libraries:
# celltrackR


# object based stylized betweenness (OBSB)
get_object_based_stylized_betweenness <- function(context,lambda=1){
   n_objects <- nrow(context)
   stylized_betweenness <- array(as.logical(FALSE),rep(n_objects,3))
   for( k in seq_len(n_objects)){
       intent_1 <- context[k,]
	   print(k)
	for( m in seq_len(n_objects)){
	    intent_2 <- pmin(intent_1 , context[m,])#pmin(intent_1,context[m,])
		extent_2 <- oofos:::compute_phi(intent_2,context)
		size_extent_2 <- sum(extent_2)
		for( l in seq_len(n_objects)){
		   intent_3 <- pmin(intent_2 , context[l,])#pmin(intent_2, context[l,])
		   extent_3 <- oofos:::compute_phi(intent_3,context)
		   size_extent_3 <- sum(extent_3)
		   stylized_betweenness[k,l,m] <-  1- (size_extent_3 - size_extent_2)/(lambda* n_objects+(1-lambda)* size_extent_3)
		}
	}
  }
return(stylized_betweenness)}
get_obsb <- get_object_based_stylized_betweenness



# object based stylized betweenness parallelized (OBSB)
get_object_based_stylized_betweenness_par <- function(context,lambda=1,n_tasks=6){
  library(doParallel)
  library(tcltk)
  library(doSNOW)
  pb <- tkProgressBar(max=nrow(context))
  progress <- function(k) setTkProgressBar(pb, k)
  opts <- list(progress=progress)

  cl <- makeSOCKcluster(n_tasks)
  registerDoSNOW(cl)

  #cl <- makeCluster(n_tasks,outfile="")
  #registerDoParallel(cl)
  n_objects <- nrow(context)
  #stylized_betweenness <- array(as.logical(FALSE),rep(n_objects,3))
  foreach( k = (1:n_objects), .options.snow=opts)%dopar%{
    stylized_betweenness <- array(as.logical(FALSE),rep(n_objects,2))
    intent_1 <- context[k,]

    for( m in seq_len(n_objects)){
      intent_2 <- pmin(intent_1 , context[m,])#pmin(intent_1,context[m,])
      extent_2 <- oofos:::compute_phi(intent_2,context)
      size_extent_2 <- sum(extent_2)
      for( l in seq_len(n_objects)){
        intent_3 <- pmin(intent_2 , context[l,])#pmin(intent_2, context[l,])
        extent_3 <- oofos:::compute_phi(intent_3,context)
        size_extent_3 <- sum(extent_3)
        stylized_betweenness[l,m] <-  1- (size_extent_3 - size_extent_2)/(lambda* n_objects+(1-lambda)* size_extent_3)
      }
    }
    return(stylized_betweenness)}
  stopCluster(cl)
  }

get_obsb_par <- get_object_based_stylized_betweenness_par


# attribute based stylized betweenness (ABSB)
get_attribute_based_stylized_betweenness <- function(context,lambda=1){
   n_objects <- nrow(context)
   n_attributes <- ncol(context)
   stylized_betweenness <- array(as.logical(FALSE),rep(n_objects,3))
   context <- t(context)
   dims <- dim(context)
   context <- as.logical(context)
   dim(context) <- dims
   for( k in seq_len(n_objects)){
       intent_1 <- context[,k]
	   print(k)
	for( m in seq_len(n_objects)){
	    intent_2 <- intent_1 & context[,m]##pmin(intent_1,context[m,])
		#extent_2 <- oofos:::compute_phi(intent_2,context)
		size_intent_2 <- sum(intent_2)
		for( l in seq_len(n_objects)){
		   intent_3 <- intent_2 & context[,l]#pmin(intent_2, context[l,])
		   #extent_3 <- oofos:::compute_phi(intent_3,context)
		   size_intent_3 <- sum(intent_3)
		   stylized_betweenness[k,l,m] <-  1- (size_intent_2 - size_intent_3)/(lambda* n_attributes+(1-lambda)* size_intent_2)
		}
	}
  }
return(stylized_betweenness)}
get_absb <- get_attribute_based_stylized_betweenness

# geometry based stylized betweenness (ABSB)
get_geometry_based_stylized_betweenness <- function(X){
  n_objects <- nrow(X)
  stylized_betweenness <- array(0,rep(n_objects,3))
  indexs <- expand.grid(seq_len(n_objects),seq_len(n_objects))
  for(l in (1:n_objects)){
    print(l)
	difference_vectors <- X
	for(k in (1:nrow(X))){difference_vectors[k,] <- X[k,]-X[l,]}
	stylized_betweenness[,l,] <- celltrackR::vecAngle(difference_vectors[indexs$Var1,],difference_vectors[indexs$Var2,])

  }

  for(l in (1:n_objects)){
    diag(stylized_betweenness[l,,]) <- 180
	diag(stylized_betweenness[,,l]) <- 180
  }
return(stylized_betweenness/180)}

get_gbsb <- get_geometry_based_stylized_betweenness


##
## veralteter Code

compute_attribute_counting_betweenness <- function(context) {
  n_rows <- nrow(context)
  n_cols <- ncol(context)
  col_means <- colMeans(context)

  betweenness <- array(0, rep(n_cols, 3))
  pb <- utils::txtProgressBar(min = 0, max = n_cols, initial = 0)
  for (k in seq_len(n_cols)) {
    utils::setTxtProgressBar(pb, k)
    for (m in seq_len(n_cols)) {
      #temp <- rep(0, n_cols)
      #i <- which(context[,k ] == 1 & context[,m ] == 1)
      mask <- context[,k] & context[,l]
      #temp[i] <- col_means[i]
      for (l in seq_len(n_cols)) {
        betweenness[k, l, m] <- sum(context[mask,l])#max(temp - context[l, ])
      }
    }
  }

  close(pb)
  return(1 - betweenness)
}
