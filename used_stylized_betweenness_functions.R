# needed libraries:
# celltrackR


# object based stylized betweenness (OBSB)
get_object_based_stylized_betweenness <- function(context){
   n_objects <- nrow(context)
   stylized_betweenness <- array(as.logical(FALSE),rep(n_objects,3))
   for( k in seq_len(n_objects)){
       intent_1 <- context[k,]
	   print(k)
	for( m in seq_len(n_objects)){
	    intent_2 <- intent_1 & context[m,]#pmin(intent_1,context[m,])
		extent_2 <- oofos:::compute_phi(intent_2,context)
		size_extent_2 <- sum(extent_2)
		for( l in seq_len(n_objects)){
		   intent_3 <- intent_2 & context[l,]#pmin(intent_2, context[l,])
		   extent_3 <- oofos:::compute_phi(intent_3,context)
		   size_extent_3 <- sum(extent_3)
		   stylized_betweenness[k,l,m] <- size_extent_2/size_extent_3
		}
	}
  }
return(stylized_betweenness)}
get_obsb <- get_object_based_stylized_betweenness

# angle based stylized betweenness (ABSB)
get_angle_based_stylized_betweenness <- function(X){
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

get_absb <- get_angle_based_stylized_betweenness
