gene_filter <- function(data, x_percent=6 ){
  # gene filter as described under https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/
  # Note that also genes with sdtandard deviation of epsression zero are also deleted
  i <- NULL
  for(k in seq_len(ncol(data))){
    if( mean(data[,k]>2) < x_percent | mean(x[,k]>0)> 1-x_percent |  sd(x[,k])==0){i <- c(i,k)}
  }
return(data[,-i])}

scaling <- function(data){
  # function that scales every column of the data matrix data by substracting the mean and dividing by the standard deviation afterwards
  for(k in seq_len(ncol(data))){
    if(sd(data[,k])==0){print("Warning: standard deviation of 0 observed")}
    data[,k] <- (data[,k]-mean(data[,k]))/sd(data[,k])
   }
return(data)}
