#' Function to calculate a mutual exclusion p-value in a gene set.
#' It keeps the number of alterations in a gene fixed, and permutates
#' the samples to calculate an empirical p-value. The minimum p-value
#' it will return is 1/n.
#' 
#' @param mutmat matrix, genes as rows, samples as columns
#' @param n integer specifying the number of permutations
me.test.permutateSamples <- function(mutmat, n=10000){
  
  # convert to binary
  bin.mutmat <- mutmat > 0
  
  # k number of samples explained
  k <- get.weight(bin.mutmat)
  
  # get random distribution of samples explained
  K <- numeric()
  for (i in 1:n){
    perm <- apply(bin.mutmat,1,sample)
    perm <- t(perm)
    W <- get.weight(perm)
    K <- c(K,W)
  }
  
  # if only one data point, pvalue cannot be calculated
  if(length(unique(K)) == 1){
  	NA
  } else {
  	# calculate distribution and return empirical p
    # add 0.1 to avoid overly optimistic pvalues when k falls in one of the extremes
  	empirical.test <- ecdf(K+0.1)
    p <- 1 - empirical.test(k)
	
	  # if p is 0, return the minimum possible p-value that makes sense
	  min.p <- 1/n
	
	  max(p,min.p)
  }
}