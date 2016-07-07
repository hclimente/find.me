#' Function to calculate a mutual exclusion p-value in a gene set.
#' It keeps the number of alterations in a gene fixed, and permutates
#' the samples.
#' 
#' @param mutmat matrix, genes as rows, samples as columns
#' @param n integer specifying the number of permutations
me.test.permutateSamples <- function(mutmat, n=1000){
  
  # convert to binary
  bin.mutmat <- mutmat > 0
  
  # k number of samples explained
  k <- sum(colSums(bin.mutmat) > 0)
  
  # get random distribution of samples explained
  K <- numeric()
  for (i in 1:n){
    perm <- apply(bin.mutmat,1,sample)
    K <- c(K,sum(rowSums(perm) > 0))
  }
  
  # calculate distribution and return empirical p
  t <- ecdf(K)
  p <- 1 - t(k)
  
  # if p is 0, return the minimum possible p-value that makes sense
  min.p <- 1/(n+1)
  
  max(p,min.p)
}