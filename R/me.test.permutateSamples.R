#' Function to calculate a mutual exclusion p-value in a gene set.
#' It keeps the number of alterations in a gene fixed, and permutates
#' the samples.
#' 
#' @param mutmat matrix, genes as rows, samples as columns
#' @param n integer specifying the number of permutations
me.test.permutateSamples <- function(mutmat, n=10000){
  
  # convert to binary
  bin.mutmat <- mutmat > 0
  
  # k number of samples explained
  k <- sum(colSums(bin.mutmat) > 0)
  
  # get random distribution of samples explained
  K <- numeric()
  for (i in 1:n){
    perm <- apply(bin.mutmat,1,sample)
    perm <- t(perm)
    W <- get.weight(perm)
    K <- c(K,W)
  }
  
  # calculate distribution and return empirical p
  empirical.test <- ecdf(K)
  p <- 1 - empirical.test(k)
  
  # if p is 0, return the minimum possible p-value that makes sense
  min.p <- 1/(n+1)
  
  max(p,min.p)
}