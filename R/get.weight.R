#' Weight function for  a binary alteration matrix of genes and patients.
#' As described in http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3266044/
#' 
#' @param bin.mutmat logical matrix, genes as rows, samples as columns
get.weight <- function(bin.mutmat){
  # number of patients with each gene altered
  g <- rowSums(bin.mutmat)
  # number of patients with at least one alteration
  M <- sum(colSums(bin.mutmat) > 0)
  
  W <- 2*M - sum(g)
  
  return(W)
}