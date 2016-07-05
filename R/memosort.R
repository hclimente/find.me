#' Function to sort gene mutation matrix
#' from https://gist.github.com/dakl/5974ac1d78dac88c1c78
#' 
#' @param mutmat matrix, genes as rows, samples as columns
#' @param sortGenes boolean wheather or not to sort genes (rows)
memoSort <- function(mutmat, sortGenes=TRUE) {
  if(sortGenes){
    geneOrder <- sort(rowSums(mutmat), decreasing=TRUE, index.return=TRUE)$ix;
  } else {
    geneOrder <- 1:nrow(mutmat)
  }
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(20+length(x)-i);        
        break
      }
    }
    score <- score + sum(x*(length(x):1))
    return(score);
  }
  scores <- apply(mutmat[geneOrder, ], 2, scoreCol);
  sampleOrder <- order(scores, decreasing = TRUE)# sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(mutmat[geneOrder, sampleOrder]);
}
