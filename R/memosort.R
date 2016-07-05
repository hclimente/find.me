#' Function to sort gene mutation matrix
#' from https://gist.github.com/dakl/5974ac1d78dac88c1c78
#' 
#' @param mutmat matrix, genes as rows, samples as columns
#' @param sortGenes boolean wheather or not to sort genes (rows)
memoSort <- function(mutmat, sortGenes=TRUE) {
  
  # order rows
  if(sortGenes){
    geneOrder <- sort(rowSums(mutmat), decreasing=TRUE, index.return=TRUE)$ix;
  } else {
    geneOrder <- 1:nrow(mutmat)
  }
  
  # order columns
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
  
  mutmat <- mutmat[geneOrder, sampleOrder]
  
  firstRow <- mutmat[1,]
  
  new.order <- c()
  for(i in unique(as.numeric(firstRow))){
    sub.mutmat <- mutmat[2:nrow(mutmat),firstRow==i]
    if (is.vector(sub.mutmat))
      new.order <- c(new.order, colnames(mutmat)[firstRow==i])
    else if (is.data.frame(sub.mutmat) & nrow(sub.mutmat)==2) 
      new.order <- c(new.order, colnames(mutmat)[firstRow==i])
    else {
      sub.mutmat <- memoSort(sub.mutmat, sortGenes = FALSE)
      new.order <- c(new.order, colnames(sub.mutmat))
    }
      
  }
  
  mutmat <- mutmat[,new.order]
  
  return(mutmat)
}