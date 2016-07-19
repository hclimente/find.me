#' Function to calculate a mutual exclusion p-value in a gene set.
#' It performs a Fisher's test for each gene agains the aggregation
#' of the rest. After adjusting the p-values, it returns the maximum one.
#' 
#' @param mutmat matrix, genes as rows, samples as columns
me.test.fisher <- function(mutmat){
  
  # convert to binary
  bin.mutmat <- mutmat > 0

  p <- c()
  for (i in 1:nrow(bin.mutmat)){
    this.gene <- bin.mutmat[i,]
    if (nrow(bin.mutmat) > 2) {
      agg.genes <- colSums(bin.mutmat[-i,]) > 0
    } else {
      agg.genes <- bin.mutmat[-i,]
    }
    
    # double intersection
    a <- sum(this.gene & agg.genes)
    # single cases
    b <- sum(this.gene) - a
    c <- sum(agg.genes) - a
    # rest
    d <- length(this.gene) - (a+b+c)
          
    p <- c(p,fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2), alternative="less")$p.value)
  }
      
  padj <- p.adjust(p)
  max(padj)

}