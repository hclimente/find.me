#' Function to increment gene matrix
#' @param M data frame, colnames=patients, rownames=genes
#' @param events data frame. Needs columns patient and gene
#' @param inc number by which to increment
incrementMatrix <- function(M, events, scoring){
  for(k in 1:nrow(events)){
    alteration <- events$alteration[k]
    increment <- as.numeric(scoring[alteration])
    gene <- as.character(events$gene[k])
    patient <- as.character(events$patient[k])
    idx.g <- which(rownames(M)==gene)
    idx.s <- which(colnames(M)==patient)
    if(!is.na(gene)){ 
      M[idx.g,idx.s] <- M[idx.g,idx.s] + increment
    }
  }
  M
}

