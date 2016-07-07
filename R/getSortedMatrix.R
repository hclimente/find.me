getSortedMatrix <- function(M){
  # convert from wide to long format
  all <- melt(M, varnames = c("gene", "patient"), value.name = "alteration")
  
  genes <- na.omit(unique(as.character(all$gene)))
  patients <- na.omit(unique( as.character(all$patient) ))
  
  # create data structure for each alteration in a list
  alterations <- list()
  for( a in names(keys) ){
    df <- all[ grep(pattern = keys[[a]], all$alteration) ,]
    if(nrow(df) > 0 ){
      df$alteration <- a
      alterations[[a]] <- df
    }
  }
  
  alterations <- do.call("rbind",alterations)
  
  ## create numerical mutation matrix
  mutmat <- as.data.frame(matrix(0, nc=length(patients), nr=length(genes)))
  colnames(mutmat) <- patients
  rownames(mutmat) <- genes
  
  # from https://github.com/gideonite/WIP/blob/gh-pages/oncoprint/MemoSort.js
  # // sorting order : amplification, deletion, mutation, splicing, mrna, rppa
  # // mutation > 0
  # // amp > del > 0
  scoringMatrix <- c("amp" = 128, "del" = 64, "somatic" = 32, "splicing" = 25, 
                     "germline" = 16, "up" = 8, "down" = 4)
  mutmat <- incrementMatrix(M=mutmat, events=alterations, scoring=scoringMatrix)
  mutmat <- memoSort(mutmat, sortGenes = sortGenes)

  return(mutmat)
}