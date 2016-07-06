#' function to create an oncoprint plot
#'
#' @param M mutation matrix
#' @param keys list with the following elements: splicing, somatic, germline, amp, del, upreg, downreg
#' @param sortGenes boolean whether or not to sort the genes, default TRUE
oncoprint <- function(M, keys=list(somatic="MUT", germline="GERMLINE", amp="AMP", 
                      del="HOMDEL", upreg="UP", downreg="DOWN", splicing="SPLICING"), 
                      sortGenes=TRUE){
  
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

  # order alterations and patients based on matrix
  alterations$gene <- factor(alterations$gene, levels=rev(rownames(mutmat)))
  alterations$patient <- factor(alterations$patient, levels=colnames(mutmat))
  
  # create background ie cases with no alteration found
  background <- as.data.frame(matrix(0, nc=length(patients), nr=length(genes)))
  colnames(background) <- patients
  background$gene <- factor(genes, levels=rev(rownames(mutmat)))
  
  background.m <- melt(background, id.vars = "gene", variable.name =  "patient")
  background.m$patient <- factor(background.m$patient, levels=colnames(mutmat))
  
  plot.params <- data.frame(alteration=c("amp","del","somatic","splicing",
                                         "germline","upreg","downreg"), 
                            size=c(2,2,1,1,1,2,2),
                            width=c(.9,.9,.9,.9,.9,.9,.9),
                            height=c(.9,.9,.4,.95,.4,.9,.9))
  alterations <- merge(alterations,plot.params)
  
  plot.colors <- c("amp" = "firebrick", "del" = "blue", "up" = NA, "down" = NA,
                   "splicing" = "forestgreen", "germline" = "purple", "somatic" = "#36454F")
    
  ggplot() + 
    geom_tile(data=background.m, aes(x=patient, y=gene), fill="#DCDCDC", colour="white", size=1.1) + 
    geom_tile(data=alterations, aes(x=patient, y=gene, fill=alteration, size=size, width=width, height=height),alpha=0.5) +
    scale_fill_manual(values=plot.colors) +
    theme_minimal() + 
    labs(x="Sample", y="Gene") +
    theme(axis.text.x=element_text(angle=90,size=9), legend.position="none")
  
}