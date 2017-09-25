#' function to create an oncoprint plot
#'
#' @param A dataframe containing the alterations
#' @param keys list with the following elements: splicing, somatic, 
#' germline, amp, del, upreg, downreg
#' @param sortGenes boolean whether or not to sort the genes, default TRUE
#' @importFrom ggplot2 ggplot labs scale_alpha_manual scale_fill_manual theme theme_minimal geom_tile element_text aes 
#' @importFrom tidyr gather
#' @export
oncoprint <- function(A, 
                      keys = list(somatic = "MUT", germline = "GERMLINE", amp = "AMP", 
                                  del = "HOMDEL", upreg = "UP", downreg = "DOWN", 
                                  splicing = "SPLICING"), 
                      sortGenes = TRUE){
  
  M <- long2wide(A)
  M.sorted <- getSortedMatrix(M,keys,sortGenes)
  mutmat <- M.sorted$mutmat
  alterations <- M.sorted$alterations
  genes <- M.sorted$genes
  patients <- M.sorted$patients

  # order alterations and patients based on matrix
  alterations$gene <- factor(alterations$gene, levels=rev(rownames(mutmat)))
  alterations$patient <- factor(alterations$patient, levels=colnames(mutmat))
  
  # create background ie cases with no alteration found
  background <- as.data.frame(matrix(0, ncol=length(patients), nrow=length(genes)))
  colnames(background) <- patients
  background$gene <- factor(genes, levels=rev(rownames(mutmat)))
  
  background.m <- gather(background, gene, patient)
  background.m$patient <- factor(background.m$patient, levels=colnames(mutmat))
  
  plot.params <- data.frame(alteration = c("amp","del","somatic","splicing",
                                           "germline","upreg","downreg"), 
                            size=c(2,2,1,1,1,2,2),
                            width=c(.9,.9,.9,.9,.9,.9,.9),
                            height=c(.9,.9,.4,.95,.4,.9,.9))
  alterations <- merge(alterations,plot.params)
  
  plot.fill <- c("amp" = "firebrick", "del" = "blue", "upreg" = NA, "downreg" = NA,
                 "splicing" = "forestgreen", "germline" = "purple", "somatic" = "#36454F")
  plot.alpha <- c("amp" = 0.6, "del" = 0.6, "upreg" = 0.6, "downreg" = 0.6,
                  "splicing" = 0.6, "germline" = 1, "somatic" = 1)
    
  ggplot() + 
    geom_tile(data = background.m, aes(x = patient, y = gene), 
              fill = "#DCDCDC", colour = "white", size = 1.1) + 
    geom_tile(data = alterations, aes(x = patient, y = gene, fill = alteration, 
                                      size = size, width = width, height = height,
                                      alpha = alteration)) +
    scale_fill_manual(values = plot.fill) +
    scale_alpha_manual(values = plot.alpha) +
    theme_minimal() + 
    labs(x = "Sample", y = "Gene") +
    theme(axis.text.x = element_text(angle = 90, size = 9), legend.position = "none")
  
}

#' @importFrom dplyr %>% group_by mutate select
long2wide <- function(alterations) {
  
  colnames(alterations) <- c("sample", "gene", "alteration")
  
  alterations <- alterations %>%
    group_by(sample, gene) %>%
    mutate(alterations = paste(alteration, collapse = ";")) %>%
    select(sample, gene, alterations) %>%
    unique
  
  wide <- spread(alterations, sample, alterations) %>%
    as.matrix
  rownames(wide) <- wide[,"gene"]
  wide <- wide[,-1]
  
  return(wide)
  
}

#' Function to order the alteration matrix to evidence the mutual exclusion pattern.
#' 
#' @param M mutation matrix
#' @param keys list with the following elements: splicing, somatic, 
#' germline, amp, del, upreg, downreg
#' @param sortGenes boolean whether or not to sort the genes, default TRUE
#' @importFrom stats na.omit
#' @importFrom reshape2 melt
getSortedMatrix <- function(M, keys=list(somatic="MUT", germline="GERMLINE", amp="AMP", 
                                         del="HOMDEL", upreg="UP", downreg="DOWN", 
                                         splicing="SPLICING"), sortGenes=TRUE){
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
  mutmat <- as.data.frame(matrix(0, ncol=length(patients), nrow=length(genes)))
  colnames(mutmat) <- patients
  rownames(mutmat) <- genes
  
  # from https://github.com/gideonite/WIP/blob/gh-pages/oncoprint/MemoSort.js
  # // sorting order : amplification, deletion, mutation, splicing, mrna, rppa
  # // mutation > 0
  # // amp > del > 0
  scoringMatrix <- c("amp" = 128, "del" = 64, "somatic" = 32, "splicing" = 25, 
                     "germline" = 16, "up" = 8, "downreg" = 4)
  mutmat <- incrementMatrix(M=mutmat, events=alterations, scoring=scoringMatrix)
  mutmat <- memoSort(mutmat, sortGenes = sortGenes)
  
  list(mutmat=mutmat,alterations=alterations,genes=genes,patients=patients)
}

#' Function to increment gene matrix
#' @param M data frame, colnames=patients, rownames=genes
#' @param events data frame. Needs columns patient and gene
#' @param scoring number by which to increment
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