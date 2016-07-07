#' function to create an oncoprint plot
#'
#' @param M mutation matrix
#' @param keys list with the following elements: splicing, somatic, 
#' germline, amp, del, upreg, downreg
#' @param sortGenes boolean whether or not to sort the genes, default TRUE
oncoprint <- function(M, keys=list(somatic="MUT", germline="GERMLINE", amp="AMP", 
                      del="HOMDEL", upreg="UP", downreg="DOWN", splicing="SPLICING"), 
                      sortGenes=TRUE){
  
  mutmat <- getSortedMatrix(M)

  # order alterations and patients based on matrix
  alterations$gene <- factor(alterations$gene, levels=rev(rownames(mutmat)))
  alterations$patient <- factor(alterations$patient, levels=colnames(mutmat))
  
  # create background ie cases with no alteration found
  background <- as.data.frame(matrix(0, nc=length(patients), nr=length(genes)))
  colnames(background) <- patients
  background$gene <- factor(genes, levels=rev(rownames(mutmat)))
  
  background.m <- melt(background, id.vars = "gene", variable.name =  "patient")
  background.m$patient <- factor(background.m$patient, levels=colnames(mutmat))
  
  plot.params <- data.frame(alteration = c("amp","del","somatic","splicing",
                                           "germline","upreg","downreg"), 
                            size=c(2,2,1,1,1,2,2),
                            width=c(.9,.9,.9,.9,.9,.9,.9),
                            height=c(.9,.9,.4,.95,.4,.9,.9))
  alterations <- merge(alterations,plot.params)
  
  plot.fill <- c("amp" = "firebrick", "del" = "blue", "up" = NA, "down" = NA,
                 "splicing" = "forestgreen", "germline" = "purple", "somatic" = "#36454F")
  plot.alpha <- c("amp" = 0.6, "del" = 0.6, "up" = 0.6, "down" = 0.6,
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