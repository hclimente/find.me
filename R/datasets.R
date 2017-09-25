#' Description of the TCGA BRCA.
#'
#' @name tcga_brca
#' @docType data
#' @description Dataframe containing the mutational status of BRCA1 and BRCA2 over a subset of patients.
#' @format A dataframe with 3 columns:
#' \describe{
#'   \item{gene}{Gene with the alteration (BRCA1 or BRCA2).}
#'   \item{sample}{Sample id.}
#'   \item{alteration}{Alteration detected in that gene (AMP for amplification, MUT for somatic mutation, HOMDEL for homozygous deletion).}
#' }
NULL