#########################################
#
# Package: 
#
# File: getImp.R
# Contains: getImp
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
##########################################

#' Impute the marker anf filter out.
#'
#' @description
#' The function does the imputation by the mean of the marker and filter by missing data.
#'
#' @param Markers matrix with markers information for all candidate parents, coded as 0,1,2.
#' @param thresh minor allele frequency threshold. Default is 0.05
#' @param IM individual missing threshold.Default is 0.80
#' @param MM marker missing threshold. Default is 0.80
#'
#' @return A dosage matrix with no missing data.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the data
#' data(lines_geno)
#' 
#' # 2. Imputation
#' Markers <- getImp(lines_geno)
#'
#' }
#'
#' @export


getImp <- function(Markers, thresh = 0.05,  IM = 0.80, MM = 0.80) {
  if (!("matrix" %in% class(Markers))) {
    stop("Argument 'Markers' is not a matrix.\n")
  }
  
  SI = missing_data(H_SNP = Markers, IM = IM, MM = MM)
  
  SF = maf_filter(H_SNP = SI, thresh = thresh)
  
  SD = ImpMarker(SF)
  
  return(SD)
}





missing_data <- function(H_SNP,IM,MM){
  #Remove individuals with values higher than the threshold
  individual.missing <- apply(H_SNP,1,function(x){
    return(length(which(is.na(x)))/ncol(H_SNP))
  })
  
  #Remove markers with values higher than the threshold
  marker.missing <- apply(H_SNP,2,function(x)
  {return(length(which(is.na(x)))/nrow(H_SNP))
  })
  
  #return the individuals
  filtered <- H_SNP[which(individual.missing<IM),which(marker.missing<MM)]
  return(filtered)
}


####>>>------------ Minor allele frequency
maf_filter<-function(H_SNP,thresh){
  freq<-colMeans(H_SNP, na.rm=T)/2
  maf<-freq
  maf[which(maf > 0.5)] <- 1-maf[which(maf > 0.5)]
  
  snps1<-H_SNP[,which(maf>thresh)]
  
  return(snps1)
  
}

# Imputation
ImpMarker = function(Markers_filtered){
  MarkersFilt <- apply(Markers_filtered, 2, FUN = function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm = TRUE), ina)))
  return(MarkersFilt)
}

