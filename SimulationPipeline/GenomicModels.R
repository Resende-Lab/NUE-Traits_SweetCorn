
#' `runBayesB`
#' Function to run a BayesB multivarite using BGLR in populations from AlphaSimR
#'
#'
#'
#' @param trainPop population where to train the model
#' @param targetPop Target population from AlphaSimR
#' @param nTraits number of traits in the simulation
#' 
#' @return estimated SNP effects for the individuals in the population
#' 
#' @import BGLR
#' 
#' @author Marco A Peixoto 


runBayesB <- function(trainPop, targetPop, nTraits = 2){
  
library(BGLR)

# 1.1 Phenotypic data
y1 = as.matrix(pheno(trainPop))
rownames(y1) = c(trainPop@id)

y2 = matrix(ncol = nTraits, nrow = length(targetPop@id))
rownames(y2) = targetPop@id

y = rbind(y1,y2)



# 1.2 SNP panel from base population
M1 = pullSnpGeno(trainPop) #Access SNP genotype data
rownames(M1) = c(trainPop@id)


# 1.3 SNP panel from target population
M2 = pullSnpGeno(targetPop) #Access SNP genotype data
rownames(M2) = targetPop@id

Markers = rbind(M1,M2)



####>>>----------- 2. The Model
# ETA
ETA_K.Linear = list(list(model='SpikeSlab', X=Markers))


# Model
gmodel = BGLR::Multitrait(y = y,
                         ETA=ETA_K.Linear,
                         nIter = 10000,
                         burnIn = 1000,
                         thin = 5)


###>>>----------- 5. Pushing data back to AlphaSimR code
addEff = as.matrix(cbind(gmodel$ETA[[1]]$b))



}


#' `runMtMGBLUP`
#' Function to run a BayesB multivarite using BGLR in populations from AlphaSimR
#'
#'
#'
#' @param trainPop population where to train the model
#' @param targetPop Target population from AlphaSimR
#' @param nTraits number of traits in the simulation
#' 
#' @return estimated SNP effects for the individuals in the population
#' 
#' @import BGLR
#' 
#' @author Marco A Peixoto 


runMtMGBLUP <- function(trainPop, targetPop, nTraits = 2){
  
  library(BGLR)
  
  # 1.1 Phenotypic data
  y1 = as.matrix(pheno(trainPop))
  rownames(y1) = c(trainPop@id)
  
  y2 = matrix(ncol = nTraits, nrow = length(targetPop@id))
  rownames(y2) = targetPop@id
  
  y = rbind(y1,y2)
  
  # 1.2 SNP panel from base population
  M1 = pullSnpGeno(trainPop) #Access SNP genotype data
  rownames(M1) = c(trainPop@id)
  
  
  # 1.3 SNP panel from target population
  M2 = pullSnpGeno(targetPop) #Access SNP genotype data
  rownames(M2) = targetPop@id
  
  Markers = rbind(M1,M2)
  
  # 1.4 Rel mat
  
  relMat = AGHmatrix::Gmatrix(Markers)
  
  
  ####>>>----------- 2. The Model
  # ETA
  ETA = list(list(model='RKHS', K=relMat))
  
  
  # Model
  gmodel = BGLR::Multitrait(y = y,
                            ETA=ETA,
                            nIter = 10000,
                            burnIn = 1000,
                            thin = 5)
  
  
  ###>>>----------- 5. Pushing data back to AlphaSimR code
  indEbv = as.matrix(gmodel$ETAHat)[((trainPop@nInd+1):(trainPop@nInd + targetPop@nInd)),]
 
  
 return(indEbv)
  
}


#' `runGBLUP`
#' Function to run a BayesB multivarite using BGLR in populations from AlphaSimR
#'
#'
#'
#' @param trainPop population where to train the model
#' @param targetPop Target population from AlphaSimR
#' 
#' @return estimated SNP effects for the individuals in the population
#' 
#' @import BGLR
#' 
#' @author Marco A Peixoto 


runGBLUP <- function(trainPop, targetPop){
  
  library(BGLR)
  
  # 1.1 Phenotypic data
  y1 = as.matrix(pheno(trainPop)[,1])
  rownames(y1) = c(trainPop@id)
  
  y2 = matrix(ncol = 1, nrow = length(targetPop@id))
  rownames(y2) = targetPop@id
  
  y = rbind(y1,y2)
  
  # 1.2 SNP panel from base population
  M1 = pullSnpGeno(trainPop) #Access SNP genotype data
  rownames(M1) = c(trainPop@id)
  
  
  # 1.3 SNP panel from target population
  M2 = pullSnpGeno(targetPop) #Access SNP genotype data
  rownames(M2) = targetPop@id
  
  Markers = rbind(M1,M2)
  
  # 1.4 Rel mat
  
  relMat = AGHmatrix::Gmatrix(Markers)
  
  
  ####>>>----------- 2. The Model
  # ETA
  ETA = list(list(model='RKHS', K=relMat))
  
  
  # Model
  gmodel = BGLR::BGLR(y = y,
                      ETA=ETA,
                      nIter = 10000,
                      burnIn = 1000,
                      thin = 5)
  
  
  ###>>>----------- 5. Pushing data back to AlphaSimR code
  indEbv = as.matrix(gmodel$yHat[((trainPop@nInd+1):(trainPop@nInd + targetPop@nInd))])
  
  
  return(indEbv)
  
}
