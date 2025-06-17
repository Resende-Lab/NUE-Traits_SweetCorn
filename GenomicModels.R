########################################################################################
##### Models for RGS paper - Finegan et al. 2025
##### Dr. Marco Peixoto
##### 12/10/2024
########################################################################################



#######--------------------------------x
### Packages and Environment
#######--------------------------------
rm(list=ls())

library(AGHmatrix)
library(BGLR)
library(tidyverse)

source("1.Functions.R")


#######--------------------------------
### Loading datasets for lines
#######--------------------------------

### Phenotypic data for lines
PhenoLines = read.table("data/LinesBLUEs", h=T)
PhenoLines = PhenoLines[order(PhenoLines$vcfID),]


### markers ~ 200k
load() # load the marker data for the SNPs

MarkersLines = SNPdata


### Relationship matrix for lines
GLines = AGHmatrix::Gmatrix(as.matrix(MarkersLines),
                            missingValue = "NA",
                            thresh.missing = 0.8,
                            maf = 0.05)

### Filtering and ordering
GMat = GLines[rownames(GLines) %in% PhenoLines$vcfID,
              colnames(GLines) %in% PhenoLines$vcfID]

GMat = GMat[order(rownames(GMat)),
            order(colnames(GMat))]


### Imputation and filtering of missing markers
LinesMarkers = getImp(Markers = as.matrix(MarkersLines),
                      MM = 0.8,
                      IM = 0.8,
                      thresh = 0.05)

### Filtering and ordering
MarkLines = LinesMarkers[rownames(LinesMarkers) %in% PhenoLines$vcfID,]
MarkLines = MarkLines[order(rownames(MarkLines)),]

### Check if they are in the same order
stopifnot(all(tolower(PhenoLines$vcfID) == tolower(colnames(GMat))))
stopifnot(all(tolower(PhenoLines$vcfID) == tolower(rownames(MarkLines))))

#######--------------------------------
### Loading datasets for Hybrids
#######--------------------------------

##>>-----  SNPs
load()# load the marker data for the hybrids

##>>-----  BLUEs
PhenoHybrid = read.table("data/HybridBLUEs", h=T, sep="\t")

PhenoHybrid = PhenoHybrid[,-3] # BLUES
PhenoHybrid = PhenoHybrid[order(PhenoHybrid$vcfID),]
colnames(PhenoHybrid) <- c( 'vcfID',     'R3.LN')



## ---------- 1.0 Parents
ParentsID = c( "K120-1-1-1-1-1-2-1",      "UFI9EC",                  "I31-1-1-1-2-1-1"  ,       "wh07166R",
               "I94-1-1-1-1-1-1-1" ,      "wuh09191i" ,              "wh09141R" ,               "I237-1-2-1-1-1-1-1" ,
               "Wh01001"   ,              "wuh12039i" ,              "wh12004",                 "wh10046R" ,
               "wuh12027i" ,              "L25-1-1-1-1-1" ,          "I112-1-1-2-1-1-1-" ,      "I58-1-2-2-1-1-1" ,
               "K115-1-1-1-1-1-1-1-1B-1", "wuh09153i" ,              "wh10006R"  ,              "N38-2-1-1-1-1-1-1",
               "K164-1-1-1-1-1"  ,        "N99-1-1-1-1-1-1-1" ,      "I136-1-1-1-1-1",          "L48-1-1-1-1-1",
               "L108-1-1-1-1-1-1" ,                               "wh08045",                 "M31-1-1-1-1-1-1",
               "I152-1-1-1-1-2-1-1",      "L125-1-1-1" ,             "K127-1-1-1-2-1-1-1-1",    "K43-1-1-1-1-1-1-1"  ,
               "K169-1-1-1-2-2" ,         "wh12005" ,                "wh10173V",                "wh10035R",
               "wh08114"  ,               "wh10057R" ,               "K15-1(X9)",               "K113-2-1-2-1-1"  ,
               "wh14101V" ,               "wh14209N",                "UFI2EDbt1-1",             "wh13041",
               "wh07061b"    ,            "wh11005",                 "wh14203N",                "wh09061"  ,
               "wh11001",                 "Wh14201N" )



# As this is C--Series, one side came from Wisconsin and the other one came from UF, so, lets use this info

ParentsWis <- ParentsID[grepl("^[wW]", ParentsID)]
ParentsFL <- ParentsID[!grepl("^[wW]", ParentsID)]


#######--------------------------------
### 3. Bayesian parameters
#######--------------------------------

nReps = 10
nFolds = 5


#######--------------------------------
### 3. Single trait - CV1 for lines
#######--------------------------------


#--- 1. Building  the ETA for BGLR
ETA_GBLUP = list(list(K = GMat, model='RKHS'))
ETA_BayesB = list(list(X = as.matrix(MarkLines), model='BayesB'))



#--- 2. Output for the estimated values

#Data for the loop
traits = colnames(PhenoLines)[-1]


#--- 2. Output for the estimated values
yHatCV <- list(GBLUPstm = data.frame("R1.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                     "R3.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                     "R6.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                     stringsAsFactors=FALSE),
               BayesBstm = data.frame("R1.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                      "R3.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                      "R6.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                      stringsAsFactors=FALSE))


#--- 3. Model
y = as.matrix(PhenoLines[,-1])

predMod = data.frame()


#Data for the loop
traits = colnames(PhenoLines)[-1]
Mod = c("GBLUPstm", "BayesBstm")


for (j in 1:length(traits)){
  
  #Data for the loop
  Y1 = y[,j]
  
  # Time for the loop
  set.seed(0928761)
  for(Rep in 1:nReps){
    
    folds=sample(1:nFolds,size=length(Y1),replace=T)
    
    for(i in 1:max(folds)){
      tst=which(folds==i)
      yNA=Y1
      yNA[tst]=NA
      
      # model
      fmGblup=BGLR(y=yNA,
                   ETA=ETA_GBLUP,
                   nIter = 10000,
                   burnIn = 1000,
                   verbose = TRUE)
      
      # Predicted values
      yHatCV$GBLUPstm[tst,paste0(traits[j])]=fmGblup$yHat[fmGblup$whichNa]
      
      
      # model
      fmBayes=BGLR(y=yNA,
                   ETA=ETA_BayesB,
                   nIter = 10000,
                   burnIn = 1000,
                   verbose = TRUE)
      
      # Predicted values
      yHatCV$BayesBstm[tst,paste0(traits[j])]=fmBayes$yHat[fmBayes$whichNa]
      
    }
    
    for(model in Mod){
      
      yHatCV_trait <- yHatCV[[model]][, traits[j]]  # Predicted values for the current trait
      head(yHatCV_trait, 20)
      y_trait <- y[, traits[j]]  # Observed values for the current trait
      head(y_trait, 20)
      
      # Calculate metrics for the current trait
      predMod <- rbind(predMod,
                       data.frame(Trait = traits[j],
                                  Model = model,
                                  Repetition = Rep,
                                  Cor = cor(yHatCV_trait[tst], y_trait[tst], use = "complete"),
                                  RMSE = sqrt(mean((yHatCV_trait[tst] - y_trait[tst])^2, na.rm = TRUE))))
    } # Closing models
  } # closing Reps
} # Closing traits

colnames(predMod) <- NULL

write.table(predMod, file = "1.CV1_Lines.txt", row.names = F, quote = F, append = T)


#######--------------------------------
### 4. Multitrait - CV1 for lines
#######--------------------------------


#--- 1. Building  the ETA for GBLUP/BayesB

ETA_GBLUP = list(list(K = GMat, model='RKHS'))
ETA_BayesB = list(list(X = as.matrix(MarkLines), model='SpikeSlab'))


#--- 2. Output for the estimated values
yHatCV <- list(GBLUPmtm = data.frame("R1.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                     "R3.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                     "R6.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                     stringsAsFactors=FALSE),
               BayesBmtm = data.frame("R1.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                      "R3.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                      "R6.LN"=as.numeric(rep(NA, dim(PhenoLines)[1])),
                                      stringsAsFactors=FALSE))


predMod = data.frame()


#Data for the loop
traits = colnames(PhenoLines)[-1]
Mod = c("GBLUPmtm", "BayesBmtm")

#--- 3. Model
y = as.matrix(PhenoLines[,-1])



# Time for the loop
set.seed(88451)
for(Rep in 1:nReps){
  
  folds=sample(1:nFolds,size=dim(y)[1],replace=T)
  
  for(i in 1:max(folds)){
    tst=which(folds==i)
    yNA=y
    yNA[tst,]=NA
    
    # model GBLUP
    fm=Multitrait(y=yNA,
                  ETA=ETA_GBLUP,
                  nIter = 10000,
                  burnIn = 1000,
                  verbose = TRUE)
    
    # Predicted values
    yHatCV$GBLUPmtm[tst,]=fm$ETAHat[tst,]
    
    # model BayesB
    fmB=Multitrait(y=yNA,
                   ETA=ETA_BayesB,
                   nIter = 10000,
                   burnIn = 1000,
                   verbose = TRUE)
    
    # Predicted values
    yHatCV$BayesBmtm[tst,]=fmB$ETAHat[tst,]
    
    
  }
  
  for (trait in traits) {
    for(model in Mod){
      
      yHatCV_trait <- yHatCV[[model]][, trait]  # Predicted values for the current trait
      y_trait <- y[, trait]  # Observed values for the current trait
      
      # Calculate metrics for the current trait
      predMod <- rbind(predMod,
                       data.frame(Trait = trait,
                                  Model = model,
                                  Repetitions = Rep,
                                  Cor = cor(yHatCV_trait, y_trait, use = "complete"),
                                  RMSE = sqrt(mean((yHatCV_trait - y_trait)^2, na.rm = TRUE)) ))
    } # Closing models
  } #Closing traits
  
}


colnames(predMod) <- NULL

write.table(predMod, file = "1.CV1_Lines.txt", row.names = F, quote = F, append = T)



#######--------------------------------
### 4. Single trait for hybrids
#######--------------------------------


#######--------------------------------
### Combining the markers from the hybrids and from the lines and creating a G matrix
#######--------------------------------

##>>----- 1. Phenotypes
PhenoLines$Id = 1
PhenoHybrid$Id = 2

dfFinal = rbind(PhenoLines[,c(1,3, 5)], PhenoHybrid)

##>>----- 2. Marker data
SNPdata_lines = SNPdata[rownames(SNPdata) %in% PhenoLines$vcfID,]
SNPdata_lines = SNPdata_lines[order(rownames(SNPdata_lines)),]

SNPMathybrid = SNPMathybrid[order(rownames(SNPMathybrid)),]


##>>-----  Check if they are in the same order
stopifnot( all(tolower(colnames(SNPMathybrid)) == tolower(colnames(SNPdata_lines))))


MkrComb = rbind(SNPdata_lines, SNPMathybrid)

##>>----- Imputation and filtering of missing markers
MkrComb_filt = getImp(Markers = as.matrix(MkrComb),
                      MM = 0.8,
                      IM = 0.8,
                      thresh = 0.05)

MkrComb_filt_Dom = MkrComb_filt
MkrComb_filt_Dom[MkrComb_filt_Dom == 2] <- 0

##>>----- Generating the matrix
GMatBoth <- Gmatrix(as.matrix(MkrComb),
                    missingValue = "NA",
                    thresh.missing = .8,
                    method = "VanRaden",
                    maf = 0.05)

GMatBoth_Dom <- Gmatrix(as.matrix(MkrComb),
                        missingValue = "NA",
                        thresh.missing = .8,
                        method = "Vitezica",
                        maf = 0.05)

##>>-----  Check if they are in the same order
all(tolower(dfFinal$vcfID) == tolower(colnames(GMatBoth)))


#######--------------------------------
### Scenario i. Using all parents in the matrix (lines that were parents of the hybrids)
#######--------------------------------

##>>----- Traits
traits = 'R3.LN'
predMod = data.frame()

##>>----- Building  the ETA for BGLR

# Incidence for hybrids
dfFinal$Id = factor(dfFinal$Id)
Z_L=model.matrix(~0+Id,data=dfFinal)

# ETA
ETA_GBLUP = list(fixed=list(model='FIXED',X=Z_L),
                 Rand = list(K = GMatBoth, model='RKHS'))

ETA_BayesB = list(fixed=list(model='FIXED',X=Z_L),
                  Rand = list(X = MkrComb_filt, model='BayesB'))

# ETA
ETA_GBLUP_Dom = list(fixed=list(model='FIXED',X=Z_L),
                     A = list(K = GMatBoth, model='RKHS'),
                     D = list(K = GMatBoth_Dom, model = "RKHS"))

ETA_BayesB_Dom = list(fixed=list(model='FIXED',X=Z_L),
                      A = list(X = MkrComb_filt, model='BayesB'),
                      D = list(X = MkrComb_filt_Dom, model='BayesB'))


##>>-----  Missing data for the hybrids
y = dfFinal[,2]
yNA = y
yNA[547:611] <- NA

##>>-----  model GLUP additive
fmGblup=BGLR(y=yNA,
             ETA=ETA_GBLUP,
             nIter = 10000,
             burnIn = 1000,
             verbose = FALSE)

# Predicted values
predMod = data.frame(Model = "GBLUP",
                     Scenario = "All",
                     Cor=cor(fmGblup$yHat[547:611], y[547:611], use = "complete"),
                     RMSE = sqrt(mean((fmGblup$yHat[547:611] - y[547:611])^2, na.rm=TRUE)),
                     Trait = "R3.LN")

colnames(predMod) <- NULL

write.table(predMod, file = "3.outMod/1.CV0_Hybrids.txt", row.names = F, quote = F, append = T)

##>>-----  model GLUP additive
fmGblup_Dom=BGLR(y=yNA,
                 ETA=ETA_GBLUP_Dom,
                 nIter = 10000,
                 burnIn = 1000,
                 verbose = FALSE)

# Predicted values
predMod = data.frame(Model = "GBLUP_Dom",
                     Scenario = "All",
                     Cor=cor(fmGblup_Dom$yHat[547:611], y[547:611], use = "complete"),
                     RMSE = sqrt(mean((fmGblup_Dom$yHat[547:611] - y[547:611])^2, na.rm=TRUE)),
                     Trait = "R3.LN")

colnames(predMod) <- NULL

write.table(predMod, file = "1.CV0_Hybrids.txt", row.names = F, quote = F, append = T)


##>>-----  model Bayes B with additive markers only
fmBayes=BGLR(y=yNA,
             ETA=ETA_BayesB,
             nIter = 10000,
             burnIn = 1000,
             verbose = FALSE)

# Predicted values
predMod = data.frame(Model = "BayesB",
                     Scenario = "All",
                     Cor=cor(fmBayes$yHat[547:611], y[547:611], use = "complete"),
                     RMSE = sqrt(mean((fmBayes$yHat[547:611] - y[547:611])^2, na.rm=TRUE)),
                     Trait = "R3.LN")

colnames(predMod) <- NULL

write.table(predMod, file = "1.CV0_Hybrids.txt", row.names = F, quote = F, append = T)


##>>-----  model Bayes B with additive + dominance coded markers
fmBayes_Dom=BGLR(y=yNA,
                 ETA=ETA_BayesB_Dom,
                 nIter = 10000,
                 burnIn = 1000,
                 verbose = FALSE)

# Predicted values
predMod = data.frame(Model = "BayesB_Dom",
                     Scenario = "All",
                     Cor=cor(fmBayes_Dom$yHat[547:611], y[547:611], use = "complete"),
                     RMSE = sqrt(mean((fmBayes_Dom$yHat[547:611] - y[547:611])^2, na.rm=TRUE)),
                     Trait = "R3.LN")

colnames(predMod) <- NULL

write.table(predMod, file = "1.CV0_Hybrids.txt", row.names = F, quote = F, append = T)








