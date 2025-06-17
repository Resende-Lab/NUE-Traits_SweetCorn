####################################################################################
###################### AlphaSimR code for sweet corn breeding
###################### This is the first script of the simulation
####################################################################################

##>>>----------- Setting the environment
rm(list = ls())

require('AlphaSimR')

setwd("/orange/mresende/marcopxt/RGS/5.Simulations/2.SweetCorn/")

options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
rep <- as.numeric(args[1])

#rep=1
##>>>----------- Setting the scenarios
source("GlobalParameters.R")

# Initialize variables for results
Accuracy =
MeanA_pop =
VarA_pop =
HBMean =
inbFalconer =
nPar =
matrix(NA, 35, 2)

###------------------------------------------------------------------------
###                Setting up the base population
###------------------------------------------------------------------------

###>>>>----------------- Creating parents and Starting the pipeline
# Create initial parents and set testers and hybrid parents
source("CreatParents.R")

# Fill breeding pipeline with unique individuals from initial parents
source("FillPipeline.R")

###>>>>----------------- p-values for Genotype-by-environment effects
b.f=burninYears+futureYears

# p-values for GxY effects
Pgy = 0.5 + rnorm(b.f,0,0.03)


# p-values for GxYxE effects
Pgye1 = 0.9 + rnorm(b.f,0,0.03)
Pgye2 = 0.1 + rnorm(b.f,0,0.03)
Pgye = c(Pgye1,Pgye2)

Pgye = sample(Pgye,b.f, replace=F)

startTrainPop = 12

# Weigths for private program
b = c(1,1)

###################################################################################
############# Burn-in Phase
###################################################################################

for(year in 1:burninYears){
  cat("Working on Year:",year,"\n")
  pgy = Pgy[year]
  pgye = Pgye[year]

  #Private program
  source("UpdateParents_private.R") #Pick new parents and testers
  source("Advance_cycle_private.R") #Advances private yield trials by a year
  source("UpdateTesters.R") #Pick new testers

  # Public
  source("UpdateParents_DH.R") #Pick new parents
  source("Advance_cycle.DH.R") #Advances yield trials by a year
  source("UpdateResults.BS.R") #Track summary data
  source("Writing_records.R")


}

#Save burn-in to load later use
save.image(paste0("BURNIN_",rep,".RData"))


##########################################################################################
############# Scenario 1 - DH scenario
##########################################################################################
# 3.0 Loading the scenarios
load(paste0("BURNIN_",rep,".RData"))


# 3.1 Looping through the years
cat("Working on Scenario 2\n")
for(year in (burninYears+1):(burninYears+futureYears)){
  cat("Working on Year:",year,"\n")
  pgy = Pgy[year]
  pgye = Pgye[year]

  #Private program
  source("UpdateParents_private.R") #Pick new parents and testers
  source("Advance_cycle_private.R") #Advances private yield trials by a year
  source("UpdateTesters.R") #Pick new testers

  # Publica program
  source("UpdateParents_DH.R") #Pick new parents
  source("Advance_cycle.DH.R") #Advances yield trials by a year
  source("UpdateResults.BS.R") #Track summary data


 }

# 3.2 Recording results
output = data.frame(rep=rep(rep, 35),
                     scenario=rep("Conv", 35),
                     AccT1 = Accuracy[,1],
                     AccT2 = Accuracy[,2],
                     MeanAT1 = MeanA_pop[,1],
                     MeanAT2 = MeanA_pop[,2],
                     VarAT1 = VarA_pop[,1],
                     VarAT2 = VarA_pop[,2],
                     HybridMT1 = HBMean[,1],
                     HybridMT2 = HBMean[,2],
                     inFT1 = inbFalconer[,1],
                     inFT2 = inbFalconer[,2],
                     nParT1 = nPar[,1],
                     nParT2 = nPar[,2],
                     stringsAsFactors=FALSE)

# 3.3 Saving the results
saveRDS(output, paste0("Scenario_PS",rep,".rds"))


##########################################################################################
############# Scenario 2 - DHGS scenario
##########################################################################################
# 4.0 Loading the scenarios
load(paste0("BURNIN_",rep,".RData"))

# 4.1 Looping through the years
cat("Working on Scenario 2\n")
for(year in (burninYears+1):(burninYears+futureYears)){
  cat("Working on Year:",year,"\n")
  pgy = Pgy[year]
  pgye = Pgye[year]

  #Private program
  source("UpdateParents_private.R") #Pick new parents and testers
  source("Advance_cycle_private.R") #Advances private yield trials by a year
  source("UpdateTesters.R") #Pick new testers

  source("GenomicModels.R")


  if(year==(burninYears+1)){
    #--- GS model

    PublicInbredYT1@ebv = runGBLUP(trainPop = trainPop, targetPop = PublicInbredYT1)
    PublicInbredYT1Sel = selectInd(PublicInbredYT1, nInd=800, use = 'ebv') # Fall

    PublicInbredYT2 = selectInd(PublicInbredYT1Sel, 150, use = "ebv") # Winter

    PublicInbredYT3 = selectInd(PublicInbredYT2, 20, use="ebv") # Winter

    
    Parents_public_update@ebv = runGBLUP(trainPop = trainPop, targetPop = Parents_public_update)
    }

  source("UpdateParents_GS.R") #Pick new parents
  source("Advance_cycle.DHGS.R") #Advances yield trials by a year
  source("UpdateResults.GS.R") #Track summary data
  source("Writing_records.R")

  cat(table(trainPop@fixEff), '\n')
  cat(names(table(trainPop@fixEff)), '\n')

}

# 4.2 Recording results
output2 = data.frame(rep=rep(rep, 35),
                     scenario=rep("GS", 35),
                     AccT1 = Accuracy[,1],
                     AccT2 = Accuracy[,2],
                     MeanAT1 = MeanA_pop[,1],
                     MeanAT2 = MeanA_pop[,2],
                     VarAT1 = VarA_pop[,1],
                     VarAT2 = VarA_pop[,2],
                     HybridMT1 = HBMean[,1],
                     HybridMT2 = HBMean[,2],
                     inFT1 = inbFalconer[,1],
                     inFT2 = inbFalconer[,2],
                     nParT1 = nPar[,1],
                     nParT2 = nPar[,2],
                     stringsAsFactors=FALSE)

# 4.3 Saving the results
saveRDS(output2, paste0("Scenario_GS",rep,".rds"))


##########################################################################################
############# Scenario 3 - OCS scenario with MPV
##########################################################################################
# 5.0 Loading the scenarios
load(paste0("BURNIN_",rep,".RData"))

# 5.1 Looping through the years
cat("Working on Scenario 3\n")
for(year in (burninYears+1):(burninYears+futureYears)){
  cat("Working on Year:",year,"\n")
  pgy = Pgy[year]
  pgye = Pgye[year]

  #Private program
  source("UpdateParents_private.R") #Pick new parents and testers
  source("Advance_cycle_private.R") #Advances private yield trials by a year
  source("UpdateTesters.R") #Pick new testers
  source("GenomicModels.R") # Genomic model
 
    if(year==(burninYears+1)){
      #--- GS model

      PublicInbredYT1@ebv = runGBLUP(trainPop = trainPop, targetPop = PublicInbredYT1)
      PublicInbredYT1Sel = selectInd(PublicInbredYT1, nInd=800, use = 'ebv') # Fall

      PublicInbredYT2 = selectInd(PublicInbredYT1Sel, 150, use = "ebv") # Winter

      PublicInbredYT3 = selectInd(PublicInbredYT2, 20, use="ebv") # Winter


      Parents_public_update@ebv = runGBLUP(trainPop = trainPop, targetPop = Parents_public_update)
      }
  

  source("UpdateParents_OCS.R") #Pick new parents
  source("Advance_cycle_OCS.R") #Advances yield trials by a year
  source("UpdateResults_Mate.R") #Track summary data
  source("Writing_records.R")

  cat(table(trainPop@fixEff), '\n')
  cat(names(table(trainPop@fixEff)), '\n')

}

# 5.2 Recording results
output3 = data.frame(rep=rep(rep, 35),
                    scenario=rep("OCS", 35),
                    AccT1 = Accuracy[,1],
                    AccT2 = Accuracy[,2],
                    MeanAT1 = MeanA_pop[,1],
                    MeanAT2 = MeanA_pop[,2],
                    VarAT1 = VarA_pop[,1],
                    VarAT2 = VarA_pop[,2],
                    HybridMT1 = HBMean[,1],
                    HybridMT2 = HBMean[,2],
                    inFT1 = inbFalconer[,1],
                    inFT2 = inbFalconer[,2],
                    nParT1 = nPar[,1],
                    nParT2 = nPar[,2],
                    stringsAsFactors=FALSE)

# 5.3 Saving the results
saveRDS(output3, paste0("Scenario_OCSrep",rep,".rds"))



# #####------------------------------ the end ------------------#################





