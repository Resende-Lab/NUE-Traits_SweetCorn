#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data

########-----------------  Year 4  -------------########
# fall - off-season environment
PublicHybridYT3 = hybridCross(Tester3_private, PublicInbredYT3) # Fall

# Winter - target environment
PublicHybridYT3 = setPheno(PublicHybridYT3,varE=VarE, reps=repTC3, p=pgy) # Winter
PublicInbredYT3@pheno = as.matrix(data.frame(Trait1 = calcGCA(PublicHybridYT3)$GCAm[,2],
                                             Trait2 = calcGCA(PublicHybridYT3)$GCAm[,3]))

PublicHybridYT4 = selectInd(PublicHybridYT3, 2, use="pheno") # Winter

HybridMean = PublicHybridYT4
DHTrain3 = PublicInbredYT3

########-----------------  Year 3  -------------########
# Fall - off-season environment
PublicHybridYT2 = hybridCross(Tester2_private, PublicInbredYT2) #Fall

# Winter - target environment
PublicHybridYT2 = setPheno(PublicHybridYT2, varE=VarE, reps=repTC2, p=pgy) # Winter
PublicInbredYT2@pheno = as.matrix(data.frame(Trait1 = calcGCA(PublicHybridYT2)$GCAm[,2],
                                             Trait2 = calcGCA(PublicHybridYT2)$GCAm[,3]))

PublicInbredYT2@ebv = runGBLUP(trainPop = trainPop, targetPop = PublicInbredYT2)
PublicInbredYT3 = selectInd(PublicInbredYT2, 20, use="ebv") # Winter

DHTrain2 = PublicInbredYT2

########-----------------  Year 2  -------------########
# Fall - off-season environment
PublicHybridYT1 = hybridCross(Tester1_private, PublicInbredYT1Sel) # Fall

# Winter - target environment
PublicHybridYT1 = setPheno(PublicHybridYT1, varE=VarE, reps=repTC1, p=pgy) # Winter
PublicInbredYT1Sel@pheno = as.matrix(data.frame(Trait1 = calcGCA(PublicHybridYT1)$GCAm[,2],
                                                Trait2 = calcGCA(PublicHybridYT1)$GCAm[,3]))

PublicInbredYT1Sel@ebv = runGBLUP(trainPop = trainPop, targetPop = PublicInbredYT1Sel)
PublicInbredYT2 = selectInd(PublicInbredYT1Sel, 150, use = "ebv") # Winter

DHTrain1 = PublicInbredYT1Sel


########-----------------  Year 1  -------------########
source("MatingAllocation.R")

# Spring
F1_public = makeCross(Parents2OCS, crossPlan = OCSPlan)


# Fall and Winter = Make DH
PublicDH = makeDH(F1_public, nDH = pub_DH, keepParents = FALSE) 
PublicInbredYT1 = PublicDH

#--- GS model
PublicInbredYT1@ebv = runGBLUP(trainPop = trainPop, targetPop = PublicInbredYT1)
PublicInbredYT1Sel = selectInd(PublicInbredYT1, nInd=800, use = 'ebv') # Fall

