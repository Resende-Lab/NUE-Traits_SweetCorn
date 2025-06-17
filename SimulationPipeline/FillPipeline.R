######-------------------------------------------------------------
###### Fillpipeline
######-------------------------------------------------------------

library(AlphaSimR)

###>>>>----------------- p-values for GxY effects
b.f=burninYears+futureYears

# p-values for GxY effects
Pgy = 0.5 + rnorm(b.f,0,0.03)


# p-values for GxYxE effects
Pgye1 = 0.9 + rnorm(b.f,0,0.03)
Pgye2 = 0.1 + rnorm(b.f,0,0.03)
Pgye = c(Pgye1,Pgye2)

Pgye = sample(Pgye,b.f, replace=F)

#From normal

P = runif(b.f)

###>>>------------- Running the Fillpipeline
for(year in 1:7){
  cat("FillPipeline year:",year,"of 7\n")
  
  #Year 1
  ########-----------------  Year 1  -------------########
  # Spring
  F1_public = randCross(Parents_public_update, nCrosses_public, ignoreSexes = TRUE)
  
  # Fall and Winter = Make DH
  PublicDH = makeDH(F1_public, nDH = pub_DH, keepParents = FALSE) 
  PublicInbredYT1 = PublicDH
  
  
  ####>>>>----Private program
  F1_private = randCross(Parents_private, nCrosses_private, ignoreSexes = TRUE)
  F2_private = self(F1_private, nProgeny = 10, parents = NULL, keepParents = FALSE)
  F2_private = setPheno(F2_private,reps=1)
  F2_private.sel = selectWithinFam(F2_private, nInd = 6, selectTop = TRUE)
  #Year 2
  if(year<7){
    pgy = Pgy[year]                                                      
    pgye = Pgye[year]
    p=P[year]
    
    ########-----------------  Year 2  -------------########
    # Fall
    PublicInbredYT1 = setPheno(PublicInbredYT1,reps=1,p=pgye, varE = VarE) # Fall
    PublicInbredYT1Sel = selectInd(PublicInbredYT1, nInd=800) # Fall
    PublicHybridYT1 = hybridCross(Tester1_private, PublicInbredYT1Sel) # Fall
    
    # Winter - target environment
    PublicHybridYT1 = setPheno(PublicHybridYT1,varE=VarE, reps=repTC1, p=pgy) # Winter
    PublicInbredYT1Sel@pheno = as.matrix(calcGCA(PublicHybridYT1)$GCAm[,2]) # Winter
    PublicInbredYT2 = selectInd(selectWithinFam(PublicInbredYT1Sel, 4, use = "pheno"),  # Winter
                                200, use = "pheno")
    
    ####>>>>----Private program
    F3_private = self(F2_private.sel, nProgeny = 8, parents = NULL, keepParents = FALSE)
    F3_private = setPheno(F3_private,reps=1, p=p)
    F3_private.sel = selectWithinFam(F3_private, nInd = 2, selectTop = TRUE)
  }
  #Year 3
  if(year<6){
    pgy = Pgy[year]                                                      
    pgye = Pgye[year]
    p=P[year]
    
    ########-----------------  Year 3  -------------########
    # Fall
    PublicHybridYT2 = hybridCross(Tester2_private, PublicInbredYT2) #Fall
    
    # Winter - target environment
    PublicHybridYT2 = setPheno(PublicHybridYT2,varE=VarE, reps=repTC2,p=pgy) # Winter
    PublicInbredYT2@pheno = as.matrix(calcGCA(PublicHybridYT2)$GCAm[,2]) # Winter
    PublicInbredYT3 = selectInd(PublicInbredYT2, 100, use="pheno") # Winter
    PublicInbredYT2Phen = PublicInbredYT2 # Parental Selection
    
    ####>>>>----Private program
    F4_private = self(F3_private.sel, nProgeny = 6, parents = NULL, keepParents = FALSE)
    F4_private = setPheno(F4_private,reps=2, p=p)
    F4_private.sel = selectWithinFam(F4_private, nInd = 2, selectTop = TRUE)
    F4_private.sel = selectInd(F4_private.sel, nInd = length(F4_private.sel@id) * 0.2, selectTop = TRUE)
  }
  #Year 4
  if(year<5){
    pgy = Pgy[year]                                                      
    pgye = Pgye[year]
    p=P[year]
    
    ########-----------------  Year 4  -------------########
    # fall
    PublicHybridYT3 = hybridCross(Tester3_private,PublicInbredYT3) # Fall
    
    # Winter - target environment
    PublicHybridYT3 = setPheno(PublicHybridYT3,varE=VarE, reps=repTC3, p=pgy) # Winter
    PublicInbredYT3@pheno = as.matrix(calcGCA(PublicHybridYT3)$GCAm[,2]) # Winter
    PublicHybridYT4 = selectInd(PublicHybridYT3, 10, use="pheno") # Winter
    
    
    ####>>>>----Private program
    F5_private = self(F4_private.sel, nProgeny = 6, parents = NULL, keepParents = FALSE)
    F5_private = setPheno(F5_private,reps=2, p=p)
    F5_private.sel = selectWithinFam(F5_private, nInd = 1, selectTop = TRUE)
    F5_private.sel = selectInd(F5_private.sel, nInd = length(F5_private.sel@id) * 0.15, selectTop = TRUE)
  }
  #Year 5
  if(year<4){
    pgy = Pgy[year]                                                      
    pgye = Pgye[year]
    p=P[year]
    
    ########-----------------  Year 5  -------------########
    PublicHybridYT4 = setPheno(PublicHybridYT4,varE=VarE,reps=25,p=pgy)
    PublicHybridYT5 = selectInd(PublicHybridYT4, 5, use="pheno")
    
    
    ####>>>>----Private program
    F6_private = self(F5_private.sel, nProgeny = 6, parents = NULL, keepParents = FALSE)
    F6_private = setPheno(F6_private,reps=2, p=p)
    F6_private.sel = selectWithinFam(F6_private, nInd = 1, selectTop = TRUE)
    F6_private.sel = selectInd(F6_private.sel, nInd = length(F6_private.sel@id) * 0.5, selectTop = TRUE)
  }
  #Year 6
  if(year<3){
    pgy = Pgy[year]                                                      
    pgye = Pgye[year]
    p=P[year]
    
    ########-----------------  Year 6  -------------########
    PublicHybridYT5 = setPheno(PublicHybridYT5,varE=VarE,reps=30,p=pgy)
    Release_HybridYT5 = selectInd(PublicHybridYT5, 2, use="pheno")
    
    
    ####>>>>----Private program
    F7_private = self(F6_private.sel, nProgeny = 6, parents = NULL, keepParents = FALSE)
    F7_private = setPheno(F7_private,reps=2, p=p)
    F7_private.sel = selectWithinFam(F7_private, nInd = 1, selectTop = TRUE)
    F7_private.sel = selectInd(F7_private.sel, nInd = 15, selectTop = TRUE)
  }
  #Year 7
  if(year<2){
    pgy = Pgy[year]                                                      
    pgye = Pgye[year]
    p=P[year]
    ####>>>>----Private program
    F8_private = self(F7_private.sel, nProgeny = 6, parents = NULL, keepParents = FALSE)
    F8_private = setPheno(F8_private,reps=2, p=p)
    F8_private.sel = selectWithinFam(F8_private, nInd = 1, selectTop = TRUE)
    F8_private.sel = selectInd(F8_private.sel, nInd = 5, selectTop = TRUE)
  }
}