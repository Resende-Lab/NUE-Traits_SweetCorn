#Track performance 

#####-----------------------------------------------
#Track population, hybrid and parents performances
#####-----------------------------------------------

###>>>----- Population

popname = Parents2OCS[Parents2OCS@id %in% unique(c(OCSPlan[,1],OCSPlan[,2]))]

genPa = genParam(popname)

###>>>-------- 1. Paramenters
Accuracy[year,]  = diag(cor(gv(popname), ebv(popname), use = "pairwise.complete.obs"))
MeanA_pop[year,] = apply(genPa$gv_a, 2, mean)
VarA_pop[year,]  = diag(genPa$varA)
HBMean[year,] = meanG(HybridMean)
nPar[year,] = popname@nInd

#####>>>>----- 3. Based on Falconer equation

# Inbreeding falconer
inbCoef = function(Pop){
  
  Markers = pullSnpGeno(Pop)
  p = colMeans(Markers)/2
  
  f = sum(p^2)/ncol(Markers)
  
  return(f)
}

inbFalconer[year,] = inbCoef(popname)

