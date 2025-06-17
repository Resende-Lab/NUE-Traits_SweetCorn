#####-----------------------------------------------
#Track population, hybrid and parents performances
#####-----------------------------------------------

###>>>----- Population
popname = Parents_public_update
genPa = genParam(popname)

###>>>-------- 1. Paramenters
Accuracy[year,]  = diag(cor(gv(popname), pheno(popname), use = "pairwise.complete.obs"))
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




