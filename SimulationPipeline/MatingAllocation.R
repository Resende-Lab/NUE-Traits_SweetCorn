################################################################
###########----------- Mating program -------------
################################################################
require("SimpleMating")
require("AGHmatrix")

#------------------------------------------------------------------------------
####----------- 1. EBVs
#------------------------------------------------------------------------------
# 1.1 Access SNP genotype data
Markers = pullSnpGeno(Parents2OCS)

# 1.2. Ebvs for MPV

Crit <- data.frame(Id = Parents2OCS@id,
                   Criterion = Parents2OCS@ebv)

rownames(Crit) <- NULL


#------------------------------------------------------------------------------
####----------- 2. Relationship matrix
#------------------------------------------------------------------------------
# 2.1 Potential parents - SNPs
relMat = Gmatrix(Markers)

#------------------------------------------------------------------------------
####----------- 3. mid-parental values
#------------------------------------------------------------------------------
# 3.1 Plan
CrossPlan = planCross(TargetPop = Parents2OCS@id,
                 MateDesign = 'half')

# 3.2 Single trait mid parental value
MTM_mpv <- getMPV(MatePlan = CrossPlan,
                 Criterion = Crit,
                 K = relMat)


#------------------------------------------------------------------------------
####----------- 4. Criterion
#------------------------------------------------------------------------------
# 4.1 Crosses selected
maxGainPlan = selectCrosses(data = MTM_mpv, 
                            n.cross = nCrosses_public,
                            max.cross = 2,
                            min.cross = 1,
                            culling.pairwise.k = -0.2,
                            max.cross.to.search = 1e+10) 


# 4.2 Mating plan
opt.planM1 = maxGainPlan[[2]]
OCSPlan = as.matrix(opt.planM1[,c(1,2)])


