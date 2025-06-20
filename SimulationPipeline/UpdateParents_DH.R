###---------
# Choosing the best lines to update the crossing blocks
###---------

PublicInbredYT2Phen@pheno = as.matrix(data.frame(Trait1 = aggregate(as.matrix(calcGCA(PublicHybridYT2)$GCAm[,2])~as.matrix(calcGCA(PublicHybridYT2)$GCAm[,1]), FUN=mean)[,2],
                                                 Trait2 = aggregate(as.matrix(calcGCA(PublicHybridYT2)$GCAm[,3])~as.matrix(calcGCA(PublicHybridYT2)$GCAm[,1]), FUN=mean)[,2]))

new_parents = selectInd(PublicInbredYT2Phen, nInd= 37)
old_parents = selectInd(Parents_public_update, nInd= 13) 

Parents_public_update = c(new_parents, old_parents)

