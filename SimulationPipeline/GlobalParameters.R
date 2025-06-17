###########################################################
#### Parameters of the base population
#########################################################


burninYears = 15
futureYears = 20

nParents = 200
nParents_public = 50
nParents_private = 100
nQtl = 300
nSnp = 500

MeanG = c(4.04,3.04)
VarG =  c(0.5, 0.5) 
  
ddVar = c(0.3, 0.3) 
ddMean = c(0.92 , 0.92)

# Dhs
pub_DH = 25

VarE=c(0.9, 0.7)
VarGE=c(0.9*3, 0.7*3)


# Additive effect
corA = matrix(c(1,	0.48,
                0.48,	1), nc=2, byrow=T)


#Dominance degree
corDD = matrix(c(1,	0,
                 0,	1), nc=2, byrow=T)

#GxE effect
corGE = matrix(c(1,	0,
                 0,	1), nc=2, byrow=T)




