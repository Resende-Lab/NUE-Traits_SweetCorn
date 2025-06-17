########################################################################################
##### Models for RGS paper - Finegan et al. 2025
##### Dr. Marco Peixoto
##### 12/10/2024
########################################################################################



#######--------------------------------
### 1. Packages and Environment
#######--------------------------------
rm(list=ls())

library(AGHmatrix)
library(BGLR)
library(SimpleMating)

source("../1.Functions.R")


#######--------------------------------
### 2. Loading datasets for lines
#######--------------------------------


### markers ~ 200k
load("../data/SNPmat200k.RData", verbose = T)

### Phenotypic data for lines
PhenoLines = read.table("../data/LinesBLUEs", h=T)
PhenoLines = PhenoLines[order(PhenoLines$vcfID),]

### Filtering and ordering
SNPdata = SNPdata[order(rownames(SNPdata)),]


# Ordering the lines and predicting all individuals
PhenoNa = data.frame(vcfID = rownames(SNPdata))

PhenoNa = plyr::join(PhenoNa, PhenoLines)
dim(PhenoNa)


#######--------------------------------
### 3. Imputation and filtering of missing markers
#######--------------------------------

### Relationship matrix for lines
relMat = AGHmatrix::Gmatrix(as.matrix(SNPdata),
                            missingValue = "NA",
                            thresh.missing = 0.8,
                            maf = 0.05)


all(tolower(PhenoNa$vcfID) == tolower(rownames(relMat)))


#######--------------------------------
### 4. Predicting individuals
#######--------------------------------

## ---------- 1. Building  the ETA for BGLR
ETAL = list(list(K = relMat, model="RKHS"))

## ---------- 2. Model
y = as.matrix(PhenoNa[,c(2,3)]) # targeting trait R1.LN and R3.LN


fm = BGLR::Multitrait(y=y,
                      ETA=ETAL,
                      nIter = 10000,
                      burnIn = 1000,
                      verbose = FALSE)

## ---------- 3. Outcomes
IndPred = data.frame(GenID = PhenoNa$vcfID,
                     "R3.LN" = fm$ETAHat[,2])



#######--------------------------------
### 5. Prediction of crosses and optimization
#######--------------------------------

## ---------- 1. Crossing plan
Plan = SimpleMating::planCross(TargetPop = PhenoNa$vcfID)


## ---------- 2. Crosses
crossPlanAll = SimpleMating::getMPV(MatePlan = Plan,
                                    Criterion = IndPred,
                                    K = relMat)

# ---------- 3. Optimization
MatePlan = SimpleMating::selectCrosses(data = crossPlanAll,
                                       n.cross = 100,
                                       max.cross = 2,
                                       min.cross = 1,
                                       culling.pairwise.k = -0.02,
                                       TargetPool = NULL)


# # --------- 4. Saving

MatePlan[[2]]

#######--------------------------------
### 6. Prediction of crosses for hybrids created in the target population
#######--------------------------------

## ---------- 1. Load the pedigree
Ped = read.table("../data/PedData", h=T)
PlanHybrid = Ped[,c(2,3)]
allPar = unique(c(PlanHybrid$Parent1, PlanHybrid$Parent2))


## ---------- 2. Crosses
crossPlan = SimpleMating::getMPV(MatePlan = PlanHybrid,
                                 Criterion = IndPred,
                                 K = relMat)


## ---------- 3. Plotting
FinalHybrids = read.table("../data/HybridBLUEs", sep = ",", h=T)[,-c(2, 4,5)]
head(FinalHybrids)

## ----------  4. Loading pedigree

Ped = read.table("../data/PedData", h=T)
Ped$CrossID = paste0(Ped$Parent1, "_", Ped$Parent2)
colnames(Ped)[1] <- 'vcfID'

# Combining
OriData = plyr::join(FinalHybrids, Ped)

head(OriData)

# 3. Combination with predicted data
CrossPred = plyr::join(OriData, crossPlan)





###>>>>>---------------
###> Plotting Figure 5 from the paper
###>

df = crossPlanAll
df$ID = paste0(df$Parent1, "_", df$Parent2)
selCrosses <- MatePlan[[2]]
selCrosses$ID <- paste0(selCrosses$Parent1, "_", selCrosses$Parent2)
df$Sel <- ifelse(df$ID %in% selCrosses$ID, "S", 
                 "NS")
dfSel <- df[order(df$Sel, decreasing = FALSE), ]
relMax <- max(dfSel$K)
relMin <- min(dfSel$K)
CritMax <- max(dfSel$Y)
CritMin <- min(dfSel$Y)
K = dfSel$K
Y = dfSel$Y
Sel = dfSel$Sel


a = ggplot(dfSel, aes(K, Y)) + 
  geom_point(aes(colour = factor(Sel)), size = 3) + 
  scale_color_manual(values = c("#ACA4E2", "wheat")) + 
  geom_vline(xintercept = -0.02, linetype = "dashed", color = "grey30", linewidth = 0.8) + 
  
  theme(panel.grid = element_line(color = "grey90", size = 0.35, linetype = 2),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = c(0.85, 0.9), 
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black", size = 0.1, linetype = "solid"),
        axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid")) +
  
  scale_x_continuous("Relatedness", limits = c(relMin - 0.05, relMax + 0.05)) + 
  scale_y_continuous("Mid-Parental values", limits = c(CritMin - 0.25, CritMax + 0.25))



# Create the plot




b = ggplot(CrossPred, aes(x = Hy_R3.LN_BLUE, y = Y)) +
  geom_point(size = 3, alpha = 1, color =  "#ACA4E2") + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="grey30")+
  labs(
    title = " ",
    x = "BLUEs",
    y = " "
  ) +
  theme(
    panel.grid = element_line(color = "grey90", size = 0.35, linetype = 2),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid"),
  ) +
  ylim(2.75, 3.4) +
  xlim(3.25,4.75)







jpeg("1.CrossPred3.jpeg", width = 18, height = 9, units = "cm", res = 300)
ggpubr::ggarrange(a,b,  # list of plots
                  labels = c('A|', 'B|'), # labels
                  font.label = list(size = 16, 
                                    color = "black", family = "sans"),
                  label.x = 0,           # Horizontal position of labels
                  label.y = 1,            # Vertical position of labels
                  common.legend = F, # COMMON LEGEND
                  #legend = "bottom", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  heights = 1,
                  nrow = 1,
                  ncol = 2)  
dev.off()




tiff("1.CrossPred3.tiff", width = 18, height = 9, units = "cm", res = 300)
ggpubr::ggarrange(a,b,  # list of plots
                  labels = c('A|', 'B|'), # labels
                  font.label = list(size = 16, 
                                    color = "black", family = "sans"),
                  label.x = 0,           # Horizontal position of labels
                  label.y = 1,            # Vertical position of labels
                  common.legend = F, # COMMON LEGEND
                  #legend = "bottom", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  heights = 1,
                  nrow = 1,
                  ncol = 2)  
dev.off()



