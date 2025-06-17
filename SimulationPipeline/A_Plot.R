################################################################################
##### Plotting considering different scenarios and se
################################################################################

##>>---------------------------------------------------------------
##> 1. Environment and packages
##>>---------------------------------------------------------------

###>---- Environment
rm(list=ls())

require(plyr)
library('ggplot2')
require(AlphaSimR)
library(RColorBrewer)
library(data.table)



setwd("/orange/mresende/marcopxt/RGS/5.Simulations/2.SweetCorn/") #*** change work D

##>>---------------------------------------------------------------
##> 2. Loop through the repetitions to read the datals
##>>---------------------------------------------------------------

# files
files = list.files("/orange/mresende/marcopxt/RGS/5.Simulations/2.SweetCorn/", pattern = 'Scenario_')

rawData = list()

# Loop
for(i in 1:length(files)){
  
  temp = data.frame(Year=1:35, readRDS(files[i]),
                    stringsAsFactors=FALSE)
  
  rawData[[i]] = temp
}

rawData = rbindlist(rawData) #collapse all Scenarios into one data frame

# Changing to better plot the results
rawData$Year = rawData$Year-15

df3<-rawData



##>>---------------------------------------------------------------
##> 3. Plotting the trait mean
##>>---------------------------------------------------------------

### Summarize the values of reps
temp = ddply(df3,c("Year","scenario"), summarize,
             mean = mean(MeanAT1),
             se = sd(MeanAT1 )/sqrt(length(rep)),
             Trait = "Trait 1")

temp2 = ddply(df3,c("Year","scenario"), summarize,
               mean = mean(MeanAT2),
               se = sd(MeanAT2)/sqrt(length(rep)),
               Trait = "Trait 2")


deltaDf = rbind(temp, temp2)


mean = ggplot(deltaDf,aes(x=Year,y=mean,group=scenario,color=scenario,fill=scenario,
                linetype = scenario))+
  facet_wrap(~Trait,ncol=2) +
  geom_line(size = 1.5)+
  geom_ribbon(aes(x=Year,ymin=mean-se,ymax=mean+se), alpha=0.2, linetype=0, show.legend = F, fill="grey55")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(colour = "grey85", linewidth = 0.3, linetype = "dotted"),
        panel.grid.major = element_line(colour = "grey85", linewidth = 0.3, linetype = "dotted"),
        legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black", family = "Helvetica", face = 'bold'),
        axis.text.y.left = element_text(size = 12, colour = "black", family = "Helvetica", face = 'bold'),
        axis.text.y.right = element_blank(),
        axis.title = element_text(size = 18, face = "bold", family = "Helvetica"),
        axis.line.x = element_line(colour = "black", linewidth = 0.5),
        panel.spacing = unit(1, "lines"),  # adjust this spacing to reduce gap
        axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
        axis.line.y.right = element_blank() ,
        strip.text = element_text(size = 18, face = "bold", family = "Helvetica") )+
    scale_linetype_manual(values = c("solid", "solid",'solid', "solid" ,'solid', 'solid'))+
    scale_color_manual(values = c("#E495A5",'#ACA4E2',"orange", 'grey90'))+
    scale_x_continuous("Year",limits=c(0,20))+
    scale_y_continuous("Genetic mean",expand=c(0,0),limits=c(NA,NA), sec.axis = dup_axis(name = element_blank())) 




##>>---------------------------------------------------------------
##> 4. Plotting the variance
##>>---------------------------------------------------------------
colnames(df3)

### Summarize the values of reps
temp = ddply(df3,c("Year","scenario"), summarize,
             mean = mean(VarAT1),
             se = sd(VarAT1 )/sqrt(length(rep)),
             Trait = "Trait 1")

temp2 = ddply(df3,c("Year","scenario"), summarize,
              mean = mean(VarAT2),
              se = sd(VarAT2)/sqrt(length(rep)),
              Trait = "Trait 2")



deltaDfvar = rbind(temp, temp2)


var = ggplot(deltaDfvar,aes(x=Year,y=mean,group=scenario,color=scenario,fill=scenario,
                          linetype = scenario))+
  facet_wrap(~Trait,ncol=2) +
  geom_line(size = 1.5)+
  geom_ribbon(aes(x=Year,ymin=mean-se,ymax=mean+se), alpha=0.2, linetype=0, show.legend = F, fill="grey55")+
  theme_minimal() +
  theme(panel.grid.minor = element_line(colour = "grey85", linewidth = 0.3, linetype = "dotted"),
        panel.grid.major = element_line(colour = "grey85", linewidth = 0.3, linetype = "dotted"),
        legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black", family = "Helvetica", face = 'bold'),
        axis.text.y.left = element_text(size = 12, colour = "black", family = "Helvetica", face = 'bold'),
        axis.text.y.right = element_blank(),
        axis.title = element_text(size = 18, face = "bold", family = "Helvetica"),
        axis.line.x = element_line(colour = "black", linewidth = 0.5),
        panel.spacing = unit(1, "lines"),  # adjust this spacing to reduce gap
        axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
        axis.line.y.right = element_blank() ,
        strip.text = element_text(size = 18, face = "bold", family = "Helvetica") )+
  scale_linetype_manual(values = c("solid", "solid",'solid', "solid" ,'solid', 'solid'))+
  scale_color_manual(values = c("#E495A5",'#ACA4E2',"orange", 'grey90'))+
  scale_x_continuous("Year",limits=c(0,20))+
  scale_y_continuous("Genetic variance",expand=c(0,0),limits=c(NA,NA), sec.axis = dup_axis(name = element_blank())) 



##>>---------------------------------------------------------------
##> 5. Combine and save
##>>---------------------------------------------------------------



# Plot
jpeg("1.PlotMeanVar.jpeg", width = 18, height = 6, units = 'in', res = 350)
ggpubr::ggarrange(mean, var, # list of plots
                  labels = c('A|', 'B|'), # labels
                  font.label = list(size = 28, 
                                    color = "black", face = "bold", family = NULL),
                  common.legend = T, # COMMON LEGEND
                  legend = "bottom", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  heights = 2,
                  nrow = 1,
                  ncol = 2)  # number of rows

dev.off()



# Plot
tiff("1.PlotMeanVar.tiff", width = 18, height = 6, units = 'in', res = 350)
ggpubr::ggarrange(mean, var, # list of plots
                  labels = c('A|', 'B|'), # labels
                  font.label = list(size = 28, 
                                    color = "black", face = "bold", family = NULL),
                  common.legend = T, # COMMON LEGEND
                  legend = "bottom", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  heights = 2,
                  nrow = 1,
                  ncol = 2)  # number of rows

dev.off()



