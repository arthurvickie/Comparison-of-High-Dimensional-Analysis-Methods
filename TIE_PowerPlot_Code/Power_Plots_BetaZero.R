##################################################################
## Make Plots
library(ggplot2)
library(cowplot)

#inputs when score is true:
inputs = c("NAT2", "1000", 20, "small", 25, FALSE)
inputs = c("NAT2", "1000", 10, "medium", 25, TRUE)
inputs = c("NAT2", "1000", 10, "large", 5, FALSE)

inputs = c("NAT2", "500", 20, "small", 25, FALSE)
inputs = c("NAT2", "500", 20, "medium", 25, TRUE)
inputs = c("NAT2", "500", 10, "large", 5, FALSE)

# inputs = c("CHI3L2", "1000", 20, "small", 25, FALSE)
# inputs = c("CHI3L2", "1000", 10, "medium", 25, TRUE)
# inputs = c("CHI3L2", "1000", 10, "large", 5, FALSE)

inputs = c("ASAH1", "1000", 20, "small", 15, FALSE)
# inputs = c("ASAH1", "1000", 10, "medium", 25, TRUE)
# inputs = c("ASAH1", "1000", 10, "large", 5, FALSE)

gene = inputs[1]
ss = inputs[2]
prev = as.numeric(inputs[3])
orSize = inputs[4]
percentAssoc = as.numeric(inputs[5])
lowLD = as.logical(inputs[6])

if(lowLD == TRUE){
  ld = "LowLD"
} else {
  ld = "HighLD"
}

scores = c("AMS", "BinMM", "IBS", "Incomp")

######## Binary Outcome ################################################################################################################
# ### LHT ###############################################################################################################################
# setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_LHT_Power_BetaZero"))
# 
# allDataLHT = list()
# 
# for(jj in 1:length(scores)){
#   LHTFile = read.csv(paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_LHT__BetaZero.csv"))
#   
#   #pull the row for the true score
#   LHTFile.trueScore = LHTFile[jj,]
#   
#   #remove column 1
#   LHTFile.trueScore = LHTFile.trueScore[,2:ncol(LHTFile.trueScore)]
#   
#   #need to stack by test
#   library(ramify)
#   LHTFile.trueScore.mat = resize(LHTFile.trueScore, nrow = 4, ncol = 3, byrow = TRUE)
#   LHTFile.trueScore.dat = as.data.frame(LHTFile.trueScore.mat)
#   
#   #Add prevalence
#   LHTFile.trueScore.dat$Prev = prev
#   
#   #add score type
#   LHTFile.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
#   
#   #add OR 
#   LHTFile.trueScore.dat$OR = orSize
#   
#   #LRT values
#   LHTFile.trueScore.dat.LRT = LHTFile.trueScore.dat[,c(1,4:6)]
#   #Wald values
#   LHTFile.trueScore.dat.Wald = LHTFile.trueScore.dat[,c(2,4:6)]
#   #Score values
#   LHTFile.trueScore.dat.Score = LHTFile.trueScore.dat[,c(3:6)]
#   
#   #add method
#   LHTFile.trueScore.dat.LRT$method = "LHT \n LRT"
#   LHTFile.trueScore.dat.Wald$method = "LHT \n Wald"
#   LHTFile.trueScore.dat.Score$method = "LHT \n Score"
#   
#   #rename Score column to Binary Y
#   colnames(LHTFile.trueScore.dat.LRT)[1] = "Binary Y" 
#   colnames(LHTFile.trueScore.dat.Wald)[1] = "Binary Y" 
#   colnames(LHTFile.trueScore.dat.Score)[1] = "Binary Y"
#   
#   #stack everything
#   LHTFile.trueScore.all = rbind(LHTFile.trueScore.dat.LRT, LHTFile.trueScore.dat.Wald, LHTFile.trueScore.dat.Score)
#   
#   #add score type
#   LHTFile.trueScore.all$True_Score = scores[jj]
#   #add ld
#   LHTFile.trueScore.all$ld = ld
#   #add % assoc
#   LHTFile.trueScore.all$PercentAssoc = percentAssoc
#   
#   #all methods give same results, so just use Score test
#   allDataLHT[[jj]] = LHTFile.trueScore.all
# }
# 
# allDataLHT.mat = do.call(rbind, allDataLHT)
# allDataLHT.mat.subset = allDataLHT.mat[allDataLHT.mat$method != "LHT \n LRT",]
### GLM ###############################################################################################################################
# setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_GLM_Power_BetaZero"))
# 
# allDataGLM = list()
# 
# for(jj in 1:length(scores)){
#   File = read.csv(paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_GLM_BetaZero.csv"))
#   
#   #pull the row for the true score
#   File.trueScore = File[File$X == scores[jj],]
#   
#   #remove column 1
#   File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
#   
#   #need to stack by test
#   library(ramify)
#   File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
#   File.trueScore.dat = as.data.frame(File.trueScore.mat)
#   
#   #Add prevalence
#   File.trueScore.dat$Prev = prev
#   
#   #add score type
#   File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
#   
#   #add OR 
#   File.trueScore.dat$OR = orSize
#   
#   #add method
#   File.trueScore.dat$method = "GLM"
#   
#   #rename Score column to Binary Y
#   colnames(File.trueScore.dat)[1] = "Binary Y" 
# 
#   #add score type
#   File.trueScore.dat$True_Score = scores[jj]
#   #add ld
#   File.trueScore.dat$ld = ld
#   #add % assoc
#   File.trueScore.dat$PercentAssoc = percentAssoc
#   
#   #all methods give same results, so just use Score test
#   allDataGLM[[jj]] = File.trueScore.dat
# }
# 
# allDataGLM.mat = do.call(rbind, allDataGLM)
### Desparse Lasso ###############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_DesparseLasso_Power_BetaZero"))

allDataDL = list()

for(jj in 1:length(scores)){
  File = read.csv(paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DL_BetaZero.csv"))
  
  #pull the row for the true score
  File.trueScore = File[jj,]
  
  #remove column 1
  File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
  
  #need to stack by test
  library(ramify)
  File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
  File.trueScore.dat = as.data.frame(File.trueScore.mat)
  
  #Add prevalence
  File.trueScore.dat$Prev = prev
  
  #add score type
  File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
  
  #add OR 
  File.trueScore.dat$OR = orSize
  
  #add method
  File.trueScore.dat$method = "Desparse \n Lasso"
  
  #rename Score column to Binary Y
  colnames(File.trueScore.dat)[1] = "Binary Y" 
  
  #add score type
  File.trueScore.dat$True_Score = scores[jj]
  #add ld
  File.trueScore.dat$ld = ld
  #add % assoc
  File.trueScore.dat$PercentAssoc = percentAssoc
  
  #all methods give same results, so just use Score test
  allDataDL[[jj]] = File.trueScore.dat
}

allDataDL.mat = do.call(rbind, allDataDL)
### AISPU ###############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_AISPU_Power_BetaZero"))

allDataAISPU = list()

for(jj in 1:length(scores)){
  File = read.csv(paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_AISPU_BetaZero.csv"))
  
  #pull the row for the true score
  File.trueScore = File[jj,]
  
  #remove column 1
  File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
  
  #need to stack by test
  library(ramify)
  File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
  File.trueScore.dat = as.data.frame(File.trueScore.mat)
  
  #Add prevalence
  File.trueScore.dat$Prev = prev
  
  #add score type
  File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
  
  #add OR 
  File.trueScore.dat$OR = orSize
  
  #add method
  File.trueScore.dat$method = "aiSPU"
  
  #rename Score column to Binary Y
  colnames(File.trueScore.dat)[1] = "Binary Y" 
  
  #add score type
  File.trueScore.dat$True_Score = scores[jj]
  #add ld
  File.trueScore.dat$ld = ld
  #add % assoc
  File.trueScore.dat$PercentAssoc = percentAssoc
  
  #all methods give same results, so just use Score test
  allDataAISPU[[jj]] = File.trueScore.dat
}

allDataAISPU.mat = do.call(rbind, allDataAISPU)
# ### Decorr Score ###############################################################################################################################
# setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_DecorrScore_Power_BetaZero"))
# 
# allDataDS = list()
# 
# for(jj in 1:length(scores)){
#   File = read.csv(paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DecorrScore_BetaZero.csv"))
#   
#   #pull the row for the true score
#   File.trueScore = File[jj,]
#   
#   #remove column 1
#   File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
#   
#   #need to stack by test
#   library(ramify)
#   File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
#   File.trueScore.dat = as.data.frame(File.trueScore.mat)
#   
#   #Add prevalence
#   File.trueScore.dat$Prev = prev
#   
#   #add score type
#   File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
#   
#   #add OR 
#   File.trueScore.dat$OR = orSize
#   
#   #add method
#   File.trueScore.dat$method = "Decorr. \n Score"
#   
#   #rename Score column to Binary Y
#   colnames(File.trueScore.dat)[1] = "Binary Y" 
#   
#   #add score type
#   File.trueScore.dat$True_Score = scores[jj]
#   #add ld
#   File.trueScore.dat$ld = ld
#   #add % assoc
#   File.trueScore.dat$PercentAssoc = percentAssoc
#   
#   #all methods give same results, so just use Score test
#   allDataDS[[jj]] = File.trueScore.dat
# }
# 
# allDataDS.mat = do.call(rbind, allDataDS)
# 
######## Make Plots  ################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/Plot_Code/Power_BetaZero"))

#Combine The data
allData = rbind(allDataDL.mat, allDataAISPU.mat)
#break up by true score
allData.AMS = allData[allData$True_Score == "AMS",]
allData.BinMM = allData[allData$True_Score == "BinMM",]
allData.IBS = allData[allData$True_Score == "IBS",]
allData.Incomp = allData[allData$True_Score == "Incomp",]

if(orSize == "small"){
  ORSizeCap = "Small"
} else if(orSize == "medium"){
  ORSizeCap = "Medium"
} else if(orSize == "large"){
  ORSizeCap = "Large"
}

if(ld == "LowLD"){
  ldSpaced = "Low LD"
} else if(ld == "HighLD"){
  ldSpaced = "High LD"
}

powerPlot.AMS = ggplot(allData.AMS, aes(fill=Score, y=`Binary Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: AMS")) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") +
  theme(plot.title = element_text(hjust = 0.5, size=16), plot.subtitle = element_text(hjust = 0.5, size=16), legend.position = "none", axis.title = element_text(size = 16), axis.text = element_text(size = 14),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

powerPlot.BinMM = ggplot(allData.BinMM, aes(fill=Score, y=`Binary Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: Binary Mismatch")) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") +
  theme(plot.title = element_text(hjust = 0.5, size=16), legend.position = "none", legend.title = element_text(size = 14), legend.text = element_text(size=12), axis.title = element_text(size = 16), axis.text = element_text(size = 14),axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

powerPlot.IBS = ggplot(allData.IBS, aes(fill=Score, y=`Binary Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: IBS"), subtitle = paste0(gene,", ",ss," Pairs, Outcome Prev. ",prev,"\n ",ORSizeCap," OR, ",percentAssoc,"% Assoc. SNPs, ",ldSpaced)) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), plot.subtitle = element_text(hjust = 0.5, size=16), legend.position = "none", axis.title = element_text(size = 16), axis.text = element_text(size = 14))

powerPlot.Incomp = ggplot(allData.Incomp, aes(fill=Score, y=`Binary Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: Incompatibility")) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), legend.position = "bottom", legend.title = element_text(size = 14), legend.text = element_text(size=12), axis.title = element_text(size = 16), axis.text = element_text(size = 14),axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

subtitle = get_plot_component(powerPlot.IBS, "subtitle") 
legend = get_legend(powerPlot.Incomp)
# and replot suppressing the legend
powerPlot.IBS = powerPlot.IBS + labs(subtitle = NULL)
powerPlot.Incomp = powerPlot.Incomp + theme(legend.position='none')

gridPlot = plot_grid(powerPlot.AMS, powerPlot.BinMM, powerPlot.IBS, powerPlot.Incomp, labels = "auto")
# Now plots are aligned vertically with the legend to the right

plotWlegend = plot_grid(subtitle, legend, ncol = 2)
finalPlot = plot_grid(gridPlot, plotWlegend, ncol = 1, rel_heights = c(1, .1))

#save plot
ggsave2(paste0(gene," ",ss,"Pairs_BinaryY_Prev ",prev,"_Power_",orSize,"OR_",ld,"_",percentAssoc,"SNPsAssoc_AllScores_BetaZero.png"), plot = finalPlot, width = 15, height = 14, units = "in")

######## Continuous Outcome ################################################################################################################
# inputs = c("NAT2", "1000", 20, "small", 25, TRUE) #for continuous
# inputs = c("CHI3L2", "1000", 20, "small", 25, TRUE) #for continuous
inputs = c("ASAH1", "1000", 20, "small", 25, TRUE) #for continuous

inputs = c("NAT2", "500", 20, "small", 25, TRUE) #for continuous
inputs = c("NAT2", "500", 20, "small", 15, TRUE) #for continuous

gene = inputs[1]
ss = inputs[2]
prev = as.numeric(inputs[3])
orSize = inputs[4]
percentAssoc = as.numeric(inputs[5])
lowLD = as.logical(inputs[6])

if(lowLD == TRUE){
  ld = "LowLD"
} else {
  ld = "HighLD"
}

scores = c("AMS", "BinMM", "IBS", "Incomp")

### LHT ###############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_LHT_Power_BetaZero"))

allDataLHTCont = list()

for(jj in 1:length(scores)){
  LHTFile = read.csv(paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_LHT_BetaZero.csv"))
  
  #pull the row for the true score
  LHTFile.trueScore = LHTFile[jj,]
  
  #remove column 1
  LHTFile.trueScore = LHTFile.trueScore[,2:ncol(LHTFile.trueScore)]
  
  #need to stack by test
  library(ramify)
  LHTFile.trueScore.mat = resize(LHTFile.trueScore, nrow = 4, ncol = 3, byrow = TRUE)
  LHTFile.trueScore.dat = as.data.frame(LHTFile.trueScore.mat)
  
  #Add prevalence
  LHTFile.trueScore.dat$Prev = prev
  
  #add score type
  LHTFile.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
  
  #add OR 
  LHTFile.trueScore.dat$OR = orSize
  
  #LRT values
  LHTFile.trueScore.dat.LRT = LHTFile.trueScore.dat[,c(1,4:6)]
  #Wald values
  LHTFile.trueScore.dat.Wald = LHTFile.trueScore.dat[,c(2,4:6)]
  #Score values
  LHTFile.trueScore.dat.Score = LHTFile.trueScore.dat[,c(3:6)]
  
  #add method
  LHTFile.trueScore.dat.LRT$method = "LHT \n LRT"
  LHTFile.trueScore.dat.Wald$method = "LHT \n Wald"
  LHTFile.trueScore.dat.Score$method = "LHT \n Score"
  
  #rename Score column to Binary Y
  colnames(LHTFile.trueScore.dat.LRT)[1] = "Cont. Y" 
  colnames(LHTFile.trueScore.dat.Wald)[1] = "Cont. Y" 
  colnames(LHTFile.trueScore.dat.Score)[1] = "Cont. Y"
  
  #stack everything
  LHTFile.trueScore.all = rbind(LHTFile.trueScore.dat.LRT, LHTFile.trueScore.dat.Wald, LHTFile.trueScore.dat.Score)
  
  #add score type
  LHTFile.trueScore.all$True_Score = scores[jj]
  #add ld
  LHTFile.trueScore.all$ld = ld
  #add % assoc
  LHTFile.trueScore.all$PercentAssoc = percentAssoc
  
  #all methods give same results, so just use Score test
  allDataLHTCont[[jj]] = LHTFile.trueScore.all
}

allDataLHTCont.mat = do.call(rbind, allDataLHTCont)

# ### GLM ###############################################################################################################################
# setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_GLM_Power_BetaZero"))
# 
# allDataGLMCont = list()
# 
# for(jj in 1:length(scores)){
#   File = read.csv(paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_GLM_BetaZero.csv"))
#   
#   #pull the row for the true score
#   File.trueScore = File[jj,]
#   
#   #remove column 1
#   File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
#   
#   #need to stack by test
#   library(ramify)
#   File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
#   File.trueScore.dat = as.data.frame(File.trueScore.mat)
#   
#   #Add prevalence
#   File.trueScore.dat$Prev = prev
#   
#   #add score type
#   File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
#   
#   #add OR 
#   File.trueScore.dat$OR = orSize
#   
#   #add method
#   File.trueScore.dat$method = "GLM"
#   
#   #rename Score column to Binary Y
#   colnames(File.trueScore.dat)[1] = "Cont. Y" 
#   
#   #add score type
#   File.trueScore.dat$True_Score = scores[jj]
#   #add ld
#   File.trueScore.dat$ld = ld
#   #add % assoc
#   File.trueScore.dat$PercentAssoc = percentAssoc
#   
#   #all methods give same results, so just use Score test
#   allDataGLMCont[[jj]] = File.trueScore.dat
# }
# 
# allDataGLMCont.mat = do.call(rbind, allDataGLMCont)
### Desparse Lasso ###############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_DesparseLasso_Power_BetaZero"))

allDataDLCont = list()

for(jj in 1:length(scores)){
  File = read.csv(paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DL_BetaZero.csv"))
  
  #pull the row for the true score
  File.trueScore = File[jj,]
  
  #remove column 1
  File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
  
  #need to stack by test
  library(ramify)
  File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
  File.trueScore.dat = as.data.frame(File.trueScore.mat)
  
  #Add prevalence
  File.trueScore.dat$Prev = prev
  
  #add score type
  File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
  
  #add OR 
  File.trueScore.dat$OR = orSize
  
  #add method
  File.trueScore.dat$method = "Desparse \n Lasso"
  
  #rename Score column to Binary Y
  colnames(File.trueScore.dat)[1] = "Cont. Y" 
  
  #add score type
  File.trueScore.dat$True_Score = scores[jj]
  #add ld
  File.trueScore.dat$ld = ld
  #add % assoc
  File.trueScore.dat$PercentAssoc = percentAssoc
  
  #all methods give same results, so just use Score test
  allDataDLCont[[jj]] = File.trueScore.dat
}

allDataDLCont.mat = do.call(rbind, allDataDLCont)
### AISPU ###############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_AISPU_Power_BetaZero"))

allDataAISPUCont = list()

for(jj in 1:length(scores)){
  File = read.csv(paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_AISPU_BetaZero.csv"))
  
  #pull the row for the true score
  File.trueScore = File[jj,]
  
  #remove column 1
  File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
  
  #need to stack by test
  library(ramify)
  File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
  File.trueScore.dat = as.data.frame(File.trueScore.mat)
  
  #Add prevalence
  File.trueScore.dat$Prev = prev
  
  #add score type
  File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
  
  #add OR 
  File.trueScore.dat$OR = orSize
  
  #add method
  File.trueScore.dat$method = "aiSPU"
  
  #rename Score column to Binary Y
  colnames(File.trueScore.dat)[1] = "Cont. Y" 
  
  #add score type
  File.trueScore.dat$True_Score = scores[jj]
  #add ld
  File.trueScore.dat$ld = ld
  #add % assoc
  File.trueScore.dat$PercentAssoc = percentAssoc
  
  #all methods give same results, so just use Score test
  allDataAISPUCont[[jj]] = File.trueScore.dat
}

allDataAISPUCont.mat = do.call(rbind, allDataAISPUCont)
### Decorr Score ###############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_DecorrScore_Power_BetaZero"))

allDataDSCont = list()

for(jj in 1:length(scores)){
  File = read.csv(paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DecorrScore_BetaZero.csv"))
  
  #pull the row for the true score
  File.trueScore = File[jj,]
  
  #remove column 1
  File.trueScore = File.trueScore[,2:ncol(File.trueScore)]
  
  #need to stack by test
  library(ramify)
  File.trueScore.mat = resize(File.trueScore, nrow = 4, ncol = 1, byrow = TRUE)
  File.trueScore.dat = as.data.frame(File.trueScore.mat)
  
  #Add prevalence
  File.trueScore.dat$Prev = prev
  
  #add score type
  File.trueScore.dat$Score = c("IBS", "Incomp", "AMS", "Bin MM")
  
  #add OR 
  File.trueScore.dat$OR = orSize
  
  #add method
  File.trueScore.dat$method = "Decorr. \n Score"
  
  #rename Score column to Binary Y
  colnames(File.trueScore.dat)[1] = "Cont. Y" 
  
  #add score type
  File.trueScore.dat$True_Score = scores[jj]
  #add ld
  File.trueScore.dat$ld = ld
  #add % assoc
  File.trueScore.dat$PercentAssoc = percentAssoc
  
  #all methods give same results, so just use Score test
  allDataDSCont[[jj]] = File.trueScore.dat
}

allDataDSCont.mat = do.call(rbind, allDataDSCont)
######## Make Plots  ################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/Plot_Code/Power_BetaZero"))

if(orSize == "small"){
  ORSizeCap = "Small"
} else if(orSize == "medium"){
  ORSizeCap = "Medium"
} else if(orSize == "large"){
  ORSizeCap = "Large"
}

if(ld == "LowLD"){
  ldSpaced = "Low LD"
} else if(ld == "HighLD"){
  ldSpaced = "High LD"
}

#Combine The data
allDataCont = rbind(allDataDLCont.mat, allDataAISPUCont.mat)
#break up by true score
allDataCont.AMS = allDataCont[allDataCont$True_Score == "AMS",]
allDataCont.BinMM = allDataCont[allDataCont$True_Score == "BinMM",]
allDataCont.IBS = allDataCont[allDataCont$True_Score == "IBS",]
allDataCont.Incomp = allDataCont[allDataCont$True_Score == "Incomp",]

powerPlot.AMS = ggplot(allDataCont.AMS, aes(fill=Score, y=`Cont. Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: AMS")) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), plot.subtitle = element_text(hjust = 0.5, size=16), legend.position = "none", axis.title = element_text(size = 16), axis.text = element_text(size = 14),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

powerPlot.BinMM = ggplot(allDataCont.BinMM, aes(fill=Score, y=`Cont. Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: Binary Mismatch")) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), legend.position = "none", legend.title = element_text(size = 14), legend.text = element_text(size=12), axis.title = element_text(size = 16), axis.text = element_text(size = 14),axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

powerPlot.IBS = ggplot(allDataCont.IBS, aes(fill=Score, y=`Cont. Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: IBS"), subtitle = paste0(gene," ",ss," Pairs, Continuous Outcome \n",ORSizeCap," OR, ",percentAssoc,"% Assoc. SNPs,",ldSpaced)) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), plot.subtitle = element_text(hjust = 0.5, size=16), legend.position = "none", axis.title = element_text(size = 16), axis.text = element_text(size = 14))

powerPlot.Incomp = ggplot(allDataCont.Incomp, aes(fill=Score, y=`Cont. Y`, x=method)) + 
  geom_bar(position="dodge", color="black", stat="identity") + ggtitle(paste0("True Score: Incompatibility")) + geom_hline(yintercept=0.65, color='blue') +
  geom_hline(yintercept=0.8, color='red') + xlab("Method Used") + ylab("Est. Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Score") + scale_fill_brewer(palette="Spectral") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), legend.position = "bottom", legend.title = element_text(size = 14), legend.text = element_text(size=12), axis.title = element_text(size = 16), axis.text = element_text(size = 14),axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

subtitle = get_plot_component(powerPlot.IBS, "subtitle") 
legend = get_legend(powerPlot.Incomp)
# and replot suppressing the legend
powerPlot.IBS = powerPlot.IBS + labs(subtitle = NULL)
powerPlot.Incomp = powerPlot.Incomp + theme(legend.position='none')

gridPlot = plot_grid(powerPlot.AMS, powerPlot.BinMM, powerPlot.IBS, powerPlot.Incomp, labels = "auto")
# Now plots are aligned vertically with the legend to the right

plotWlegend = plot_grid(subtitle, legend, ncol = 2)
finalPlot = plot_grid(gridPlot, plotWlegend, ncol = 1, rel_heights = c(1, .1))

#save plot
ggsave2(paste0(gene," ",ss,"Pairs_ContinuousY_Power_",orSize,"OR_",ld,"_",percentAssoc,"SNPsAssoc_AllScores_BetaZero.png"), plot = finalPlot, width = 15, height = 14, units = "in")

