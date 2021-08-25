#################################################################
## Gamma = 0 Only Testing, When Beta = 0
# gene = "NAT2"
# gene = "CHI3L2"
gene = "ASAH1"
ss = 500
# ss = 1000
Prev = c(20,10,5)

path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_TIE_GammaBetaZero/",gene,"_",ss,"_All_TIE_GammaBetaZero")
setwd(path)
#################################################################
## Read in LHT data
setwd(paste0(path,"/",gene,"_",ss,"_LHT_TIE_GammaBetaZero"))

LHT.file = read.csv(paste0(gene,"_",ss,"Pairs_LHT_TIE_3KSims_GammaBetaZero.csv"))

## Read in GLM data
setwd(paste0(path,"/",gene,"_",ss,"_GLM_TIE_GammaBetaZero"))

GLM.file = read.csv(paste0(gene,"_",ss,"Pairs_GLM_TIE_3000Sims_GammaBetaZero.csv"))

## Read in Desparse Lasso data
setwd(paste0(path,"/",gene,"_",ss,"_DesparseLasso_TIE_GammaBetaZero"))

DL.file = read.csv(paste0(gene,"_",ss,"Pairs_DesparseLasso_TIE_3000Sims_GammaBetaZero.csv"))

## Read in Decorrelated Score data
setwd(paste0(path,"/",gene,"_",ss,"_DecorrScore_TIE_GammaBetaZero"))

DS.file = read.csv(paste0(gene,"_",ss,"Pairs_DecorrScore_TIE_3000Sims_GammaBetaZero.csv"))

## Read in aiSPU data 
setwd(paste0(path,"/",gene,"_",ss,"_AISPU_TIE_GammaBetaZero"))

aispu.file = read.csv(paste0(gene,"_",ss,"Pairs_AISPU_TIE_3000Sims_scadPenalty_GammaBetaZero.csv"))
#################################################################
#reorganize LHT data
TestList = rep(c("TL","TW","TS"),8)
PrevsList = c(rep(c("20%","10%","5%","Cont."), each=6))
covsList = c(rep(c("No Covs", "Covs"),times=4, each = 3))
LHT.file = LHT.file[,2:ncol(LHT.file)]

LHT.file.reorg = cbind(PrevsList, covsList, TestList, LHT.file)

#pull LHT file by test type
TL.data = LHT.file.reorg[LHT.file.reorg$TestList == "TL",]
TW.data = LHT.file.reorg[LHT.file.reorg$TestList == "TW",]
TS.data = LHT.file.reorg[LHT.file.reorg$TestList == "TS",]

#need to stack by Score used
TL.dataIBS = TL.data[,1:4]
TL.dataIBS$Score = "IBS"
TL.dataIncomp = TL.data[,c(1:3,5)]
TL.dataIncomp$Score = "Incomp."
TL.dataAMS = TL.data[,c(1:3,6)]
TL.dataAMS$Score = "AMS"
TL.dataBinMM = TL.data[,c(1:3,7)]
TL.dataBinMM$Score = "Binary Mismatch"
colnames(TL.dataIBS)[4] = colnames(TL.dataIncomp)[4] = colnames(TL.dataAMS)[4] = colnames(TL.dataBinMM)[4] = "TL"
TL.stacked = rbind(TL.dataIBS, TL.dataIncomp, TL.dataAMS, TL.dataBinMM)

TS.dataIBS = TS.data[,1:4]
TS.dataIncomp = TS.data[,c(1:3,5)]
TS.dataAMS = TS.data[,c(1:3,6)]
TS.dataBinMM = TS.data[,c(1:3,7)]
TS.dataIBS$Score = "IBS"
TS.dataIncomp$Score = "Incomp."
TS.dataAMS$Score = "AMS"
TS.dataBinMM$Score = "Binary Mismatch"
colnames(TS.dataIBS)[4] = colnames(TS.dataIncomp)[4] = colnames(TS.dataAMS)[4] = colnames(TS.dataBinMM)[4] = "TS"
TS.stacked = rbind(TS.dataIBS, TS.dataIncomp, TS.dataAMS, TS.dataBinMM)

TW.dataIBS = TW.data[,1:4]
TW.dataIncomp = TW.data[,c(1:3,5)]
TW.dataAMS = TW.data[,c(1:3,6)]
TW.dataBinMM = TW.data[,c(1:3,7)]
TW.dataIBS$Score = "IBS"
TW.dataIncomp$Score = "Incomp."
TW.dataAMS$Score = "AMS"
TW.dataBinMM$Score = "Binary Mismatch"
colnames(TW.dataIBS)[4] = colnames(TW.dataIncomp)[4] = colnames(TW.dataAMS)[4] = colnames(TW.dataBinMM)[4] = "TW"
TW.stacked = rbind(TW.dataIBS, TW.dataIncomp, TW.dataAMS, TW.dataBinMM)

#reorganize glm data
PrevsListShort = c(rep(c("20%","10%","5%","Cont."), each=2))
covsListShort = c(rep(c("No Covs", "Covs"),times=4))
GLM.file = GLM.file[,2:ncol(GLM.file)]
GLM.file.reorg = cbind(PrevsListShort, covsListShort, GLM.file)

#need to stack by Score used
GLM.file.reorg.IBS = GLM.file.reorg[,1:3]
GLM.file.reorg.Incomp = GLM.file.reorg[,c(1:2,4)]
GLM.file.reorg.AMS = GLM.file.reorg[,c(1:2,5)]
GLM.file.reorg.BinMM = GLM.file.reorg[,c(1:2,6)]
GLM.file.reorg.IBS$Score = "IBS"
GLM.file.reorg.Incomp$Score = "Incomp."
GLM.file.reorg.AMS$Score = "AMS"
GLM.file.reorg.BinMM$Score = "Binary Mismatch"
colnames(GLM.file.reorg.IBS)[3] = colnames(GLM.file.reorg.Incomp)[3] = colnames(GLM.file.reorg.AMS)[3] = colnames(GLM.file.reorg.BinMM)[3] = "GLM"
GLM.stacked = rbind(GLM.file.reorg.IBS, GLM.file.reorg.Incomp, GLM.file.reorg.AMS, GLM.file.reorg.BinMM)

#reorganize desparse Lasso data
DL.file = DL.file[,2:ncol(DL.file)]
DL.file.reorg = cbind(PrevsListShort, covsListShort, DL.file)

#need to stack by Score used
DL.file.reorg.IBS = DL.file.reorg[,1:3]
DL.file.reorg.Incomp = DL.file.reorg[,c(1:2,4)]
DL.file.reorg.AMS = DL.file.reorg[,c(1:2,5)]
DL.file.reorg.BinMM = DL.file.reorg[,c(1:2,6)]
DL.file.reorg.IBS$Score = "IBS"
DL.file.reorg.Incomp$Score = "Incomp."
DL.file.reorg.AMS$Score = "AMS"
DL.file.reorg.BinMM$Score = "Binary Mismatch"
colnames(DL.file.reorg.IBS)[3] = colnames(DL.file.reorg.Incomp)[3] = colnames(DL.file.reorg.AMS)[3] = colnames(DL.file.reorg.BinMM)[3] = "Desparse Lasso"
DL.stacked = rbind(DL.file.reorg.IBS, DL.file.reorg.Incomp, DL.file.reorg.AMS, DL.file.reorg.BinMM)

#reorganize decorr score data
DS.file = DS.file[,2:ncol(DS.file)]
DS.file.reorg = cbind(PrevsListShort, covsListShort, DS.file)

#need to stack by Score used
DS.file.IBS = DS.file.reorg[,1:3]
DS.file.Incomp = DS.file.reorg[,c(1:2,4)]
DS.file.AMS = DS.file.reorg[,c(1:2,5)]
DS.file.BinMM = DS.file.reorg[,c(1:2,6)]
DS.file.IBS$Score = "IBS"
DS.file.Incomp$Score = "Incomp."
DS.file.AMS$Score = "AMS"
DS.file.BinMM$Score = "Binary Mismatch"
colnames(DS.file.IBS)[3] = colnames(DS.file.Incomp)[3] = colnames(DS.file.AMS)[3] = colnames(DS.file.BinMM)[3] = "Decorr. Score" 
DS.stacked = rbind(DS.file.IBS, DS.file.Incomp, DS.file.AMS, DS.file.BinMM)

#reorganize the aispu data
aispu.file = aispu.file[,2:ncol(aispu.file)]
aispu.file.reorg = cbind(PrevsListShort, covsListShort, aispu.file)

#need to stack by Score used
aispu.file.IBS = aispu.file.reorg[,1:3]
aispu.file.Incomp = aispu.file.reorg[,c(1:2,4)]
aispu.file.AMS = aispu.file.reorg[,c(1:2,5)]
aispu.file.BinMM = aispu.file.reorg[,c(1:2,6)]
aispu.file.IBS$Score = "IBS"
aispu.file.Incomp$Score = "Incomp."
aispu.file.AMS$Score = "AMS"
aispu.file.BinMM$Score = "Binary Mismatch"
colnames(aispu.file.IBS)[3] = colnames(aispu.file.Incomp)[3] = colnames(aispu.file.AMS)[3] = colnames(aispu.file.BinMM)[3] = "aiSPU"
aispu.stacked = rbind(aispu.file.IBS, aispu.file.Incomp, aispu.file.AMS, aispu.file.BinMM)

##########################################################
#Combine all the data
allData = cbind(TL.stacked$PrevsList, TL.stacked$Score, TL.stacked$covsList, round(TL.stacked$TL,2), round(TW.stacked$TW,2),
                round(TS.stacked$TS, 2), round(GLM.stacked$GLM,2), round(DL.stacked$`Desparse Lasso`,4), round(DS.stacked$`Decorr. Score`,2), round(aispu.stacked$aiSPU, 2))
allData.df = as.data.frame(allData)
colnames(allData.df)[1] = "Outcome Prevalence"
colnames(allData.df)[2] = "Score Used"
colnames(allData.df)[3] = "Covariates, Y/N"
colnames(allData.df)[4] = "TL"
colnames(allData.df)[5] = "TW"
colnames(allData.df)[6] = "TS"
colnames(allData.df)[7] = "GLM"
colnames(allData.df)[8] = "Desparse Lasso"
colnames(allData.df)[9] = "Decorr. Score"
colnames(allData.df)[10] = "aiSPU"

allData.df$`Outcome Prevalence` = factor(allData.df$`Outcome Prevalence`, levels = c("20%", "10%", "5%", "Cont."))

allData.ordered = allData.df[order(allData.df$`Outcome Prevalence`, decreasing = FALSE),]

### Without Covs, No covs#############
allData.allCovs = allData.df[allData.df$`Covariates, Y/N` == "Covs",]
allData.allCovs = allData.allCovs[,c(1:2,4:10)]
allData.allCovs.ordered = allData.allCovs[order(allData.allCovs$`Outcome Prevalence`),]

library(flextable)
library(dplyr)

ft = flextable(allData.ordered) %>%
merge_v(j = c("Outcome Prevalence","Score Used")) %>%
colformat_double(j= c(4:7,9:10), digits = 2) %>%
colformat_double(j = "Desparse Lasso", digits = 4) %>%
theme_vanilla() %>%
fontsize(part = "all", size = 8) %>%
align(j = c("Outcome Prevalence", "Score Used", "Covariates, Y/N", "TL", "TW", "TS", "GLM", "Desparse Lasso", "Decorr. Score", "aiSPU"), align = "center", part = "all" )

setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_TIE_GammaBetaZero/",gene,"_",ss,"_All_TIE_GammaBetaZero"))
save_as_docx(ft, path = paste0(gene,"_",ss,"_All_TIE_GammaBetaZero_Table.docx"))

ftallcovs = flextable(allData.allCovs.ordered) %>%
merge_v(j = c("Outcome Prevalence")) %>%
colformat_double(j= c(3:7,8:9), digits = 2) %>%
colformat_double(j = "Desparse Lasso", digits = 4) %>%
theme_vanilla() %>%
fontsize(part = "all", size = 8) %>%
align(j = c("Outcome Prevalence", "Score Used", "TL", "TW", "TS", "GLM", "Desparse Lasso", "Decorr. Score", "aiSPU"), align = "center", part = "all" )
save_as_docx(ftallcovs, path = paste0(gene,"_",ss,"_All_TIE_GammaBetaZero_Table_AllCovsIncluded.docx"))

allData.allCovs.ordered.subset = allData.allCovs.ordered[allData.allCovs.ordered$`Outcome Prevalence` == "20%" | allData.allCovs.ordered$`Outcome Prevalence`=="Cont.",]
allData.allCovs.ordered.subset[1,6] = "0.40"
allData.allCovs.ordered.subset[4,6] = "0.40"
allData.allCovs.ordered.subset$TS = as.numeric(allData.allCovs.ordered.subset$TS)
allData.allCovs.ordered.subset$TW = as.numeric(allData.allCovs.ordered.subset$TW)
allData.allCovs.ordered.subset$TL = as.numeric(allData.allCovs.ordered.subset$TL)
allData.allCovs.ordered.subset$GLM = as.numeric(allData.allCovs.ordered.subset$GLM)
allData.allCovs.ordered.subset$`Desparse Lasso` = as.numeric(allData.allCovs.ordered.subset$`Desparse Lasso`)
allData.allCovs.ordered.subset$`Decorr. Score` = as.numeric(allData.allCovs.ordered.subset$`Decorr. Score`)
allData.allCovs.ordered.subset$aiSPU = as.numeric(allData.allCovs.ordered.subset$aiSPU)

ftallcovs = flextable(allData.allCovs.ordered.subset) %>%
  merge_v(j = c("Outcome Prevalence")) %>%
  colformat_double(j= c(3:9), digits = 2) %>%
  # colformat_double(j = "Desparse Lasso", digits = 2) %>%
  theme_vanilla() %>%
  fontsize(part = "all", size = 8) %>%
  align(j = c("Outcome Prevalence", "Score Used", "TL", "TW", "TS", "GLM", "Desparse Lasso", "Decorr. Score", "aiSPU"), align = "center", part = "all" ) %>%
  bold(i = c(1:4), j = 3, bold = TRUE, part = "body") %>%
  bold(i = c(5:8), j = 5, bold = TRUE, part = "body") %>%
  bold(i = c(1:4,6:8), j = 6, bold = TRUE, part = "body") %>%
  bold(i = c(1:8), j = 8, bold = TRUE, part = "body")

save_as_docx(ftallcovs, path = paste0(gene,"_",ss,"_All_TIE_GammaBetaZero_Table_AllCovsIncluded_Prev20Cont.docx"))


setwd("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Thesis_Defense")
save_as_image(ftallcovs, path = paste0("PaperIII_TIETable1.png"))
#################################################################
## Gamma = 0 Only Testing, When Beta is associated
# gene = "NAT2"
# gene = "CHI3L2"
gene = "ASAH1"
# ss = 500
ss = 1000
Prev = c(20,10,5)

path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_TIE_GammaZero_BetaAssoc/",gene,"_",ss,"_All_TIE_GammaZero_BetaAssoc")
setwd(path)
#################################################################
## Read in LHT data
setwd(paste0(path,"/",gene,"_",ss,"_LHT_TIE_GammaZero_BetaAssoc"))

LHT.file = read.csv(paste0(gene,"_",ss,"Pairs_LHT_TIE_3KSims_GammaZero_BetaAssoc.csv"))

## Read in GLM data
setwd(paste0(path,"/",gene,"_",ss,"_GLM_TIE_GammaZero_BetaAssoc"))

GLM.file = read.csv(paste0(gene,"_",ss,"Pairs_GLM_TIE_3000Sims_GammaZero_BetaAssoc.csv"))

## Read in Desparse Lasso data
setwd(paste0(path,"/",gene,"_",ss,"_DesparseLasso_TIE_GammaZero_BetaAssoc"))

DL.file = read.csv(paste0(gene,"_",ss,"Pairs_DesparseLasso_TIE_3000Sims_GammaZero_BetaAssoc.csv"))

## Read in Decorrelated Score data
setwd(paste0(path,"/",gene,"_",ss,"_DecorrScore_TIE_GammaZero_BetaAssoc"))

DS.file = read.csv(paste0(gene,"_",ss,"Pairs_DecorrScore_TIE_3000Sims_GammaZero_BetaAssoc.csv"))

## Read in aiSPU data 
setwd(paste0(path,"/",gene,"_",ss,"_AISPU_TIE_GammaZero_BetaAssoc"))

aispu.file = read.csv(paste0(gene,"_",ss,"Pairs_AISPU_TIE_3000Sims_scadPenalty_GammaZero_BetaAssoc.csv"))

###
#reorganize LHT data
TestList = rep(c("TL","TW","TS"),12)
PrevsList = c(rep(c("20%","10%","5%","Cont."), each=9))
ORList = c(rep(c("Small OR", "Medium OR", "Large OR"), times=4, each = 3))
LHT.file = LHT.file[,2:ncol(LHT.file)]

LHT.file.reorg = cbind(PrevsList, TestList, ORList, LHT.file)

#pull LHT file by test type
TL.data = LHT.file.reorg[LHT.file.reorg$TestList == "TL",]
TW.data = LHT.file.reorg[LHT.file.reorg$TestList == "TW",]
TS.data = LHT.file.reorg[LHT.file.reorg$TestList == "TS",]

#need to stack by % Assoc SNP and LD
TL.data5Low = TL.data[,1:4]
TL.data5Low$PercentAssoc = "5% R SNPs"
TL.data5Low$LD = "Low LD"
TL.data5High = TL.data[,c(1:3,5)]
TL.data5High$PercentAssoc = "5% R SNPs"
TL.data5High$LD = "High LD"
TL.data15Low = TL.data[,c(1:3,6)]
TL.data15Low$PercentAssoc = "15% R SNPs"
TL.data15Low$LD = "Low LD"
TL.data15High = TL.data[,c(1:3,7)]
TL.data15High$PercentAssoc = "15% R SNPs"
TL.data15High$LD = "High LD"
TL.data25Low = TL.data[,c(1:3,8)]
TL.data25Low$PercentAssoc = "25% R SNPs"
TL.data25Low$LD = "Low LD"
TL.data25High = TL.data[,c(1:3,9)]
TL.data25High$PercentAssoc = "25% R SNPs"
TL.data25High$LD = "High LD"
colnames(TL.data5Low)[4] = colnames(TL.data5High)[4] = colnames(TL.data15Low)[4] = colnames(TL.data15High)[4] = colnames(TL.data25Low)[4] = colnames(TL.data25High)[4] = "TL"
TL.stacked = rbind(TL.data5Low, TL.data5High, TL.data15Low, TL.data15High, TL.data25Low, TL.data25High)

TS.data5Low = TS.data[,1:4]
TS.data5Low$PercentAssoc = "5% R SNPs"
TS.data5Low$LD = "Low LD"
TS.data5High = TS.data[,c(1:3,5)]
TS.data5High$PercentAssoc = "5% R SNPs"
TS.data5High$LD = "High LD"
TS.data15Low = TS.data[,c(1:3,6)]
TS.data15Low$PercentAssoc = "15% R SNPs"
TS.data15Low$LD = "Low LD"
TS.data15High = TS.data[,c(1:3,7)]
TS.data15High$PercentAssoc = "15% R SNPs"
TS.data15High$LD = "High LD"
TS.data25Low = TS.data[,c(1:3,8)]
TS.data25Low$PercentAssoc = "25% R SNPs"
TS.data25Low$LD = "Low LD"
TS.data25High = TS.data[,c(1:3,9)]
TS.data25High$PercentAssoc = "25% R SNPs"
TS.data25High$LD = "High LD"
colnames(TS.data5Low)[4] = colnames(TS.data5High)[4] = colnames(TS.data15Low)[4] = colnames(TS.data15High)[4] = colnames(TS.data25Low)[4] = colnames(TS.data25High)[4] = "TS"
TS.stacked = rbind(TS.data5Low, TS.data5High, TS.data15Low, TS.data15High, TS.data25Low, TS.data25High)

TW.data5Low = TW.data[,1:4]
TW.data5Low$PercentAssoc = "5% R SNPs"
TW.data5Low$LD = "Low LD"
TW.data5High = TW.data[,c(1:3,5)]
TW.data5High$PercentAssoc = "5% R SNPs"
TW.data5High$LD = "High LD"
TW.data15Low = TW.data[,c(1:3,6)]
TW.data15Low$PercentAssoc = "15% R SNPs"
TW.data15Low$LD = "Low LD"
TW.data15High = TW.data[,c(1:3,7)]
TW.data15High$PercentAssoc = "15% R SNPs"
TW.data15High$LD = "High LD"
TW.data25Low = TW.data[,c(1:3,8)]
TW.data25Low$PercentAssoc = "25% R SNPs"
TW.data25Low$LD = "Low LD"
TW.data25High = TW.data[,c(1:3,9)]
TW.data25High$PercentAssoc = "25% R SNPs"
TW.data25High$LD = "High LD"
colnames(TW.data5Low)[4] = colnames(TW.data5High)[4] = colnames(TW.data15Low)[4] = colnames(TW.data15High)[4] = colnames(TW.data25Low)[4] = colnames(TW.data25High)[4] = "TW"
TW.stacked = rbind(TW.data5Low, TW.data5High, TW.data15Low, TW.data15High, TW.data25Low, TW.data25High)

#reorganize glm data
PrevsListShort = c(rep(c("20%","10%","5%","Cont."), each=3))
ORListShort = c(rep(c("Small OR", "Medium OR", "Large OR"), times=4))
GLM.file = GLM.file[,2:ncol(GLM.file)]
GLM.file.reorg = cbind(PrevsListShort, ORListShort, GLM.file)

#need to stack by Score used
GLM.file.reorg.5Low = GLM.file.reorg[,1:3]
GLM.file.reorg.5High = GLM.file.reorg[,c(1:2,4)]
GLM.file.reorg.15Low = GLM.file.reorg[,c(1:2,5)]
GLM.file.reorg.15High = GLM.file.reorg[,c(1:2,6)]
GLM.file.reorg.25Low = GLM.file.reorg[,c(1:2,7)]
GLM.file.reorg.25High = GLM.file.reorg[,c(1:2,8)]

GLM.file.reorg.5High$PercentAssoc = GLM.file.reorg.5Low$PercentAssoc = "5% R SNPs"
GLM.file.reorg.15Low$PercentAssoc = GLM.file.reorg.15High$PercentAssoc = "15% R SNPs"
GLM.file.reorg.25Low$PercentAssoc = GLM.file.reorg.25High$PercentAssoc = "25% R SNPs"
GLM.file.reorg.5Low$LD = GLM.file.reorg.15Low$LD = GLM.file.reorg.25Low$LD = "Low LD"
GLM.file.reorg.5High$LD = GLM.file.reorg.15High$LD = GLM.file.reorg.25High$LD = "High LD"

colnames(GLM.file.reorg.5Low)[3] = colnames(GLM.file.reorg.5High)[3] = colnames(GLM.file.reorg.15Low)[3] = colnames(GLM.file.reorg.15High)[3] = colnames(GLM.file.reorg.25Low)[3] = colnames(GLM.file.reorg.25High)[3] = "GLM"
GLM.stacked = rbind(GLM.file.reorg.5Low, GLM.file.reorg.5High, GLM.file.reorg.15Low, GLM.file.reorg.15High, GLM.file.reorg.25Low, GLM.file.reorg.25High)

#reorganize desparse lasso data
PrevsListShort = c(rep(c("20%","10%","5%","Cont."), each=3))
ORListShort = c(rep(c("Small OR", "Medium OR", "Large OR"), times=4))
DL.file = DL.file[,2:ncol(DL.file)]
DL.file.reorg = cbind(PrevsListShort, ORListShort, DL.file)

#need to stack by Score used
DL.file.reorg.5Low = DL.file.reorg[,1:3]
DL.file.reorg.5High = DL.file.reorg[,c(1:2,4)]
DL.file.reorg.15Low = DL.file.reorg[,c(1:2,5)]
DL.file.reorg.15High = DL.file.reorg[,c(1:2,6)]
DL.file.reorg.25Low = DL.file.reorg[,c(1:2,7)]
DL.file.reorg.25High = DL.file.reorg[,c(1:2,8)]

DL.file.reorg.5High$PercentAssoc = DL.file.reorg.5Low$PercentAssoc = "5% R SNPs"
DL.file.reorg.15Low$PercentAssoc = DL.file.reorg.15High$PercentAssoc = "15% % SNPs"
DL.file.reorg.25Low$PercentAssoc = DL.file.reorg.25High$PercentAssoc = "25% % SNPs"
DL.file.reorg.5Low$LD = DL.file.reorg.15Low$LD = DL.file.reorg.25Low$LD = "Low LD"
DL.file.reorg.5High$LD = DL.file.reorg.15High$LD = DL.file.reorg.25High$LD = "High LD"

colnames(DL.file.reorg.5Low)[3] = colnames(DL.file.reorg.5High)[3] = colnames(DL.file.reorg.15Low)[3] = colnames(DL.file.reorg.15High)[3] = colnames(DL.file.reorg.25Low)[3] = colnames(DL.file.reorg.25High)[3] = "DL"
DL.stacked = rbind(DL.file.reorg.5Low, DL.file.reorg.5High, DL.file.reorg.15Low, DL.file.reorg.15High, DL.file.reorg.25Low, DL.file.reorg.25High)

#reorganize decorr score data
PrevsListShort = c(rep(c("20%","10%","5%","Cont."), each=3))
ORListShort = c(rep(c("Small OR", "Medium OR", "Large OR"), times=4))
DS.file = DS.file[,2:ncol(DS.file)]
DS.file.reorg = cbind(PrevsListShort, ORListShort, DS.file)

#need to stack by Score used
DS.file.reorg.5Low = DS.file.reorg[,1:3]
DS.file.reorg.5High = DS.file.reorg[,c(1:2,4)]
DS.file.reorg.15Low = DS.file.reorg[,c(1:2,5)]
DS.file.reorg.15High = DS.file.reorg[,c(1:2,6)]
DS.file.reorg.25Low = DS.file.reorg[,c(1:2,7)]
DS.file.reorg.25High = DS.file.reorg[,c(1:2,8)]

DS.file.reorg.5High$PercentAssoc = DS.file.reorg.5Low$PercentAssoc = "5% R SNPs"
DS.file.reorg.15Low$PercentAssoc = DS.file.reorg.15High$PercentAssoc = "15% % SNPs"
DS.file.reorg.25Low$PercentAssoc = DS.file.reorg.25High$PercentAssoc = "25% % SNPs"
DS.file.reorg.5Low$LD = DS.file.reorg.15Low$LD = DS.file.reorg.25Low$LD = "Low LD"
DS.file.reorg.5High$LD = DS.file.reorg.15High$LD = DS.file.reorg.25High$LD = "High LD"

colnames(DS.file.reorg.5Low)[3] = colnames(DS.file.reorg.5High)[3] = colnames(DS.file.reorg.15Low)[3] = colnames(DS.file.reorg.15High)[3] = colnames(DS.file.reorg.25Low)[3] = colnames(DS.file.reorg.25High)[3] = "DS"
DS.stacked = rbind(DS.file.reorg.5Low, DS.file.reorg.5High, DS.file.reorg.15Low, DS.file.reorg.15High, DS.file.reorg.25Low, DS.file.reorg.25High)

#reorganize aispu data
PrevsListShort = c(rep(c("20%","10%","5%","Cont."), each=3))
ORListShort = c(rep(c("Small OR", "Medium OR", "Large OR"), times=4))
aispu.file = aispu.file[,2:ncol(aispu.file)]
aispu.file.reorg = cbind(PrevsListShort, ORListShort, aispu.file)

#need to stack by Score used
aispu.file.reorg.5Low = aispu.file.reorg[,1:3]
aispu.file.reorg.5High = aispu.file.reorg[,c(1:2,4)]
aispu.file.reorg.15Low = aispu.file.reorg[,c(1:2,5)]
aispu.file.reorg.15High = aispu.file.reorg[,c(1:2,6)]
aispu.file.reorg.25Low = aispu.file.reorg[,c(1:2,7)]
aispu.file.reorg.25High = aispu.file.reorg[,c(1:2,8)]

aispu.file.reorg.5High$PercentAssoc = aispu.file.reorg.5Low$PercentAssoc = "5% R SNPs"
aispu.file.reorg.15Low$PercentAssoc = aispu.file.reorg.15High$PercentAssoc = "15% % SNPs"
aispu.file.reorg.25Low$PercentAssoc = aispu.file.reorg.25High$PercentAssoc = "25% % SNPs"
aispu.file.reorg.5Low$LD = aispu.file.reorg.15Low$LD = aispu.file.reorg.25Low$LD = "Low LD"
aispu.file.reorg.5High$LD = aispu.file.reorg.15High$LD = aispu.file.reorg.25High$LD = "High LD"

colnames(aispu.file.reorg.5Low)[3] = colnames(aispu.file.reorg.5High)[3] = colnames(aispu.file.reorg.15Low)[3] = colnames(aispu.file.reorg.15High)[3] = colnames(aispu.file.reorg.25Low)[3] = colnames(aispu.file.reorg.25High)[3] = "aispu"
aispu.stacked = rbind(aispu.file.reorg.5Low, aispu.file.reorg.5High, aispu.file.reorg.15Low, aispu.file.reorg.15High, aispu.file.reorg.25Low, aispu.file.reorg.25High)

#Combine all the data
allData = cbind(TL.stacked$PrevsList, TL.stacked$ORList, TL.stacked$PercentAssoc, TL.stacked$LD, round(TL.stacked$TL,2), round(TW.stacked$TW,2),
                round(TS.stacked$TS, 2), round(GLM.stacked$GLM,2), round(DL.stacked$DL,4), round(DS.stacked$DS,2), round(aispu.stacked$aispu, 2))
allData.df = as.data.frame(allData)
colnames(allData.df)[1] = "Outcome Prevalence"
colnames(allData.df)[2] = "OR Size"
colnames(allData.df)[3] = "Percent R SNPs Assoc."
colnames(allData.df)[4] = "LD"
colnames(allData.df)[5] = "TL"
colnames(allData.df)[6] = "TW"
colnames(allData.df)[7] = "TS"
colnames(allData.df)[8] = "GLM"
colnames(allData.df)[9] = "Desparse Lasso"
colnames(allData.df)[10] = "Decorr. Score"
colnames(allData.df)[11] = "aiSPU"
allData.df$`Outcome Prevalence` = factor(allData.df$`Outcome Prevalence`, levels = c("Cont.", "5%", "10%",  "20%"))

allData.ordered = allData.df[order(allData.df$`Outcome Prevalence`, allData.df$`OR Size`, decreasing = TRUE),]

library(flextable)
library(dplyr)

ft = flextable(allData.ordered) %>%
  merge_v(j = c("Outcome Prevalence","OR Size", "Percent R SNPs Assoc.")) %>%
  colformat_double(j= c(5:8,10:11), digits = 2) %>%
  colformat_double(j = "Desparse Lasso", digits = 4) %>%
  theme_vanilla() %>%
  fontsize(part = "all", size = 8) %>%
  align(j = c("Outcome Prevalence", "OR Size", "Percent R SNPs Assoc.", "LD", "TL", "TW", "TS", "GLM", "Desparse Lasso", "Decorr. Score", "aiSPU"), align = "center", part = "all" )

plot(as.factor(DL.stacked$PrevsListShort), DL.stacked$DL)
plot(as.factor(DL.stacked$ORListShort), DL.stacked$DL)
plot(as.factor(DL.stacked$PercentAssoc), DL.stacked$DL)
plot(as.factor(DL.stacked$LD), DL.stacked$DL)

allData.ordered.subset = allData.ordered[allData.ordered$`Outcome Prevalence` == "20%",]
allData.ordered.subset$TS = as.numeric(allData.ordered.subset$TS)
allData.ordered.subset$TW = as.numeric(allData.ordered.subset$TW)
allData.ordered.subset$TL = as.numeric(allData.ordered.subset$TL)
allData.ordered.subset$GLM = as.numeric(allData.ordered.subset$GLM)
allData.ordered.subset$`Desparse Lasso` = as.numeric(allData.ordered.subset$`Desparse Lasso`)
allData.ordered.subset$`Decorr. Score` = as.numeric(allData.ordered.subset$`Decorr. Score`)
allData.ordered.subset$aiSPU = as.numeric(allData.ordered.subset$aiSPU)

ftallcovs = flextable(allData.ordered.subset) %>%
  merge_v(j = c("Outcome Prevalence","OR Size", "Percent R SNPs Assoc.")) %>%
  colformat_double(j= c(5:11), digits = 2) %>%
  colformat_double(j= c(9), digits = 3) %>%
  theme_vanilla() %>%
  fontsize(part = "all", size = 8) %>%
  align(j = c("Outcome Prevalence", "OR Size", "Percent R SNPs Assoc.", "LD", "TL", "TW", "TS", "GLM", "Desparse Lasso", "Decorr. Score", "aiSPU"), align = "center", part = "all" ) %>%
  bold(i = c(1:18), j = c(5,8), bold = TRUE, part = "body") %>%
  bold(i = c(2:18), j = 6, bold = TRUE, part = "body") %>%
  bold(i = c(2:6, 8:18), j = 7, bold = TRUE, part = "body")


save_as_docx(ftallcovs, path = paste0(gene,"_",ss,"_All_TIE_GammaBetaZero_Table_AllCovsIncluded_Prev20Cont.docx"))


setwd("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Thesis_Defense")
save_as_image(ftallcovs, path = paste0("PaperIII_TIETable2-binary.png"))

allData.ordered.subset = allData.ordered[allData.ordered$`Outcome Prevalence` == "Cont.",]
allData.ordered.subset$TS = as.numeric(allData.ordered.subset$TS)
allData.ordered.subset$TW = as.numeric(allData.ordered.subset$TW)
allData.ordered.subset$TL = as.numeric(allData.ordered.subset$TL)
allData.ordered.subset$GLM = as.numeric(allData.ordered.subset$GLM)
allData.ordered.subset$`Desparse Lasso` = as.numeric(allData.ordered.subset$`Desparse Lasso`)
allData.ordered.subset$`Decorr. Score` = as.numeric(allData.ordered.subset$`Decorr. Score`)
allData.ordered.subset$aiSPU = as.numeric(allData.ordered.subset$aiSPU)

ftallcovs = flextable(allData.ordered.subset) %>%
  merge_v(j = c("Outcome Prevalence","OR Size", "Percent R SNPs Assoc.")) %>%
  colformat_double(j= c(5:11), digits = 2) %>%
  colformat_double(j= c(9), digits = 2) %>%
  theme_vanilla() %>%
  fontsize(part = "all", size = 8) %>%
  align(j = c("Outcome Prevalence", "OR Size", "Percent R SNPs Assoc.", "LD", "TL", "TW", "TS", "GLM", "Desparse Lasso", "Decorr. Score", "aiSPU"), align = "center", part = "all" ) %>%
  bold(i = c(1:6,8:12,14:18), j = c(5), bold = TRUE, part = "body") %>%
  bold(i = c(1:18), j = 8, bold = TRUE, part = "body") %>%
  bold(i = c(2:6,8:12,14:18), j = 6, bold = TRUE, part = "body") %>%
  bold(i = c(6,10,12,14,16,18), j = 11, bold = TRUE, part = "body") %>%
  bold(i = c(1:12,14:18), j = 7, bold = TRUE, part = "body")
setwd("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Thesis_Defense")
save_as_image(ftallcovs, path = paste0("PaperIII_TIETable2-cont.png"))

