#inputs when score is true:
# inputs = c("NAT2", "1000", 20, "small", 25, FALSE)
# inputs = c("NAT2", "1000", 10, "medium", 25, TRUE)
# inputs = c("NAT2", "1000", 10, "large", 5, FALSE)

inputs = c("NAT2", "500", 20, "small", 25, FALSE)
inputs = c("NAT2", "500", 20, "medium", 25, TRUE)
inputs = c("NAT2", "500", 10, "large", 5, FALSE)

# inputs = c("CHI3L2", "1000", 20, "small", 25, FALSE)
# inputs = c("CHI3L2", "1000", 10, "medium", 25, TRUE)
# inputs = c("CHI3L2", "1000", 10, "large", 5, FALSE)

inputs = c("ASAH1", "1000", 20, "small", 15, FALSE)
inputs = c("ASAH1", "1000", 10, "medium", 15, TRUE)
inputs = c("ASAH1", "1000", 10, "large", 5, FALSE)

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

alpha = 0.05

#### LHT ##########################################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_Power_BetaZero/",gene,"_",ss,"_LHT_Power_BetaZero"))

# scores = c("AMS.gene", "BinMM.gene", "IBS.gene", "Incomp.gene") #for old code
scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat = matrix(NA, nrow = length(scores), ncol = 12)

for(ii in 1:length(scores)){
  # filenames = c(paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",percentAssoc,"SNPsAssoc_GammaOnly_OR",orSize,"_",ld,"_Sim1to1000_Stats_ForceCovFit.csv"),
  #               paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",percentAssoc,"SNPsAssoc_GammaOnly_OR",orSize,"_",ld,"_Sim1001to2000_Stats_ForceCovFit.csv"),
  #               paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",percentAssoc,"SNPsAssoc_GammaOnly_OR",orSize,"_",ld,"_Sim2001to3000_Stats_ForceCovFit.csv"))
  filenames = c(paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1to1000_Stats.csv"),
                paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1001to2000_Stats.csv"),
                paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim2001to3000_Stats.csv"))
  files = lapply(filenames, read.csv) 
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  # invalidTL = which(files.mat$V1 < 0 | files.mat$V4 < 0 | files.mat$V7 < 0 | files.mat$V10 < 0)
  # invalidTW = which(files.mat$V2 < 0 | files.mat$V5 < 0 | files.mat$V8 < 0 | files.mat$V11 < 0)
  # invalidTS = which(files.mat$V3 < 0 | files.mat$V6 < 0 | files.mat$V9 < 0 | files.mat$V12 < 0)
  invalidTL = which(files.mat$IBS < 0 | files.mat$Incomp < 0 | files.mat$AMS < 0 | files.mat$BinMM < 0)
  invalidTW = which(files.mat$IBS.1 < 0 | files.mat$Incomp.1 < 0 | files.mat$AMS.1< 0 | files.mat$BinMM.1 < 0)
  invalidTS = which(files.mat$IBS.2 < 0 | files.mat$Incomp.2 < 0 | files.mat$AMS.2 < 0 | files.mat$BinMM.2 < 0)
  
  # pValfiles = c(paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",percentAssoc,"SNPsAssoc_GammaOnly_OR",orSize,"_",ld,"_Sim1to1000_PValues_ForceCovFit.csv"),
  #               paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",percentAssoc,"SNPsAssoc_GammaOnly_OR",orSize,"_",ld,"_Sim1001to2000_PValues_ForceCovFit.csv"),
  #               paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",percentAssoc,"SNPsAssoc_GammaOnly_OR",orSize,"_",ld,"_Sim2001to3000_PValues_ForceCovFit.csv"))
  pValfiles = c(paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1to1000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1001to2000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_LinHypTest_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim2001to3000_PValues.csv"))
  
  allFiles = lapply(pValfiles, read.csv)
  allFiles.mat = rbind(allFiles[[1]], allFiles[[2]], allFiles[[3]])
  #remove column 1
  allFiles.mat = allFiles.mat[,2:ncol(allFiles.mat)]
  
  if(length(invalidTL) > 0){
    # allFiles.mat.validTL.IBS = allFiles.mat$V1[-invalidTL]
    # allFiles.mat.validTL.Incomp = allFiles.mat$V4[-invalidTL]
    # allFiles.mat.validTL.AMS = allFiles.mat$V7[-invalidTL]
    # allFiles.mat.validTL.BinMM = allFiles.mat$V10[-invalidTL]
    allFiles.mat.validTL.IBS = allFiles.mat$IBS[-invalidTL]
    allFiles.mat.validTL.Incomp = allFiles.mat$Incomp[-invalidTL]
    allFiles.mat.validTL.AMS = allFiles.mat$AMS[-invalidTL]
    allFiles.mat.validTL.BinMM = allFiles.mat$BinMM[-invalidTL]
  } else {
    # allFiles.mat.validTL.IBS = allFiles.mat$V1
    # allFiles.mat.validTL.Incomp = allFiles.mat$V4
    # allFiles.mat.validTL.AMS = allFiles.mat$V7
    # allFiles.mat.validTL.BinMM = allFiles.mat$V10
    allFiles.mat.validTL.IBS = allFiles.mat$IBS
    allFiles.mat.validTL.Incomp = allFiles.mat$Incomp
    allFiles.mat.validTL.AMS = allFiles.mat$AMS
    allFiles.mat.validTL.BinMM = allFiles.mat$BinMM
  }
  
  if(length(invalidTW) > 0){
    # allFiles.mat.validTW.IBS = allFiles.mat$V2[-invalidTW]
    # allFiles.mat.validTW.Incomp = allFiles.mat$V5[-invalidTW]
    # allFiles.mat.validTW.AMS = allFiles.mat$V8[-invalidTW]
    # allFiles.mat.validTW.BinMM = allFiles.mat$V11[-invalidTW]
    allFiles.mat.validTW.IBS = allFiles.mat$IBS.1[-invalidTW]
    allFiles.mat.validTW.Incomp = allFiles.mat$Incomp.1[-invalidTW]
    allFiles.mat.validTW.AMS = allFiles.mat$AMS.1[-invalidTW]
    allFiles.mat.validTW.BinMM = allFiles.mat$BinMM.1[-invalidTW]
  } else {
    # allFiles.mat.validTW.IBS = allFiles.mat$V2
    # allFiles.mat.validTW.Incomp = allFiles.mat$V5
    # allFiles.mat.validTW.AMS = allFiles.mat$V8
    # allFiles.mat.validTW.BinMM = allFiles.mat$V11
    allFiles.mat.validTW.IBS = allFiles.mat$IBS.1
    allFiles.mat.validTW.Incomp = allFiles.mat$Incomp.1
    allFiles.mat.validTW.AMS = allFiles.mat$AMS.1
    allFiles.mat.validTW.BinMM = allFiles.mat$BinMM.1
  }
  
  if(length(invalidTS) > 0){
    # allFiles.mat.validTS.IBS = allFiles.mat$V3[-invalidTS]
    # allFiles.mat.validTS.Incomp = allFiles.mat$V6[-invalidTS]
    # allFiles.mat.validTS.AMS = allFiles.mat$V9[-invalidTS]
    # allFiles.mat.validTS.BinMM = allFiles.mat$V12[-invalidTS]
    allFiles.mat.validTS.IBS = allFiles.mat$IBS.2[-invalidTS]
    allFiles.mat.validTS.Incomp = allFiles.mat$Incomp.2[-invalidTS]
    allFiles.mat.validTS.AMS = allFiles.mat$AMS.2[-invalidTS]
    allFiles.mat.validTS.BinMM = allFiles.mat$BinMM.2[-invalidTS]
  } else {
    # allFiles.mat.validTS.IBS = allFiles.mat$V3
    # allFiles.mat.validTS.Incomp = allFiles.mat$V6
    # allFiles.mat.validTS.AMS = allFiles.mat$V9
    # allFiles.mat.validTS.BinMM = allFiles.mat$V12
    allFiles.mat.validTS.IBS = allFiles.mat$IBS.2
    allFiles.mat.validTS.Incomp = allFiles.mat$Incomp.2
    allFiles.mat.validTS.AMS = allFiles.mat$AMS.2
    allFiles.mat.validTS.BinMM = allFiles.mat$BinMM.2
  }
  
  PowerrateTL.IBS = mean(allFiles.mat.validTL.IBS > 0)
  PowerrateTL.Incomp = mean(allFiles.mat.validTL.Incomp > 0)
  PowerrateTL.AMS = mean(allFiles.mat.validTL.AMS > 0)
  PowerrateTL.BinMM = mean(allFiles.mat.validTL.BinMM > 0)
  
  PowerrateTW.IBS = mean(allFiles.mat.validTW.IBS > 0)
  PowerrateTW.Incomp = mean(allFiles.mat.validTW.Incomp > 0)
  PowerrateTW.AMS = mean(allFiles.mat.validTW.AMS > 0)
  PowerrateTW.BinMM = mean(allFiles.mat.validTW.BinMM > 0)
  
  PowerrateTS.IBS = mean(allFiles.mat.validTS.IBS > 0)
  PowerrateTS.Incomp = mean(allFiles.mat.validTS.Incomp > 0)
  PowerrateTS.AMS = mean(allFiles.mat.validTS.AMS > 0)
  PowerrateTS.BinMM = mean(allFiles.mat.validTS.BinMM > 0)
  
  PowerRate.IBS = cbind(PowerrateTL.IBS, PowerrateTW.IBS, PowerrateTS.IBS)
  PowerRate.Incomp = cbind(PowerrateTL.Incomp, PowerrateTW.Incomp, PowerrateTS.Incomp)
  PowerRate.AMS = cbind(PowerrateTL.AMS, PowerrateTW.AMS, PowerrateTS.AMS)
  PowerRate.BinMM = cbind(PowerrateTL.BinMM, PowerrateTW.BinMM, PowerrateTS.BinMM)
  
  powerRates.mat[ii,1:3] = PowerRate.IBS
  powerRates.mat[ii,4:6] = PowerRate.Incomp
  powerRates.mat[ii,7:9] = PowerRate.AMS
  powerRates.mat[ii,10:12] = PowerRate.BinMM
}

rownames(powerRates.mat) = scores
colnames(powerRates.mat) = c("LRT, IBS", "Wald, IBS", "Score, IBS", "LRT, Incomp", "Wald, Incomp", "Score, Incomp",
                             "LRT, AMS", "Wald, AMS", "Score, AMS", "LRT, Bin MM", "Wald, Bin MM", "Score, Bin MM")

write.csv(powerRates.mat, paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_LHT__BetaZero.csv"))

### GLM ############################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_GLM_Power_BetaZero"))

scores = c("IBS", "Incomp", "AMS", "BinMM")

powerRatesGLM.mat = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_CatPhenos_Prev",prev,"_GLM_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs_Sim1to1000_StatsAndPValues.csv"),
                paste0("Power_CatPhenos_Prev",prev,"_GLM_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs_Sim1001to2000_StatsAndPValues.csv"),
                paste0("Power_CatPhenos_Prev",prev,"_GLM_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs_Sim2001to3000_StatsAndPValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  Powerrate.IBS = mean(files.mat$IBS.P.Val <= alpha)
  Powerrate.Incomp = mean(files.mat$Incomp.P.Val <= alpha)
  Powerrate.AMS = mean(files.mat$AMS.P.Val <= alpha)
  Powerrate.BinMM = mean(files.mat$BinMM.P.Val <= alpha)
  
  powerRatesGLM.mat[ii,1] = Powerrate.IBS
  powerRatesGLM.mat[ii,2] = Powerrate.Incomp
  powerRatesGLM.mat[ii,3] = Powerrate.AMS
  powerRatesGLM.mat[ii,4] = Powerrate.BinMM
}

rownames(powerRatesGLM.mat) = scores
colnames(powerRatesGLM.mat) = c("GLM, IBS","GLM, Incomp","GLM, AMS", "GLM, Bin MM")

write.csv(powerRatesGLM.mat, paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_GLM_BetaZero.csv"))

### Desparse Lasso #################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_DesparseLasso_Power_BetaZero"))

scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_CatY_Prev",prev,"_DesparseLasso_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1to1000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_DesparseLasso_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1001to2000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_DesparseLasso_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim2001to3000_PValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  PowerRate.IBS = mean(files.mat$IBS <= alpha)
  PowerRate.Incomp = mean(files.mat$Incomp <= alpha)
  PowerRate.AMS = mean(files.mat$AMS <= alpha)
  PowerRate.BinMM = mean(files.mat$BinMM <= alpha)
  
  powerRates.mat[ii,1] = PowerRate.IBS
  powerRates.mat[ii,2] = PowerRate.Incomp
  powerRates.mat[ii,3] = PowerRate.AMS
  powerRates.mat[ii,4] = PowerRate.BinMM
}

rownames(powerRates.mat) = scores
colnames(powerRates.mat) = c("DL, IBS","DL, Incomp","DL, AMS", "DL, Bin MM")

write.csv(powerRates.mat, paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DL_BetaZero.csv"))

### AISPU ###################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_AISPU_Power_BetaZero"))

# scores = c("AMS.gene", "BinMM.gene", "IBS.gene", "Incomp.gene") #for old code
scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat = matrix(NA, nrow = length(scores), ncol = 4)
powerRates.mat_inf = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_CatY_Prev",prev,"_AISPU_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs__Sim1to1000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_AISPU_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs__Sim1001to2000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_AISPU_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs__Sim2001to3000_PValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]

  files.mat = na.omit(files.mat)
  
  Powerrate.IBS = mean(files.mat$IBS <= alpha)
  Powerrate.Incomp = mean(files.mat$Incomp <= alpha)
  Powerrate.AMS = mean(files.mat$AMS <= alpha)
  Powerrate.BinMM = mean(files.mat$BinMM <= alpha)
  
  Powerrate_inf.IBS = mean(files.mat$IBS.1 <= alpha)
  Powerrate_inf.Incomp = mean(files.mat$Incomp.1 <= alpha)
  Powerrate_inf.AMS = mean(files.mat$AMS.1 <= alpha)
  Powerrate_inf.BinMM = mean(files.mat$BinMM.1 <= alpha)
  
  powerRates.mat[ii,1] = Powerrate.IBS
  powerRates.mat[ii,2] = Powerrate.Incomp
  powerRates.mat[ii,3] = Powerrate.AMS
  powerRates.mat[ii,4] = Powerrate.BinMM
  
  powerRates.mat_inf[ii,1] = Powerrate_inf.IBS
  powerRates.mat_inf[ii,2] = Powerrate_inf.Incomp
  powerRates.mat_inf[ii,3] = Powerrate_inf.AMS
  powerRates.mat_inf[ii,4] = Powerrate_inf.BinMM
}

#check that they are equal
powerRates.mat == powerRates.mat_inf

rownames(powerRates.mat_inf) = scores
colnames(powerRates.mat_inf) = c("AISPU, IBS", "AISPU, Incomp", "AISPU, AMS", "AISPU, Bin MM")

write.csv(powerRates.mat_inf, paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_AISPU_BetaZero.csv"))

### Decorr Score #############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_DecorrScore_Power_BetaZero"))

scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_CatY_Prev",prev,"_DecorrScore_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1to1000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_DecorrScore_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1001to2000_PValues.csv"),
                paste0("Power_CatY_Prev",prev,"_DecorrScore_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim2001to3000_PValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  Powerrate.IBS = mean(files.mat$IBS <= alpha)
  Powerrate.Incomp = mean(files.mat$Incomp <= alpha)
  Powerrate.AMS = mean(files.mat$AMS <= alpha)
  Powerrate.BinMM = mean(files.mat$BinMM <= alpha)
  
  powerRates.mat[ii,1] = Powerrate.IBS
  powerRates.mat[ii,2] = Powerrate.Incomp
  powerRates.mat[ii,3] = Powerrate.AMS
  powerRates.mat[ii,4] = Powerrate.BinMM
}

rownames(powerRates.mat) = scores
colnames(powerRates.mat) = c("DS, IBS","DS, Incomp","DS, AMS", "DS, Bin MM")

write.csv(powerRates.mat, paste0(gene,"_",ss,"_PowerResults_Prev",prev,"_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DecorrScore_BetaZero.csv"))

##### For Continuous ##################################################################################################
# inputs = c("NAT2", "1000", 20, "small", 25, TRUE) #for continuous
# inputs = c("CHI3L2", "1000", 20, "small", 25, TRUE) #for continuous
inputs = c("ASAH1", "1000", 20, "small", 15, TRUE) #for continuous

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

#### LHT ##########################################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_LHT_Power_BetaZero"))

scores = c("AMS.gene", "BinMM.gene", "IBS.gene", "Incomp.gene") #for old code
# scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat.cont = matrix(NA, nrow = length(scores), ncol = 12)

for(jj in 1:length(scores)){
  filenames = c(paste0("Power_ContY_LinHypTest_Score",scores[jj],"_GammaOnly_",percentAssoc,"SNPsAssoc_JointTesting_OR",orSize,"_",ld,"_Sim1to1000_Stats_ForceCovs.csv"),
                paste0("Power_ContY_LinHypTest_Score",scores[jj],"_GammaOnly_",percentAssoc,"SNPsAssoc_JointTesting_OR",orSize,"_",ld,"_Sim1001to2000_Stats_ForceCovs.csv"),
                paste0("Power_ContY_LinHypTest_Score",scores[jj],"_GammaOnly_",percentAssoc,"SNPsAssoc_JointTesting_OR",orSize,"_",ld,"_Sim2001to3000_Stats_ForceCovs.csv"))
  files = lapply(filenames, read.csv) 
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  invalidTL = which(files.mat$V1 < 0 | files.mat$V4 < 0 | files.mat$V7 < 0 | files.mat$V10 < 0)
  invalidTW = which(files.mat$V2 < 0 | files.mat$V5 < 0 | files.mat$V8 < 0 | files.mat$V11 < 0)
  invalidTS = which(files.mat$V3 < 0 | files.mat$V6 < 0 | files.mat$V9 < 0 | files.mat$V12 < 0)
  
  pValfiles = c(paste0("Power_ContY_LinHypTest_Score",scores[jj],"_GammaOnly_",percentAssoc,"SNPsAssoc_JointTesting_OR",orSize,"_",ld,"_Sim1to1000_PValues_ForceCovs.csv"),
                paste0("Power_ContY_LinHypTest_Score",scores[jj],"_GammaOnly_",percentAssoc,"SNPsAssoc_JointTesting_OR",orSize,"_",ld,"_Sim1001to2000_PValues_ForceCovs.csv"),
                paste0("Power_ContY_LinHypTest_Score",scores[jj],"_GammaOnly_",percentAssoc,"SNPsAssoc_JointTesting_OR",orSize,"_",ld,"_Sim2001to3000_PValues_ForceCovs.csv"))
  allFiles = lapply(pValfiles, read.csv)
  allFiles.mat = rbind(allFiles[[1]], allFiles[[2]], allFiles[[3]])
  #remove column 1
  allFiles.mat = allFiles.mat[,2:ncol(allFiles.mat)]
  
  if(length(invalidTL) > 0){
    allFiles.mat.validTL.IBS = allFiles.mat$V1[-invalidTL]
    allFiles.mat.validTL.Incomp = allFiles.mat$V4[-invalidTL]
    allFiles.mat.validTL.AMS = allFiles.mat$V7[-invalidTL]
    allFiles.mat.validTL.BinMM = allFiles.mat$V10[-invalidTL]
  } else {
    allFiles.mat.validTL.IBS = allFiles.mat$V1
    allFiles.mat.validTL.Incomp = allFiles.mat$V4
    allFiles.mat.validTL.AMS = allFiles.mat$V7
    allFiles.mat.validTL.BinMM = allFiles.mat$V10
  }
  
  if(length(invalidTW) > 0){
    allFiles.mat.validTW.IBS = allFiles.mat$V2[-invalidTW]
    allFiles.mat.validTW.Incomp = allFiles.mat$V5[-invalidTW]
    allFiles.mat.validTW.AMS = allFiles.mat$V8[-invalidTW]
    allFiles.mat.validTW.BinMM = allFiles.mat$V11[-invalidTW]
  } else {
    allFiles.mat.validTW.IBS = allFiles.mat$V2
    allFiles.mat.validTW.Incomp = allFiles.mat$V5
    allFiles.mat.validTW.AMS = allFiles.mat$V8
    allFiles.mat.validTW.BinMM = allFiles.mat$V11
  }
  
  if(length(invalidTS) > 0){
    allFiles.mat.validTS.IBS = allFiles.mat$V3[-invalidTS]
    allFiles.mat.validTS.Incomp = allFiles.mat$V6[-invalidTS]
    allFiles.mat.validTS.AMS = allFiles.mat$V9[-invalidTS]
    allFiles.mat.validTS.BinMM = allFiles.mat$V12[-invalidTS]
  } else {
    allFiles.mat.validTS.IBS = allFiles.mat$V3
    allFiles.mat.validTS.Incomp = allFiles.mat$V6
    allFiles.mat.validTS.AMS = allFiles.mat$V9
    allFiles.mat.validTS.BinMM = allFiles.mat$V12
  }
  
  PowerrateTL.IBS = mean(allFiles.mat.validTL.IBS > 0)
  PowerrateTL.Incomp = mean(allFiles.mat.validTL.Incomp > 0)
  PowerrateTL.AMS = mean(allFiles.mat.validTL.AMS > 0)
  PowerrateTL.BinMM = mean(allFiles.mat.validTL.BinMM > 0)
  
  PowerrateTW.IBS = mean(allFiles.mat.validTW.IBS > 0)
  PowerrateTW.Incomp = mean(allFiles.mat.validTW.Incomp > 0)
  PowerrateTW.AMS = mean(allFiles.mat.validTW.AMS > 0)
  PowerrateTW.BinMM = mean(allFiles.mat.validTW.BinMM > 0)
  
  PowerrateTS.IBS = mean(allFiles.mat.validTS.IBS > 0)
  PowerrateTS.Incomp = mean(allFiles.mat.validTS.Incomp > 0)
  PowerrateTS.AMS = mean(allFiles.mat.validTS.AMS > 0)
  PowerrateTS.BinMM = mean(allFiles.mat.validTS.BinMM > 0)
  
  PowerRate.IBS = cbind(PowerrateTL.IBS, PowerrateTW.IBS, PowerrateTS.IBS)
  PowerRate.Incomp = cbind(PowerrateTL.Incomp, PowerrateTW.Incomp, PowerrateTS.Incomp)
  PowerRate.AMS = cbind(PowerrateTL.AMS, PowerrateTW.AMS, PowerrateTS.AMS)
  PowerRate.BinMM = cbind(PowerrateTL.BinMM, PowerrateTW.BinMM, PowerrateTS.BinMM)
  
  powerRates.mat.cont[jj,1:3] = PowerRate.IBS
  powerRates.mat.cont[jj,4:6] = PowerRate.Incomp
  powerRates.mat.cont[jj,7:9] = PowerRate.AMS
  powerRates.mat.cont[jj,10:12] = PowerRate.BinMM
} 

rownames(powerRates.mat.cont) = scores
colnames(powerRates.mat.cont) = c("LRT, IBS", "Wald, IBS", "Score, IBS", "LRT, Incomp", "Wald, Incomp", "Score, Incomp",
                                  "LRT, AMS", "Wald, AMS", "Score, AMS", "LRT, Bin MM", "Wald, Bin MM", "Score, Bin MM")

write.csv(powerRates.mat.cont, paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_LHT_BetaZero.csv"))

### GLM ############################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_GLM_Power_BetaZero"))

scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRatesGLM.mat.cont = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_ContPhenos_Phenos_GLM_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs_Sim1to1000_StatsAndPValues.csv"),
                paste0("Power_ContPhenos_Phenos_GLM_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs_Sim1001to2000_StatsAndPValues.csv"),
                paste0("Power_ContPhenos_Phenos_GLM_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs_Sim2001to3000_StatsAndPValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  Powerrate.IBS = mean(files.mat$IBS.P.Val <= 0.05)
  Powerrate.Incomp = mean(files.mat$Incomp.P.Val <= 0.05)
  Powerrate.AMS = mean(files.mat$AMS.P.Val <= 0.05)
  Powerrate.BinMM = mean(files.mat$BinMM.P.Val <= 0.05)
  
  powerRatesGLM.mat.cont[ii,1] = Powerrate.IBS
  powerRatesGLM.mat.cont[ii,2] = Powerrate.Incomp
  powerRatesGLM.mat.cont[ii,3] = Powerrate.AMS
  powerRatesGLM.mat.cont[ii,4] = Powerrate.BinMM
}

rownames(powerRatesGLM.mat.cont) = scores
colnames(powerRatesGLM.mat.cont) = c("GLM, IBS","GLM, Incomp","GLM, AMS", "GLM, Bin MM")

write.csv(powerRatesGLM.mat.cont, paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_GLM_BetaZero.csv"))

### Desparse Lasso #################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_DesparseLasso_Power_BetaZero"))

scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_ContY_DesparseLasso_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1to1000_PValues.csv"),
                paste0("Power_ContY_DesparseLasso_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1001to2000_PValues.csv"),
                paste0("Power_ContY_DesparseLasso_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim2001to3000_PValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  PowerRate.IBS = mean(files.mat$IBS <= 0.05)
  PowerRate.Incomp = mean(files.mat$Incomp <= 0.05)
  PowerRate.AMS = mean(files.mat$AMS <= 0.05)
  PowerRate.BinMM = mean(files.mat$BinMM <= 0.05)
  
  powerRates.mat[ii,1] = PowerRate.IBS
  powerRates.mat[ii,2] = PowerRate.Incomp
  powerRates.mat[ii,3] = PowerRate.AMS
  powerRates.mat[ii,4] = PowerRate.BinMM
}

rownames(powerRates.mat) = scores
colnames(powerRates.mat) = c("DL, IBS","DL, Incomp","DL, AMS", "DL, Bin MM")

write.csv(powerRates.mat, paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DL_BetaZero.csv"))

### AISPU ###################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_AISPU_Power_BetaZero"))

# scores = c("AMS.gene", "BinMM.gene", "IBS.gene", "Incomp.gene") #for old code
scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat = matrix(NA, nrow = length(scores), ncol = 4)
powerRates.mat_inf = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_ContY_AISPU_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs__Sim1to1000_PValues.csv"),
                paste0("Power_ContY_AISPU_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs__Sim1001to2000_PValues.csv"),
                paste0("Power_ContY_AISPU_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"_GammaSig__WCovs__Sim2001to3000_PValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  Powerrate.IBS = mean(files.mat$IBS <= 0.05)
  Powerrate.Incomp = mean(files.mat$Incomp <= 0.05)
  Powerrate.AMS = mean(files.mat$AMS <= 0.05)
  Powerrate.BinMM = mean(files.mat$BinMM <= 0.05)
  
  Powerrate_inf.IBS = mean(files.mat$IBS.1 <= 0.05)
  Powerrate_inf.Incomp = mean(files.mat$Incomp.1 <= 0.05)
  Powerrate_inf.AMS = mean(files.mat$AMS.1 <= 0.05)
  Powerrate_inf.BinMM = mean(files.mat$BinMM.1 <= 0.05)
  
  powerRates.mat[ii,1] = Powerrate.IBS
  powerRates.mat[ii,2] = Powerrate.Incomp
  powerRates.mat[ii,3] = Powerrate.AMS
  powerRates.mat[ii,4] = Powerrate.BinMM
  
  powerRates.mat_inf[ii,1] = Powerrate_inf.IBS
  powerRates.mat_inf[ii,2] = Powerrate_inf.Incomp
  powerRates.mat_inf[ii,3] = Powerrate_inf.AMS
  powerRates.mat_inf[ii,4] = Powerrate_inf.BinMM
}

#check that they are equal
powerRates.mat == powerRates.mat_inf

rownames(powerRates.mat_inf) = scores
colnames(powerRates.mat_inf) = c("AISPU, IBS", "AISPU, Incomp", "AISPU, AMS", "AISPU, Bin MM")

write.csv(powerRates.mat_inf, paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_AISPU_BetaZero.csv"))

### Decorr Score #############################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_DecorrScore_Power_BetaZero"))

scores = c("AMS", "BinMM", "IBS", "Incomp")

powerRates.mat = matrix(NA, nrow = length(scores), ncol = 4)

for(ii in 1:length(scores)){
  filenames = c(paste0("Power_ContY_DecorrScore_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1to1000_PValues.csv"),
                paste0("Power_ContY_DecorrScore_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim1001to2000_PValues.csv"),
                paste0("Power_ContY_DecorrScore_Score",scores[ii],"_",orSize,"OR_",percentAssoc,"Assoc_LowLD",lowLD,"GammaSig_WCovs__Sim2001to3000_PValues.csv"))
  files = lapply(filenames, read.csv)
  files.mat = rbind(files[[1]], files[[2]], files[[3]])
  #remove column 1
  files.mat = files.mat[,2:ncol(files.mat)]
  
  Powerrate.IBS = mean(files.mat$IBS <= 0.05)
  Powerrate.Incomp = mean(files.mat$Incomp <= 0.05)
  Powerrate.AMS = mean(files.mat$AMS <= 0.05)
  Powerrate.BinMM = mean(files.mat$BinMM <= 0.05)
  
  powerRates.mat[ii,1] = Powerrate.IBS
  powerRates.mat[ii,2] = Powerrate.Incomp
  powerRates.mat[ii,3] = Powerrate.AMS
  powerRates.mat[ii,4] = Powerrate.BinMM
}

rownames(powerRates.mat) = scores
colnames(powerRates.mat) = c("DS, IBS","DS, Incomp","DS, AMS", "DS, Bin MM")

write.csv(powerRates.mat, paste0(gene,"_",ss,"_PowerResults_Cont_ORSize",orSize,"_",percentAssoc,"Assoc_",ld,"_DecorrScore_BetaZero.csv"))
