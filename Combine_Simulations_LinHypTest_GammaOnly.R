#############################################
## Combining the 5K Simulation Results
#############################################
### For GLM ###################################################
gene = "ASAH1"
# gene = "CHI3L2"
# gene = "NAT2"
# ss = 500
ss = 1000
nSims = 3000

Prev = c(20, 10, 5)
FitScore = c("IBS", "Incomp", "AMS", "BinMM")
covariates = c("NoCovs", "WCovs")

setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_GLM_TIE_GammaBetaZero"))
powerList = list()
powerByScore = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(ii in 1:length(Prev)){
  for(jj in 1:length(FitScore)){
    files_NoCovs = read.csv(paste0("TIE_CatPhenos_Prev",Prev[ii],"_GLM_Score",FitScore[jj],"__",covariates[1],"_Sim1to3000_StatsAndPValues.csv"))
    files_WCovs = read.csv(paste0("TIE_CatPhenos_Prev",Prev[ii],"_GLM_Score",FitScore[jj],"__",covariates[2],"_Sim1to5000_StatsAndPValues.csv"))
    
    files_NoCovs = files_NoCovs[1:nSims,2:3]
    files_WCovs = files_WCovs[1:nSims,2:3]
    
    power_NoCovs = sum(files_NoCovs$V1 < 0.05)/dim(files_NoCovs)[1]
    power_WCovs = sum(files_WCovs$V1 < 0.05)/dim(files_WCovs)[1]  
    
    powerByScore[1,jj] = power_NoCovs
    powerByScore[2,jj] = power_WCovs
  }
  powerList[[ii]] = powerByScore
}

powerByScoreCont = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(jj in 1:length(FitScore)){
  files_NoCovs_Cont = read.csv(paste0("TIE_ContPhenos_Phenos_GLM_Score",FitScore[jj],"__",covariates[1],"_Sim1to3000_StatsAndPValues.csv"))
  files_WCovs_Cont = read.csv(paste0("TIE_ContPhenos_Phenos_GLM_Score",FitScore[jj],"__",covariates[2],"_Sim1to5000_StatsAndPValues.csv"))
  
  files_NoCovs_Cont = files_NoCovs_Cont[1:nSims,2:3]
  files_WCovs_Cont = files_WCovs_Cont[1:nSims,2:3]
  
  power_NoCovs_Cont = sum(files_NoCovs_Cont$V1 < 0.05)/dim(files_NoCovs_Cont)[1]
  power_WCovs_Cont = sum(files_WCovs_Cont$V1 < 0.05)/dim(files_WCovs_Cont)[1]
  
  powerByScoreCont[1,jj] = power_NoCovs_Cont
  powerByScoreCont[2,jj] = power_WCovs_Cont
}  

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByScoreCont)
rownames(PowerData) = c("Prev 20 - NoCovs", "Prev 20 - WCovs", "Prev 10 - NoCovs", "Prev 10 - WCovs",
                        "Prev 5 - NoCovs", "Prev 5 - WCovs", "Cont. - NoCovs", "Cont. - WCovs")
colnames(PowerData) = c("IBS Fit", "Incomp. Fit", "AMS Fit", "Binary Mismatch Fit")
write.csv(PowerData, paste0(gene,"_",ss,"Pairs_GLM_TIE_",nSims,"Sims_GammaBetaZero.csv"))

#For LHT #################################################################
# gene = "ASAH1"
gene = "CHI3L2"
# gene = "NAT2"
ss = 500
# ss = 1000
nSims = 5000

Prev = c(20, 10, 5)
FitScore = c("IBS", "Incomp", "AMS", "BinMM")
covariates = c("NoCovs", "WCovs")

setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_LHT_TIE_GammaBetaZero"))
powerList = list()
powerByScore = matrix(data = NA, nrow = length(covariates)*3, ncol = length(FitScore))

for(ii in 1:length(Prev)){
  for(jj in 1:length(FitScore)){
    #read in stats
    statsFilesNoCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_StatsCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_StatsCV_ForceCovFit.csv"))
    
    statsFilesWCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_StatsCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_StatsCV_ForceCovFit.csv"))
    
    statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
    statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
    
    statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
    statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
    
    statFiles_NoCovs.all = statFiles_NoCovs.all[,2:ncol(statFiles_NoCovs.all)]
    statFiles_WCovs.all = statFiles_WCovs.all[,2:ncol(statFiles_WCovs.all)]
    
    # determine which stats were bad
    if(length(which(statFiles_NoCovs.all$V1 < 0)) > 0){
      badTests_TL.NoCovs = which(statFiles_NoCovs.all$V1 < 0)
    } else {
      badTests_TL.NoCovs = c(5001)
    }
    
    if(length(which(statFiles_WCovs.all$V1 < 0)) > 0){
      badTests_TL.WCovs = which(statFiles_WCovs.all$V1 < 0)
    } else {
      badTests_TL.WCovs = c(5001)
    }
    
    ##
    if(length(which(statFiles_NoCovs.all$V2 < 0)) > 0){
      badTests_TW.NoCovs = which(statFiles_NoCovs.all$V2 < 0)
    } else {
      badTests_TW.NoCovs = c(5001)
    }

    if(length(which(statFiles_WCovs.all$V2 < 0)) > 0){
      badTests_TW.WCovs = which(statFiles_WCovs.all$V2 < 0)
    } else {
      badTests_TW.WCovs = c(5001)
    }
    
    ##    
    if(length(which(statFiles_NoCovs.all$V3 < 0)) > 0){
      badTests_TS.NoCovs = which(statFiles_NoCovs.all$V3 < 0)
    } else {
      badTests_TS.NoCovs = c(5001)
    }
    
    if(length(which(statFiles_WCovs.all$V3 < 0)) > 0){
      badTests_TS.WCovs = which(statFiles_WCovs.all$V3 < 0)
    } else {
      badTests_TS.WCovs = c(5001)
    }

    #read in p-values
    pvalsFilesNoCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
    
    pvalsFilesWCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
    
    pvalsFiles_NoCovs = lapply(pvalsFilesNoCovs, read.csv)
    pvalsFiles_WCovs = lapply(pvalsFilesWCovs, read.csv)
    
    pvalsFiles_NoCovs.all = rbind(pvalsFiles_NoCovs[[1]][1:1000,], pvalsFiles_NoCovs[[2]][1:1000,], pvalsFiles_NoCovs[[3]][1:1000,])
    pvalsFiles_WCovs.all = rbind(pvalsFiles_WCovs[[1]][1:1000,], pvalsFiles_WCovs[[2]][1:1000,], pvalsFiles_WCovs[[3]][1:1000,])

    pvalsFiles_NoCovs.all = pvalsFiles_NoCovs.all[,2:ncol(pvalsFiles_NoCovs.all)]
    pvalsFiles_WCovs.all = pvalsFiles_WCovs.all[,2:ncol(pvalsFiles_WCovs.all)]
    
    TL_goodPvals.NoCovs = pvalsFiles_NoCovs.all$V1[-c(badTests_TL.NoCovs)]
    TW_goodPvals.NoCovs = pvalsFiles_NoCovs.all$V2[-c(badTests_TW.NoCovs)]
    TS_goodPvals.NoCovs = pvalsFiles_NoCovs.all$V3[-c(badTests_TS.NoCovs)]
    
    TL_goodPvals.WCovs = pvalsFiles_WCovs.all$V1[-c(badTests_TL.WCovs)]
    TW_goodPvals.WCovs = pvalsFiles_WCovs.all$V2[-c(badTests_TW.WCovs)]
    TS_goodPvals.WCovs = pvalsFiles_WCovs.all$V3[-c(badTests_TS.WCovs)]
    
    TL_TIE.NoCovs = mean(TL_goodPvals.NoCovs)*nSims
    TW_TIE.NoCovs = mean(TW_goodPvals.NoCovs)*nSims
    TS_TIE.NoCovs = mean(TS_goodPvals.NoCovs)*nSims

    TIE.NoCovs = rbind(round(TL_TIE.NoCovs, 2), round(TW_TIE.NoCovs,2), round(TS_TIE.NoCovs,2))
    
    TL_TIE.WCovs = mean(TL_goodPvals.WCovs)*nSims
    TW_TIE.WCovs = mean(TW_goodPvals.WCovs)*nSims
    TS_TIE.WCovs = mean(TS_goodPvals.WCovs)*nSims
    
    TIE.WCovs = rbind(round(TL_TIE.WCovs,2), round(TW_TIE.WCovs,2), round(TS_TIE.WCovs,2))
    
    All.TIE = rbind(TIE.NoCovs, TIE.WCovs)
    powerByScore[,jj] = All.TIE
  }
  powerList[[ii]] = powerByScore
}

powerByScoreCont = matrix(data = NA, nrow = length(covariates)*3, ncol = length(FitScore))

for(jj in 1:length(FitScore)){
  #read in stats
  statsFilesNoCovs = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_StatsCV_ForceCovFit.csv"),
                       paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                       paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_StatsCV_ForceCovFit.csv"))
  
  statsFilesWCovs = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_StatsCV_ForceCovFit.csv"),
                      paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                      paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_StatsCV_ForceCovFit.csv"))
  
  statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
  statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
  
  statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
  statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
  
  statFiles_NoCovs.all = statFiles_NoCovs.all[,2:ncol(statFiles_NoCovs.all)]
  statFiles_WCovs.all = statFiles_WCovs.all[,2:ncol(statFiles_WCovs.all)]
  
  # determine which stats were bad
  if(length(which(statFiles_NoCovs.all$V1 < 0)) > 0){
    badTests_TL.NoCovs = which(statFiles_NoCovs.all$V1 < 0)
  } else {
    badTests_TL.NoCovs = c(5001)
  }
  
  if(length(which(statFiles_WCovs.all$V1 < 0)) > 0){
    badTests_TL.WCovs = which(statFiles_WCovs.all$V1 < 0)
  } else {
    badTests_TL.WCovs = c(5001)
  }
  
  ##
  if(length(which(statFiles_NoCovs.all$V2 < 0)) > 0){
    badTests_TW.NoCovs = which(statFiles_NoCovs.all$V2 < 0)
  } else {
    badTests_TW.NoCovs = c(5001)
  }
  
  if(length(which(statFiles_WCovs.all$V2 < 0)) > 0){
    badTests_TW.WCovs = which(statFiles_WCovs.all$V2 < 0)
  } else {
    badTests_TW.WCovs = c(5001)
  }
  
  ##    
  if(length(which(statFiles_NoCovs.all$V3 < 0)) > 0){
    badTests_TS.NoCovs = which(statFiles_NoCovs.all$V3 < 0)
  } else {
    badTests_TS.NoCovs = c(5001)
  }
  
  if(length(which(statFiles_WCovs.all$V3 < 0)) > 0){
    badTests_TS.WCovs = which(statFiles_WCovs.all$V3 < 0)
  } else {
    badTests_TS.WCovs = c(5001)
  }
  
  #read in p-values
  pvalsFilesNoCovs = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                       paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                       paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
  
  pvalsFilesWCovs = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                      paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                      paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
  
  pvalsFiles_NoCovs = lapply(pvalsFilesNoCovs, read.csv)
  pvalsFiles_WCovs = lapply(pvalsFilesWCovs, read.csv)
  
  pvalsFiles_NoCovs.all = rbind(pvalsFiles_NoCovs[[1]][1:1000,], pvalsFiles_NoCovs[[2]][1:1000,], pvalsFiles_NoCovs[[3]][1:1000,])
  pvalsFiles_WCovs.all = rbind(pvalsFiles_WCovs[[1]][1:1000,], pvalsFiles_WCovs[[2]][1:1000,], pvalsFiles_WCovs[[3]][1:1000,])
  
  pvalsFiles_NoCovs.all = pvalsFiles_NoCovs.all[,2:ncol(pvalsFiles_NoCovs.all)]
  pvalsFiles_WCovs.all = pvalsFiles_WCovs.all[,2:ncol(pvalsFiles_WCovs.all)]
  
  TL_goodPvals.NoCovs = pvalsFiles_NoCovs.all$V1[-c(badTests_TL.NoCovs)]
  TW_goodPvals.NoCovs = pvalsFiles_NoCovs.all$V2[-c(badTests_TW.NoCovs)]
  TS_goodPvals.NoCovs = pvalsFiles_NoCovs.all$V3[-c(badTests_TS.NoCovs)]
  
  TL_goodPvals.WCovs = pvalsFiles_WCovs.all$V1[-c(badTests_TL.WCovs)]
  TW_goodPvals.WCovs = pvalsFiles_WCovs.all$V2[-c(badTests_TW.WCovs)]
  TS_goodPvals.WCovs = pvalsFiles_WCovs.all$V3[-c(badTests_TS.WCovs)]
  
  TL_TIE.NoCovs = mean(TL_goodPvals.NoCovs)*nSims
  TW_TIE.NoCovs = mean(TW_goodPvals.NoCovs)*nSims
  TS_TIE.NoCovs = mean(TS_goodPvals.NoCovs)*nSims
  
  TIE.NoCovs = rbind(round(TL_TIE.NoCovs, 2), round(TW_TIE.NoCovs,2), round(TS_TIE.NoCovs,2))
  
  TL_TIE.WCovs = mean(TL_goodPvals.WCovs)*nSims
  TW_TIE.WCovs = mean(TW_goodPvals.WCovs)*nSims
  TS_TIE.WCovs = mean(TS_goodPvals.WCovs)*nSims
  
  TIE.WCovs = rbind(round(TL_TIE.WCovs,2), round(TW_TIE.WCovs,2), round(TS_TIE.WCovs,2))
  
  All.TIE = rbind(TIE.NoCovs, TIE.WCovs)
  powerByScoreCont[,jj] = All.TIE
}  

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByScoreCont)
rownames(PowerData) = c("TL - Prev 20 - NoCovs", "TW - Prev 20 - NoCovs", "TS - Prev 20 - NoCovs", 
                        "TL - Prev 20 - WCovs", "TW - Prev 20 - WCovs", "TS - Prev 20 - WCovs",
                        "TL - Prev 10 - NoCovs", "TW - Prev 10 - NoCovs", "TS - Prev 10 - NoCovs",
                        "TL - Prev 10 - WCovs", "TW - Prev 10 - WCovs", "TS - Prev 10 - WCovs",
                        "TL - Prev 5 - NoCovs", "TW - Prev 5 - NoCovs",  "TS - Prev 5 - NoCovs", 
                        "TL - Prev 5 - WCovs", "TW - Prev 5 - WCovs", "TS - Prev 5 - WCovs", 
                        "TL - Cont. - NoCovs", "TW - Cont. - NoCovs", "TS - Cont. - NoCovs", 
                        "TL - Cont. - WCovs", "TW - Cont. - WCovs", "TS - Cont. - WCovs")
colnames(PowerData) = c("IBS Fit", "Incomp. Fit", "AMS Fit", "Binary Mismatch Fit")
write.csv(PowerData, paste0(gene,"_",ss,"Pairs_LHT_TIE_3KSims_GammaBetaZero.csv"))
## For Desparsified Lasso ########################################################################
# gene = "ASAH1"
# gene = "CHI3L2"
gene = "NAT2"
# ss = 500
ss = 1000
nSims = 3000

Prev = c(20, 10, 5)
FitScore = c("IBS", "Incomp", "AMS", "BinMM")
covariates = c("NoCovs", "WCovs")

setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/All_TIE_GammaBetaZero/",gene,"_",ss,"_All_TIE_GammaBetaZero/",gene,"_",ss,"_DesparseLasso_TIE_GammaBetaZero"))
powerList = list()
powerByScore = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(ii in 1:length(Prev)){
  for(jj in 1:length(FitScore)){
    statsFilesNoCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValues.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValues.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValues.csv"))
    
    statsFilesWCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValues.csv"),
                        paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValues.csv"),
                        paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValues.csv"))
    
    statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
    statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
    
    statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
    statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
    
    power_NoCovs = sum(statFiles_NoCovs.all$V1 < 0.05)/dim(statFiles_NoCovs.all)[1]
    power_WCovs = sum(statFiles_WCovs.all$V1 < 0.05)/dim(statFiles_WCovs.all)[1]  
    
    powerByScore[1,jj] = power_NoCovs
    powerByScore[2,jj] = power_WCovs
  }
  powerList[[ii]] = powerByScore
}

powerByScoreCont = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(jj in 1:length(FitScore)){
  statsFilesNoCovs = c(paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValues.csv"),
                       paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValues.csv"),
                       paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValues.csv"))
  
  statsFilesWCovs = c(paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValues.csv"),
                      paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValues.csv"),
                      paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValues.csv"))
  
  statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
  statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
  
  statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
  statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
  
  power_NoCovs_Cont = sum(statFiles_NoCovs.all$V1 < 0.05)/dim(statFiles_NoCovs.all)[1]
  power_WCovs_Cont = sum(statFiles_WCovs.all$V1 < 0.05)/dim(statFiles_WCovs.all)[1]
  
  powerByScoreCont[1,jj] = power_NoCovs_Cont
  powerByScoreCont[2,jj] = power_WCovs_Cont
}  

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByScoreCont)
rownames(PowerData) = c("Prev 20 - NoCovs", "Prev 20 - WCovs", "Prev 10 - NoCovs", "Prev 10 - WCovs",
                        "Prev 5 - NoCovs", "Prev 5 - WCovs", "Cont. - NoCovs", "Cont. - WCovs")
colnames(PowerData) = c("IBS Fit", "Incomp. Fit", "AMS Fit", "Binary Mismatch Fit")
PowerData = round(PowerData, digits = 4)
write.csv(PowerData, paste0(gene,"_",ss,"Pairs_DesparseLasso_TIE_",nSims,"Sims_GammaBetaZero.csv"))

## For AISPU ########################################################################
gene = "ASAH1"
# gene = "CHI3L2"
# gene = "NAT2"
ss = 500
# ss = 1000
nSims = 3000
#what penalty was used?
# penalty = "tlp"
penalty = "scad"

Prev = c(20, 10, 5)
FitScore = c("IBS", "Incomp", "AMS", "BinMM")
covariates = c("NoCovs", "WCovs")

setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_AISPU_TIE_GammaBetaZero"))
powerList = list()
powerList_inf = list()
powerByScore = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))
powerByScore_inf = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(ii in 1:length(Prev)){
  for(jj in 1:length(FitScore)){
    statsFilesNoCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValues.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValues.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValues.csv"))
    
    statsFilesWCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValues.csv"),
                        paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValues.csv"),
                        paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValues.csv"))
    
    statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
    statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
    
    statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
    statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
    
    statFiles_NoCovs.all = na.omit(statFiles_NoCovs.all)
    statFiles_WCovs.all = na.omit(statFiles_WCovs.all)
    
    power_NoCovs = sum(statFiles_NoCovs.all$V1 < 0.05)/dim(statFiles_NoCovs.all)[1]
    power_WCovs = sum(statFiles_WCovs.all$V1 < 0.05)/dim(statFiles_WCovs.all)[1]  
    
    power_NoCovs_inf = sum(statFiles_NoCovs.all$V2 < 0.05)/dim(statFiles_NoCovs.all)[1]
    power_WCovs_inf = sum(statFiles_WCovs.all$V2 < 0.05)/dim(statFiles_WCovs.all)[1]  
    
    powerByScore[1,jj] = power_NoCovs
    powerByScore[2,jj] = power_WCovs
    
    powerByScore_inf[1,jj] = power_NoCovs_inf
    powerByScore_inf[2,jj] = power_WCovs_inf
  }
  powerList[[ii]] = powerByScore
  powerList_inf[[ii]] = powerByScore_inf
}

powerByScoreCont = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))
powerByScoreCont_inf = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(jj in 1:length(FitScore)){
  statsFilesNoCovs = c(paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValues.csv"),
                       paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValues.csv"),
                       paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValues.csv"))
  
  statsFilesWCovs = c(paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValues.csv"),
                      paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValues.csv"),
                      paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValues.csv"))
  
  statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
  statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
  
  statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
  statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
  
  power_NoCovs_Cont = sum(statFiles_NoCovs.all$V1 < 0.05)/dim(statFiles_NoCovs.all)[1]
  power_WCovs_Cont = sum(statFiles_WCovs.all$V1 < 0.05)/dim(statFiles_WCovs.all)[1]
  
  power_NoCovs_Cont_inf = sum(statFiles_NoCovs.all$V2 < 0.05)/dim(statFiles_NoCovs.all)[1]
  power_WCovs_Cont_inf = sum(statFiles_WCovs.all$V2 < 0.05)/dim(statFiles_WCovs.all)[1]
  
  powerByScoreCont[1,jj] = power_NoCovs_Cont
  powerByScoreCont[2,jj] = power_WCovs_Cont
  
  powerByScoreCont_inf[1,jj] = power_NoCovs_Cont_inf
  powerByScoreCont_inf[2,jj] = power_WCovs_Cont_inf
}  

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByScoreCont)
PowerData_inf = rbind(powerList_inf[[1]], powerList_inf[[2]], powerList_inf[[3]], powerByScoreCont_inf)

#see if they equal:
PowerData == PowerData_inf

rownames(PowerData) = c("Prev 20 - NoCovs", "Prev 20 - WCovs", "Prev 10 - NoCovs", "Prev 10 - WCovs",
                        "Prev 5 - NoCovs", "Prev 5 - WCovs", "Cont. - NoCovs", "Cont. - WCovs")
colnames(PowerData) = c("IBS Fit", "Incomp. Fit", "AMS Fit", "Binary Mismatch Fit")
PowerData = round(PowerData, digits = 2)
write.csv(PowerData, paste0(gene,"_",ss,"Pairs_AISPU_TIE_",nSims,"Sims_",penalty,"Penalty_GammaBetaZero.csv"))

## Decorrelated Score Test ########################################################################
# gene = "ASAH1"
gene = "CHI3L2"
# gene = "NAT2"
# ss = 500
ss = 1000
nSims = 3000

Prev = c(20, 10, 5)
FitScore = c("IBS", "Incomp", "AMS", "BinMM")
covariates = c("NoCovs", "WCovs")

setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_DecorrScore_TIE_GammaBetaZero"))
powerList = list()
powerByScore = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(ii in 1:length(Prev)){
  for(jj in 1:length(FitScore)){
    statsFilesNoCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValues.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValues.csv"),
                         paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValues.csv"))
    
    statsFilesWCovs = c(paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValues.csv"),
                        paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValues.csv"),
                        paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValues.csv"))
    
    statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
    statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
    
    statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
    statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
    
    statFiles_NoCovs.all = na.omit(statFiles_NoCovs.all)
    statFiles_WCovs.all = na.omit(statFiles_WCovs.all)
    
    power_NoCovs = sum(statFiles_NoCovs.all$V1 < 0.05)/dim(statFiles_NoCovs.all)[1]
    power_WCovs = sum(statFiles_WCovs.all$V1 < 0.05)/dim(statFiles_WCovs.all)[1]  
    
    powerByScore[1,jj] = power_NoCovs
    powerByScore[2,jj] = power_WCovs
  }
  powerList[[ii]] = powerByScore
}

powerByScoreCont = matrix(data = NA, nrow = length(covariates), ncol = length(FitScore))

for(jj in 1:length(FitScore)){
  statsFilesNoCovs = c(paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1to1000_PValues.csv"),
                       paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim1001to2000_PValues.csv"),
                       paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[1],"_Sim2001to3000_PValues.csv"))
  
  statsFilesWCovs = c(paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1to1000_PValues.csv"),
                      paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim1001to2000_PValues.csv"),
                      paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_GammaOnly_",covariates[2],"_Sim2001to3000_PValues.csv"))
  
  statFiles_NoCovs = lapply(statsFilesNoCovs, read.csv)
  statFiles_WCovs = lapply(statsFilesWCovs, read.csv)
  
  statFiles_NoCovs.all = rbind(statFiles_NoCovs[[1]][1:1000,], statFiles_NoCovs[[2]][1:1000,], statFiles_NoCovs[[3]][1:1000,])
  statFiles_WCovs.all = rbind(statFiles_WCovs[[1]][1:1000,], statFiles_WCovs[[2]][1:1000,], statFiles_WCovs[[3]][1:1000,])
  
  power_NoCovs_Cont = sum(statFiles_NoCovs.all$V1 < 0.05)/dim(statFiles_NoCovs.all)[1]
  power_WCovs_Cont = sum(statFiles_WCovs.all$V1 < 0.05)/dim(statFiles_WCovs.all)[1]
  
  powerByScoreCont[1,jj] = power_NoCovs_Cont
  powerByScoreCont[2,jj] = power_WCovs_Cont
}  

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByScoreCont)
rownames(PowerData) = c("Prev 20 - NoCovs", "Prev 20 - WCovs", "Prev 10 - NoCovs", "Prev 10 - WCovs",
                        "Prev 5 - NoCovs", "Prev 5 - WCovs", "Cont. - NoCovs", "Cont. - WCovs")
colnames(PowerData) = c("IBS Fit", "Incomp. Fit", "AMS Fit", "Binary Mismatch Fit")
PowerData = round(PowerData, digits = 2)
write.csv(PowerData, paste0(gene,"_",ss,"Pairs_DecorrScore_TIE_",nSims,"Sims_GammaBetaZero.csv"))

