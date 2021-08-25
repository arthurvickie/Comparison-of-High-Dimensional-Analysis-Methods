#############################################
## Combining the 5K Simulation Results
#############################################
####  For GLM ##################################################################################################################
## Gamma = 0 Only Testing, When Beta is Sig
# gene = "NAT2"
# gene = "CHI3L2"
gene = "ASAH1"
# ss = 500
ss = 1000

nSims = 3000

Prev = c(20, 10, 5)
FitScore = c("IBS")

ORSize = c("small", "medium", "large")
percentageAssoc = c(5,15,25)
LowLD = c(TRUE, FALSE)

path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_GLM_TIE_GammaZero_BetaAssoc")
setwd(path)

powerList = list()
powerByPercentAssoc = list()
powerByLD = matrix(data = NA, nrow = length(ORSize), ncol = length(LowLD))

for(jj in 1:length(FitScore)){
  for(ii in 1:length(Prev)){
    for(kk in 1:length(percentageAssoc)){
      for(aa in 1:length(LowLD)){     
        files_small = read.csv(paste0("TIE_CatPhenos_Prev",Prev[ii],"_GLM_Score",FitScore[jj],"_",ORSize[1],"OR_",percentageAssoc[kk],"Assoc_LowLD",LowLD[aa],"_BetaSig__WCovs_Sim1to3000_StatsAndPValues.csv"))
        files_med = read.csv(paste0("TIE_CatPhenos_Prev",Prev[ii],"_GLM_Score",FitScore[jj],"_",ORSize[2],"OR_",percentageAssoc[kk],"Assoc_LowLD",LowLD[aa],"_BetaSig__WCovs_Sim1to3000_StatsAndPValues.csv"))
        files_large = read.csv(paste0("TIE_CatPhenos_Prev",Prev[ii],"_GLM_Score",FitScore[jj],"_",ORSize[3],"OR_",percentageAssoc[kk],"Assoc_LowLD",LowLD[aa],"_BetaSig__WCovs_Sim1to3000_StatsAndPValues.csv"))
        
        files_small = files_small[1:nSims,2:3]
        files_med = files_med[1:nSims,2:3]
        files_large = files_large[1:nSims,2:3]
        
        tie_files_small = sum(files_small$V1 < 0.05)/dim(files_small)[1]
        tie_files_med = sum(files_med$V1 < 0.05)/dim(files_med)[1]
        tie_files_large = sum(files_large$V1 < 0.05)/dim(files_large)[1]
        
        powerByLD[1,aa] = tie_files_small
        powerByLD[2,aa] = tie_files_med
        powerByLD[3,aa] = tie_files_large
      }
      powerByPercentAssoc[[kk]] =  powerByLD
    }
    powerByPercentAssoc.mat = do.call(cbind, powerByPercentAssoc)
    powerList[[ii]] = powerByPercentAssoc.mat
  }
}

powerByPercentAssocCont = list()
powerByLDCont = matrix(data = NA, nrow = length(ORSize), ncol = length(LowLD))

for(jj in 1:length(FitScore)){
  for(kk in 1:length(percentageAssoc)){
    for(aa in 1:length(LowLD)){
      files_small_Cont = read.csv(paste0("TIE_ContPhenos_Phenos_GLM_Score",FitScore[jj],"_",ORSize[1],"OR_",percentageAssoc[kk],"Assoc_LowLD",LowLD[aa],"_BetaSig__WCovs_Sim1to3000_StatsAndPValues.csv"))
      files_medium_Cont = read.csv(paste0("TIE_ContPhenos_Phenos_GLM_Score",FitScore[jj],"_",ORSize[2],"OR_",percentageAssoc[kk],"Assoc_LowLD",LowLD[aa],"_BetaSig__WCovs_Sim1to3000_StatsAndPValues.csv"))
      files_large_Cont = read.csv(paste0("TIE_ContPhenos_Phenos_GLM_Score",FitScore[jj],"_",ORSize[3],"OR_",percentageAssoc[kk],"Assoc_LowLD",LowLD[aa],"_BetaSig__WCovs_Sim1to3000_StatsAndPValues.csv"))
      
      files_small_Cont = files_small_Cont[1:nSims,2:3]
      files_medium_Cont = files_medium_Cont[1:nSims,2:3]
      files_large_Cont = files_large_Cont[1:nSims,2:3]
      
      power_small_Cont = sum(files_small_Cont$V1 < 0.05)/dim(files_small_Cont)[1]
      power_med_Cont = sum(files_medium_Cont$V1 < 0.05)/dim(files_medium_Cont)[1]
      power_large_Cont = sum(files_large_Cont$V1 < 0.05)/dim(files_large_Cont)[1]
      
      powerByLDCont[1,aa] = power_small_Cont
      powerByLDCont[2,aa] = power_med_Cont
      powerByLDCont[3,aa] = power_large_Cont
    }
    powerByPercentAssocCont[[kk]] = powerByLD
  } 
  powerByPercentAssocCont.mat = do.call(cbind, powerByPercentAssocCont)
}

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByPercentAssocCont.mat)
rownames(PowerData) = c("Prev 20 - Small OR", "Prev 20 - Medium OR", "Prev 20 - Large OR", 
                        "Prev 10 - Small OR", "Prev 10 - Medium OR", "Prev 10 - Large OR", 
                        "Prev 5 - Small OR", "Prev 5 - Medium OR", "Prev 5 - Large OR", 
                        "Cont. - Small OR", "Cont. - Medium OR", "Cont. - Large OR")
colnames(PowerData) = c("5% Assoc - Low LD", "5% Assoc - High LD",
                        "15% Assoc - Low LD", "15% Assoc - High LD",
                        "25% Assoc - Low LD", "25% Assoc - High LD")
write.csv(PowerData, paste0(gene,"_",ss,"Pairs_GLM_TIE_",nSims,"Sims_GammaZero_BetaAssoc.csv"))

# For LHT ##################################################################################################################
## Gamma = 0 Only Testing, When Beta is Sig
# gene = "NAT2"
# gene = "CHI3L2"
gene = "ASAH1"
# ss = 500
ss = 1000

nSims = 3000

Prev = c(20, 10, 5)
FitScore = c("IBS")

OR = c("small", "medium", "large")
assocSNP = c(5,15,25)
lowLD = c(TRUE, FALSE)

path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_LHT_TIE_GammaZero_BetaAssoc")
setwd(path)

powerList = list()
powerByPercentAssoc = list()
powerByLD = matrix(data = NA, nrow = length(OR)*3, ncol = length(lowLD))

for(jj in 1:length(FitScore)){
  for(ii in 1:length(Prev)){
    for(kk in 1:length(assocSNP)){
      for(aa in 1:length(lowLD)){   
        #read in stats  
        statsFilesSmall = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[1],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1to1000_StatsCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[1],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[1],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim2001to3000_StatsCV_ForceCovFit.csv"))

        statsFilesMedium = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[2],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1to1000_StatsCV_ForceCovFit.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[2],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[2],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim2001to3000_StatsCV_ForceCovFit.csv"))

        statsFilesLarge = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[3],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1to1000_StatsCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[3],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[3],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim2001to3000_StatsCV_ForceCovFit.csv"))

        # statsFilesSmall = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_StatsCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_StatsCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_StatsCV_ForceCovFit.csv"))
        # 
        # statsFilesMedium = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_StatsCV_ForceCovFit.csv"),
        #                      paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_StatsCV_ForceCovFit.csv"),
        #                      paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_StatsCV_ForceCovFit.csv"))
        # 
        # statsFilesLarge = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_StatsCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_StatsCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_StatsCV_ForceCovFit.csv"))

        statFiles_Small = lapply(statsFilesSmall, read.csv)
        statFiles_Med = lapply(statsFilesMedium, read.csv)
        statFiles_Large = lapply(statsFilesLarge, read.csv)
        
        statFiles_Small.all = rbind(statFiles_Small[[1]][1:1000,], statFiles_Small[[2]][1:1000,], statFiles_Small[[3]][1:1000,])
        statFiles_Medium.all = rbind(statFiles_Med[[1]][1:1000,], statFiles_Med[[2]][1:1000,], statFiles_Med[[3]][1:1000,])
        statFiles_Large.all = rbind(statFiles_Large[[1]][1:1000,], statFiles_Large[[2]][1:1000,], statFiles_Large[[3]][1:1000,])
        
        statFiles_Small.all = statFiles_Small.all[,2:ncol(statFiles_Small.all)]
        statFiles_Medium.all = statFiles_Medium.all[,2:ncol(statFiles_Medium.all)]
        statFiles_Large.all = statFiles_Large.all[,2:ncol(statFiles_Large.all)]
        
        # determine which stats were bad
        if(length(which(statFiles_Small.all$V1 < 0)) > 0){
          badTests_TL.small = which(statFiles_Small.all$V1 < 0)
        } else {
          badTests_TL.small = c(5001)
        }
        if(length(which(statFiles_Medium.all$V1 < 0)) > 0){
          badTests_TL.medium = which(statFiles_Medium.all$V1 < 0)
        } else {
          badTests_TL.medium = c(5001)
        }
        if(length(which(statFiles_Large.all$V1 < 0)) > 0){
          badTests_TL.large = which(statFiles_Large.all$V1 < 0)
        } else {
          badTests_TL.large = c(5001)
        }
        
        ##
        if(length(which(statFiles_Small.all$V2 < 0)) > 0){
          badTests_TW.small = which(statFiles_Small.all$V2 < 0)
        } else {
          badTests_TW.small = c(5001)
        }
        if(length(which(statFiles_Medium.all$V2 < 0)) > 0){
          badTests_TW.medium = which(statFiles_Medium.all$V2 < 0)
        } else {
          badTests_TW.medium = c(5001)
        }
        if(length(which(statFiles_Large.all$V2 < 0)) > 0){
          badTests_TW.large = which(statFiles_Large.all$V2 < 0)
        } else {
          badTests_TW.large = c(5001)
        }
        
        ##    
        if(length(which(statFiles_Small.all$V3 < 0)) > 0){
          badTests_TS.small = which(statFiles_Small.all$V3 < 0)
        } else {
          badTests_TS.small = c(5001)
        }
        if(length(which(statFiles_Medium.all$V3 < 0)) > 0){
          badTests_TS.medium = which(statFiles_Medium.all$V3 < 0)
        } else {
          badTests_TS.medium = c(5001)
        }
        if(length(which(statFiles_Large.all$V3 < 0)) > 0){
          badTests_TS.large = which(statFiles_Large.all$V3 < 0)
        } else {
          badTests_TS.large = c(5001)
        }
        
        #read in p-values
        pvalsFilesSmall = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[1],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[1],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[1],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim2001to3000_PValuesCV_ForceCovFit.csv"))

        pvalsFilesMedium = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[2],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[2],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[2],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim2001to3000_PValuesCV_ForceCovFit.csv"))

        pvalsFilesLarge = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[3],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[3],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_GammaOnly_BetaSig_",OR[3],"OR_",assocSNP[kk],"AssocSNPs_",lowLD[aa],"_Sim2001to3000_PValuesCV_ForceCovFit.csv"))

        # pvalsFilesSmall = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValuesCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
        # 
        # pvalsFilesMedium = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValuesCV_ForceCovFit.csv"),
        #                      paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
        #                      paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
        # 
        # pvalsFilesLarge = c(paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValuesCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
        #                     paste0("TIE_CatY_Prev",Prev[ii],"_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValuesCV_ForceCovFit.csv"))

        pvalsFiles_small = lapply(pvalsFilesSmall, read.csv)
        pvalsFiles_medium = lapply(pvalsFilesMedium, read.csv)
        pvalsFiles_large = lapply(pvalsFilesLarge, read.csv)
        
        pvalsFiles_small.all = rbind(pvalsFiles_small[[1]][1:1000,], pvalsFiles_small[[2]][1:1000,], pvalsFiles_small[[3]][1:1000,])
        pvalsFiles_medium.all = rbind(pvalsFiles_medium[[1]][1:1000,], pvalsFiles_medium[[2]][1:1000,], pvalsFiles_medium[[3]][1:1000,])
        pvalsFiles_large.all = rbind(pvalsFiles_large[[1]][1:1000,], pvalsFiles_large[[2]][1:1000,], pvalsFiles_large[[3]][1:1000,])
        
        pvalsFiles_small.all = pvalsFiles_small.all[,2:ncol(pvalsFiles_small.all)]
        pvalsFiles_medium.all = pvalsFiles_medium.all[,2:ncol(pvalsFiles_medium.all)]
        pvalsFiles_large.all = pvalsFiles_large.all[,2:ncol(pvalsFiles_large.all)]
        
        TL_goodPvals.Small = pvalsFiles_small.all$V1[-c(badTests_TL.small)]
        TW_goodPvals.Small = pvalsFiles_small.all$V2[-c(badTests_TW.small)]
        TS_goodPvals.Small = pvalsFiles_small.all$V3[-c(badTests_TS.small)]
        
        TL_goodPvals.Medium = pvalsFiles_medium.all$V1[-c(badTests_TL.medium)]
        TW_goodPvals.Medium = pvalsFiles_medium.all$V2[-c(badTests_TW.medium)]
        TS_goodPvals.Medium = pvalsFiles_medium.all$V3[-c(badTests_TS.medium)]
        
        TL_goodPvals.Large = pvalsFiles_large.all$V1[-c(badTests_TL.large)]
        TW_goodPvals.Large = pvalsFiles_large.all$V2[-c(badTests_TW.large)]
        TS_goodPvals.Large = pvalsFiles_large.all$V3[-c(badTests_TS.large)]
        
        TL_TIE.small = sum(TL_goodPvals.Small > 0)/length(TL_goodPvals.Small)
        TW_TIE.small = sum(TW_goodPvals.Small > 0)/length(TW_goodPvals.Small)
        TS_TIE.small = sum(TS_goodPvals.Small > 0)/length(TS_goodPvals.Small)
        
        TIE.Small = rbind(round(TL_TIE.small, 2), round(TW_TIE.small,2), round(TS_TIE.small,2))
        
        TL_TIE.Medium = sum(TL_goodPvals.Medium > 0)/length(TL_goodPvals.Medium)
        TW_TIE.Medium = sum(TW_goodPvals.Medium > 0)/length(TW_goodPvals.Medium)
        TS_TIE.Medium = sum(TS_goodPvals.Medium > 0)/length(TS_goodPvals.Medium)
        
        TIE.Medium = rbind(round(TL_TIE.Medium, 2), round(TW_TIE.Medium,2), round(TS_TIE.Medium,2))
        
        TL_TIE.Large = sum(TL_goodPvals.Large > 0)/length(TL_goodPvals.Large)
        TW_TIE.Large = sum(TW_goodPvals.Large > 0)/length(TW_goodPvals.Large)
        TS_TIE.Large = sum(TS_goodPvals.Large > 0)/length(TS_goodPvals.Large)
        
        TIE.Large = rbind(round(TL_TIE.Large,2), round(TW_TIE.Large,2), round(TS_TIE.Large,2))
        
        powerByLD[1:3,aa] = TIE.Small
        powerByLD[4:6,aa] = TIE.Medium
        powerByLD[7:9,aa] = TIE.Large
      }
      powerByPercentAssoc[[kk]] =  powerByLD
    }
    powerByPercentAssoc.mat = do.call(cbind, powerByPercentAssoc)
    powerList[[ii]] = powerByPercentAssoc.mat
  }
}

powerByPercentAssocCont = list()
powerByLDCont = matrix(data = NA, nrow = length(OR)*3, ncol = length(lowLD))

for(jj in 1:length(FitScore)){
  for(kk in 1:length(assocSNP)){
    for(aa in 1:length(lowLD)){
      #read in stats                
      statsFilesSmall = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_StatsCV_ForceCovFit.csv"),
                           paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                           paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_StatsCV_ForceCovFit.csv"))
      
      statsFilesMedium = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_StatsCV_ForceCovFit.csv"),
                          paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                          paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_StatsCV_ForceCovFit.csv"))
      
      statsFilesLarge = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_StatsCV_ForceCovFit.csv"),
                           paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_StatsCV_ForceCovFit.csv"),
                           paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_StatsCV_ForceCovFit.csv"))
      
      statFiles_small = lapply(statsFilesSmall, read.csv)
      statFiles_medium = lapply(statsFilesMedium, read.csv)
      statFiles_large = lapply(statsFilesLarge, read.csv)
      
      statFiles_small.all = rbind(statFiles_small[[1]][1:1000,], statFiles_small[[2]][1:1000,], statFiles_small[[3]][1:1000,])
      statFiles_medium.all = rbind(statFiles_medium[[1]][1:1000,], statFiles_medium[[2]][1:1000,], statFiles_medium[[3]][1:1000,])
      statFiles_large.all = rbind(statFiles_large[[1]][1:1000,], statFiles_large[[2]][1:1000,], statFiles_large[[3]][1:1000,])
      
      statFiles_small.all = statFiles_small.all[,2:ncol(statFiles_small.all)]
      statFiles_medium.all = statFiles_medium.all[,2:ncol(statFiles_medium.all)]
      statFiles_large.all = statFiles_large.all[,2:ncol(statFiles_large.all)]
      
      # determine which stats were bad
      if(length(which(statFiles_small.all$V1 < 0)) > 0){
        badTests_TL.small = which(statFiles_small.all$V1 < 0)
      } else {
        badTests_TL.small = c(5001)
      }
      
      if(length(which(statFiles_medium.all$V1 < 0)) > 0){
        badTests_TL.medium = which(statFiles_medium.all$V1 < 0)
      } else {
        badTests_TL.medium = c(5001)
      }
      
      if(length(which(statFiles_large.all$V1 < 0)) > 0){
        badTests_TL.large = which(statFiles_large.all$V1 < 0)
      } else {
        badTests_TL.large = c(5001)
      }
      
      ##
      if(length(which(statFiles_small.all$V2 < 0)) > 0){
        badTests_TW.small = which(statFiles_small.all$V2 < 0)
      } else {
        badTests_TW.small = c(5001)
      }
      
      if(length(which(statFiles_medium.all$V2 < 0)) > 0){
        badTests_TW.medium = which(statFiles_medium.all$V2 < 0)
      } else {
        badTests_TW.medium = c(5001)
      }
      
      if(length(which(statFiles_large.all$V2 < 0)) > 0){
        badTests_TW.large = which(statFiles_large.all$V2 < 0)
      } else {
        badTests_TW.large = c(5001)
      }
      
      ##    
      if(length(which(statFiles_small.all$V3 < 0)) > 0){
        badTests_TS.small = which(statFiles_small.all$V3 < 0)
      } else {
        badTests_TS.small = c(5001)
      }
      
      if(length(which(statFiles_medium.all$V3 < 0)) > 0){
        badTests_TS.medium = which(statFiles_medium.all$V3 < 0)
      } else {
        badTests_TS.medium = c(5001)
      }
      
      if(length(which(statFiles_large.all$V3 < 0)) > 0){
        badTests_TS.large = which(statFiles_large.all$V3 < 0)
      } else {
        badTests_TS.large = c(5001)
      }
      
      #read in p-values
      pvalsFilesSmall = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                          paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                          paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
      
      pvalsFilesMedium = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                           paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                           paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
      
      pvalsFilesLarge = c(paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValuesCV_ForceCovFit.csv"),
                          paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValuesCV_ForceCovFit.csv"),
                          paste0("TIE_ContY_Prev_LinHypTest_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValuesCV_ForceCovFit.csv"))
      
      pvalsFiles_small = lapply(pvalsFilesSmall, read.csv)
      pvalsFiles_medium = lapply(pvalsFilesMedium, read.csv)
      pvalsFiles_large = lapply(pvalsFilesLarge, read.csv)
      
      pvalsFiles_small.all = rbind(pvalsFiles_small[[1]][1:1000,], pvalsFiles_small[[2]][1:1000,], pvalsFiles_small[[3]][1:1000,])
      pvalsFiles_medium.all = rbind(pvalsFiles_medium[[1]][1:1000,], pvalsFiles_medium[[2]][1:1000,], pvalsFiles_medium[[3]][1:1000,])
      pvalsFiles_large.all = rbind(pvalsFiles_large[[1]][1:1000,], pvalsFiles_large[[2]][1:1000,], pvalsFiles_large[[3]][1:1000,])
      
      pvalsFiles_small.all = pvalsFiles_small.all[,2:ncol(pvalsFiles_small.all)]
      pvalsFiles_medium.all = pvalsFiles_medium.all[,2:ncol(pvalsFiles_medium.all)]
      pvalsFiles_large.all = pvalsFiles_large.all[,2:ncol(pvalsFiles_large.all)]
      
      TL_goodPvals.small = pvalsFiles_small.all$V1[-c(badTests_TL.small)]
      TW_goodPvals.small = pvalsFiles_small.all$V2[-c(badTests_TW.small)]
      TS_goodPvals.small = pvalsFiles_small.all$V3[-c(badTests_TS.small)]
      
      TL_goodPvals.medium = pvalsFiles_medium.all$V1[-c(badTests_TL.medium)]
      TW_goodPvals.medium = pvalsFiles_medium.all$V2[-c(badTests_TW.medium)]
      TS_goodPvals.medium = pvalsFiles_medium.all$V3[-c(badTests_TS.medium)]
      
      TL_goodPvals.large = pvalsFiles_large.all$V1[-c(badTests_TL.large)]
      TW_goodPvals.large = pvalsFiles_large.all$V2[-c(badTests_TW.large)]
      TS_goodPvals.large = pvalsFiles_large.all$V3[-c(badTests_TS.large)]
      
      TL_TIE.small = sum(TL_goodPvals.small > 0)/length(TL_goodPvals.small)
      TW_TIE.small = sum(TW_goodPvals.small > 0)/length(TW_goodPvals.small)
      TS_TIE.small = sum(TS_goodPvals.small > 0)/length(TS_goodPvals.small)
      
      TIE.small = rbind(round(TL_TIE.small, 2), round(TW_TIE.small,2), round(TS_TIE.small,2))
      
      TL_TIE.medium = sum(TL_goodPvals.medium > 0)/length(TL_goodPvals.medium)
      TW_TIE.medium = sum(TW_goodPvals.medium > 0)/length(TW_goodPvals.medium)
      TS_TIE.medium = sum(TS_goodPvals.medium > 0)/length(TS_goodPvals.medium)
      
      TIE.medium = rbind(round(TL_TIE.medium, 2), round(TW_TIE.medium,2), round(TS_TIE.medium,2))
      
      TL_TIE.large = sum(TL_goodPvals.large > 0)/length(TL_goodPvals.large)
      TW_TIE.large = sum(TW_goodPvals.large > 0)/length(TW_goodPvals.large)
      TS_TIE.large = sum(TS_goodPvals.large > 0)/length(TS_goodPvals.large)
      
      TIE.large = rbind(round(TL_TIE.large,2), round(TW_TIE.large,2), round(TS_TIE.large,2))
      
      powerByLDCont[1:3,aa] = TIE.small
      powerByLDCont[4:6,aa] = TIE.medium
      powerByLDCont[7:9,aa] = TIE.large
    } 
    powerByPercentAssocCont[[kk]] =  powerByLDCont
  }
}

powerByPercentAssocCont.mat = do.call(cbind, powerByPercentAssocCont)

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByPercentAssocCont.mat)
rownames(PowerData) = c("TL - Prev 20 - Small OR", "TW - Prev 20 - Small OR", "TS - Prev 20 - Small OR", 
                        "TL - Prev 20 - Medium OR", "TW - Prev 20 - Medium OR", "TS - Prev 20 - Medium OR",
                        "TL - Prev 20 - Large OR", "TW - Prev 20 - Large OR", "TS - Prev 20 - Large OR",
                        "TL - Prev 10 - Small OR", "TW - Prev 10 - Small OR", "TS - Prev 10 - Small OR",
                        "TL - Prev 10 - Medium OR", "TW - Prev 10 - Medium OR", "TS - Prev 10 - Medium OR",
                        "TL - Prev 10 - Large OR", "TW - Prev 10 - Large OR", "TS - Prev 10 - Large OR",
                        "TL - Prev 5 - Small OR", "TW - Prev 5 - Small OR",  "TS - Prev 5 - Small OR", 
                        "TL - Prev 5 - Medium OR", "TW - Prev 5 - Medium OR",  "TS - Prev 5 - Medium OR", 
                        "TL - Prev 5 - Large OR", "TW - Prev 5 - Large OR", "TS - Prev 5 - Large OR", 
                        "TL - Cont. - Small OR", "TW - Cont. - Small OR", "TS - Cont. - Small OR", 
                        "TL - Cont. - Medium OR", "TW - Cont. - Medium OR", "TS - Cont. - Medium OR",
                        "TL - Cont. - Large OR", "TW - Cont. - Large OR", "TS - Cont. - Large OR")

colnames(PowerData) = c("5% Assoc SNPs - Low LD", "5% Assoc SNPs - High LD", "15% Assoc SNPs - Low LD", "15% Assoc SNPs - High LD",
                        "25% Assoc SNPs - Low LD", "25% Assoc SNPs - High LD")
write.csv(PowerData, paste0(gene,"_",ss,"Pairs_LHT_TIE_3KSims_GammaZero_BetaAssoc.csv"))
## For Desparsified Lasso ########################################################################
## Gamma = 0 Only Testing, When Beta is Sig
# gene = "NAT2"
# gene = "CHI3L2"
gene = "ASAH1"
# ss = 500
ss = 1000

nSims = 3000

Prev = c(20, 10, 5)
FitScore = c("IBS")

OR = c("small", "medium", "large")
assocSNP = c(5,15,25)
lowLD = c(TRUE, FALSE)

path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_DesparseLasso_TIE_GammaZero_BetaAssoc")
setwd(path)

powerList = list()
powerByPercentAssoc = list()
powerByLD = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))

for(jj in 1:length(FitScore)){ 
  for(ii in 1:length(Prev)){
    for(kk in 1:length(assocSNP)){
      for(aa in 1:length(lowLD)){
        statsFilesSmall = c(paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statsFilesMedium = c(paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statsFilesLarge = c(paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                            paste0("TIE_CatY_Prev",Prev[ii],"_DesparseLasso_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statFiles_Small = lapply(statsFilesSmall, read.csv)
        statFiles_Medium = lapply(statsFilesMedium, read.csv)
        statFiles_Large = lapply(statsFilesLarge, read.csv)
        
        statFiles_Small.all = rbind(statFiles_Small[[1]][1:1000,], statFiles_Small[[2]][1:1000,], statFiles_Small[[3]][1:1000,])
        statFiles_Medium.all = rbind(statFiles_Medium[[1]][1:1000,], statFiles_Medium[[2]][1:1000,], statFiles_Medium[[3]][1:1000,])
        statFiles_Large.all = rbind(statFiles_Large[[1]][1:1000,], statFiles_Large[[2]][1:1000,], statFiles_Large[[3]][1:1000,])
        
        tie_small = sum(statFiles_Small.all$V1 < 0.05)/dim(statFiles_Small.all)[1]
        tie_medium = sum(statFiles_Medium.all$V1 < 0.05)/dim(statFiles_Medium.all)[1]  
        tie_large = sum(statFiles_Large.all$V1 < 0.05)/dim(statFiles_Large.all)[1]
        
        powerByLD[1,aa] = tie_small
        powerByLD[2,aa] = tie_medium
        powerByLD[3,aa] = tie_large
      } 
      powerByPercentAssoc[[kk]] = powerByLD
    }
    powerByPercentAssoc.mat = do.call(cbind, powerByPercentAssoc)
    powerList[[ii]] = powerByPercentAssoc.mat
  }
}

powerByPercentAssocCont = list()
powerByLDCont = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))

for(jj in 1:length(FitScore)){
  for(kk in 1:length(assocSNP)){
    for(aa in 1:length(lowLD)){
      statsFilesSmall = c(paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statsFilesMedium = c(paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                          paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                          paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statsFilesLarge = c(paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                          paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                          paste0("TIE_ContY_Prev_DesparseLasso_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statFiles_Small = lapply(statsFilesSmall, read.csv)
      statFiles_Medium = lapply(statsFilesMedium, read.csv)
      statFiles_Large = lapply(statsFilesLarge, read.csv)
      
      statFiles_Small.all = rbind(statFiles_Small[[1]][1:1000,], statFiles_Small[[2]][1:1000,], statFiles_Small[[3]][1:1000,])
      statFiles_Medium.all = rbind(statFiles_Medium[[1]][1:1000,], statFiles_Medium[[2]][1:1000,], statFiles_Medium[[3]][1:1000,])
      statFiles_Large.all = rbind(statFiles_Large[[1]][1:1000,], statFiles_Large[[2]][1:1000,], statFiles_Large[[3]][1:1000,])
      
      TIE_Small_Cont = sum(statFiles_Small.all$V1 < 0.05)/dim(statFiles_Small.all)[1]
      TIE_Medium_Cont = sum(statFiles_Medium.all$V1 < 0.05)/dim(statFiles_Medium.all)[1]
      TIE_Large_Cont = sum(statFiles_Large.all$V1 < 0.05)/dim(statFiles_Large.all)[1]
      
      powerByLDCont[1,aa] = TIE_Small_Cont
      powerByLDCont[2,aa] = TIE_Medium_Cont
      powerByLDCont[3,aa] = TIE_Large_Cont
    }
    powerByPercentAssocCont[[kk]] = powerByLDCont
  }
}

powerByPercentAssocCont.mat = do.call(cbind, powerByPercentAssocCont)

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByPercentAssocCont.mat)
rownames(PowerData) = c("Prev 20 - Small OR","Prev 20 - Medium OR", "Prev 20 - Large OR", 
                        "Prev 10 - Small OR", "Prev 10 - Medium OR", "Prev 10 - Large OR", 
                        "Prev 5 - Small OR", "Prev 5 - Medium OR", "Prev 5 - Large OR", 
                        "Cont. - Small OR", "Cont. - Medium OR", "Cont. - Large OR")

colnames(PowerData) = c("5% Assoc SNPs - Low LD", "5% Assoc SNPs - High LD", "15% Assoc SNPs - Low LD", "15% Assoc SNPs - High LD",
                        "25% Assoc SNPs - Low LD", "25% Assoc SNPs - High LD")

write.csv(PowerData, paste0(gene,"_",ss,"Pairs_DesparseLasso_TIE_",nSims,"Sims_GammaZero_BetaAssoc.csv"))

## For AISPU ########################################################################
## Gamma = 0 Only Testing, When Beta is Sig
# gene = "NAT2"
# gene = "CHI3L2"
gene = "ASAH1"
# ss = 500
ss = 1000

nSims = 3000

penalty = "scad"

Prev = c(20, 10, 5)
FitScore = c("IBS")

OR = c("small", "medium", "large")
assocSNP = c(5,15,25)
lowLD = c(TRUE, FALSE)

path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_AISPU_TIE_GammaZero_BetaAssoc")
setwd(path)

powerList = list()
powerList_inf = list()
powerByPercentAssoc = list()
powerByPercentAssoc_inf = list()
powerByLD = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))
powerByLD_inf = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))

for(jj in 1:length(FitScore)){ 
  for(ii in 1:length(Prev)){
    for(kk in 1:length(assocSNP)){
      for(aa in 1:length(lowLD)){
        statsFilesSmall = c(paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statsFilesMedium = c(paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statsFilesLarge = c(paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_AISPU_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statFiles_Small = lapply(statsFilesSmall, read.csv)
        statFiles_Medium = lapply(statsFilesMedium, read.csv)
        statFiles_Large = lapply(statsFilesLarge, read.csv)
        
        statFiles_Small.all = rbind(statFiles_Small[[1]][1:1000,], statFiles_Small[[2]][1:1000,], statFiles_Small[[3]][1:1000,])
        statFiles_Medium.all = rbind(statFiles_Medium[[1]][1:1000,], statFiles_Medium[[2]][1:1000,], statFiles_Medium[[3]][1:1000,])
        statFiles_Large.all = rbind(statFiles_Large[[1]][1:1000,], statFiles_Large[[2]][1:1000,], statFiles_Large[[3]][1:1000,])
        
        statFiles_Small.all = na.omit(statFiles_Small.all)
        statFiles_Medium.all = na.omit(statFiles_Medium.all)
        statFiles_Large.all = na.omit(statFiles_Large.all)
        
        power_Small = sum(statFiles_Small.all$V1 < 0.05)/dim(statFiles_Small.all)[1]
        power_Medium = sum(statFiles_Medium.all$V1 < 0.05)/dim(statFiles_Medium.all)[1]  
        power_Large = sum(statFiles_Large.all$V1 < 0.05)/dim(statFiles_Large.all)[1]  
        
        power_Small_inf = sum(statFiles_Small.all$V2 < 0.05)/dim(statFiles_Small.all)[1]
        power_Medium_inf = sum(statFiles_Medium.all$V2 < 0.05)/dim(statFiles_Medium.all)[1]  
        power_Large_inf = sum(statFiles_Large.all$V2 < 0.05)/dim(statFiles_Large.all)[1]  
        
        powerByLD[1,aa] = power_Small
        powerByLD[2,aa] = power_Medium
        powerByLD[3,aa] = power_Large
        
        powerByLD_inf[1,aa] = power_Small_inf
        powerByLD_inf[2,aa] = power_Medium_inf
        powerByLD_inf[3,aa] = power_Large_inf
      }
      powerByPercentAssoc[[kk]] = powerByLD
      powerByPercentAssoc_inf[[kk]] = powerByLD_inf
    }
    powerByPercentAssoc.mat = do.call(cbind, powerByPercentAssoc)
    powerByPercentAssoc_inf.mat = do.call(cbind, powerByPercentAssoc_inf)
    powerList[[ii]] = powerByPercentAssoc.mat
    powerList_inf[[ii]] = powerByPercentAssoc_inf.mat
  }
}

powerByPercentAssocCont = list()
powerByPercentAssocCont_inf = list()
powerByLDCont = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))
powerByLDCont_inf = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))

for(jj in 1:length(FitScore)){
  for(kk in 1:length(assocSNP)){
    for(aa in 1:length(lowLD)){
      statsFilesSmall = c(paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                          paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                          paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statsFilesMedium = c(paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                           paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                           paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statsFilesLarge = c(paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                          paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                          paste0("TIE_ContY_Prev_AISPU_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statFiles_Small = lapply(statsFilesSmall, read.csv)
      statFiles_Medium = lapply(statsFilesMedium, read.csv)
      statFiles_Large = lapply(statsFilesLarge, read.csv)
      
      statFiles_Small.all = rbind(statFiles_Small[[1]][1:1000,], statFiles_Small[[2]][1:1000,], statFiles_Small[[3]][1:1000,])
      statFiles_Medium.all = rbind(statFiles_Medium[[1]][1:1000,], statFiles_Medium[[2]][1:1000,], statFiles_Medium[[3]][1:1000,])
      statFiles_Large.all = rbind(statFiles_Large[[1]][1:1000,], statFiles_Large[[2]][1:1000,], statFiles_Large[[3]][1:1000,])
      
      TIE_Small_Cont = sum(statFiles_Small.all$V1 < 0.05)/dim(statFiles_Small.all)[1]
      TIE_Medium_Cont = sum(statFiles_Medium.all$V1 < 0.05)/dim(statFiles_Medium.all)[1]
      TIE_Large_Cont = sum(statFiles_Large.all$V1 < 0.05)/dim(statFiles_Large.all)[1]
      
      TIE_Small_Cont_inf = sum(statFiles_Small.all$V2 < 0.05)/dim(statFiles_Small.all)[1]
      TIE_Medium_Cont_inf = sum(statFiles_Medium.all$V2 < 0.05)/dim(statFiles_Medium.all)[1]
      TIE_Large_Cont_inf = sum(statFiles_Large.all$V2 < 0.05)/dim(statFiles_Large.all)[1]
      
      powerByLDCont[1,aa] = TIE_Small_Cont
      powerByLDCont[2,aa] = TIE_Medium_Cont
      powerByLDCont[3,aa] = TIE_Large_Cont
      
      powerByLDCont_inf[1,aa] = TIE_Small_Cont_inf
      powerByLDCont_inf[2,aa] = TIE_Medium_Cont_inf
      powerByLDCont_inf[3,aa] = TIE_Large_Cont_inf
    } 
    powerByPercentAssocCont[[kk]] = powerByLDCont
    powerByPercentAssocCont_inf[[kk]] = powerByLDCont_inf
  } 
}

powerByPercentAssocCont.mat = do.call(cbind, powerByPercentAssocCont)
powerByPercentAssocCont_inf.mat = do.call(cbind, powerByPercentAssocCont_inf)

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByPercentAssocCont.mat)
PowerData_inf = rbind(powerList_inf[[1]], powerList_inf[[2]], powerList_inf[[3]], powerByPercentAssocCont_inf.mat)

#see if they equal:
PowerData == PowerData_inf

rownames(PowerData) = c("Prev 20 - Small OR","Prev 20 - Medium OR", "Prev 20 - Large OR", 
                        "Prev 10 - Small OR", "Prev 10 - Medium OR", "Prev 10 - Large OR", 
                        "Prev 5 - Small OR", "Prev 5 - Medium OR", "Prev 5 - Large OR", 
                        "Cont. - Small OR", "Cont. - Medium OR", "Cont. - Large OR")

colnames(PowerData) = c("5% Assoc SNPs - Low LD", "5% Assoc SNPs - High LD", "15% Assoc SNPs - Low LD", "15% Assoc SNPs - High LD",
                        "25% Assoc SNPs - Low LD", "25% Assoc SNPs - High LD")

write.csv(PowerData, paste0(gene,"_",ss,"Pairs_AISPU_TIE_",nSims,"Sims_",penalty,"Penalty_GammaZero_BetaAssoc.csv"))

## For Decorrelated Score ########################################################################
## Gamma = 0 Only Testing, When Beta is Sig
gene = "NAT2"
# gene = "CHI3L2"
# gene = "ASAH1"
# ss = 500
ss = 1000

nSims = 3000

Prev = c(20, 10, 5)
FitScore = c("IBS")

OR = c("small", "medium", "large")
assocSNP = c(5,15,25)
lowLD = c(TRUE, FALSE)

path = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_III/",gene,"_",ss,"_DecorrScore_TIE_GammaZero_BetaAssoc")
setwd(path)

powerList = list()
powerByPercentAssoc = list()
powerByLD = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))

for(jj in 1:length(FitScore)){ 
  for(ii in 1:length(Prev)){
    for(kk in 1:length(assocSNP)){
      for(aa in 1:length(lowLD)){    
        statsFilesSmall = c(paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statsFilesMedium = c(paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statsFilesLarge = c(paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                             paste0("TIE_CatY_Prev",Prev[ii],"_DecorrScore_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
        
        statFiles_Small = lapply(statsFilesSmall, read.csv)
        statFiles_Medium = lapply(statsFilesMedium, read.csv)
        statFiles_Large = lapply(statsFilesLarge, read.csv)
        
        statFiles_Small.all = rbind(statFiles_Small[[1]][1:1000,], statFiles_Small[[2]][1:1000,], statFiles_Small[[3]][1:1000,])
        statFiles_Medium.all = rbind(statFiles_Medium[[1]][1:1000,], statFiles_Medium[[2]][1:1000,], statFiles_Medium[[3]][1:1000,])
        statFiles_Large.all = rbind(statFiles_Large[[1]][1:1000,], statFiles_Large[[2]][1:1000,], statFiles_Large[[3]][1:1000,])
        
        statFiles_Small.all = na.omit(statFiles_Small.all)
        statFiles_Medium.all = na.omit(statFiles_Medium.all)
        statFiles_Large.all = na.omit(statFiles_Large.all)
        
        TIE_Small = mean(statFiles_Small.all$V1 < 0.05)
        TIE_Medium = mean(statFiles_Medium.all$V1 < 0.05)  
        TIE_Large = mean(statFiles_Large.all$V1 < 0.05)  
        
        powerByLD[1,aa] = TIE_Small
        powerByLD[2,aa] = TIE_Medium
        powerByLD[3,aa] = TIE_Large
      }
      powerByPercentAssoc[[kk]] = powerByLD
    }
    powerByPercentAssoc.mat = do.call(cbind, powerByPercentAssoc)
    powerList[[ii]] = powerByPercentAssoc.mat
  }
}

powerByPercentAssocCont = list()
powerByLDCont = matrix(data = NA, nrow = length(OR), ncol = length(lowLD))

for(jj in 1:length(FitScore)){
  for(kk in 1:length(assocSNP)){
    for(aa in 1:length(lowLD)){
      statsFilesSmall = c(paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[1],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statsFilesMedium = c(paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[2],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statsFilesLarge = c(paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1to1000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim1001to2000_PValues.csv"),
                           paste0("TIE_ContY_Prev_DecorrScore_Score",FitScore[jj],"_",OR[3],"OR_",assocSNP[kk],"Assoc_LowLD",lowLD[aa],"_GammaOnly_BetaSig__WCovs_Sim2001to3000_PValues.csv"))
      
      statFiles_Small = lapply(statsFilesSmall, read.csv)
      statFiles_Medium = lapply(statsFilesMedium, read.csv)
      statFiles_Large = lapply(statsFilesLarge, read.csv)
      
      statFiles_Small.all = rbind(statFiles_Small[[1]][1:1000,], statFiles_Small[[2]][1:1000,], statFiles_Small[[3]][1:1000,])
      statFiles_Medium.all = rbind(statFiles_Medium[[1]][1:1000,], statFiles_Medium[[2]][1:1000,], statFiles_Medium[[3]][1:1000,])
      statFiles_Large.all = rbind(statFiles_Large[[1]][1:1000,], statFiles_Large[[2]][1:1000,], statFiles_Large[[3]][1:1000,])
      
      TIE_Small_Cont = sum(statFiles_Small.all$V1 < 0.05)/dim(statFiles_Small.all)[1]
      TIE_Medium_Cont = sum(statFiles_Medium.all$V1 < 0.05)/dim(statFiles_Medium.all)[1]
      TIE_Large_Cont = sum(statFiles_Large.all$V1 < 0.05)/dim(statFiles_Large.all)[1]
      
      powerByLDCont[1,aa] = TIE_Small_Cont
      powerByLDCont[2,aa] = TIE_Medium_Cont
      powerByLDCont[3,aa] = TIE_Large_Cont
    }
    powerByPercentAssocCont[[kk]] = powerByLDCont
  }
}

powerByPercentAssocCont.mat = do.call(cbind, powerByPercentAssocCont)

PowerData = rbind(powerList[[1]], powerList[[2]], powerList[[3]], powerByPercentAssocCont.mat)
rownames(PowerData) = c("Prev 20 - Small OR","Prev 20 - Medium OR", "Prev 20 - Large OR", 
                        "Prev 10 - Small OR", "Prev 10 - Medium OR", "Prev 10 - Large OR", 
                        "Prev 5 - Small OR", "Prev 5 - Medium OR", "Prev 5 - Large OR", 
                        "Cont. - Small OR", "Cont. - Medium OR", "Cont. - Large OR")

colnames(PowerData) = c("5% Assoc SNPs - Low LD", "5% Assoc SNPs - High LD", "15% Assoc SNPs - Low LD", "15% Assoc SNPs - High LD",
                        "25% Assoc SNPs - Low LD", "25% Assoc SNPs - High LD")

write.csv(PowerData, paste0(gene,"_",ss,"Pairs_DecorrScore_TIE_",nSims,"Sims_GammaZero_BetaAssoc.csv"))
