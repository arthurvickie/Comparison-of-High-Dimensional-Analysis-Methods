##########################################################################################################################################
## Shi, Song, Chen, and Li Paper Pipelines
##########################################################################################################################################
# rely on code from 2019 Annals of Statistics Paper: Linear Hyp Testing for High Dim GLMs
##########################################################################################################################################
#################           Gamma = 0 Testing Only         ###############################################################################
#################### TIE Pipeline with AnnalsOfStats Paper
#################### Case 1: Testing Gamma = 0 when Beta = 0

#TIE for GLM comparison (Gamma and Beta = 0)
RunTIEPipelineLinHypTest_GammaTest_GLM = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, start, stop){
  #function to run whole TIE pipeline, Calculates LRT test stat for Lin Hyp test method and whether p-value <0.05
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #start - simulation number to start at
  #stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files

  library(parallel)
  suppressMessages(library(lmtest))

  #define stopping point
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")

  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")

  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()

  statsAndPVals = mclapply(myList, function(ii){
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
	  # and Cat then Cont phenos in list form
	  allData = CalcNullPhenotypeData(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs)
    
	  output = RunGLM.LRT(allData, FitScore, Cat = Cat,covs = covs)
	  
    print(paste0("Simulation ",ii," is complete."))

    finalOutput[[ii]] = output
  }, mc.cores=4)

  reorgData = OrganizeOutput(statsAndPVals)
  
  statsAndPValsAll.mat = reorgData[[1]]
  nulValuesAll.mat = reorgData[[2]]
  altValuesAll.mat = reorgData[[3]]
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_Score",FitScore,"_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatPhenos_NullBetasSummary_Prev",YPrev*100,"_GLM_Score",FitScore,"_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatPhenos_AltBetasSummary_Prev",YPrev*100,"_GLM_Score",FitScore,"_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  } else {
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/TIE_ContPhenos_Phenos_GLM_Score",FitScore,"_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContPhenos_NullBetasSummary_GLM_Score",FitScore,"_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContPhenos_AltBetasSummary_GLM_Score",FitScore,"_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  }
  
}

#TIE for LHT (Gamma and Beta = 0)
RunTIEPipelineLinHypTest_GammaTest = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  if(Cat == TRUE){
    source("/home/vlynn/Paper_III_Sims/Scripts/Logistic_ADMM0.r")
  } else {
    source("/home/vlynn/Paper_III_Sims/Scripts/Linear_ADMM0.r")
  }
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcNullPhenotypeData(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs)
    
    if(Cat == TRUE){
      results = RunLHT_Cat(allData = allData, FitScore = FitScore, covs = covs)
    } else {
      results = RunLHT_Cont(allData = allData, FitScore = FitScore, covs = covs)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  pValsCV = matrix(nrow = n, ncol = 3)
  ScoresCV = matrix(nrow = n, ncol = 3)
  BetasCV = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pValsCV[ll,] = unlist(statsAndPVals[[ll]]$pv)
    ScoresCV[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    BetasCV[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values, and Summary stats
    #write out the Stats and p values
    write.csv(pValsCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValuesCV_ForceCovFit.csv"))
    write.csv(ScoresCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_StatsCV_ForceCovFit.csv"))
    write.csv(BetasCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_BetasCV_ForceCovFit.csv"))
  } else {
    #write out the Stats and p values, and Summary stats
    write.csv(pValsCV, file = paste0(path,"/TIE_ContY_Prev_LinHypTest_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValuesCV_ForceCovFit.csv"))
    write.csv(ScoresCV, file = paste0(path,"/TIE_ContY_Prev_LinHypTest_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_StatsCV_ForceCovFit.csv"))
    write.csv(BetasCV, file = paste0(path,"/TIE_ContY_Prev_LinHypTest_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_BetasCV_ForceCovFit.csv"))
  }
}

#TIE for Desparsified Lasso (Gamma and Beta = 0)
RunTIEPipelineDesparseLasso_GammaTest = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, start, stop){
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  # start - simulation number to start at
  # stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcNullPhenotypeData(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs)

    results = RunDesparseLasso(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 1)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]])
  }
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the p values
    write.csv(pVals, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_DesparseLasso_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  } else {
    #write out the p values
    write.csv(pVals, file = paste0(path,"/TIE_ContY_Prev_DesparseLasso_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  }
}

#TIE for AISPU (Gamma and Beta = 0)
RunTIEPipelineAISPU_GammaTest = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, seed, start, stop){
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #seed - seed that was input for reproducibility
  #start - simulation number to start at
  #stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  library(rlist)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcNullPhenotypeData(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs)
    
    results = RunAISPU(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs, seed = seed)
    
    if(length(results) != 3){
      #then there was an error, so set results to something
      results = list("NA", "NA", "NA")
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  #save the list in case it errors:
   if(Cat == TRUE){
     #write out the p values
     list.save(statsAndPVals, paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,".rdata"))
  } else {
    #write out the p values
    list.save(statsAndPVals, paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,".rdata"))
  }
  
  for(ll in 1:n){
		lengthBetas = max(length(statsAndPVals[[ll]][[3]]))
	}
  
  Scores = matrix(nrow = n, ncol = 2)
  pVals = matrix(nrow = n, ncol = 2)
  betas = matrix(nrow = n, ncol=lengthBetas)
  
  for(ll in 1:n){
    Scores[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals[ll,] = unlist(statsAndPVals[[ll]][[2]])
    betas[ll,1:length(unlist(statsAndPVals[[ll]][[3]]))] = unlist(statsAndPVals[[ll]][[3]])
  }
  
  if(Cat == TRUE){
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(betas, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_Betas.csv"))
  } else {
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(betas, file = paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_Betas.csv"))
  }
}

#TIE for Decorrelated Score Test (Gamma and Beta = 0)
RunTIEPipelineDecorrScore_GammaTest = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, start, stop){
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #start - simulation number to start at
  #stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcNullPhenotypeData(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs)
    
    results = RunDecorrScore(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  Scores = matrix(nrow = n, ncol = 1)
  pVals = matrix(nrow = n, ncol = 1)

  for(ll in 1:n){
    Scores[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals[ll,] = unlist(statsAndPVals[[ll]][[2]])
  }
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  } else {
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_ContY_Prev_DecorrScore_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_ContY_Prev_DecorrScore_Score",FitScore,"_GammaOnly",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  }
}

#################### Case 2: Testing Gamma = 0 when Beta not = 0
#TIE for GLM comparison (Gamma and Beta neq 0)
RunTIEPipelineLinHypTest_GammaTest_GLM_BetaSig = function(chr, gene, numPairs, YPrev, ORSize, percentageAssoc, LowLD, covs, FitScore, Cat, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #ORSize = small, medium, or large for OR of 1.25, 1.50, 2.00
  #percentageAssoc = 5, 15, or 25% of RSNPs associated with outcome
  #LowLD = T or F, are the associated SNPs in low or high LD with the rest of the SNPs in the region?
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #start - simulation number to start at
  #stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  suppressMessages(library(lmtest))
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")

  myList = lapply(start:numSims, rep, times = 1)
  nullValues =  altValues = statsAndPValsmat = finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    
    allData = CalcAltPhenotypeData_RSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSize = ORSize,
                                         LowLD = LowLD, percentAssoc = percentageAssoc, TrueScore = FitScore)
    
    output = RunGLM.LRT(allData, FitScore, Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = output
  }, mc.cores=4)
  
  reorgData = OrganizeOutput(statsAndPVals)
  
  statsAndPValsAll.mat = reorgData[[1]]
  nulValuesAll.mat = reorgData[[2]]
  altValuesAll.mat = reorgData[[3]]
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values, and Summary stats  , 
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/TIE_CatPhenos_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_CatPhenos_NullBetasSummary_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_CatPhenos_AltBetasSummary_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  } else {
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/TIE_ContPhenos_Phenos_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValuesAll.mat, file = paste0(path,"/TIE_ContPhenos_NullBetasSummary_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValuesAll.mat, file = paste0(path,"/TIE_ContPhenos_AltBetasSummary_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  }
}

#TIE for LHT (Gamma and Beta neq 0)
RunTIEPipelineLinHypTest_GammaTest_BetaSig = function(chr, gene, numPairs, YPrev, ORSize, percentageAssoc, LowLD, covs, FitScore, Cat, start, stop){
  #function to run whole TIE pipeline, Calculates Score test stat for Lin Hyp test method and whether p-value <0.05
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #ORSize = small, medium, or large for OR of 1.25, 1.50, 2.00
  #percentageAssoc = 5, 15, or 25% of RSNPs associated with outcome
  #LowLD = T or F, are the associated SNPs in low or high LD with the rest of the SNPs in the region?
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #start - simulation number to start at
  #stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  if(Cat == TRUE){
    source("/home/vlynn/Paper_III_Sims/Scripts/Logistic_ADMM0.r")
  } else {
    source("/home/vlynn/Paper_III_Sims/Scripts/Linear_ADMM0.r")
  }
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcAltPhenotypeData_RSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSize = ORSize, LowLD = LowLD, percentAssoc = percentageAssoc, TrueScore = FitScore)
     
	if(Cat == TRUE){
      results = RunLHT_Cat(allData = allData, FitScore = FitScore, covs = covs)
    } else {
      results = RunLHT_Cont(allData = allData, FitScore = FitScore, covs = covs)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  pValsCV = matrix(nrow = n, ncol = 3)
  ScoresCV = matrix(nrow = n, ncol = 3)
  BetasCV = matrix(nrow = n, ncol = 2)
  
  for(ll in 1:n){
    pValsCV[ll,] = unlist(statsAndPVals[[ll]]$pv)
    ScoresCV[ll,] = unlist(statsAndPVals[[ll]]$TScores)
    BetasCV[ll,] = unlist(statsAndPVals[[ll]]$beta)
  }
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values, and Summary stats
    #write out the Stats and p values
    write.csv(pValsCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValuesCV_ForceCovFit.csv"))
    write.csv(ScoresCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsCV_ForceCovFit.csv"))
    write.csv(BetasCV, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_BetasCV_ForceCovFit.csv"))
  } else {
    #write out the Stats and p values, and Summary stats
    write.csv(pValsCV, file = paste0(path,"/TIE_ContY_Prev_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValuesCV_ForceCovFit.csv"))
    write.csv(ScoresCV, file = paste0(path,"/TIE_ContY_Prev_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_StatsCV_ForceCovFit.csv"))
    write.csv(BetasCV, file = paste0(path,"/TIE_ContY_Prev_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_BetasCV_ForceCovFit.csv"))
  }
}

#TIE for Desparsified Lasso (Gamma and Beta neq 0)
RunTIEPipelineDesparseLasso_GammaTest_BetaSig = function(chr, gene, numPairs, YPrev, ORSize, percentageAssoc, LowLD, covs, FitScore, Cat, start, stop){
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #ORSize = small, medium, or large for OR of 1.25, 1.50, 2.00
  #percentageAssoc = 5, 15, or 25% of RSNPs associated with outcome
  #LowLD = T or F, are the associated SNPs in low or high LD with the rest of the SNPs in the region?
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #start - simulation number to start at
  #stop - simulation number to stop at#Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcAltPhenotypeData_RSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSize = ORSize, LowLD = LowLD, percentAssoc = percentageAssoc, TrueScore = FitScore)
    
    results = RunDesparseLasso(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals = matrix(nrow = n, ncol = 1)
  
  for(ll in 1:n){
    pVals[ll,] = unlist(statsAndPVals[[ll]])
  }
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the p values
    write.csv(pVals, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_DesparseLasso_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  } else {
    #write out the p values
    write.csv(pVals, file = paste0(path,"/TIE_ContY_Prev_DesparseLasso_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  }
}

#TIE for AISPU (Gamma and Beta neq 0)
RunTIEPipelineAISPU_GammaTest_BetaSig = function(chr, gene, numPairs, YPrev, ORSize, percentageAssoc, LowLD, covs, FitScore, Cat, seed, start, stop){
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #ORSize = small, medium, or large for OR of 1.25, 1.50, 2.00
  #percentageAssoc = 5, 15, or 25% of RSNPs associated with outcome
  #LowLD = T or F, are the associated SNPs in low or high LD with the rest of the SNPs in the region?
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #start - simulation number to start at
  #stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  library(rlist)

  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcAltPhenotypeData_RSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSize = ORSize,LowLD = LowLD, percentAssoc = percentageAssoc, TrueScore = FitScore)
    
    results = RunAISPU(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs, seed = seed)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  for(ll in 1:n){
		lengthBetas = max(length(statsAndPVals[[ll]][[3]]))
	}
  
  Scores = matrix(nrow = n, ncol = 2)
  pVals = matrix(nrow = n, ncol = 2)
  betas = matrix(nrow = n, ncol=lengthBetas)
   
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  #save the list in case it errors:
   if(Cat == TRUE){
     #write out the p values
     list.save(statsAndPVals, paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,".rdata"))
  } else {
    #write out the p values
    list.save(statsAndPVals, paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,".rdata"))
  }
  
  for(ll in 1:n){
    Scores[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals[ll,] = unlist(statsAndPVals[[ll]][[2]])
    betas[ll,1:length(unlist(statsAndPVals[[ll]][[3]]))] = unlist(statsAndPVals[[ll]][[3]])
  }
  
  if(Cat == TRUE){
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(betas, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_Betas.csv"))
  } else {
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(betas, file = paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_Betas.csv"))
  }
}

#TIE for Decorrelated Score Test (Gamma and Beta neq 0)
RunTIEPipelineDecorrScore_GammaTest_BetaSig = function(chr, gene, numPairs, YPrev, ORSize, percentageAssoc, LowLD, covs, FitScore, Cat, start, stop){
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #ORSize = small, medium, or large for OR of 1.25, 1.50, 2.00
  #percentageAssoc = 5, 15, or 25% of RSNPs associated with outcome
  #LowLD = T or F, are the associated SNPs in low or high LD with the rest of the SNPs in the region?
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #start - simulation number to start at
  #stop - simulation number to stop at
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes TIE values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  statsAndPVals = mclapply(myList, function(ii){
    #define matrix to hold all Stats and Pvalues
    statsAndPVals = matrix(NA, nrow = numSims, ncol = 4)
    
    #calc null phenos, and output RGenos, Gene Score, Covariates (if wanted)
    # and Cat then Cont phenos in list form
    allData = CalcAltPhenotypeData_RSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSize = ORSize,LowLD = LowLD, percentAssoc = percentageAssoc, TrueScore = FitScore)
    
    results = RunDecorrScore(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    finalOutput[[ii]] = results
  }, mc.cores=6)
  
  n = length(statsAndPVals)
  
  Scores = matrix(nrow = n, ncol = 1)
  pVals = matrix(nrow = n, ncol = 1)
  
  for(ll in 1:n){
    Scores[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals[ll,] = unlist(statsAndPVals[[ll]][[2]])
  }
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  } else {
    #write out the p values
    write.csv(Scores, file = paste0(path,"/TIE_ContY_Prev_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_Scores.csv"))
    write.csv(pVals, file = paste0(path,"/TIE_ContY_Prev_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaOnly_BetaSig_",covariates,"Sim",start,"to",numSims,"_PValues.csv"))
  }
}

##########################################################
## Power Pipelines
## Score true (gamma neq 0), beta = 0

#Power for GLM
RunPowerPipelineGLM_GammaTest_Score_BetaZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSize, start, stop, percentageAssoc, LowLD){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  suppressMessages(library(lmtest))

  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_Scores(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, percentageAssoc = percentageAssoc, LowLD = LowLD, TrueScore = FitScore)
    
    results.FitScore = RunGLM.LRT(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    results.OtherScore1 = RunGLM.LRT(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs) 
    results.OtherScore2 = RunGLM.LRT(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs)
    results.OtherScore3 = RunGLM.LRT(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=4)
  
  reorgData = OrganizeOutputPower(statsAndPVals)
  
  statsAndPValsAll.mat = reorgData[[1]]
  nulValues1.mat = reorgData[[2]]
  altValues1.mat = reorgData[[3]]
  nulValues2.mat = reorgData[[4]]
  altValues2.mat = reorgData[[5]]
  nulValues3.mat = reorgData[[6]]
  altValues3.mat = reorgData[[7]]
  nulValues4.mat = reorgData[[8]]
  altValues4.mat = reorgData[[9]]
  Score.mat = reorgData[[10]]
  
  nameColsStatPVal = c(paste0(Score.mat[[1]]," P Val"), paste0(Score.mat[[1]]," Score"),
                       paste0(Score.mat[[2]]," P Val"), paste0(Score.mat[[2]]," Score"),
                       paste0(Score.mat[[3]]," P Val"), paste0(Score.mat[[3]]," Score"),
                       paste0(Score.mat[[4]]," P Val"), paste0(Score.mat[[4]]," Score"))
  colnames(statsAndPValsAll.mat) = nameColsStatPVal
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values, and Summary stats  , 
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues1.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[1]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues2.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[2]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues3.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[3]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues4.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[4]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues1.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[1]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues2.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[2]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues3.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[3]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues4.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[4]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  } else {
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/Power_ContPhenos_Phenos_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues1.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[1]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues2.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[2]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues3.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[3]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues4.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[4]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues1.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[1]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues2.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[2]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues3.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[3]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues4.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[4]],"_GLM_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  }
}

RunPowerPipelineLinHypTest_GammaTest_Score_BetaZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSize, start, stop, percentageAssoc, LowLD){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  if(Cat == TRUE){
    source("/home/vlynn/Paper_III_Sims/Scripts/Logistic_ADMM0.r")
  } else {
    source("/home/vlynn/Paper_III_Sims/Scripts/Linear_ADMM0.r")
  }
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_Scores(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, percentageAssoc = percentageAssoc, LowLD = LowLD, TrueScore = FitScore)
    
    if(Cat == TRUE){
      results.FitScore = RunLHT_Cat(allData = allData, FitScore = FitScore, covs = covs)
      results.OtherScore1 = RunLHT_Cat(allData = allData, FitScore = OtherScores[1], covs = covs) 
      results.OtherScore2 = RunLHT_Cat(allData = allData, FitScore = OtherScores[2], covs = covs)
      results.OtherScore3 = RunLHT_Cat(allData = allData, FitScore = OtherScores[3], covs = covs)
    } else {
      results.FitScore = RunLHT_Cont(allData = allData, FitScore = FitScore, covs = covs)
      results.OtherScore1 = RunLHT_Cont(allData = allData, FitScore = OtherScores[1], covs = covs)
      results.OtherScore2 = RunLHT_Cont(allData = allData, FitScore = OtherScores[2], covs = covs)
      results.OtherScore3 = RunLHT_Cont(allData = allData, FitScore = OtherScores[3], covs = covs)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 3)
  Scores.1 = Scores.2 = Scores.3 = Scores.4 = matrix(nrow = n, ncol = 3)
  Betas.1 = Betas.2 = Betas.3 = Betas.4 = matrix(nrow = n, ncol = 2)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  for(ll in 1:n){
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    Scores.1[ll,] = unlist(statsAndPVals[[ll]][[2]])
    Betas.1[ll,] = unlist(statsAndPVals[[ll]][[3]])
    Score.1[1] = unlist(statsAndPVals[[1]][[4]])
    
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[5]])
    Scores.2[ll,] = unlist(statsAndPVals[[ll]][[6]])
    Betas.2[ll,] = unlist(statsAndPVals[[ll]][[7]])
    Score.2[1] = unlist(statsAndPVals[[1]][[8]])
    
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[9]])
    Scores.3[ll,] = unlist(statsAndPVals[[ll]][[10]])
    Betas.3[ll,] = unlist(statsAndPVals[[ll]][[11]])
    Score.3[1] = unlist(statsAndPVals[[1]][[12]])
    
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[13]])
    Scores.4[ll,] = unlist(statsAndPVals[[ll]][[14]])
    Betas.4[ll,] = unlist(statsAndPVals[[ll]][[15]])
    Score.4[1] = unlist(statsAndPVals[[1]][[16]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  Scores = cbind(Scores.1, Scores.2, Scores.3, Scores.4)
  Betas = cbind(Betas.1, Betas.2, Betas.3, Betas.4)
  
  ScoreLabelPValsScores = c(rep(Score.1,3), rep(Score.2,3), rep(Score.3,3), rep(Score.4,3))
  ScoreLabelBetas = c(rep(Score.1,2), rep(Score.2,2), rep(Score.3,2), rep(Score.4,2))
  
  colnames(pVals) = colnames(Scores) = ScoreLabelPValsScores
  colnames(Betas) = ScoreLabelBetas
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(Scores, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
    write.csv(Betas, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_Betas.csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(Scores, file = paste0(path,"/Power_ContY_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
    write.csv(Betas, file = paste0(path,"/Power_ContY_LinHypTest_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_Betas.csv"))
  }
}

RunPowerPipelineDesparseLasso_GammaTest_Score_BetaZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSize, start, stop, percentageAssoc, LowLD){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_Scores(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, percentageAssoc = percentageAssoc, LowLD = LowLD, TrueScore = FitScore)
    
    results.FitScore = RunDesparseLasso(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    results.OtherScore1 = RunDesparseLasso(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs) 
    results.OtherScore2 = RunDesparseLasso(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs)
    results.OtherScore3 = RunDesparseLasso(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 1)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  for(ll in 1:n){
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    Score.1[1] = unlist(statsAndPVals[[1]][[2]])
    
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[3]])
    Score.2[1] = unlist(statsAndPVals[[ll]][[4]])
   
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[5]])
    Score.3[1] = unlist(statsAndPVals[[1]][[6]])
    
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[7]])
    Score.4[1] = unlist(statsAndPVals[[1]][[8]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  
  ScoreLabelPVals = c(Score.1, Score.2, Score.3, Score.4)

  colnames(pVals) = ScoreLabelPVals

  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_DesparseLasso_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_DesparseLasso_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
  }
}

RunPowerPipelineAISPU_GammaTest_Score_BetaZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, seed, ORSize, start, stop, percentageAssoc, LowLD){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  library(rlist)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_Scores(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, percentageAssoc = percentageAssoc, LowLD = LowLD, TrueScore = FitScore)
    
    results.FitScore = RunAISPU(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs, seed = seed)
    results.OtherScore1 = RunAISPU(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs, seed = seed) 
    results.OtherScore2 = RunAISPU(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs, seed = seed)
    results.OtherScore3 = RunAISPU(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs, seed = seed)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  for(ll in 1:n){
    lengthBetas = max(length(statsAndPVals[[ll]][[3]]),length(statsAndPVals[[ll]][[7]]),length(statsAndPVals[[ll]][[11]]), length(statsAndPVals[[ll]][[15]]))
  }
  
  Scores.1 = Scores.2 = Scores.3 = Scores.4 = matrix(nrow = n, ncol = 2)
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 2)
  betas.1 = betas.2 = betas.3 = betas.4 = matrix(nrow = n, ncol=lengthBetas)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  #save the list in case it errors:
  if(Cat == TRUE){
    #write out the p values
    list.save(statsAndPVals, paste0(path,"/TIE_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"Sim",start,"to",numSims,".rdata"))
  } else {
    #write out the p values
    list.save(statsAndPVals, paste0(path,"/TIE_ContY_Prev_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"Sim",start,"to",numSims,".rdata"))
  }
  
  for(ll in 1:n){
    Scores.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[2]])
    betas.1[ll,1:length(unlist(statsAndPVals[[ll]][[3]]))] = unlist(statsAndPVals[[ll]][[3]])
    Score.1 = unlist(statsAndPVals[[1]][[4]])
      
    Scores.2[ll,] = unlist(statsAndPVals[[ll]][[5]])
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[6]])
    betas.2[ll,1:length(unlist(statsAndPVals[[ll]][[7]]))] = unlist(statsAndPVals[[ll]][[7]])
    Score.2 = unlist(statsAndPVals[[1]][[8]])
    
    Scores.3[ll,] = unlist(statsAndPVals[[ll]][[9]])
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[10]])
    betas.3[ll,1:length(unlist(statsAndPVals[[ll]][[11]]))] = unlist(statsAndPVals[[ll]][[11]])
    Score.3 = unlist(statsAndPVals[[1]][[12]])
    
    Scores.4[ll,] = unlist(statsAndPVals[[ll]][[13]])
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[14]])
    betas.4[ll,1:length(unlist(statsAndPVals[[ll]][[15]]))] = unlist(statsAndPVals[[ll]][[15]])
    Score.4 = unlist(statsAndPVals[[1]][[16]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  ScoresAll = cbind(Scores.1, Scores.2, Scores.3, Scores.4)
  
  ScoreLabelPValsScore = c(rep(Score.1,2), rep(Score.2,2), rep(Score.3,2), rep(Score.4,2))
  
  colnames(pVals) = ScoreLabelPValsScore
  colnames(ScoresAll) = ScoreLabelPValsScore
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoresAll, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_ScoreValues.csv"))
    write.csv(betas.1, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.1,".csv"))
    write.csv(betas.2, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.2,".csv"))
    write.csv(betas.3, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.3,".csv"))
    write.csv(betas.4, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.4,".csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoresAll, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_ScoreValues.csv"))
    write.csv(betas.1, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.1,".csv"))
    write.csv(betas.2, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.2,".csv"))
    write.csv(betas.3, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.3,".csv"))
    write.csv(betas.4, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"_GammaSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.4,".csv"))
  }
}

RunPowerPipelineDecorrScore_GammaTest_Score_BetaZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSize, start, stop, percentageAssoc, LowLD){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_Scores(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, percentageAssoc = percentageAssoc, LowLD = LowLD, TrueScore = FitScore)
    
    results.FitScore = RunDecorrScore(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    results.OtherScore1 = RunDecorrScore(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs) 
    results.OtherScore2 = RunDecorrScore(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs)
    results.OtherScore3 = RunDecorrScore(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 1)
  ScoreVal.1 = ScoreVal.2 = ScoreVal.3 = ScoreVal.4 = matrix(nrow = n, ncol = 1)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  for(ll in 1:n){
    ScoreVal.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[2]])
    Score.1[1] = unlist(statsAndPVals[[1]][[3]])
    
    ScoreVal.2[ll,] = unlist(statsAndPVals[[ll]][[4]])
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[5]])
    Score.2[1] = unlist(statsAndPVals[[ll]][[6]])
    
    ScoreVal.3[ll,] = unlist(statsAndPVals[[ll]][[7]])
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[8]])
    Score.3[1] = unlist(statsAndPVals[[1]][[9]])
    
    ScoreVal.4[ll,] = unlist(statsAndPVals[[ll]][[10]])
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[11]])
    Score.4[1] = unlist(statsAndPVals[[1]][[12]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  ScoreVals = cbind(ScoreVal.1, ScoreVal.2, ScoreVal.3, ScoreVal.4)
  
  ScoreLabelPVals = c(Score.1, Score.2, Score.3, Score.4)
  
  colnames(pVals) = colnames(ScoreVals) = ScoreLabelPVals
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoreVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoreVals, file = paste0(path,"/Power_ContY_DecorrScore_Score",FitScore,"_",ORSize,"OR_",percentageAssoc,"Assoc_LowLD",LowLD,"GammaSig",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
  }
}

##########################################################
## Score is assoc (gamma neq 0), beta associated (neq 0)
RunPowerPipelineGLM_GammaTest_Score_BetaNZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat,  ORSizeGene, ORSizeSNP, percentageAssocGene, percentageAssocSNP, LowLDGene, LowLDSNP, start, stop){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  suppressMessages(library(lmtest))
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_ScoreAndRSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSizeGene = ORSizeGene, ORSizeSNP = ORSizeSNP,
                                                 LowLDGene = LowLDGene, LowLDSNP = LowLDSNP, percentAssocGene = percentageAssocGene, percentAssocSNP = percentageAssocSNP, TrueScore = FitScore)
    
    results.FitScore = RunGLM.LRT(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    results.OtherScore1 = RunGLM.LRT(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs) 
    results.OtherScore2 = RunGLM.LRT(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs)
    results.OtherScore3 = RunGLM.LRT(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=4)
  
  reorgData = OrganizeOutputPower(statsAndPVals)
  
  statsAndPValsAll.mat = reorgData[[1]]
  nulValues1.mat = reorgData[[2]]
  altValues1.mat = reorgData[[3]]
  nulValues2.mat = reorgData[[4]]
  altValues2.mat = reorgData[[5]]
  nulValues3.mat = reorgData[[6]]
  altValues3.mat = reorgData[[7]]
  nulValues4.mat = reorgData[[8]]
  altValues4.mat = reorgData[[9]]
  Score.mat = reorgData[[10]]
  
  nameColsStatPVal = c(paste0(Score.mat[[1]]," P Val"), paste0(Score.mat[[1]]," Score"),
                       paste0(Score.mat[[2]]," P Val"), paste0(Score.mat[[2]]," Score"),
                       paste0(Score.mat[[3]]," P Val"), paste0(Score.mat[[3]]," Score"),
                       paste0(Score.mat[[4]]," P Val"), paste0(Score.mat[[4]]," Score"))
  colnames(statsAndPValsAll.mat) = nameColsStatPVal
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/Power_CatPhenos_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues1.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[1]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues2.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[2]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues3.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[3]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues4.mat, file = paste0(path,"/Power_CatPhenos_NullBetas",Score.mat[[4]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues1.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[1]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues2.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[2]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues3.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[3]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues4.mat, file = paste0(path,"/Power_CatPhenos_AltBetas",Score.mat[[4]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  } else {
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path,"/Power_ContPhenos_Phenos_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues1.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[1]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues2.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[2]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues3.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[3]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues4.mat, file = paste0(path,"/Power_ContPhenos_NullBetas",Score.mat[[4]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues1.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[1]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues2.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[2]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues3.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[3]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues4.mat, file = paste0(path,"/Power_ContPhenos_AltBetas",Score.mat[[4]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  }
}

RunPowerPipelineLinHypTest_GammaTest_Score_BetaNZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSizeGene, ORSizeSNP, percentageAssocGene, percentageAssocSNP, LowLDGene, LowLDSNP, start, stop){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score (one for SNPs, one for Score)
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssocGene = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #percentageAssocSNP = percentage of SNPs associated with outcome (5, 15, or 25)
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD (one for SNPs, one for Score)
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  if(Cat == TRUE){
    source("/home/vlynn/Paper_III_Sims/Scripts/Logistic_ADMM0.r")
  } else {
    source("/home/vlynn/Paper_III_Sims/Scripts/Linear_ADMM0.r")
  }
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_ScoreAndRSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSizeGene = ORSizeGene, ORSizeSNP = ORSizeSNP,
                                                 LowLDGene = LowLDGene, LowLDSNP = LowLDSNP, percentAssocGene = percentageAssocGene, percentAssocSNP = percentageAssocSNP, TrueScore = FitScore)
  
    if(Cat == TRUE){
      results.FitScore = RunLHT_Cat(allData = allData, FitScore = FitScore, covs = covs)
      results.OtherScore1 = RunLHT_Cat(allData = allData, FitScore = OtherScores[1], covs = covs) 
      results.OtherScore2 = RunLHT_Cat(allData = allData, FitScore = OtherScores[2], covs = covs)
      results.OtherScore3 = RunLHT_Cat(allData = allData, FitScore = OtherScores[3], covs = covs)
    } else {
      results.FitScore = RunLHT_Cont(allData = allData, FitScore = FitScore, covs = covs)
      results.OtherScore1 = RunLHT_Cont(allData = allData, FitScore = OtherScores[1], covs = covs)
      results.OtherScore2 = RunLHT_Cont(allData = allData, FitScore = OtherScores[2], covs = covs)
      results.OtherScore3 = RunLHT_Cont(allData = allData, FitScore = OtherScores[3], covs = covs)
    }
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 3)
  Scores.1 = Scores.2 = Scores.3 = Scores.4 = matrix(nrow = n, ncol = 3)
  Betas.1 = Betas.2 = Betas.3 = Betas.4 = matrix(nrow = n, ncol = 2)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  for(ll in 1:n){
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    Scores.1[ll,] = unlist(statsAndPVals[[ll]][[2]])
    Betas.1[ll,] = unlist(statsAndPVals[[ll]][[3]])
    Score.1[1] = unlist(statsAndPVals[[1]][[4]])
    
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[5]])
    Scores.2[ll,] = unlist(statsAndPVals[[ll]][[6]])
    Betas.2[ll,] = unlist(statsAndPVals[[ll]][[7]])
    Score.2[1] = unlist(statsAndPVals[[1]][[8]])
    
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[9]])
    Scores.3[ll,] = unlist(statsAndPVals[[ll]][[10]])
    Betas.3[ll,] = unlist(statsAndPVals[[ll]][[11]])
    Score.3[1] = unlist(statsAndPVals[[1]][[12]])
    
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[13]])
    Scores.4[ll,] = unlist(statsAndPVals[[ll]][[14]])
    Betas.4[ll,] = unlist(statsAndPVals[[ll]][[15]])
    Score.4[1] = unlist(statsAndPVals[[1]][[16]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  Scores = cbind(Scores.1, Scores.2, Scores.3, Scores.4)
  Betas = cbind(Betas.1, Betas.2, Betas.3, Betas.4)
  
  ScoreLabelPValsScores = c(rep(Score.1,3), rep(Score.2,3), rep(Score.3,3), rep(Score.4,3))
  ScoreLabelBetas = c(rep(Score.1,2), rep(Score.2,2), rep(Score.3,2), rep(Score.4,2))
  
  colnames(pVals) = colnames(Scores) = ScoreLabelPValsScores
  colnames(Betas) = ScoreLabelBetas
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LHT_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(Scores, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LHT_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
    write.csv(Betas, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_LHT_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas.csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_LHT_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(Scores, file = paste0(path,"/Power_ContY_LHT_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
    write.csv(Betas, file = paste0(path,"/Power_ContY_LHT_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas.csv"))
  }
}

RunPowerPipelineDesparseLasso_GammaTest_Score_BetaNZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSizeGene, ORSizeSNP, percentageAssocGene, percentageAssocSNP, LowLDGene, LowLDSNP, start, stop){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_ScoreAndRSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSizeGene = ORSizeGene, ORSizeSNP = ORSizeSNP,
                                                 LowLDGene = LowLDGene, LowLDSNP = LowLDSNP, percentAssocGene = percentageAssocGene, percentAssocSNP = percentageAssocSNP, TrueScore = FitScore)
    
    results.FitScore = RunDesparseLasso(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    results.OtherScore1 = RunDesparseLasso(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs) 
    results.OtherScore2 = RunDesparseLasso(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs)
    results.OtherScore3 = RunDesparseLasso(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 1)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  for(ll in 1:n){
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    Score.1[1] = unlist(statsAndPVals[[1]][[2]])
    
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[3]])
    Score.2[1] = unlist(statsAndPVals[[ll]][[4]])
    
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[5]])
    Score.3[1] = unlist(statsAndPVals[[1]][[6]])
    
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[7]])
    Score.4[1] = unlist(statsAndPVals[[1]][[8]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  
  ScoreLabelPVals = c(Score.1, Score.2, Score.3, Score.4)
  
  colnames(pVals) = ScoreLabelPVals
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_DesparseLasso_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_DesparseLasso_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
  }
}

RunPowerPipelineAISPU_GammaTest_Score_BetaNZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSizeGene, ORSizeSNP, percentageAssocGene, percentageAssocSNP, LowLDGene, LowLDSNP, start, stop){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  library(rlist)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_ScoreAndRSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSizeGene = ORSizeGene, ORSizeSNP = ORSizeSNP,
                                                 LowLDGene = LowLDGene, LowLDSNP = LowLDSNP, percentAssocGene = percentageAssocGene, percentAssocSNP = percentageAssocSNP, TrueScore = FitScore)
    
    results.FitScore = RunAISPU(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs, seed = seed)
    results.OtherScore1 = RunAISPU(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs, seed = seed) 
    results.OtherScore2 = RunAISPU(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs, seed = seed)
    results.OtherScore3 = RunAISPU(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs, seed = seed)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=10)
  
  n = length(statsAndPVals)
  
  for(ll in 1:n){
    lengthBetas = max(length(statsAndPVals[[ll]][[3]]),length(statsAndPVals[[ll]][[7]]),length(statsAndPVals[[ll]][[11]]), length(statsAndPVals[[ll]][[15]]))
  }
  
  Scores.1 = Scores.2 = Scores.3 = Scores.4 = matrix(nrow = n, ncol = 2)
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 2)
  betas.1 = betas.2 = betas.3 = betas.4 = matrix(nrow = n, ncol=lengthBetas)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  #save the list in case it errors:
  if(Cat == TRUE){
    #write out the p values
    list.save(statsAndPVals, paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,".rdata"))
  } else {
    #write out the p values
    list.save(statsAndPVals, paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,".rdata"))
  }
  
  for(ll in 1:n){
    Scores.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[2]])
    betas.1[ll,1:length(unlist(statsAndPVals[[ll]][[3]]))] = unlist(statsAndPVals[[ll]][[3]])
    Score.1 = unlist(statsAndPVals[[1]][[4]])
    
    Scores.2[ll,] = unlist(statsAndPVals[[ll]][[5]])
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[6]])
    betas.2[ll,1:length(unlist(statsAndPVals[[ll]][[7]]))] = unlist(statsAndPVals[[ll]][[7]])
    Score.2 = unlist(statsAndPVals[[1]][[8]])
    
    Scores.3[ll,] = unlist(statsAndPVals[[ll]][[9]])
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[10]])
    betas.3[ll,1:length(unlist(statsAndPVals[[ll]][[11]]))] = unlist(statsAndPVals[[ll]][[11]])
    Score.3 = unlist(statsAndPVals[[1]][[12]])
    
    Scores.4[ll,] = unlist(statsAndPVals[[ll]][[13]])
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[14]])
    betas.4[ll,1:length(unlist(statsAndPVals[[ll]][[15]]))] = unlist(statsAndPVals[[ll]][[3]])
    Score.4 = unlist(statsAndPVals[[1]][[16]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  ScoresAll = cbind(Scores.1, Scores.2, Scores.3, Scores.4)
  
  ScoreLabelPValsScore = c(rep(Score.1,2), rep(Score.2,2), rep(Score.3,2), rep(Score.4,2))
  
  colnames(pVals) = ScoreLabelPValsScore
  colnames(ScoresAll) = ScoreLabelPValsScore
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoresAll, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_ScoreValues.csv"))
    write.csv(betas.1, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.1,".csv"))
    write.csv(betas.2, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.2,".csv"))
    write.csv(betas.3, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.3,".csv"))
    write.csv(betas.4, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.4,".csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoresAll, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_ScoreValues.csv"))
    write.csv(betas.1, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.1,".csv"))
    write.csv(betas.2, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.2,".csv"))
    write.csv(betas.3, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.3,".csv"))
    write.csv(betas.4, file = paste0(path,"/Power_ContY_AISPU_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Betas",Score.4,".csv"))
  }
}

RunPowerPipelineDecorrScore_GammaTest_Score_BetaNZero = function(chr, gene, numPairs, YPrev, covs, FitScore, Cat, ORSizeGene, ORSizeSNP, percentageAssocGene, percentageAssocSNP, LowLDGene, LowLDSNP, start, stop){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  # Only use when Y is binary
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path = paste0("/home/vlynn/Paper_III_Sims/",gene,"_Results_",numPairs,"Pairs")
  
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  # p-value initialization
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_ScoreAndRSNPs(chr = chr, numPairs = numPairs, simNum = ii, YPrev = YPrev, gene = gene, path = path, covs = covs, ORSizeGene = ORSizeGene, ORSizeSNP = ORSizeSNP,
                                                 LowLDGene = LowLDGene, LowLDSNP = LowLDSNP, percentAssocGene = percentageAssocGene, percentAssocSNP = percentageAssocSNP, TrueScore = FitScore)
    
    results.FitScore = RunDecorrScore(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    results.OtherScore1 = RunDecorrScore(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs) 
    results.OtherScore2 = RunDecorrScore(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs)
    results.OtherScore3 = RunDecorrScore(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=4)
  
  n = length(statsAndPVals)
  
  pVals.1 = pVals.2 = pVals.3 = pVals.4 = matrix(nrow = n, ncol = 1)
  ScoreVal.1 = ScoreVal.2 = ScoreVal.3 = ScoreVal.4 = matrix(nrow = n, ncol = 1)
  Score.1 = Score.2 = Score.3 = Score.4 = matrix(nrow = 1, ncol = 1)
  
  for(ll in 1:n){
    ScoreVal.1[ll,] = unlist(statsAndPVals[[ll]][[1]])
    pVals.1[ll,] = unlist(statsAndPVals[[ll]][[2]])
    Score.1[1] = unlist(statsAndPVals[[1]][[3]])
    
    ScoreVal.2[ll,] = unlist(statsAndPVals[[ll]][[4]])
    pVals.2[ll,] = unlist(statsAndPVals[[ll]][[5]])
    Score.2[1] = unlist(statsAndPVals[[ll]][[6]])
    
    ScoreVal.3[ll,] = unlist(statsAndPVals[[ll]][[7]])
    pVals.3[ll,] = unlist(statsAndPVals[[ll]][[8]])
    Score.3[1] = unlist(statsAndPVals[[1]][[9]])
    
    ScoreVal.4[ll,] = unlist(statsAndPVals[[ll]][[10]])
    pVals.4[ll,] = unlist(statsAndPVals[[ll]][[11]])
    Score.4[1] = unlist(statsAndPVals[[1]][[12]])
  }
  
  pVals = cbind(pVals.1, pVals.2, pVals.3, pVals.4)
  ScoreVals = cbind(ScoreVal.1, ScoreVal.2, ScoreVal.3, ScoreVal.4)
  
  ScoreLabelPVals = c(Score.1, Score.2, Score.3, Score.4)
  
  colnames(pVals) = colnames(ScoreVals) = ScoreLabelPVals
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values
    write.csv(pVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoreVals, file = paste0(path,"/Power_CatY_Prev",YPrev*100,"_DecorrScore_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
  } else {
    write.csv(pVals, file = paste0(path,"/Power_ContY_DecorrScore_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_PValues.csv"))
    write.csv(ScoreVals, file = paste0(path,"/Power_ContY_DecorrScore_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"_Sim",start,"to",numSims,"_Stats.csv"))
  }
}

###############################################################
## R SNPs and Matching both associated, but not on same gene
RunPowerPipelineGLM_GammaTest_Score_BetaNZero = function(chrSNPs, chrMatch, geneSNPs, geneMatch, numPairs, YPrev, covs, FitScore, Cat,  ORSizeGene, ORSizeSNP, percentageAssocGene, percentageAssocSNP, LowLDGene, LowLDSNP, start, stop){
  #Function to determine power of linear hyp testing method when score is associated, gamma = 0 testing only
  #Inputs:
  #chr = chromosome number
  #gene = gene name, in quotes
  #numPairs = number of D/R pairs
  #YPrev = prevalence of binary outcome Y
  #covs - TRUE or FALSE to include covariates or not
  #FitScore - which score we are using to fit the LRT
  #Cat - T or F, is outcome categorical
  #ORSize = Small, Medium, or Large for what OR was used for the associated SNP/score
  # start - simulation number to start at
  # stop - simulation number to stop at
  #percentageAssoc = percentage of SNPs associated with outcome (either 5, 25, 50, 75, or 100) 
  #LowLD = True or FALSE whether the associated SNPs are in low LD or high LD
  #Outputs:
  #No direct outputs, writes scores and pvalues to csv files
  #also writes power values to csv files
  
  library(parallel)
  suppressMessages(library(lmtest))
  
  #always the same
  numSims = stop
  
  #make sure T/F values are correct
  if(covs == "TRUE"){
    covs = TRUE
  } else {
    covs = FALSE
  }
  
  if(Cat == "TRUE"){
    Cat = TRUE
  } else {
    Cat = FALSE
  }
  
  #define path to data
  #for HapGen generated data
  path1 = paste0("/home/vlynn/Paper_III_Sims/",geneSNPs,"_Results_",numPairs,"Pairs")
  path2 = paste0("/home/vlynn/Paper_III_Sims/",geneMatch,"_Results_",numPairs,"Pairs")
    
  #source the needed functions
  source("/home/vlynn/Paper_III_Sims/Scripts/ProjectIISourceFunctions_v2.R")
  
  myList = lapply(start:numSims, rep, times = 1)
  finalOutput = list()
  
  Scores = c("IBS", "Incomp", "AMS", "BinMM")
  fitScoreNum = which(FitScore == Scores)
  OtherScores = Scores[-fitScoreNum]
  
  statsAndPVals = mclapply(myList, function(ii){
    #calc alt phenos
    allData = CalcAltPhenotypeData_ScoreAndRSNPs_DiffGenes(chrSNPs=chrSNPs, chrMatch=chrMatch, numPairs = numPairs, simNum = ii, 
                                                           YPrev = YPrev, geneSNPs=geneSNPs, geneMatch=geneMatch, path1 = path1, path2=path2, covs = covs, 
                                                           ORSizeGene = ORSizeGene, ORSizeSNP = ORSizeSNP, LowLDGene = LowLDGene, LowLDSNP = LowLDSNP, 
                                                           percentAssocGene = percentageAssocGene, percentAssocSNP = percentageAssocSNP, TrueScore = FitScore)
    
    results.FitScore = RunGLM.LRT(allData = allData, FitScore = FitScore, Cat = Cat, covs = covs)
    results.OtherScore1 = RunGLM.LRT(allData = allData, FitScore = OtherScores[1], Cat = Cat, covs = covs) 
    results.OtherScore2 = RunGLM.LRT(allData = allData, FitScore = OtherScores[2], Cat = Cat, covs = covs)
    results.OtherScore3 = RunGLM.LRT(allData = allData, FitScore = OtherScores[3], Cat = Cat, covs = covs)
    
    print(paste0("Simulation ",ii," is complete."))
    
    results.FitScore$Score = FitScore
    results.OtherScore1$Score = OtherScores[1]
    results.OtherScore2$Score = OtherScores[2]
    results.OtherScore3$Score = OtherScores[3]
    
    results = c(results.FitScore, results.OtherScore1, results.OtherScore2, results.OtherScore3)
    
    finalOutput[[ii]] = results
  }, mc.cores=4)
  
  reorgData = OrganizeOutputPower(statsAndPVals)
  
  statsAndPValsAll.mat = reorgData[[1]]
  nulValues1.mat = reorgData[[2]]
  altValues1.mat = reorgData[[3]]
  nulValues2.mat = reorgData[[4]]
  altValues2.mat = reorgData[[5]]
  nulValues3.mat = reorgData[[6]]
  altValues3.mat = reorgData[[7]]
  nulValues4.mat = reorgData[[8]]
  altValues4.mat = reorgData[[9]]
  Score.mat = reorgData[[10]]
  
  nameColsStatPVal = c(paste0(Score.mat[[1]]," P Val"), paste0(Score.mat[[1]]," Score"),
                       paste0(Score.mat[[2]]," P Val"), paste0(Score.mat[[2]]," Score"),
                       paste0(Score.mat[[3]]," P Val"), paste0(Score.mat[[3]]," Score"),
                       paste0(Score.mat[[4]]," P Val"), paste0(Score.mat[[4]]," Score"))
  colnames(statsAndPValsAll.mat) = nameColsStatPVal
  
  #define values for naming conventions
  if(covs == FALSE){
    covariates = "_NoCovs_"
  } else {
    covariates = "_WCovs_"
  }
  
  if(Cat == TRUE){
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues1.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[1]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues2.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[2]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues3.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[3]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues4.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[4]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues1.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[1]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues2.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[2]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues3.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[3]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues4.mat, file = paste0(path2,"/Power_CatPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[4]],"_Prev",YPrev*100,"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  } else {
    #write out the Stats and p values, and Summary stats
    write.csv(statsAndPValsAll.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_Phenos_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues1.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[1]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues2.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[2]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues3.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[3]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(nulValues4.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_NullBetas",Score.mat[[4]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues1.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[1]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues2.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[2]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues3.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[3]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
    write.csv(altValues4.mat, file = paste0(path2,"/Power_ContPhenos_",geneSNPs,"RSNPs_AltBetas",Score.mat[[4]],"_GLM_Score",FitScore,"_",ORSizeGene,"OR_",ORSizeSNP,"ORsnp_",percentageAssocGene,"Assoc_",percentageAssocSNP,"AssocSnp_LowLDGene",LowLDGene,"_LowLDSNP",LowLDSNP,"_GBSig_",covariates,"Sim",start,"to",numSims,"_StatsAndPValues.csv"))
  }
}
