############################################################
### Type I Error Simulation code
############################################################
#to obtain the arguments from the bash files (chr, ss)
args <- commandArgs()

chr <- args[6]
numPairs <- args[7]
seed <- args[8]
gene <- args[9]
YPrev <- args[10]
startVal <- args[11]
stopVal <- args[12]
FitScore <- args[13]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
YPrev = as.integer(YPrev)/100
startVal = as.integer(startVal)
stopVal = as.integer(stopVal)

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_III_Sims/Scripts/MainPipeLineFunctionsToSource_LinHypTest.R')

############################################################
#for unweighted and unstandardized scores
#no covs, continuous outcome
covs = FALSE
Cat = FALSE
RunTIEPipelineLinHypTest_GammaTest_GLM(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, covs = covs, FitScore = FitScore, Cat = Cat, start = startVal, stop = stopVal)
#no covs, binary outcome
Cat = TRUE
RunTIEPipelineLinHypTest_GammaTest_GLM(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, covs = covs, FitScore = FitScore, Cat = Cat, start = startVal, stop = stopVal)
#covs, continuous outcome
covs = TRUE
Cat = FALSE
RunTIEPipelineLinHypTest_GammaTest_GLM(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, covs = covs, FitScore = FitScore, Cat = Cat, start = startVal, stop = stopVal)
#covs, binary outcome
Cat = TRUE
RunTIEPipelineLinHypTest_GammaTest_GLM(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, covs = covs, FitScore = FitScore, Cat = Cat, start = startVal, stop = stopVal)