############################################################
### Power Simulation code
############################################################
#edit 5/21/20 added in ability to standardize scores and weight scores

#to obtain the arguments from the bash files (chr, ss)
args <- commandArgs()

chr <- args[6]
numPairs <- args[7]
seed <- args[8]
gene <- args[9]
YPrev <- args[10]
FitScore <- args[11]
covs <- args[12]
Cat <- args[13]
startVal <- args[14]
stopVal <- args[15]
ORSize <- args[16]
percentageAssoc <- args[17]
lowLDTF <- args[18]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
YPrev = as.integer(YPrev)/100
start = as.integer(startVal)
stop = as.integer(stopVal)
percentageAssoc = as.integer(percentageAssoc)

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_III_Sims/Scripts/MainPipeLineFunctionsToSource_LinHypTest.R')

#make sure low LD is T or FALSE
if(lowLDTF == "TRUE"){
	LowLD = TRUE
} else {
	LowLD = FALSE
}

############################################################
RunPowerPipelineLinHypTest_GammaTest_Score_BetaZero(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, covs = covs, FitScore = FitScore, Cat=Cat, ORSize = ORSize, start = startVal, stop = stopVal, percentageAssoc = percentageAssoc, LowLD = LowLD)
