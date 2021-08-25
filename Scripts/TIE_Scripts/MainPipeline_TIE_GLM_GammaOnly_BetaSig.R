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
FitScore <- args[11]
startVal <- args[12]
stopVal <- args[13]
ORSize <- args[14]
percentageAssoc <- args[15]
lowLDTF <- args[16]

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
#for unweighted and unstandardized scores
#no covs, continuous outcome
#covs, continuous outcome
covs = TRUE
#Cat = FALSE
#RunTIEPipelineLinHypTest_GammaTest_GLM_BetaSig(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, ORSize = ORSize, percentageAssoc = percentageAssoc, LowLD = LowLD, covs = covs, FitScore = FitScore, Cat = Cat, start = start, stop = stop)
#covs, binary outcome
Cat = TRUE
RunTIEPipelineLinHypTest_GammaTest_GLM_BetaSig(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, ORSize = ORSize, percentageAssoc = percentageAssoc, LowLD = LowLD, covs = covs, FitScore = FitScore, Cat = Cat, start = start, stop = stop)
