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
ORSize <- args[10]
percentageAssoc <- args[11]
lowLDTF <- args[12]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
percentageAssoc = as.integer(percentageAssoc)

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts/MainPipeLineFunctionsToSource_v3.R')

#make sure low LD is T or FALSE
if(lowLDTF == "TRUE"){
	LowLD = TRUE
} else {
	LowLD = FALSE
}

#define score weights if wanted
scoreWeights = c()

#define gammas
null_Gamma = c(0)

############################################################
#for original, unstandardized and unweighted tests
RunPowerPipelineLinHypTestCont_RSNPs(chr = chr, gene = gene, numPairs = numPairs, Gamma = null_Gamma, TrueScore = "IBS.gene", ORSize = ORSize, standardizeScores = FALSE, weightedScores = FALSE, percentageAssoc = percentageAssoc, LowLD = LowLD)
