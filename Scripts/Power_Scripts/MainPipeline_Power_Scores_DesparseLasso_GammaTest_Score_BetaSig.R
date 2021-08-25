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
ORSizeGene <- args[11]
percentageAssocGene <- args[12]
LowLDGeneTF <- args[13]
FitScore <- args[14]
startVal <- args[15]
stopVal <- args[16]
ORSizeSNP <- args[17]
percentageAssocSNP <- args[18]
LowLDSNPTF <- args[19]
covs <- args[20]
Cat <- args[21]

chr = as.integer(chr)
seed = as.integer(seed)
numPairs = as.integer(numPairs)
YPrev = as.integer(YPrev)/100
percentageAssocGene = as.integer(percentageAssocGene)
startVal = as.integer(startVal)
stopVal = as.integer(stopVal)
percentageAssocSNP = as.integer(percentageAssocSNP)

#set random seed for reproducibility
set.seed(seed)
#source needed functions
source('/home/vlynn/Paper_III_Sims/Scripts/MainPipeLineFunctionsToSource_LinHypTest.R')

#make sure low LD is T or FALSE
if(LowLDGeneTF == "TRUE"){
	LowLDGene = TRUE
} else {
	LowLDGene = FALSE
}

if(LowLDSNPTF == "TRUE"){
	LowLDSNP = TRUE
} else {
	LowLDSNP = FALSE
}

############################################################
#for original, unstandardized and unweighted tests
#Low LD SNPs associated
RunPowerPipelineDesparseLasso_GammaTest_Score_BetaNZero(chr = chr, gene = gene, numPairs = numPairs, YPrev = YPrev, covs = covs, FitScore = FitScore, Cat = Cat, ORSizeGene = ORSizeGene, ORSizeSNP = ORSizeSNP, percentageAssocGene = percentageAssocGene, percentageAssocSNP = percentageAssocSNP, LowLDGene = LowLDGene, LowLDSNP = LowLDSNP, start = startVal, stop = stopVal)