#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
YPREV=${5?Error: no outcome prevalence given}
ORSizeGene=${6?Error: no OR effect size given}
percentageAssocGene=${7?Error: no % of associated SNPs given}
LowLDGeneTF=${8?Error: no true or false given for low LD}
FitScore=${9?Error: no true score given}
startVal=${10?Error: no start value given}
stopVal=${11?Error: no stop value given}
ORSizeSNP=${12?Error: no OR for R SNPs given}
percentageAssocSNP=${13?Error: No % assoc R SNPs given}
LowLDSNPTF=${14?Error: No T/F LD value SNPs given}
covs=${15?Error: No covariate T/F given}
Cat=${16?Error: No Categorical T/F given}

export CHR
export SS
export SEED
export GENE
export YPREV
export ORSizeGene
export percentageAssocGene
export LowLDGeneTF
export FitScore
export startVal
export stopVal
export ORSizeSNP
export percentageAssocSNP
export LowLDSNPTF
export covs
export Cat

path="/home/vlynn/Paper_III_Sims/Scripts/Power_Scripts"
cd $path

module load R/4.0
Rscript MainPipeline_Power_Scores_DecorrScore_GammaTest_Score_BetaSig.R $CHR $SS $SEED $GENE $YPREV $ORSizeGene $percentageAssocGene $LowLDGeneTF $FitScore $startVal $stopVal $ORSizeSNP $percentageAssocSNP $LowLDSNPTF $covs $Cat --nosave
