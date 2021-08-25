#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
YPREV=${5?Error: no outcome prevalence given}
ORSize=${6?Error: no OR effect size given}
percentageAssoc=${7?Error: no % of associated SNPs given}
lowLDTF=${8?Error: no true or false given for low LD}
startVal=${9?Error: no start value given}
stopVal=${10?Error: no stop value given}
FitScore=${11?Error: no Fit Score given}
covs=${12?Error: No T/F covs value given}
Cat=${13?Error: no T/F Cat value given}

export CHR
export SS
export SEED
export GENE
export YPREV
export ORSize
export percentageAssoc
export lowLDTF
export FitScore
export startVal
export stopVal

path="/home/vlynn/Paper_III_Sims/Scripts/Power_Scripts"
cd $path

module load R/4.0
Rscript MainPipeline_Power_Scores_DecorrScore_GammaTest_Score_BetaZero.R $CHR $SS $SEED $GENE $YPREV $FitScore $covs $Cat $startVal $stopVal $ORSize $percentageAssoc $lowLDTF --nosave 
