#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
YPREV=${5?Error: no outcome prevalence given}
SCORE=${6?Error: no score prevalence given}
START=${7?Error: no start value given}
STOP=${8?Error: no stop value given}
ORSIZE=${9?Error: no OR size value given}
PERCENTASSOC=${10?Error: no percent assoc value given}
LOWLDTF=${11?Error: no low LD T/F value given}
COVS=${12?Error: no T/F covs value given}
CAT=${13?Error: no T/F cat value given}

export CHR
export SS
export SEED
export GENE
export YPREV
export SCORE
export START
export STOP
export ORSIZE
export PERCENTASSOC
export LOWLDTF
export COVS
export CAT

path="/home/vlynn/Paper_III_Sims/Scripts/TIE_Scripts"
cd $path

module load R/4.0
Rscript MainPipeline_TIE_DecorrScore_GammaOnly_BetaSig.R $CHR $SS $SEED $GENE $YPREV $SCORE $COVS $CAT $START $STOP $ORSIZE $PERCENTASSOC $LOWLDTF --nosave
