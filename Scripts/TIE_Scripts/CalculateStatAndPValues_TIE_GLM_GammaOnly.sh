#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
YPREV=${5?Error: no outcome prevalence given}
FITSCORE=${6?Error: no score prevalence given}
START=${7?Error: no start value given}
STOP=${8?Error: no stop value given}

export CHR
export SS
export SEED
export GENE
export YPREV
export FITSCORE
export START
export STOP

path="/home/vlynn/Paper_III_Sims/Scripts/TIE_Scripts"
cd $path

module load R/4.0
Rscript MainPipeline_TIE_GLM_GammaOnly.R $CHR $SS $SEED $GENE $YPREV $START $STOP $FITSCORE --nosave
