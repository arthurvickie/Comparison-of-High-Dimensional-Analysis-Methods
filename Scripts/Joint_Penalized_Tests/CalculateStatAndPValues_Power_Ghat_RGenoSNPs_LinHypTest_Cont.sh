#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi

CHR=${1?Error: no chromosome given}
SS=${2?Error: no sample size given}
SEED=${3?Error: no random seed given}
GENE=${4?Error: no gene name given}
ORSize=${5?Error: no OR effect size given}
percentageAssoc=${6?Error: no % of associated SNPs given}
lowLDTF=${7?Error: no true or false given for low LD}

export CHR
export SS
export SEED
export GENE
export ORSize
export percentageAssoc
export lowLDTF

path="/home/vlynn/Paper_II_Sims/HapGen_Files/Scripts"
cd $path

module load R/3.5.3
Rscript MainPipeline_Power_Ghat_RGenoSNPs_LinHypTest_Cont.R $CHR $SS $SEED $GENE $ORSize $percentageAssoc $lowLDTF --nosave