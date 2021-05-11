#!/bin/bash
#SBATCH --time=12:00:00 --cpus-per-task=22 --mem=72g --mail-type=all -p standard

for iCond in $*
do
	hostname >> ../logs/condition-$iCond.txt
	matlab -nodesktop -nosplash -r "RunSimulation($iCond, 11.5); exit;" >> ../logs/condition-$iCond.txt 2>&1
done