#!/bin/bash
#SBATCH  -p standard --time=12:00:00 --cpus-per-task=22 --mem=72g

hostname >> ../logs/condition-$1.txt

matlab -nodesktop -nosplash -r "RunSimulation($1, 11.5);" >> ../logs/condition-$1.txt 2>>&1