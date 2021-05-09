#!/bin/bash
#SBATCH  -p standard --time=12:00:00 --cpus-per-task=22 --mem=72g -o out.%a.txt

echo This is job $[4*$SLURM_ARRAY_TASK_ID]
echo $1

matlab -nodesktop -nosplash -r "disp($1)"