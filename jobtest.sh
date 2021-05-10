#!/bin/bash
#SBATCH -t 00:05:00

# s=$(echo `hostname`)
# s=$(hostname)

hostname >> ../condition-$1.txt

# module load matlab/r2020b
# matlab -nodesktop -nosplash -r "disp(2*$1); disp($s)" >> ../logs/condition-$1.txt 2>>&1