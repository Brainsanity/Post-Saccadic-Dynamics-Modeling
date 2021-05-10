#!/bin/bash
#SBATCH -o job_out_test.$1.txt -t 00:05:00

s=hostname
module load matlab/r2020b
matlab -nodesktop -nosplash -r "disp($1); disp($s)"