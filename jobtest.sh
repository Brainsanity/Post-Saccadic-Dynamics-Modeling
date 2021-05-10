#!/bin/bash
#SBATCH -o job_out_test.$1.txt -t 00:05:00

hostname
s=`expr hostname`
echo $s
module load matlab/r2020b
matlab -nodesktop -nosplash -r "disp(2*$1); disp($s)"