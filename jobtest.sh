#!/bin/bash
#SBATCH -o job_out_test.$1.txt -t 00:05:00

echo This is job $SLURM_ARRAY_TASK_ID