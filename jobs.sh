#!/bin/bash
#SBATCH  -p standard --time=12:00:00 --cpus-per-task=22 --mem=72g -o ../out_01.txt

hostname

matlab -nodesktop -nosplash -r "disp($1)"