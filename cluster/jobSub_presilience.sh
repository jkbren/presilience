#!/bin/bash
#SBATCH --job-name=presilience
#SBATCH -o presilience%a.out     
#SBATCH -e presilience%a.err
#SBATCH -c 16
#SBATCH -p netsi_standard
#SBATCH --time=23:00:00
#SBATCH --mem=135GB

work=/scratch/klein.br/presilience/cluster/
cd $work

declare -a commands
readarray -t commands < batchPresil.conf # Exclude newline.
eval ${commands[$SLURM_ARRAY_TASK_ID]}
