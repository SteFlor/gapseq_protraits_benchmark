#!/bin/bash

#SBATCH --array 0-3616
#SBATCH --job-name=gap
#SBATCH --error=error_log/gap_%A_%a.err
#SBATCH --output=output_log/gap_%A_%a.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-23:59:00
#SBATCH --mem=64G


echo "starting"
source ~/.bashrc
micromamba activate gapseq_env
cd ../AllModels
files=(../Genomes/*.fna)
file="${files[$SLURM_ARRAY_TASK_ID]}"
echo $file

gapseq doall $file
