#!/bin/bash
#SBATCH --array 0-63
#SBATCH --job-name=gapseqFinalMAGs-%A
#SBATCH --error=gapseqFinalMAGs-%A.err
#SBATCH --output=gapseqFinalMAGs-%A.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-23:59:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

source ~/.bashrc
micromamba activate gapseq_env
files=(FinalMAGs/*.fna)
file="${files[$SLURM_ARRAY_TASK_ID]}"
echo $file

gapseq doall $file
