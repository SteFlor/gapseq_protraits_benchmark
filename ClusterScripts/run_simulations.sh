#!/bin/bash

# SLURM job settings
#SBATCH --array 0-19
#SBATCH --job-name=simulate
#SBATCH --output=output_log/job_output_%A_%a.out
#SBATCH --error=error_log/job_error_%A_%a.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=62G

source ~/.bashrc
micromamba activate gapseq_benchmark

FRACTION_OF_OPTIMUM=0.01
MODEL_DIR="../AllModels"
NO_MODEL_DIR="../AllGapseqOutputsButModels"
OUTDIR="../DataFramesFromCluster"
#move all the files of the folder AllModels that are not the finished xml model into a new folder called AllGapseqOutputsButModels
mv $MODEL_DIR/*[!genomic.xml] $NO_MODEL_DIR/

# Script to calculate the total number of batches
BATCH_SIZE=200  # Change this to your desired batch size
NUM_MODELS=$(ls -1q $MODEL_DIR | wc -l)
NUM_BATCHES=$((($NUM_MODELS + $BATCH_SIZE - 1) / $BATCH_SIZE))

# Run the Python script
python3 simulateFVA_FullDiet.py $SLURM_ARRAY_TASK_ID $BATCH_SIZE $MODEL_DIR $FRACTION_OF_OPTIMUM $OUTDIR

