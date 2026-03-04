#!/bin/bash
#SBATCH --job-name=gets_bbd_SD
#SBATCH --mem=8G 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-user=lukas.vashold@wu.ac.at
#SBATCH --mail-type=ALL
#SBATCH --array=1-5600%400	         # index values

echo "Running job $SLURM_JOB_NAME with ID $SLURM_JOB_ID, array task ID $SLURM_ARRAY_TASK_ID"
echo "Allocated hosts: $SLURM_JOB_NODELIST"

# Set up R environmentexport R_LIBS=~/R_LIBS
. /opt/apps/2023-04-11_lmod.bash
module load r-4.5.1-gcc-12.2.0-tfchyxd

# Run the R script with array task ID and count
srun Rscript 01_estimation_gets_bisam_comparison_SD.R $SLURM_ARRAY_TASK_ID

