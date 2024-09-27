#!/bin/bash
#SBATCH --job-name="Schisto"
#SBATCH --nodes=1
#SBATCH --error=error_%j.txt       # Error file
#SBATCH --ntasks=1                 
#SBATCH --cpus-per-task=16         # Number of CPU cores per task 
#SBATCH --time=01:00:00
#SBATCH --partition=rome

module load 2023
module load R/4.3.2-gfbf-2023a

export R_LIBS=$HOME/rpackages

Rscript 05_Run_model.R
