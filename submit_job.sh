#!/bin/bash
#SBATCH --job-name=fread_task
#SBATCH --output=fread_task_%j.out
#SBATCH --error=fread_task_%j.err
#SBATCH --mem=8G
#SBATCH --time=00:10:00

module load R/4.5.1-mkl
Rscript job_script.R
