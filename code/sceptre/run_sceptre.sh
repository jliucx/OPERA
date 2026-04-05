#!/bin/bash

#SBATCH --job-name=repologle
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=250GB
#SBATCH --time=240:00:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=out/repologle_%j.out
#SBATCH --error=err/repologle_%j.err

module load R/4.1.0
Rscript sceptre_repologle.R

echo "Job completed successfully"
