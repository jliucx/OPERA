#!/bin/bash

#SBATCH --job-name=split
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=250GB
#SBATCH --time=360:30:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=rcc_out/topic_%j.out
#SBATCH --error=rcc_err/topic_%j.err

module load R/4.4.1
Rscript fast_topic.R

echo "Job completed successfully"

