#!/bin/bash

#SBATCH --job-name=simulation
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=24:00:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=rcc_out/simulation%j.out
#SBATCH --error=rcc_err/simulation%j.err

module load R/4.4.1
