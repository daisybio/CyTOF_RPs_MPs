#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytof_RPs_MPs_DE_DNA_normalized
#SBATCH --error=logging/%x_%j.err
#SBATCH --output=logging/%x.%j.txt
#SBATCH --mem-per-cpu=70G
#SBATCH --partition=shared-cpu
#SBATCH --time=1:00:00

Rscript run_DE.R /nfs/data/Bongiovanni-KrdIsar-platelets/Cyanus_RPsMPs/data/sce_DNA2_RPs_MPs_rest_DNA.rds DNA2 