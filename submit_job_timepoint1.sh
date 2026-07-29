#!/bin/bash

#SBATCH --job-name=job_timepoint1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=90G
#SBATCH --time=1:00:00
#SBATCH --partition=exbio-cpu
#SBATCH --output=job_all_%j.out
#SBATCH --error=job_all_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

Rscript --version

Rscript job_timepoint1.R
