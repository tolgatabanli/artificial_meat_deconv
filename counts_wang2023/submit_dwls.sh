#!/bin/bash

#SBATCH --job-name=wang_dwls
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --time=5:00:00
#SBATCH --partition=exbio-cpu
#SBATCH --output=job_%x_%j.out
#SBATCH --error=job_%x_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

Rscript --version
Rscript dwls.R
