#!/bin/bash

#SBATCH --job-name=messmer_dwls_inj1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --partition=exbio-cpu
#SBATCH --output=job_%x_%j.out
#SBATCH --error=job_%x_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

Rscript --version
Rscript dwls_inj1.R
