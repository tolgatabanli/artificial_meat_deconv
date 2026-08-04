#!/bin/bash

#SBATCH --job-name=cattle_dwls
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=380G
#SBATCH --time=2-0:00:00
#SBATCH --partition=exbio-cpu
#SBATCH --output=job_%x_%j.out
#SBATCH --error=job_%x_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

Rscript --version
Rscript dwls_sig_and_deconv.R
