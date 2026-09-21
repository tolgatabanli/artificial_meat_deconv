#!/bin/bash

#SBATCH --job-name=scaden_messmer
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=120G
#SBATCH --time=3:00:00
#SBATCH --partition=exbio-gpu
#SBATCH --output=job_%x_%j.out
#SBATCH --error=job_%x_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

source /opt/mambaforge/etc/profile.d/conda.sh
conda activate r421
export R_LIBS_USER=/nfs/home/students/t.tabanli/.conda/envs/r421/lib/R/library
Rscript --version
Rscript scaden_train.R
