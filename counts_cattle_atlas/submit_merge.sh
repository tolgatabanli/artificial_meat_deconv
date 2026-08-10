#!/bin/bash

#SBATCH --job-name=cattle_merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=260G
#SBATCH --time=0:45:00
#SBATCH --partition=exbio-cpu
#SBATCH --output=job_%x_%j.out
#SBATCH --error=job_%x_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

Rscript --version
export TMPDIR=/tmp/t.tabanli/Rtmp
Rscript sample_and_merge.R
