#!/bin/bash -l
#SBATCH --job-name=subread_dabj
#SBATCH --account=y95
#SBATCH --time=6:00:00
#SBATCH --ntasks=12
#SBATCH --ntasks-per-core=1

module load bedtools/2.25.0

/usr/bin/time -v make -j 12 star_coverage
