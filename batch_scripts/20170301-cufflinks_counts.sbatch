#!/bin/bash -l
#SBATCH --job-name=subread_dabj
#SBATCH --account=y95
#SBATCH --time=6:00:00
#SBATCH --nodes=1

module load subread/1.5.1-gcc

/usr/bin/time -v make star_cufflinks_count
