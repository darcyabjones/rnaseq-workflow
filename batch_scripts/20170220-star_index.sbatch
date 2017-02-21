#!/bin/bash -l
#SBATCH --job-name=star_dabj
#SBATCH --account=y95
#SBATCH --time=24:00:00
#SBATCH --nodes=1

module load star/2.5.0a-gcc

/usr/bin/time -v make star_index
