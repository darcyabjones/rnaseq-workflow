#!/bin/bash -l
#SBATCH --job-name=star_dabj
#SBATCH --account=y95
#SBATCH --time=24:00:00
#SBATCH --nodes=1

module load star/2.5.0a-gcc
module load samtools/1.3.1-gcc

/usr/bin/time -v make star_align
