#!/bin/bash -l
#SBATCH --job-name=cufflinks_dabj
#SBATCH --account=y95
#SBATCH --time=24:00:00
#SBATCH --nodes=1

module unload PrgEnv-cray/5.2.82
module load PrgEnv-gnu/5.2.82

module load cufflinks/2.2.1-bin

/usr/bin/time -v make star_cufflinks_ass
