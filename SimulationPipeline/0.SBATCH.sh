#!/bin/bash
#SBATCH --job-name=SC_OCS.%j
#SBATCH --mail-type=END
#SBATCH --mail-user=deamorimpeixotom@ufl.edu
#SBATCH --account=mresende
#SBATCH --qos=mresende-b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --time=96:00:00
#SBATCH --output=1.tmp/alphasim.%a.array.%A.out
#SBATCH --error=2.error/alphasim.%a.array.%A.err
#SBATCH --array=1-25

module purge; module load R

INPUT=$(sed -n ${SLURM_ARRAY_TASK_ID}p INPUT.FILE.txt)

Rscript RUNME.R $INPUT ${rep} 
