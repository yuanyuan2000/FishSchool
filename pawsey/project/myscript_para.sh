#!/bin/sh

#SBATCH --account=courses0101
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:03:00

module load gcc

gcc -fopenmp project1_para.c -o project1_para -lm

export OMP_NUM_THREADS=2

srun ./project1_para


