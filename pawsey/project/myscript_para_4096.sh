#!/bin/sh

#SBATCH --account=courses0101
#SBATCH --partition=debug
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=64
#SBATCH --time=00:03:00

module load gcc

gcc -fopenmp -O3 project1_para.c -o project1_para -lm

export OMP_NUM_THREADS=4096

srun ./project1_para

