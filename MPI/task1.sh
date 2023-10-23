#!/bin/sh

#SBATCH --account=courses0101
#SBATCH --partition=debug
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --time=00:04:00

module load mpicc

mpicc -fopenmp project2_task1.c -o project2_task1 -lm

export OMP_NUM_THREADS=8

srun mpirun -np 4 ./project2_task1
