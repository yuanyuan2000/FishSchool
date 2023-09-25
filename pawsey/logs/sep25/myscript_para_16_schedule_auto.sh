#!/bin/sh

#SBATCH --account=courses0101
#SBATCH --partition=debug
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --time=00:04:00

module load gcc

gcc -fopenmp project1_para_schedule_auto.c -o project1_para_schedule_auto -lm

export OMP_NUM_THREADS=16

srun ./project1_para_schedule_auto

