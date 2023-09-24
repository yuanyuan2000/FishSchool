#!/bin/sh

#SBATCH --account=courses0101
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00

module load gcc

# Compile project1_seq.c
gcc -fopenmp project1_seq.c -o project1_seq -lm

# Run the compiled program
srun ./project1_seq

