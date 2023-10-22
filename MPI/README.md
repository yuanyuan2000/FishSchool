### project2_task1.c: just for task1
- complie: `mpicc -fopenmp project2_task1.c -o project2_task1 -lm`
- run: `mpirun -np 4 ./project2_task1`

### project2_task2.c: for task2, combine both MPI with openMP
- complie: `mpicc -fopenmp project2_task2.c -o project2_task2.out -lm`
- run: `mpirun -np 4 ./project2_task2.out 1`
(4 is the number of process, 1 is the option, only the option 1-5 work well, the option 6 has some bugs)

### project2_task2_simplify.c: for task2, only the MPI, no openMP
- complie: `mpicc -fopenmp project2_task2_simplify.c -o project2_task2_simplify.out -lm`
- run: `mpirun -np 4 ./project2_task2_simplify.out 1`