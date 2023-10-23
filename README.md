# Project 1
The code of project 1 is mainly in the folder `sequential` and `parallel`

## sequential
- `project1_seq.c` is a basic sequential C program for this simulation, it can be complied by `gcc -o project1_seq project1_seq.c -lm`.
- The result will be output in the `results1.txt` and we can use the `plot_results.py` to generate the plots

## parallel
- `project1_para.c` is a multi-threaded code using OpenMP for this problem, but it still not work (can compile and run, but much more slower than the basic sequential C program). It can be compiled by `gcc -o project1_para project1_para.c -lm -fopenmp`
- The result will be output in the `results1.txt` and we can use the `plot_results.py` to generate the plots

# Project 2
The code of project 2 is mainly in the folder `MPI`

## delivery 1
The code of delivery 1 is `project2_task1.c`:
- complie: `mpicc -fopenmp project2_task1.c -o project2_task1 -lm`
- run: `mpirun -np 4 ./project2_task1`

## delivery 2
The code of delivery 2 is `project2_task2.c`, it combines both MPI with openMP.
- complie: `mpicc -fopenmp project2_task2.c -o project2_task2.out -lm`
- run: `mpirun -np 4 ./project2_task2.out 1`
(4 is the number of process, 1 is the option)

And also have another simply version is `project2_task2_simplify.c`, which only uses the MPI without openMP.
- complie: `mpicc -fopenmp project2_task2_simplify.c -o project2_task2_simplify.out -lm`
- run: `mpirun -np 4 ./project2_task2_simplify.out 1`