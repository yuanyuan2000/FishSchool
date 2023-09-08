## sequential
- `project1_seq.c` is a basic sequential C program for this simulation, it can be complied by `gcc -o project1_seq project1_seq.c -lm`.
- The result will be output in the `results1.txt` and we can use the `plot_results.py` to generate the plots

## parallel
- `project1_para.c` is a multi-threaded code using OpenMP for this problem, but it still not work (can compile and run, but much more slower than the basic sequential C program). It can be compiled by `gcc -o project1_para project1_para.c -lm -fopenmp`
- The result will be output in the `results1.txt` and we can use the `plot_results.py` to generate the plots