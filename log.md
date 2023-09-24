feifan@setonix-02:~/project> sbatch myscript_para_2.sh 
Submitted batch job 4614293
feifan@setonix-02:~/project> sbatch myscript_para_4.sh 
Submitted batch job 4614294
feifan@setonix-02:~/project> sbatch myscript_para_64.sh 
Submitted batch job 4614295
feifan@setonix-02:~/project> sbatch myscript_para_16.sh 
Submitted batch job 4614296
feifan@setonix-02:~/project> ls
myscript_para_16.sh  project1_para		   slurm-4614261.out
myscript_para_2.sh   project1_para.c		   slurm-4614266.out
myscript_para_4.sh   project1_seq		   slurm-4614287.out
myscript_para_64.sh  project1_seq.c		   slurm-4614293.out
myscript_para.sh     results1_seq.txt
myscript_seq.sh      singleCore_4194304Fishes.out
feifan@setonix-02:~/project> ls
myscript_para_16.sh  project1_para		   slurm-4614261.out
myscript_para_2.sh   project1_para.c		   slurm-4614266.out
myscript_para_4.sh   project1_seq		   slurm-4614287.out
myscript_para_64.sh  project1_seq.c		   slurm-4614293.out
myscript_para.sh     results1_seq.txt
myscript_seq.sh      singleCore_4194304Fishes.out
feifan@setonix-02:~/project> ls
myscript_para_16.sh  project1_para		   slurm-4614261.out
myscript_para_2.sh   project1_para.c		   slurm-4614266.out
myscript_para_4.sh   project1_seq		   slurm-4614287.out
myscript_para_64.sh  project1_seq.c		   slurm-4614293.out
myscript_para.sh     results1_seq.txt		   slurm-4614295.out
myscript_seq.sh      singleCore_4194304Fishes.out
feifan@setonix-02:~/project> pwd
/home/feifan/project
feifan@setonix-02:~/project> mv /home/feifan/myscript_para_128.sh /home/feifan/project/
feifan@setonix-02:~/project> ls
myscript_para_128.sh  myscript_seq.sh	singleCore_4194304Fishes.out
myscript_para_16.sh   project1_para	slurm-4614261.out
myscript_para_2.sh    project1_para.c	slurm-4614266.out
myscript_para_4.sh    project1_seq	slurm-4614287.out
myscript_para_64.sh   project1_seq.c	slurm-4614293.out
myscript_para.sh      results1_seq.txt	slurm-4614295.out
feifan@setonix-02:~/project> sbatch myscript_para_128.sh 
Submitted batch job 4614299
feifan@setonix-02:~/project> ls
myscript_para_128.sh  project1_para		    slurm-4614266.out
myscript_para_16.sh   project1_para.c		    slurm-4614287.out
myscript_para_2.sh    project1_seq		    slurm-4614293.out
myscript_para_4.sh    project1_seq.c		    slurm-4614294.out
myscript_para_64.sh   results1_seq.txt		    slurm-4614295.out
myscript_para.sh      singleCore_4194304Fishes.out
myscript_seq.sh       slurm-4614261.out
feifan@setonix-02:~/project> 

feifan@setonix-02:~/project> sbatch myscript_para_16.sh 
Submitted batch job 4614310
feifan@setonix-02:~/project> ls
myscript_para_128.sh  project1_para		    slurm-4614266.out
myscript_para_16.sh   project1_para.c		    slurm-4614287.out
myscript_para_2.sh    project1_seq		    slurm-4614293.out
myscript_para_4.sh    project1_seq.c		    slurm-4614294.out
myscript_para_64.sh   results1_seq.txt		    slurm-4614295.out
myscript_para.sh      singleCore_4194304Fishes.out  slurm-4614296.out
myscript_seq.sh       slurm-4614261.out		    slurm-4614299.out
feifan@setonix-02:~/project> 

feifan@setonix-02:~/project> sbatch myscript_para_16.sh 
Submitted batch job 4614341
