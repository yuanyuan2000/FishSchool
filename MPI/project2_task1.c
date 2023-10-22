#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#define NUM_FISH 262144
#define MAX_FISH_X 100
#define MAX_FISH_Y 100
#define MIN_FISH_X -MAX_FISH_X
#define MIN_FISH_Y -MAX_FISH_Y
#define INIT_WEIGHT 10
#define SIMULATION_STEPS 500
#define MAX_THREAD 8
#define FISH_PER_THREAD 150000

struct Fish
{
    float x;          // Current x position
    float y;          // Current y position
    float plan_x;     // Planned x position
    float plan_y;     // Planned y position
    float weight;     // Current weight of the fish
    float delta_f;    // Change in fish's position
};

// Function for phase one of the base configuration
float phaseOneBase(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

    // Parallel loop for calculating plan positions and delta_f for each fish
#pragma omp parallel for private(seed) reduction(max : max_delta_f)
    for (int i = 0; i < NUM_FISH; i++)
    {
        seed = 123456789 + omp_get_thread_num(); // Initialize seed for each thread
        fishes[i].plan_x = fishes[i].x + ((float)rand_r(&seed) / RAND_MAX) * 0.2 - 0.1;
        fishes[i].plan_y = fishes[i].y + ((float)rand_r(&seed) / RAND_MAX) * 0.2 - 0.1;

        // Ensure planned positions are within limits
        if (fishes[i].plan_x > MAX_FISH_X)
            fishes[i].plan_x = MAX_FISH_X;
        if (fishes[i].plan_x < MIN_FISH_X)
            fishes[i].plan_x = MIN_FISH_X;
        if (fishes[i].plan_y > MAX_FISH_Y)
            fishes[i].plan_y = MAX_FISH_Y;
        if (fishes[i].plan_y < MIN_FISH_Y)
            fishes[i].plan_y = MIN_FISH_Y;

        // Calculate delta_f for the fish
        fishes[i].delta_f = sqrt(fishes[i].plan_x * fishes[i].plan_x + fishes[i].plan_y * fishes[i].plan_y) - sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);
        if (fishes[i].delta_f > max_delta_f)
        {
            max_delta_f = fishes[i].delta_f;
        }
    }
    return max_delta_f;
}

// Function for phase two of the base configuration
float phaseTwoBase(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

    // Parallel loop for updating fish positions and calculating distances
#pragma omp parallel for reduction(+ : total_distance_weighted, total_distance)
    for (int i = 0; i < NUM_FISH; i++)
    {
        float planWeight = fishes[i].weight + fishes[i].delta_f / max_delta_f;
        
        // Ensure weight does not exceed 2 times the initial weight
        if (planWeight > 2 * INIT_WEIGHT)
        {
            fishes[i].weight = 2 * INIT_WEIGHT;
        }
        
        // Update fish's position and weight if planWeight is greater
        if (planWeight > fishes[i].weight)
        {
            fishes[i].x = fishes[i].plan_x;
            fishes[i].y = fishes[i].plan_y;
            fishes[i].weight = planWeight;
        }
        
        // Calculate current distance for the fish
        float curr_distance = sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);
        total_distance_weighted += curr_distance * fishes[i].weight;
        total_distance += curr_distance;
    }

    return total_distance_weighted / total_distance;
}

// Function to initialize the fish positions and weights
void initFish(struct Fish *fishes)
{
    // Parallel loop for initializing fish
#pragma omp parallel for
    for (int i = 0; i < NUM_FISH; i++)
    {
        unsigned int seed = omp_get_thread_num(); // Get unique seed for each thread
        fishes[i].x = (rand_r(&seed) % (2 * MAX_FISH_X + 1)) - MAX_FISH_X;
        fishes[i].y = (rand_r(&seed) % (2 * MAX_FISH_X + 1)) - MAX_FISH_Y;
        fishes[i].weight = INIT_WEIGHT; // Set initial weight for each fish
    }
}

void writeFishDataToFile(struct Fish *fishes, const char *filename)
{
    FILE *fp = fopen(filename, "w+");
    for (int i = 0; i < NUM_FISH; i++)
    {
        fprintf(fp, "%f %f\n", fishes[i].x, fishes[i].y);
    }
    fclose(fp);
}

int main()
{
    int process_id, number_of_processes;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    struct Fish *fishes = (struct Fish *)malloc(sizeof(struct Fish) * NUM_FISH);

    // Generate fish data only in master process
    if (process_id == 0)
    {
        omp_set_num_threads(MAX_THREAD);
        srand(time(0));
        initFish(fishes);
        writeFishDataToFile(fishes, "out1.txt");
    }

    // Distribute fish data equally among all processes
    int fishes_per_process = NUM_FISH / number_of_processes;
    struct Fish *local_fishes = (struct Fish *)malloc(sizeof(struct Fish) * fishes_per_process);
    MPI_Scatter(fishes, fishes_per_process * sizeof(struct Fish), MPI_BYTE, local_fishes, fishes_per_process * sizeof(struct Fish), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Send back the data to the master
    MPI_Gather(local_fishes, fishes_per_process * sizeof(struct Fish), MPI_BYTE, fishes, fishes_per_process * sizeof(struct Fish), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Write to file after master gets back all the data
    if (process_id == 0)
    {
        writeFishDataToFile(fishes, "out2.txt");
    }

    free(local_fishes);
    free(fishes);
    MPI_Finalize();
    return 0;
}

