#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h> // Include the OpenMP header

#define NUM_FISH 4194304
#define MAX_FISH_X 100
#define MAX_FISH_Y 100
#define MIN_FISH_X -MAX_FISH_X
#define MIN_FISH_Y -MAX_FISH_Y
#define INIT_WEIGHT 10
#define SIMULATION_STEPS 500
#define MAX_THREAD 16
struct Fish
{
    float x;
    float y;
    float plan_x;
    float plan_y;
    float weight;
    float delta_f;
};

int main()
{
    omp_set_num_threads(MAX_THREAD); // Set the number of threads

#pragma omp parallel for ordered
    for (int i = 0; i < 10; ++i)
    {
#pragma omp ordered
        {
            printf("Thread %d reporting with iteration %d!\n", omp_get_thread_num(), i);
        }
    }

    //    FILE *fp;
    //    fp = fopen("results1.txt", "w");
    srand(time(0));

    struct Fish *fishes = (struct Fish *)malloc(sizeof(struct Fish) * NUM_FISH);

// Init fish position and weight
#pragma omp parallel for
    for (int i = 0; i < NUM_FISH; i++)
    {
        unsigned int seed = omp_get_thread_num(); // Get unique seed for each thread
        fishes[i].x = (rand_r(&seed) % (2 * MAX_FISH_X + 1)) - MAX_FISH_X;
        fishes[i].y = (rand_r(&seed) % (2 * MAX_FISH_X + 1)) - MAX_FISH_Y;
        fishes[i].weight = INIT_WEIGHT;
    }

    double start = omp_get_wtime(); // Start the timer

    // Simulation
    for (int t = 0; t < SIMULATION_STEPS; t++)
    {
        float max_delta_f = -99999.0;
        unsigned int seed; // Declare seed outside the loop

// Get the fish planned position and calculate maximum δ(f)
#pragma omp parallel for private(seed) reduction(max : max_delta_f)
        for (int i = 0; i < NUM_FISH; i++)
        {
            seed = 123456789 + omp_get_thread_num(); // Initialize seed for each thread
            fishes[i].plan_x = fishes[i].x + ((float)rand_r(&seed) / RAND_MAX) * 0.2 - 0.1;
            fishes[i].plan_y = fishes[i].y + ((float)rand_r(&seed) / RAND_MAX) * 0.2 - 0.1;
            if (fishes[i].plan_x > MAX_FISH_X)
                fishes[i].plan_x = MAX_FISH_X;
            if (fishes[i].plan_x < MIN_FISH_X)
                fishes[i].plan_x = MIN_FISH_X;
            if (fishes[i].plan_y > MAX_FISH_Y)
                fishes[i].plan_y = MAX_FISH_Y;
            if (fishes[i].plan_y < MIN_FISH_Y)
                fishes[i].plan_y = MIN_FISH_Y;

            fishes[i].delta_f = sqrt(fishes[i].plan_x * fishes[i].plan_x + fishes[i].plan_y * fishes[i].plan_y) - sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);
            if (fishes[i].delta_f > max_delta_f)
            {
                max_delta_f = fishes[i].delta_f;
            }
        }

#pragma omp barrier // Wait for all threads to finish

        float total_distance_weighted = 0.0;
        float total_distance = 0.0;

// Update the fish position weight if it can eat and swim
#pragma omp parallel for reduction(+ : total_distance_weighted, total_distance)
        for (int i = 0; i < NUM_FISH; i++)
        {
            float planWeight = fishes[i].weight + fishes[i].delta_f / max_delta_f;
            if (planWeight > 2 * INIT_WEIGHT)
            {
                fishes[i].weight = 2 * INIT_WEIGHT;
            }
            if (planWeight > fishes[i].weight)
            {
                fishes[i].x = fishes[i].plan_x;
                fishes[i].y = fishes[i].plan_y;
                fishes[i].weight = planWeight;
            }
            float curr_distance = sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);
            total_distance_weighted += curr_distance * fishes[i].weight;
            total_distance += curr_distance;
        }

#pragma omp barrier // Wait for all threads to finish

        float barycentre = total_distance_weighted / total_distance;
        printf("Bari %d: %.5f\n", t, barycentre);

        double end = omp_get_wtime();

        double cpu_time_used = ((double)(end - start));
        //	fprintf(fp, "%d, %f\n", t, cpu_time_used); // Changed from SIMULATION_STEPS to t
        printf("Simulation step. Time taken: %f seconds.\n", cpu_time_used);
    }

    free(fishes); // Free the allocated memory
    return 0;
}
