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
#define FISH_PER_THREAD (NUM_FISH / MAX_THREAD)

struct Fish
{
    float x;
    float y;
    float plan_x;
    float plan_y;
    float weight;
    float delta_f;
};

float phaseOneBase(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

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
    return max_delta_f;
}

float phaseTwoBase(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

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

    return total_distance_weighted / total_distance;
}

float phaseOneStaticOne(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel for schedule(static, 1) private(seed) reduction(max : max_delta_f)
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
    return max_delta_f;
}

float phaseOneStaticFPT(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel for schedule(static, FISH_PER_THREAD) private(seed) reduction(max : max_delta_f)
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
    return max_delta_f;
}

float phaseOneDynamicOne(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel for schedule(dynamic, 1) private(seed) reduction(max : max_delta_f)
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
    return max_delta_f;
}

float phaseOneDynamicFPT(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel for schedule(dynamic, FISH_PER_THREAD) private(seed) reduction(max : max_delta_f)
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
    return max_delta_f;
}

float phaseOneGuided(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel for schedule(guided) private(seed) reduction(max : max_delta_f)
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
    return max_delta_f;
}

float phaseTwoStaticOne(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

#pragma omp parallel for schedule(static, 1) reduction(+ : total_distance_weighted, total_distance)
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

    return total_distance_weighted / total_distance;
}

float phaseTwoStaticFTP(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

#pragma omp parallel for schedule(static, FISH_PER_THREAD) reduction(+ : total_distance_weighted, total_distance)
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

    return total_distance_weighted / total_distance;
}

float phaseTwoDynamicOne(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+ : total_distance_weighted, total_distance)
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

    return total_distance_weighted / total_distance;
}

float phaseTwoDynamicFTP(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

#pragma omp parallel for schedule(dynamic, FISH_PER_THREAD) reduction(+ : total_distance_weighted, total_distance)
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

    return total_distance_weighted / total_distance;
}

float phaseTwoGuided(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

#pragma omp parallel for schedule(guided) reduction(+ : total_distance_weighted, total_distance)
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

    return total_distance_weighted / total_distance;
}

void initFish(struct Fish *fishes)
{
// Init fish position and weight
#pragma omp parallel for
    for (int i = 0; i < NUM_FISH; i++)
    {
        unsigned int seed = omp_get_thread_num(); // Get unique seed for each thread
        fishes[i].x = (rand_r(&seed) % (2 * MAX_FISH_X + 1)) - MAX_FISH_X;
        fishes[i].y = (rand_r(&seed) % (2 * MAX_FISH_X + 1)) - MAX_FISH_Y;
        fishes[i].weight = INIT_WEIGHT;
    }
}

int main(int argc, char *argv[])
{
    int option = atoi(argv[1]);
    omp_set_num_threads(MAX_THREAD); // Set the number of threads

    double start, end;

    srand(time(0));

    struct Fish *fishes = (struct Fish *)malloc(sizeof(struct Fish) * NUM_FISH);

    switch (option)
    {
    case 1:
        initFish(fishes);
        // Simulation for base configuration
        start = omp_get_wtime();
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            float barycentreBase = phaseTwoBase(fishes, phaseOneBase(fishes));
        }
        end = omp_get_wtime();
        printf("Base configuration simulation step. Time taken: %f seconds.\n", end - start);
        break;
    case 2:
        initFish(fishes);
        // Simulation for static one configuration
        start = omp_get_wtime();
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            float barycentreStaticOne = phaseTwoStaticOne(fishes, phaseOneStaticOne(fishes));
        }
        end = omp_get_wtime();
        printf("Static one configuration simulation step. Time taken: %f seconds.\n", end - start);
        break;
    case 3:
        initFish(fishes);
        // Simulation for static FTP configuration
        start = omp_get_wtime();
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            float barycentreStaticFPT = phaseTwoStaticFTP(fishes, phaseOneStaticFPT(fishes));
        }
        end = omp_get_wtime();
        printf("Static FTP configuration simulation step. Time taken: %f seconds.\n", end - start);
        break;
    case 4:
        initFish(fishes);
        // Simulation for dynamic FTP configuration
        start = omp_get_wtime();
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            float barycentreDynamicFTP = phaseTwoDynamicFTP(fishes, phaseOneDynamicFPT(fishes));
        }
        end = omp_get_wtime();
        printf("Dynamic FTP configuration simulation step. Time taken: %f seconds.\n", end - start);
        break;
    case 5:
        initFish(fishes);
        // Simulation for guided configuration
        start = omp_get_wtime();
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            float barycentreGuided = phaseTwoGuided(fishes, phaseOneGuided(fishes));
        }
        end = omp_get_wtime();
        printf("Guided configuration simulation step. Time taken: %f seconds.\n", end - start);
        break;

    default:
        break;
    }

    free(fishes); // Free the allocated memory
    return 0;
}
