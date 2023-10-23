#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#define NUM_FISH 4194304
#define MAX_FISH_X 100
#define MAX_FISH_Y 100
#define MIN_FISH_X -MAX_FISH_X
#define MIN_FISH_Y -MAX_FISH_Y
#define INIT_WEIGHT 10
#define SIMULATION_STEPS 500
#define MAX_THREAD 1
#define FISH_PER_THREAD 10000

struct Fish
{
    float x;
    float y;
    float plan_x;
    float plan_y;
    float weight;
    float delta_f;
};

static int fish_per_process;  // Define the number of fish in each process as a static variable

float phaseOneBase(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel for private(seed) reduction(max : max_delta_f)
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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

float phaseOneBaseOPT(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel private(seed) reduction(max : max_delta_f)
    {
        seed = 123456789 + omp_get_thread_num(); // Initialize seed for each thread
        struct Fish *fishesCurrThread = (struct Fish *)malloc(sizeof(struct Fish) * FISH_PER_THREAD);
        memcpy(fishesCurrThread, fishes + omp_get_thread_num() * FISH_PER_THREAD, sizeof(struct Fish) * FISH_PER_THREAD);

#pragma omp for
        for (int i = 0; i < MAX_THREAD; i++)
        {
            for (int j = 0; j < FISH_PER_THREAD; j++)
            {
                fishesCurrThread[j].plan_x = fishesCurrThread[j].x + ((float)rand_r(&seed) / RAND_MAX) * 0.2 - 0.1;
                fishesCurrThread[j].plan_y = fishesCurrThread[j].y + ((float)rand_r(&seed) / RAND_MAX) * 0.2 - 0.1;
                if (fishesCurrThread[j].plan_x > MAX_FISH_X)
                    fishesCurrThread[j].plan_x = MAX_FISH_X;
                if (fishesCurrThread[j].plan_x < MIN_FISH_X)
                    fishesCurrThread[j].plan_x = MIN_FISH_X;
                if (fishesCurrThread[j].plan_y > MAX_FISH_Y)
                    fishesCurrThread[j].plan_y = MAX_FISH_Y;
                if (fishesCurrThread[j].plan_y < MIN_FISH_Y)
                    fishesCurrThread[j].plan_y = MIN_FISH_Y;

                fishesCurrThread[j].delta_f = sqrt(fishesCurrThread[j].plan_x * fishesCurrThread[j].plan_x + fishesCurrThread[j].plan_y * fishesCurrThread[j].plan_y) - sqrt(fishesCurrThread[j].x * fishesCurrThread[j].x + fishesCurrThread[j].y * fishesCurrThread[j].y);
                if (fishesCurrThread[j].delta_f > max_delta_f)
                {
                    max_delta_f = fishesCurrThread[j].delta_f;
                }
            }
        }
        free(fishesCurrThread);
    }
    return max_delta_f;
}

float phaseTwoBaseOPT(struct Fish *fishes, float max_delta_f)
{
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

#pragma omp parallel reduction(+ : total_distance_weighted, total_distance)
    {
        struct Fish *fishesCurrThread = (struct Fish *)malloc(sizeof(struct Fish) * FISH_PER_THREAD);
        memcpy(fishesCurrThread, fishes + omp_get_thread_num() * FISH_PER_THREAD, sizeof(struct Fish) * FISH_PER_THREAD);

#pragma omp for
        for (int i = 0; i < MAX_THREAD; i++)
        {
            for (int j = 0; j < FISH_PER_THREAD; j++)
            {
                float planWeight = fishesCurrThread[j].weight + fishesCurrThread[j].delta_f / max_delta_f;
                if (planWeight > 2 * INIT_WEIGHT)
                {
                    fishesCurrThread[j].weight = 2 * INIT_WEIGHT;
                }
                if (planWeight > fishesCurrThread[j].weight)
                {
                    fishesCurrThread[j].x = fishesCurrThread[j].plan_x;
                    fishesCurrThread[j].y = fishesCurrThread[j].plan_y;
                    fishesCurrThread[j].weight = planWeight;
                }
                float curr_distance = sqrt(fishesCurrThread[j].x * fishesCurrThread[j].x + fishesCurrThread[j].y * fishesCurrThread[j].y);
                total_distance_weighted += curr_distance * fishesCurrThread[j].weight;
                total_distance += curr_distance;
            }
        }
        free(fishesCurrThread);
    }
    return total_distance_weighted / total_distance;
}

float phaseOneStaticOne(struct Fish *fishes)
{
    float max_delta_f = -99999.0;
    unsigned int seed;

#pragma omp parallel for schedule(static, 1) private(seed) reduction(max : max_delta_f)
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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
    for (int i = 0; i < fish_per_process; i++)
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

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
        printf("MPI does not provide the required thread support!\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    omp_set_num_threads(MAX_THREAD); // Set the number of threads

    int process_id, number_of_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    fish_per_process = NUM_FISH / number_of_processes;

    srand(time(0) + process_id);

    struct Fish *fishes = NULL;
    if (process_id == 0) {
        fishes = (struct Fish *)malloc(sizeof(struct Fish) * NUM_FISH);
        initFish(fishes);
    }

    struct Fish *local_fishes = (struct Fish *)malloc(sizeof(struct Fish) * fish_per_process);

    MPI_Scatter(fishes, fish_per_process, MPI_BYTE, local_fishes, fish_per_process, MPI_BYTE, 0, MPI_COMM_WORLD);

    double start_time, end_time;
    if (process_id == 0) {
        start_time = MPI_Wtime();
    }

    float local_max_delta_f, global_max_delta_f;
    switch (option)
    {
    case 1:
        // Simulation for base configuration
        if (process_id == 0) printf("Base configuration simulation.\n");
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            local_max_delta_f = phaseOneBase(local_fishes);
            MPI_Allreduce(&local_max_delta_f, &global_max_delta_f, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            phaseTwoBase(local_fishes, global_max_delta_f);
        }
        break;
    case 2:
        // Simulation for static one configuration
        if (process_id == 0) printf("Static one configuration simulation.\n");
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            local_max_delta_f = phaseOneStaticOne(local_fishes);
            MPI_Allreduce(&local_max_delta_f, &global_max_delta_f, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            phaseTwoStaticOne(local_fishes, global_max_delta_f);
        }
        break;
    case 3:
        // Simulation for static FTP configuration
        if (process_id == 0) printf("Static FTP configuration simulation.\n");
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            local_max_delta_f = phaseOneStaticFPT(local_fishes);
            MPI_Allreduce(&local_max_delta_f, &global_max_delta_f, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            phaseTwoStaticFTP(local_fishes, global_max_delta_f);
        }
        break;
    case 4:
        // Simulation for dynamic FTP configuration
        if (process_id == 0) printf("Dynamic FTP configuration.\n");
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            local_max_delta_f = phaseOneDynamicFPT(local_fishes);
            MPI_Allreduce(&local_max_delta_f, &global_max_delta_f, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            phaseTwoDynamicFTP(local_fishes, global_max_delta_f);
        }
        break;
    case 5:
        // Simulation for guided configuration
        if (process_id == 0) printf("Guided configuration simulation.\n");
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            local_max_delta_f = phaseOneGuided(local_fishes);
            MPI_Allreduce(&local_max_delta_f, &global_max_delta_f, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            phaseTwoGuided(local_fishes, global_max_delta_f);
        }     
        break;
    case 6:
        // Simulation for base configuration
        if (process_id == 0) printf("Base OPT configuration simulation.\n");
        for (int t = 0; t < SIMULATION_STEPS; t++)
        {
            local_max_delta_f = phaseOneBaseOPT(local_fishes);
            MPI_Allreduce(&local_max_delta_f, &global_max_delta_f, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            phaseTwoBaseOPT(local_fishes, global_max_delta_f);
        }
        break;
    default:
        break;
    }

    if (process_id == 0) {
        end_time = MPI_Wtime();
        printf("Time taken: %f seconds.\n", end_time - start_time);
    }

    free(local_fishes);
    if (process_id == 0) {
        free(fishes);
    }

    MPI_Finalize();
    return 0;
}
