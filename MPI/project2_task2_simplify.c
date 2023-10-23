#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define NUM_FISH 4194304
#define MAX_FISH_X 100
#define MAX_FISH_Y 100
#define MIN_FISH_X -MAX_FISH_X
#define MIN_FISH_Y -MAX_FISH_Y
#define INIT_WEIGHT 10
#define SIMULATION_STEPS 500

struct Fish {
    float x;
    float y;
    float plan_x;
    float plan_y;
    float weight;
    float delta_f;
};

static int fish_per_process;  // Define the number of fish in each process as a static variable

float phaseOneBase(struct Fish *fishes) {
    float max_delta_f = -99999.0;

    for (int i = 0; i < fish_per_process; i++) {
        fishes[i].plan_x = fishes[i].x + ((float)rand() / RAND_MAX) * 0.2 - 0.1;
        fishes[i].plan_y = fishes[i].y + ((float)rand() / RAND_MAX) * 0.2 - 0.1;
        if (fishes[i].plan_x > MAX_FISH_X) fishes[i].plan_x = MAX_FISH_X;
        if (fishes[i].plan_x < MIN_FISH_X) fishes[i].plan_x = MIN_FISH_X;
        if (fishes[i].plan_y > MAX_FISH_Y) fishes[i].plan_y = MAX_FISH_Y;
        if (fishes[i].plan_y < MIN_FISH_Y) fishes[i].plan_y = MIN_FISH_Y;

        fishes[i].delta_f = sqrt(fishes[i].plan_x * fishes[i].plan_x + fishes[i].plan_y * fishes[i].plan_y) - sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);
        if (fishes[i].delta_f > max_delta_f) {
            max_delta_f = fishes[i].delta_f;
        }
    }
    return max_delta_f;
}

float phaseTwoBase(struct Fish *fishes, float max_delta_f) {
    float total_distance_weighted = 0.0;
    float total_distance = 0.0;

    for (int i = 0; i < fish_per_process; i++) {
        float planWeight = fishes[i].weight + fishes[i].delta_f / max_delta_f;
        if (planWeight > 2 * INIT_WEIGHT) {
            fishes[i].weight = 2 * INIT_WEIGHT;
        }
        if (planWeight > fishes[i].weight) {
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

void initFish(struct Fish *fishes) {
    for (int i = 0; i < NUM_FISH; i++) {
        fishes[i].x = (rand() % (2 * MAX_FISH_X + 1)) - MAX_FISH_X;
        fishes[i].y = (rand() % (2 * MAX_FISH_X + 1)) - MAX_FISH_Y;
        fishes[i].weight = INIT_WEIGHT;
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
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

    for (int t = 0; t < SIMULATION_STEPS; t++) {
        local_max_delta_f = phaseOneBase(local_fishes);
        MPI_Allreduce(&local_max_delta_f, &global_max_delta_f, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
        phaseTwoBase(local_fishes, global_max_delta_f);
    }

    if (process_id == 0) {
        end_time = MPI_Wtime();
        printf("Simulation time: %f seconds\n", end_time - start_time);
    }

    free(local_fishes);
    if (process_id == 0) {
        free(fishes);
    }

    MPI_Finalize();
    return 0;
}
