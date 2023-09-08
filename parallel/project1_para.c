#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h> // include OpenMP header

#define NUM_FISH 100000
#define LAKE_SIZE 1000
#define INIT_WEIGHT 10
#define NUM_THREADS 4

struct Fish {
    float x;
    float y;
    float old_x;
    float old_y;
    float weight;
};

float calculate_distance(float x, float y) {
    return sqrt(x * x + y * y);
}

void swim(struct Fish* fish) {
    fish->old_x = fish->x;
    fish->old_y = fish->y;
    fish->x += ((float) rand() / RAND_MAX) * 0.1;
    fish->y += ((float) rand() / RAND_MAX) * 0.1;
}

int main() {
    FILE *fp;
    fp = fopen("results1.txt", "w");

    // Setting the number of threads
    omp_set_num_threads(NUM_THREADS);

    for(int SIMULATION_STEPS = 10; SIMULATION_STEPS <= 100; SIMULATION_STEPS += 10) {
        struct Fish fishes[NUM_FISH];
        float barycentre;
        int t, i;

        srand(time(0));

        for(i = 0; i < NUM_FISH; i++) {
            fishes[i].x = (float) rand() / RAND_MAX * LAKE_SIZE - (LAKE_SIZE / 2);
            fishes[i].y = (float) rand() / RAND_MAX * LAKE_SIZE - (LAKE_SIZE / 2);
            fishes[i].weight = INIT_WEIGHT;
        }

        double start = omp_get_wtime();

        for(t = 0; t < SIMULATION_STEPS; t++) {
            float total_distance_weighted = 0;
            float total_distance = 0;
            float max_delta_f = -1;

            // Parallel loop for updating fish position and calculating maximum δ(f)
            #pragma omp parallel for schedule(dynamic) reduction(max:max_delta_f) reduction(+:total_distance_weighted,total_distance)
            for(i = 0; i < NUM_FISH; i++) {
                float old_distance = calculate_distance(fishes[i].x, fishes[i].y);
                swim(&fishes[i]);
                float new_distance = calculate_distance(fishes[i].x, fishes[i].y);
                float delta_f = new_distance - old_distance;

                if(delta_f > max_delta_f) {
                    max_delta_f = delta_f;
                }
                
                total_distance_weighted += new_distance * fishes[i].weight;
                total_distance += new_distance;
            }

            // Parallel loop for updating fish weight
            #pragma omp parallel for
            for(i = 0; i < NUM_FISH; i++) {
                float old_distance = calculate_distance(fishes[i].old_x, fishes[i].old_y);
                float new_distance = calculate_distance(fishes[i].x, fishes[i].y);
                float delta_f = new_distance - old_distance;
                
                fishes[i].weight += delta_f / max_delta_f;
                if(fishes[i].weight > 2 * INIT_WEIGHT) {
                    fishes[i].weight = 2 * INIT_WEIGHT;
                }
            }

            barycentre = total_distance_weighted / total_distance;
        }

        double end = omp_get_wtime();

        fprintf(fp, "%d, %f\n", SIMULATION_STEPS, end - start);
    }

    fclose(fp);
    return 0;
}
