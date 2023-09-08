#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NUM_FISH 1000
#define LAKE_SIZE 100
#define INIT_WEIGHT 10

struct Fish {
    float x;
    float y;
    float old_x; // For storing the old x-coordinate
    float old_y; // For storing the old y-coordinate
    float weight;
};

float calculate_distance(float x, float y) {
    return sqrt(x * x + y * y);
}

// Calculate the new position when swim
void swim(struct Fish* fish) {
    fish->old_x = fish->x;
    fish->old_y = fish->y;
    fish->x += ((float) rand() / RAND_MAX) * 0.1;
    fish->y += ((float) rand() / RAND_MAX) * 0.1;
}

int main() {
    FILE *fp;
    fp = fopen("results1.txt", "w");

    for(int SIMULATION_STEPS = 50; SIMULATION_STEPS <= 5000; SIMULATION_STEPS += 50) {
        struct Fish fishes[NUM_FISH];
        float barycentre;
        int t, i;

        srand(time(0));

        // Init fish position and weight
        for(i = 0; i < NUM_FISH; i++) {
            fishes[i].x = (float) rand() / RAND_MAX * LAKE_SIZE - (LAKE_SIZE / 2);
            fishes[i].y = (float) rand() / RAND_MAX * LAKE_SIZE - (LAKE_SIZE / 2);
            fishes[i].weight = INIT_WEIGHT;
        }

        clock_t start = clock();

        // Simulation
        for(t = 0; t < SIMULATION_STEPS; t++) {
            float total_distance_weighted = 0;
            float total_distance = 0;
            float max_delta_f = -1;

            // Update the fish position and calculate maximum δ(f)
            for(i = 0; i < NUM_FISH; i++) {
                float old_distance = calculate_distance(fishes[i].x, fishes[i].y);
                
                swim(&fishes[i]); // Assume all fish can eat food in every round
                
                float new_distance = calculate_distance(fishes[i].x, fishes[i].y);
                
                float delta_f = new_distance - old_distance;

                if(delta_f > max_delta_f) {
                    max_delta_f = delta_f;
                }
                
                total_distance_weighted += new_distance * fishes[i].weight;
                total_distance += new_distance;
            }

            // Update the fish weight
            for(i = 0; i < NUM_FISH; i++) {
                float old_distance = calculate_distance(fishes[i].old_x, fishes[i].old_y);
                float new_distance = calculate_distance(fishes[i].x, fishes[i].y);
                float delta_f = new_distance - old_distance;
                
                fishes[i].weight += delta_f / max_delta_f;
                if(fishes[i].weight > 2 * INIT_WEIGHT) {
                    fishes[i].weight = 2 * INIT_WEIGHT;
                }
            }

            // Calculate the barycentre of the fish school
            barycentre = total_distance_weighted / total_distance;
        }

        clock_t end = clock();

        double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        fprintf(fp, "%d, %f\n", SIMULATION_STEPS, cpu_time_used);
        // printf("Simulation complete. Time taken: %f seconds.\n", cpu_time_used);
    }

    return 0;
}
