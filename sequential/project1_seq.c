#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

int main() {
    FILE *fp;
    fp = fopen("results1.txt", "w");
    srand(time(0));
    
    struct Fish *fishes = (struct Fish *)malloc(sizeof(struct Fish) * NUM_FISH);

    // Init fish position and weight
    for(int i = 0; i < NUM_FISH; i++) {
        fishes[i].x = (rand() % (2 * MAX_FISH_X + 1)) - MAX_FISH_X;
        fishes[i].y = (rand() % (2 * MAX_FISH_X + 1)) - MAX_FISH_Y;
        fishes[i].weight = INIT_WEIGHT;
    }

    clock_t start = clock();

    // Simulation
    for(int t = 0; t < SIMULATION_STEPS; t++) {
        // Get the fish planned position and calculate maximum δ(f)
        float max_delta_f = -99999.0;
        for(int i = 0; i < NUM_FISH; i++) {
            // Get a radom planned position
            fishes[i].plan_x = fishes[i].x + ((float) rand() / RAND_MAX) * 0.2 - 0.1;
            fishes[i].plan_y = fishes[i].y + ((float) rand() / RAND_MAX) * 0.2 - 0.1;
            if (fishes[i].plan_x > MAX_FISH_X) fishes[i].plan_x = MAX_FISH_X;
            if (fishes[i].plan_x < MIN_FISH_X) fishes[i].plan_x = MIN_FISH_X;
            if (fishes[i].plan_y > MAX_FISH_Y) fishes[i].plan_y = MAX_FISH_Y;
            if (fishes[i].plan_y < MIN_FISH_Y) fishes[i].plan_y = MIN_FISH_Y;
            
            // Calculate the change of the objective function of if the i-fish swims
            fishes[i].delta_f = sqrt(fishes[i].plan_x * fishes[i].plan_x + fishes[i].plan_y * fishes[i].plan_y) - sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);
            if(fishes[i].delta_f > max_delta_f) {
                max_delta_f = fishes[i].delta_f;
            }
        }

        // Update the fish position weight if it can eat and swim
        float total_distance_weighted = 0.0;
        float total_distance = 0.0;
        for(int i = 0; i < NUM_FISH; i++) {
            float planWeight = fishes[i].weight + fishes[i].delta_f / max_delta_f;
            if(planWeight > 2 * INIT_WEIGHT) {
                fishes[i].weight = 2 * INIT_WEIGHT;
            }

            // if the weight increases, it can eat and swim, so update the position and weight
            if (planWeight > fishes[i].weight) {
                fishes[i].x = fishes[i].plan_x;
                fishes[i].y = fishes[i].plan_y;
                fishes[i].weight = planWeight;
            }

            float curr_distance = sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);
            total_distance_weighted += curr_distance * fishes[i].weight;
            total_distance += curr_distance;
        }

        // Calculate the barycentre of the fish school
        float barycentre = total_distance_weighted / total_distance;

        clock_t end = clock();

        printf("Bari: %.5f\n", barycentre);

        double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        fprintf(fp, "%d, %f\n", SIMULATION_STEPS, cpu_time_used);
        printf("Simulation complete. Time taken: %f seconds.\n", cpu_time_used);
    }


    free(fishes);

    return 0;
}
