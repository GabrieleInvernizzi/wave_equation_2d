#pragma once

#include <stddef.h>

#define PARSING_EXIT_SUCCESS (2)
#define PARSING_EXIT_FAILURE (1)

typedef struct {
    // File
    char *filepath;
    size_t save_period;
    unsigned int framerate;

    // Simulation
    double duration;
    double dt;
    size_t n_steps;

    double c;
    double domain_size;
    double dx;
    size_t n_cells;
    size_t tot_n_cells;
} SimConf;

SimConf get_sim_conf(int argc, char **argv, int is_master, int *err);
