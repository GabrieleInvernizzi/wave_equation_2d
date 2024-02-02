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
    double domain_width;
    double domain_height;
    double dx;
    double dy;
    size_t cols;
    size_t rows;
    size_t tot_cols;
    size_t tot_rows;
} SimConf;

SimConf get_sim_conf(int argc, char **argv, int is_master, int *err);
