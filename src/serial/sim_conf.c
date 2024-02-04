#include "sim_conf.h"

#include <argp.h>
#include <stdlib.h>

const char *argp_program_version = "wave_eq_s 1.0.0";
error_t argp_err_exit_status = PARSING_ERROR_EXIT_CODE;

static const char doc[] = "";
static const char args_doc[] = "";

enum ArgpOptionCodes {
    ARGP_FILEPATH = 'f',
    ARGP_SIZE = 's',
    ARGP_VEL = 'c',
    ARGP_DURATION = 'd',
    // Must be non-ascii
    ARGP_DX = 0x80,
    ARGP_DT = 0x81,
    ARGP_SAVE_PERIOD = 0x82
};

static struct argp_option options[] = {
    {"filepath", ARGP_FILEPATH, "PATH", 0,
     "Filepath that will be used as simulation output."},
    {"size", ARGP_SIZE, "SIZE", 0,
     "Size of the domain, since it is a square domain this is the length of "
     "the side."},
    {0, ARGP_VEL, "SPEED", 0, "Speed of the medium."},
    {"dx", ARGP_DX, "DX", 0, "dx value (dy will be the same)."},
    {"dt", ARGP_DT, "DT", 0, "dt value."},
    {"save-period", ARGP_SAVE_PERIOD, "PERIOD", 0,
     "Save to file one every PERIOD frames."},
    {"duration", ARGP_DURATION, "DURATION", 0,
     "Duration of the simulation in seconds."},
    {0}};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    SimConf *conf = state->input;

    switch (key) {
    case ARGP_FILEPATH:
        conf->filepath = arg;
        break;
    case ARGP_SAVE_PERIOD: {
        int val = atoi(arg);
        if (val < 1)
            return 1;
        conf->save_period = val;
    }

    case ARGP_SIZE: {

        double val = atof(arg);
        if (val == 0.0)
            return 1;
        conf->domain_size = val;
        break;
    }

    case ARGP_VEL: {
        double val = atof(arg);
        if (val == 0.0)
            return 1;
        conf->c = val;
        break;
    }

    case ARGP_DURATION: {
        double val = atof(arg);
        if (val == 0.0)
            return 1;
        conf->c = val;
        break;
    }

    case ARGP_DX: {
        double val = atof(arg);
        if (val == 0.0)
            return 1;
        conf->dx = val;
        break;
    }
    case ARGP_DT: {
        double val = atof(arg);
        if (val == 0.0)
            return 1;
        conf->dt = val;
        break;
    }
    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

SimConf get_sim_conf(int argc, char **argv) {
    SimConf c = {.filepath = "out.sim",
                 .save_period = 2,
                 .duration = 10,
                 .dt = 0.005,
                 .c = 1.0,
                 .domain_size = 5,
                 .dx = 0.01};

    argp_parse(&argp, argc, argv, 0, 0, &c);

    c.n_steps = (size_t)(c.duration / c.dt);
    c.framerate = 1.0 / (c.dt * c.save_period);
    c.n_cells = (size_t)(c.domain_size / c.dx);
    c.tot_n_cells = c.n_cells + 2;

    return c;
}
