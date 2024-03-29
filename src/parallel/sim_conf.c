#include "sim_conf.h"

#include <argp.h>
#include <stdio.h>
#include <stdlib.h>

const char *argp_program_version = "wave_eq_p 1.0.0";
error_t argp_err_exit_status = PARSING_EXIT_FAILURE;

static const char doc[] =
    "wave_eq_p is a parallel 2d wave equation solver.\n"
    "The number of processors must be a perfect square + 1.\n"
    "The simulation is saved in the output file as the list of raw frames (the "
    "first line is the metadata). To create a video of the simulation use the "
    "script sim2video.";
static const char args_doc[] = "";

enum ArgpOptionCodes {
    // Redef of argp default flags
    ARGP_HELP = '?',
    ARGP_VERSION = 'V',
    ARGP_USAGE = 0x90,
    // Short options
    ARGP_FILEPATH = 'f',
    ARGP_SIZE = 's',
    ARGP_VEL = 'c',
    ARGP_DURATION = 'd',
    // Must be non-ascii
    ARGP_DX = 0x80,
    ARGP_DT = 0x81,
    ARGP_SAVE_PERIOD = 0x82,
    ARGP_IGNORE_CFL = 0x83
};

static struct argp_option options[] = {
    {"help", ARGP_HELP, 0, 0, "Show this message.", -1},
    {"version", ARGP_VERSION, 0, 0, "Show the version.", -1},
    {"usage", ARGP_USAGE, 0, 0, "Show usage.", 0},

    {"filepath", ARGP_FILEPATH, "PATH", 0,
     "Filepath that will be used as simulation output (default: \"out.sim\")."},
    {"size", ARGP_SIZE, "SIZE", 0,
     "Size of the domain, since it is a square domain this is the length of "
     "the sides (default: \"5.0\")."},
    {0, ARGP_VEL, "SPEED", 0, "Speed of the medium (default: \"1.0\")."},
    {"dx", ARGP_DX, "DX", 0,
     "dx value (dy will be the same, default: \"0.01\")."},
    {"dt", ARGP_DT, "DT", 0, "dt value (default: \"0.005\")."},
    {"save-period", ARGP_SAVE_PERIOD, "PERIOD", 0,
     "Save to file one every PERIOD frames (default: \"4\")."},
    {"duration", ARGP_DURATION, "DURATION", 0,
     "Duration of the simulation in seconds (default: \"10.0\")."},
    {"ignore-cfl", ARGP_IGNORE_CFL, 0, 0,
     "Ignore CFL condition (waring: could lead to unwanted results)."},
    {0}};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    SimConf *conf = state->input;

    switch (key) {
    case ARGP_HELP:
        argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
        return PARSING_EXIT_SUCCESS;
    case ARGP_USAGE:
        argp_state_help(state, state->out_stream,
                        ARGP_HELP_STD_USAGE | ARGP_HELP_EXIT_OK);
        return PARSING_EXIT_SUCCESS;
    case ARGP_VERSION: {
        if ((state->flags & ARGP_SILENT) != ARGP_SILENT)
            fprintf(state->out_stream, "%s\n", argp_program_version);
        return PARSING_EXIT_SUCCESS;
    }

    case ARGP_FILEPATH:
        conf->filepath = arg;
        break;
    case ARGP_SAVE_PERIOD: {
        int val = atoi(arg);
        if (val < 1) {
            argp_error(state, "--save-period accepts ints grater than 0.");
            return PARSING_EXIT_FAILURE;
        }
        conf->save_period = val;
    }

    case ARGP_SIZE: {
        double val = atof(arg);
        if (val <= 0.0) {
            argp_error(state, "--size accepts floats grater than 0.");
            return PARSING_EXIT_FAILURE;
        }
        conf->domain_size = val;
        break;
    }

    case ARGP_VEL: {
        double val = atof(arg);
        if (val <= 0.0) {
            argp_error(state, "-c accepts floats grater than 0.");
            return PARSING_EXIT_FAILURE;
        }
        conf->c = val;
        break;
    }

    case ARGP_DURATION: {
        double val = atof(arg);
        if (val <= 0.0) {
            argp_error(state, "--duration accepts floats grater than 0.");
            return PARSING_EXIT_FAILURE;
        }
        conf->duration = val;
        break;
    }

    case ARGP_DX: {
        double val = atof(arg);
        if (val <= 0.0) {
            argp_error(state, "--dx accepts floats grater than 0.");
            return PARSING_EXIT_FAILURE;
        }
        conf->dx = val;
        break;
    }
    case ARGP_DT: {
        double val = atof(arg);
        if (val <= 0.0) {
            argp_error(state, "--dt accepts floats grater than 0.");
            return PARSING_EXIT_FAILURE;
        }
        conf->dt = val;
        break;
    }

    case ARGP_IGNORE_CFL:
        conf->ignore_cfl = 1;
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

SimConf get_sim_conf(int argc, char **argv, int is_master, int *err) {
    SimConf c = {.filepath = "out.sim",
                 .save_period = 4,
                 .duration = 10,
                 .dt = 0.005,
                 .c = 1.0,
                 .domain_size = 5,
                 .dx = 0.01,
                 .ignore_cfl = 0};

    (*err) = argp_parse(
        &argp, argc, argv,
        (is_master ? 0 : ARGP_SILENT) | ARGP_NO_HELP | ARGP_NO_EXIT, 0, &c);
    if ((*err) != 0)
        return c;

    c.n_steps = (size_t)(c.duration / c.dt);
    c.framerate = 1.0 / (c.dt * c.save_period);
    c.n_cells = (size_t)(c.domain_size / c.dx);
    c.tot_n_cells = c.n_cells + 2;

    return c;
}
