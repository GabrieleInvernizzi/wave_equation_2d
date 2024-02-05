#include "sim.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "arena_alloc.h"
#include "log.h"
#include "timings.h"

int sim(SimConf c, FILE *f) {
    // Init arrays
    const size_t arr_size = c.tot_n_cells * sizeof(double[c.tot_n_cells]);
    const size_t arena_size = 3 * arr_size;
    Arena arena;
    arena_init(&arena, arena_size);

    double(*u_tmp)[c.tot_n_cells] = NULL;
    double(*u0)[c.tot_n_cells] = arena_alloc(&arena, arr_size); // u(k)
    assert(u0);
    double(*u1)[c.tot_n_cells] = arena_alloc(&arena, arr_size); // u(k-1)
    assert(u1);
    double(*u2)[c.tot_n_cells] = arena_alloc(&arena, arr_size); // u(k-2)
    assert(u2);

    // Courant numbers
    double C = c.c * (c.dt / c.dx);
    double C_sq = C * C;

    // Check CFL condition
    if (C > 0.5 && (!c.ignore_cfl)) {
        LOGF("The CFL condition is not satisfied.\nC = c * (dt / dx) = %f > "
             "0.5.\nChoose the parameters so that C <= 0.5.\nExiting.",
             C);
        return 1;
    }

    LOGF("Starting simulation.\nstep: 0 / %zu.", c.n_steps);

    // set initial conditions
    for (size_t i = 0; i < c.tot_n_cells; i++) {
        for (size_t j = 0; j < c.tot_n_cells; j++) {
            u2[i][j] = 0.0;
            u1[i][j] = 0.0;
            u0[i][j] = 0.0;
        }
    }

    double t = c.dt;
    for (size_t step = 0; step < c.n_steps; step++) {
        t += c.dt;

        START_TIMER("c");
        for (size_t i = 1; i < c.tot_n_cells - 1; i++) {
            for (size_t j = 1; j < c.tot_n_cells - 1; j++) {
                u0[i][j] =
                    2 * u1[i][j] - u2[i][j] +
                    (0.5 * C_sq) *
                        (u1[i][j - 1] - 2 * u1[i][j] + u1[i][j + 1]) +
                    (0.5 * C_sq) * (u1[i - 1][j] - 2 * u1[i][j] + u1[i + 1][j]);
            }
        }

        // Forcing term
        u0[(size_t)(c.tot_n_cells / 2)][(size_t)(c.tot_n_cells / 2)] +=
            c.dt * c.dt * 200 * sin(2 * M_PI * 2 * t);

        // Reflective boundary
        for (size_t j = 0; j < c.tot_n_cells; j++) {
            u0[0][j] = u0[2][j];
            u0[c.tot_n_cells - 1][j] = u0[c.tot_n_cells - 3][j];
        }

        for (size_t i = 0; i < c.tot_n_cells; i++) {
            u0[i][0] = u0[i][2];
            u0[i][c.tot_n_cells - 1] = u0[i][c.tot_n_cells - 3];
        }

        // move the queue
        u_tmp = u2;
        u2 = u1;
        u1 = u0;
        u0 = u_tmp;

        END_TIMER;

        // Save the frame
        if (step % c.save_period == 0) {
            LOGF("step: %zu / %zu.", step, c.n_steps);
            fwrite((void *)u1, sizeof(double), c.tot_n_cells * c.tot_n_cells,
                   f);
        }
    }

    LOGF("step: %zu / %zu.\nFinished.", c.n_steps, c.n_steps);

    area_deinit(&arena);

    return 0;
}
