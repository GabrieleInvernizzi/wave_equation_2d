#include "sim.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "arena_alloc.h"
#include "log.h"

void sim(SimConf c, FILE *f) {
    // Init arrays
    const size_t arr_size = c.tot_rows * sizeof(double[c.tot_cols]);
    const size_t arena_size = 3 * arr_size;
    Arena arena;
    arena_init(&arena, arena_size);

    double(*u_tmp)[c.tot_cols] = NULL;
    double(*u0)[c.tot_cols] = arena_alloc(&arena, arr_size); // u(k)
    assert(u0);
    double(*u1)[c.tot_cols] = arena_alloc(&arena, arr_size); // u(k-1)
    assert(u1);
    double(*u2)[c.tot_cols] = arena_alloc(&arena, arr_size); // u(k-2)
    assert(u2);

    // Courant numbers
    double Cx = c.c * (c.dt / c.dx);
    double Cx_sq = Cx * Cx;
    double Cy = c.c * (c.dt / c.dy);
    double Cy_sq = Cy * Cy;

    LOGF("Starting simulation.\nstep: 0 / %zu.", c.n_steps);

    // u2 initial conditions
    for (size_t j = 0; j < c.tot_cols; j++)
        for (size_t i = 0; i < c.tot_rows; i++)
            u2[i][j] = 0.0;

    // u1 initial conditions
    for (size_t j = 1; j < c.tot_cols; j++)
        for (size_t i = 1; i < c.tot_rows; i++)
            u1[i][j] = 0.0;

    double t = c.dt;
    for (size_t step = 0; step < c.n_steps; step++) {
        t += c.dt;

        for (size_t j = 1; j < c.tot_cols - 1; j++) {
            for (size_t i = 1; i < c.tot_rows - 1; i++) {
                u0[i][j] = 2 * u1[i][j] - u2[i][j] +
                           (0.5 * Cx_sq) *
                               (u1[i][j - 1] - 2 * u1[i][j] + u1[i][j + 1]) +
                           (0.5 * Cy_sq) *
                               (u1[i - 1][j] - 2 * u1[i][j] + u1[i + 1][j]);
            }
        }

        // Forcing term
        u0[(size_t)(c.tot_rows / 1.2)][(size_t)(c.tot_cols / 1.2)] +=
            c.dt * c.dt * 200 * sin(2 * M_PI * 2 * t);

        // Reflective boundary
        for (size_t j = 0; j < c.tot_cols; j++) {
            u0[0][j] = u0[2][j];
            u0[c.tot_rows - 1][j] = u0[c.tot_rows - 3][j];
        }

        for (size_t i = 0; i < c.tot_rows; i++) {
            u0[i][0] = u0[i][2];
            u0[i][c.tot_cols - 1] = u0[i][c.tot_cols - 3];
        }

        // move the queue
        u_tmp = u2;
        u2 = u1;
        u1 = u0;
        u0 = u_tmp;

        // Save the frame
        if (step % c.save_period == 0) {
            LOGF("step: %zu / %zu.", step, c.n_steps);
            fwrite(buffer, sizeof(double), buffIndex, f);
        }
    }

    LOGF("step: %zu / %zu.\nFinished.", c.n_steps, c.n_steps);

    area_deinit(&arena);
}
