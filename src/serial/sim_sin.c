#include "sims.h"

#include <math.h>
#include <stdlib.h>

void sim_sin(SimConf c, FILE *f) {
    // Init arrays
    double(*u_tmp)[c.tot_cols] = NULL;
    double(*u0)[c.tot_cols] =
        malloc(c.tot_rows * sizeof(double[c.tot_cols])); // u(k)
    double(*u1)[c.tot_cols] =
        malloc(c.tot_rows * sizeof(double[c.tot_cols])); // u(k-1)
    double(*u2)[c.tot_cols] =
        malloc(c.tot_rows * sizeof(double[c.tot_cols])); // u(k-2)

    // Courant numbers
    double Cx = c.c * (c.dt / c.dx);
    double Cx_sq = Cx * Cx;
    double Cy = c.c * (c.dt / c.dy);
    double Cy_sq = Cy * Cy;

    // u2 initial conditions
    for (size_t j = 0; j < c.tot_cols; j++) {
        double y = j * c.dy;
        for (size_t i = 0; i < c.tot_rows; i++) {
            double x = i * c.dx;
            u2[i][j] = sin(2 * M_PI * (x + y));
        }
    }

    // u1 initial conditions
    for (size_t j = 1; j < c.tot_cols - 1; j++) {
        for (size_t i = 1; i < c.tot_rows - 1; i++) {
            u1[i][j] =
                u2[i][j] +
                (0.25 * Cx_sq) * (u2[i - 1][j] - 2 * u2[i][j] + u2[i + 1][j]) +
                (0.25 * Cy_sq) * (u2[i][j - 1] - 2 * u2[i][j] + u2[i][j + 1]);
        }
    }

    for (size_t j = 0; j < c.tot_cols; j++) {
        double y = j * c.dy;
        double x = (c.tot_rows - 1) * c.dx;
        u1[0][j] = sin(2 * M_PI * (y + c.dt));
        u1[c.tot_rows - 1][j] = sin(2 * M_PI * (x + y + c.dt));
    }

    for (size_t i = 0; i < c.tot_rows; i++) {
        double x = i * c.dx;
        double y = (c.tot_cols - 1) * c.dy;
        u1[i][0] = sin(2 * M_PI * (x + c.dt));
        u1[i][c.tot_cols - 1] = sin(2 * M_PI * (x + y + c.dt));
    }

    double t = c.dt;
    for (size_t step = 0; step < c.n_steps; step++) {
        t += c.dt;

        for (size_t j = 1; j < c.tot_cols - 1; j++) {
            for (size_t i = 1; i < c.tot_rows - 1; i++) {
                u0[i][j] = 2 * u1[i][j] - u2[i][j] +
                           (0.5 * Cx_sq) *
                               (u1[i - 1][j] - 2 * u1[i][j] + u1[i + 1][j]) +
                           (0.5 * Cy_sq) *
                               (u1[i][j - 1] - 2 * u1[i][j] + u1[i][j + 1]);
            }
        }

        // u0 boundary conditions
        for (size_t j = 0; j < c.tot_cols; j++) {
            double y = j * c.dy;
            double x = (c.tot_rows - 1) * c.dx;
            u0[0][j] = sin(2 * M_PI * (y + t));
            u0[c.tot_rows - 1][j] = sin(2 * M_PI * (x + y + t));
        }

        for (size_t i = 0; i < c.tot_rows; i++) {
            double x = i * c.dx;
            double y = (c.tot_cols - 1) * c.dy;
            u0[i][0] = sin(2 * M_PI * (x + t));
            u0[i][c.tot_cols - 1] = sin(2 * M_PI * (x + y + t));
        }

        // move the queue
        u_tmp = u2;
        u2 = u1;
        u1 = u0;
        u0 = u_tmp;

        // Save the frame
        if (step % c.save_period == 0)
            fwrite((void *)u1, sizeof(double), c.tot_rows * c.tot_cols, f);
    }

    free(u2);
    free(u1);
    free(u0);
}
