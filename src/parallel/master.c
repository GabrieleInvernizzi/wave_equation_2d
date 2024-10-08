#include "master_worker.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "log.h"

int master(SimConf c, int n_procs_world) {
    FILE *f = NULL;

    // Check for CFL condition
    if ((c.c * (c.dt / c.dx) > 0.5) && (!c.ignore_cfl)) {
        return 1;
    }

    int n_workers = n_procs_world - 1;
    int dim = (int)sqrt(n_workers);
    MPI_Status last_status;

    LOGF("M[]: started.");
    // Gather the coords
    int *recv_coords_with_master = malloc(3 * n_procs_world * sizeof(int));
    assert(recv_coords_with_master);
    {
        int dummy_coords[3] = {MASTER_RANK, -1, -1};
        MPI_Gather(&dummy_coords, 3, MPI_INT, recv_coords_with_master, 3,
                   MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
        int j = 0;
        for (int i = 0; i < n_procs_world; i++) {
            LOGF("M[]: W[%d, world: %d] coords are: [%d, %d].", i,
                   recv_coords_with_master[j], recv_coords_with_master[j + 1],
                   recv_coords_with_master[j + 2]);
            j += 3;
        }
    }
    int *recv_coords = recv_coords_with_master + 3;

    const size_t w_n_cells = c.n_cells / dim;
    const size_t w_tot_n_cells = w_n_cells + 2;

    const size_t recv_buf_size = w_tot_n_cells * sizeof(double[w_tot_n_cells]);
    double(*recv_buf)[w_tot_n_cells] = malloc(recv_buf_size);
    assert(recv_buf);

    double(*frame)[c.tot_n_cells] =
        malloc(c.tot_n_cells * sizeof(double[c.tot_n_cells]));
    assert(frame);

    f = fopen(c.filepath, "wb");
    if (f == NULL) {
        LOGF("Could not open the file \"%s\".", c.filepath);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    fprintf(f, "%zu-%zu-%zu-%zu-%u\n", sizeof(double), c.tot_n_cells, c.tot_n_cells,
            (size_t)(c.n_steps / c.save_period), c.framerate);

    for (size_t s = 0; s < c.n_steps; s += c.save_period) {
        LOGF("M[]: step: %zu / %zu.", s, c.n_steps);
        for (int w = 0; w < n_workers; w++) {
            int *w_info = recv_coords + w * 3;
            // Recv from worker
            MPI_Recv(recv_buf, recv_buf_size, MPI_DOUBLE, w_info[0],
                     SEND_MASTER_TAG, MPI_COMM_WORLD, &last_status);
            // Position it correctly
            size_t i_base = w_info[1] == 0 ? 0 : (w_info[1] * w_tot_n_cells - 3);
            size_t j_base = w_info[2] == 0 ? 0 : (w_info[2] * w_tot_n_cells - 3);
            for (size_t i = 0; i < w_tot_n_cells; i++) {
                for (size_t j = 0; j < w_tot_n_cells; j++) {
                    frame[i_base + i][j_base + j] = recv_buf[i][j];
                }
            }
        }
        // Save frame to file
        fwrite((void *)frame, sizeof(double), c.tot_n_cells * c.tot_n_cells, f);
    }

    LOGF("M[]: step: %zu / %zu.", c.n_steps, c.n_steps);

    fclose(f);

    free(recv_coords_with_master);
    free(recv_buf);
    free(frame);

    LOGF("M[]: finished.");

    return 0;
}
