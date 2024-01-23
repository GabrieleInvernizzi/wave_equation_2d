#include "master_worker.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "log.h"
#include "sim_conf.h"

int worker(MPI_Comm comm, int my_rank_world, int master_rank) {
    int my_rank;
    int n_procs;

    SimConf c = get_sim_conf();

    MPI_Comm_size(comm, &n_procs);
    MPI_Comm_rank(comm, &my_rank);

    // Check if the number of worker procs is a perfect square
    {
        int dim = (int)sqrt(n_procs);
        if (n_procs != (dim * dim)) {
            if (my_rank == 0) {
                LOGF("W[0]: worker processes must be a perfect square but "
                     "are %d.",
                     n_procs);
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 1;
            }
        }
    }

    LOGF("W[%d]: started.", my_rank);

    // Create cartesian topology
    int dims[N_DIMS] = {0};
    MPI_Dims_create(n_procs, N_DIMS, dims);
    if (my_rank == 0)
        LOGF("W[0]: worker grid = [%d, %d].", dims[0], dims[1]);

    // Create cartesian mapping
    MPI_Comm comm2d;
    int wrap_around[N_DIMS] = {0};
    int reorder = 1;

    int err =
        MPI_Cart_create(comm, N_DIMS, dims, wrap_around, reorder, &comm2d);
    if (err != 0)
        LOGF("W[%d]: Error creating the cart (%d).", my_rank, err);

    int coord[N_DIMS];
    MPI_Cart_coords(comm2d, my_rank, N_DIMS, coord);
    int my_cart_rank;
    MPI_Cart_rank(comm2d, coord, &my_cart_rank);

    LOGF("W[%d]: my_cart_rank CT[%d], my coords = (%d,%d).", my_rank,
         my_cart_rank, coord[0], coord[1]);

    // Send to master my coords
    {
        int send_gather[3] = {my_rank_world, coord[0], coord[1]};
        MPI_Gather(send_gather, 3, MPI_INT, NULL, 0, MPI_INT, master_rank,
                   MPI_COMM_WORLD);
    }

    // neigh = {top, left, down, right }
    int neigh[4] = {[0 ... 3] = MPI_PROC_NULL};
    // Get index of neighbors
    MPI_Cart_shift(comm2d, 0, 1, &neigh[0], &neigh[2]);
    MPI_Cart_shift(comm2d, 1, 1, &neigh[1], &neigh[3]);
    LOGF("W[%d, CT: %d] - neighs: [%d, %d, %d, %d].", my_rank, my_cart_rank,
         neigh[0], neigh[1], neigh[2], neigh[3]);

    // Request for async
    MPI_Request gh_cells_req[4] = {[0 ... 3] = MPI_REQUEST_NULL};
    MPI_Request send_to_mastert_req = MPI_REQUEST_NULL;
    MPI_Status last_status;

    LOGF("W[%d, CT: %d]: CT created. Starting computations.", my_rank,
         my_cart_rank);

    // Sim
    size_t cols = c.cols / dims[0];
    size_t rows = c.rows / dims[1];
    size_t tot_cols = cols + 2;
    size_t tot_rows = rows + 2;

    // Init arrays
    double(*u_tmp)[tot_cols] = NULL;
    double(*u0)[tot_cols] = malloc(tot_rows * sizeof(double[tot_cols])); // u(k)
    double(*u1)[tot_cols] =
        malloc(tot_rows * sizeof(double[tot_cols])); // u(k-1)
    double(*u2)[tot_cols] =
        malloc(tot_rows * sizeof(double[tot_cols])); // u(k-2)

    // Init ghost cells
    size_t tot_gh_cells = 2 * (rows + cols);

    double *send_buf = malloc(tot_gh_cells * sizeof(double));
    double *send_buf_sides[4] = {send_buf, send_buf + cols,
                                 send_buf + cols + rows,
                                 send_buf + 2 * cols + rows};
    double *recv_buf = malloc(tot_gh_cells * sizeof(double));
    // {top, left, down, right} as before
    double *recv_buf_sides[4] = {recv_buf, recv_buf + cols,
                                 recv_buf + cols + rows,
                                 recv_buf + 2 * cols + rows};
    double gh_cells_counts[2] = {cols, rows};

    // Forcing origin
    int f_coord = dims[0] / 2;
    size_t f_offset = (dims[0] % 2 == 0) ? 1 : (tot_rows / 2);

    // Courant numbers
    double Cx = c.c * (c.dt / c.dx);
    double Cx_sq = Cx * Cx;
    double Cy = c.c * (c.dt / c.dy);
    double Cy_sq = Cy * Cy;

    // u1 and u2 initial conditions
    for (size_t j = 0; j < cols; j++) {
        for (size_t i = 0; i < rows; i++) {
            u2[i][j] = 0.0;
            u1[i][j] = 0.0;
            u0[i][j] = 0.0;
        }
    }

    // init send and recv cells
    for (size_t i = 0; i < tot_gh_cells; i++) {
        recv_buf[i] = 0.0;
        send_buf[i] = 0.0;
    }

    double t = c.dt;
    for (size_t step = 0; step < c.n_steps; step++) {
        t += c.dt;

        // Cell calc
        for (size_t j = 1; j < tot_cols - 1; j++) {
            for (size_t i = 1; i < tot_rows - 1; i++) {
                u0[i][j] = 2 * u1[i][j] - u2[i][j] +
                           (0.5 * Cx_sq) *
                               (u1[i][j - 1] - 2 * u1[i][j] + u1[i][j + 1]) +
                           (0.5 * Cy_sq) *
                               (u1[i - 1][j] - 2 * u1[i][j] + u1[i + 1][j]);
            }
        }

        // Forcing term (only in the region where it is applied)
        if (coord[0] == f_coord && coord[1] == f_coord) {
            u0[f_offset][f_offset] += c.dt * c.dt * 200 * sin(2 * M_PI * 2 * t);
        }

        // calc top and down boundary conds
        for (size_t s = 0; s < 4; s += 2) {
            // There is no neighbor so we enforce boundary conds
            if (neigh[s] == MPI_PROC_NULL) {
                size_t i = (s == 0 ? 0 : (tot_cols - 1));
                for (size_t j = 1; j < tot_cols - 1; j++)
                    u0[i][j] = (s == 0 ? u0[2][j] : u0[i - 2][j]);
            }
        }

        // calc left and right boundary conds
        for (size_t s = 1; s < 4; s += 2) {
            // There is no neighbor so we enforce boundary conds
            if (neigh[s] == MPI_PROC_NULL) {
                size_t i = (s == 1 ? 0 : (tot_rows - 1));
                for (size_t j = 1; j < tot_rows - 1; j++)
                    u0[j][i] = (s == 1 ? u0[j][2] : u0[j][i - 2]);
            }
        }

        // Send the ghost shells
        for (size_t s = 0; s < 4; s++) {
            if (neigh[s] != MPI_PROC_NULL) {
                // Copy corresponding side in send_buf
                double *send_buf_side = send_buf_sides[s];
                double *recv_buf_side = recv_buf_sides[s];

                // Copy the outer cells to send buf
                for (size_t i = 1; i < gh_cells_counts[s % 2] + 1; i++) {
                    if ((s % 2) == 0) // top and bottom
                        send_buf_side[i] = u0[s == 0 ? 1 : (tot_rows - 2)][i];
                    else // left and right
                        send_buf_side[i] = u0[i][s == 1 ? 1 : (tot_cols - 2)];
                }

                // Send the buf
                MPI_Wait(&gh_cells_req[s], &last_status);
                MPI_Isend(send_buf_side, gh_cells_counts[s % 2], MPI_DOUBLE,
                          neigh[s], GH_CELLS_TAG, comm2d, &gh_cells_req[s]);
                // Wait and recv the corresponding side
                MPI_Recv(recv_buf_side, gh_cells_counts[s % 2], MPI_DOUBLE,
                         neigh[s], GH_CELLS_TAG, comm2d, &last_status);

                // Copy back recved gh_cells
                for (size_t i = 1; i < gh_cells_counts[s % 2] + 1; i++) {
                    if ((s % 2) == 0) // top and bottom
                        u0[s == 0 ? 0 : (tot_rows - 1)][i] = recv_buf_side[i];
                    else // left and right
                        u0[i][s == 1 ? 0 : (tot_cols - 1)] = recv_buf_side[i];
                }
            }
        }

        // move the queue
        u_tmp = u2;
        u2 = u1;
        u1 = u0;
        u0 = u_tmp;

        // Send frame (u1) to master
        if (step % c.save_period == 0) {
            MPI_Wait(&send_to_mastert_req, &last_status);
            MPI_Isend(u1, tot_rows * tot_cols, MPI_DOUBLE, master_rank,
                      SEND_MASTER_TAG, MPI_COMM_WORLD, &send_to_mastert_req);
        }
    }

    LOGF("W[%d, CT: %d]: finished.", my_rank, my_cart_rank);
    // To be sure that every communication has ended
    MPI_Wait(&send_to_mastert_req, &last_status);
    MPI_Barrier(comm2d);

    free(u2);
    free(u1);
    free(u0);
    free(send_buf);
    free(recv_buf);

    MPI_Comm_free(&comm2d);

    return 0;
}