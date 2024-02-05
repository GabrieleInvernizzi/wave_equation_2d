#include "arena_alloc.h"
#include "master_worker.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "log.h"
#include "timings.h"

static void set_init_conds(size_t tot_n_cells, double (*u0)[tot_n_cells],
                           double (*u1)[tot_n_cells],
                           double (*u2)[tot_n_cells]) {
    for (size_t i = 0; i < tot_n_cells; i++) {
        for (size_t j = 0; j < tot_n_cells; j++) {
            u2[i][j] = 0.0;
            u1[i][j] = 0.0;
            u0[i][j] = 0.0;
        }
    }
}

static inline void calc_u0_inner_cells(size_t tot_n_cells,
                                       double (*u0)[tot_n_cells],
                                       double (*u1)[tot_n_cells],
                                       double (*u2)[tot_n_cells], double C_sq) {

    // Cell calc
    for (size_t i = 1; i < tot_n_cells - 1; i++) {
        for (size_t j = 1; j < tot_n_cells - 1; j++) {
            u0[i][j] =
                2 * u1[i][j] - u2[i][j] +
                (0.5 * C_sq) * (u1[i][j - 1] - 2 * u1[i][j] + u1[i][j + 1]) +
                (0.5 * C_sq) * (u1[i - 1][j] - 2 * u1[i][j] + u1[i + 1][j]);
        }
    }
}

static inline void calc_forcing(size_t tot_n_cells, double (*u0)[tot_n_cells],
                                int coords[2], int f_coord, size_t f_offset,
                                double dt, double t) {
    // Forcing term (only in the region where it is applied)
    if (coords[0] == f_coord && coords[1] == f_coord) {
        u0[f_offset][f_offset] += dt * dt * 200 * sin(2 * M_PI * 2 * t);
    }
}

static inline void calc_boundary_conds(size_t tot_n_cells,
                                       double (*u0)[tot_n_cells],
                                       int neighs[4]) {

    // calc top and down boundary conds
    for (size_t s = 0; s < 4; s += 2) {
        // There is no neighbor so we enforce boundary conds
        if (neighs[s] == MPI_PROC_NULL) {
            size_t i = (s == 0 ? 0 : (tot_n_cells - 1));
            for (size_t j = 1; j < tot_n_cells - 1; j++)
                u0[i][j] = (s == 0 ? u0[2][j] : u0[i - 2][j]);
        }
    }

    // calc left and right boundary conds
    for (size_t s = 1; s < 4; s += 2) {
        // There is no neighbor so we enforce boundary conds
        if (neighs[s] == MPI_PROC_NULL) {
            size_t i = (s == 1 ? 0 : (tot_n_cells - 1));
            for (size_t j = 1; j < tot_n_cells - 1; j++)
                u0[j][i] = (s == 1 ? u0[j][2] : u0[j][i - 2]);
        }
    }
}

int worker(SimConf c, MPI_Comm comm, int my_rank_world, int master_rank) {
    int my_rank;
    int n_procs;

    START_TIMER("w init mpi", my_rank_world);

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
    MPI_Comm comm2d;
    int dims[N_DIMS] = {0};
    int coords[N_DIMS];
    int my_cart_rank;
    {
        MPI_Dims_create(n_procs, N_DIMS, dims);
        if (my_rank == 0)
            LOGF("W[0]: worker grid = [%d, %d].", dims[0], dims[1]);

        // Create cartesian mapping
        int wrap_around[N_DIMS] = {0};
        int reorder = 1;

        int err =
            MPI_Cart_create(comm, N_DIMS, dims, wrap_around, reorder, &comm2d);
        if (err != 0)
            LOGF("W[%d]: Error creating the cart (%d).", my_rank, err);

        MPI_Cart_coords(comm2d, my_rank, N_DIMS, coords);
        MPI_Cart_rank(comm2d, coords, &my_cart_rank);

        LOGF("W[%d]: my_cart_rank CT[%d], my coords = (%d,%d).", my_rank,
             my_cart_rank, coords[0], coords[1]);
    }
    // Send to master my coords
    {
        int send_gather[3] = {my_rank_world, coords[0], coords[1]};
        MPI_Gather(send_gather, 3, MPI_INT, NULL, 0, MPI_INT, master_rank,
                   MPI_COMM_WORLD);
    }

    // neigh = {top, left, down, right }
    int neighs[4] = {[0 ... 3] = MPI_PROC_NULL};
    // Get index of neighbors
    MPI_Cart_shift(comm2d, 0, 1, &neighs[0], &neighs[2]);
    MPI_Cart_shift(comm2d, 1, 1, &neighs[1], &neighs[3]);
    LOGF("W[%d, CT: %d] - neighs: [%d, %d, %d, %d].", my_rank, my_cart_rank,
         neighs[0], neighs[1], neighs[2], neighs[3]);

    // Request for async
    MPI_Request gh_cells_req[4] = {[0 ... 3] = MPI_REQUEST_NULL};
    MPI_Request send_to_mastert_req = MPI_REQUEST_NULL;
    MPI_Status last_status;

    LOGF("W[%d, CT: %d]: CT created. Starting computations.", my_rank,
         my_cart_rank);

    END_TIMER;

    START_TIMER("w init sim", my_rank_world);
    // Sim
    size_t n_cells = c.n_cells / dims[0];
    size_t tot_n_cells = n_cells + 2;

    // Init arena
    const size_t tot_gh_cells = 2 * (n_cells + n_cells);
    const size_t arena_size = 3 * tot_n_cells * sizeof(double[tot_n_cells]) +
                              2 * tot_gh_cells * sizeof(double);

    Arena arena;
    if (arena_init(&arena, arena_size) != 0)
        LOGF("W[%d, CT: %d]: Error allocating %zu B of memory.", my_rank,
             my_cart_rank, arena_size);

    // Init arrays
    double(*u_tmp)[tot_n_cells] = NULL;
    double(*u0)[tot_n_cells] =
        arena_alloc(&arena, tot_n_cells * sizeof(double[tot_n_cells])); // u(k)
    double(*u1)[tot_n_cells] = arena_alloc(
        &arena, tot_n_cells * sizeof(double[tot_n_cells])); // u(k-1)
    double(*u2)[tot_n_cells] = arena_alloc(
        &arena, tot_n_cells * sizeof(double[tot_n_cells])); // u(k-2)
    assert(u0);
    assert(u1);
    assert(u2);

    // Init ghost cells
    double *send_buf = arena_alloc(&arena, tot_gh_cells * sizeof(double));
    double *send_buf_sides[4] = {send_buf, send_buf + n_cells,
                                 send_buf + n_cells + n_cells,
                                 send_buf + 2 * n_cells + n_cells};
    double *recv_buf = arena_alloc(&arena, tot_gh_cells * sizeof(double));
    assert(send_buf);
    assert(recv_buf);

    // {top, left, down, right} as before
    double *recv_buf_sides[4] = {recv_buf, recv_buf + n_cells,
                                 recv_buf + n_cells + n_cells,
                                 recv_buf + 2 * n_cells + n_cells};
    double gh_cells_counts[2] = {n_cells, n_cells};
    // init send and recv cells
    for (size_t i = 0; i < tot_gh_cells; i++) {
        recv_buf[i] = 0.0;
        send_buf[i] = 0.0;
    }

    // Forcing origin
    int f_coord = dims[0] / 2;
    size_t f_offset = (dims[0] % 2 == 0) ? 1 : (tot_n_cells / 2);

    // Courant numbers
    double C = c.c * (c.dt / c.dx);
    double C_sq = C * C;

    // Check CFL condition
    if (C > 0.5 && (!c.ignore_cfl)) {
        if (my_cart_rank == 0) {
            LOGF(
                "The CFL condition is not satisfied.\nC = c * (dt / dx) = %f > "
                "0.5.\nChoose the parameters so that C <= 0.5.\nExiting.",
                C);
        }
        return 1;
    }

    set_init_conds(tot_n_cells, u0, u1, u2);

    END_TIMER;

    double t = c.dt;
    for (size_t step = 0; step < c.n_steps; step++) {
        t += c.dt;

        START_TIMER("w calc", my_rank_world);
        calc_u0_inner_cells(tot_n_cells, u0, u1, u2, C_sq);
        calc_forcing(tot_n_cells, u0, coords, f_coord, f_offset, c.dt, t);
        calc_boundary_conds(tot_n_cells, u0, neighs);
        END_TIMER;

        START_TIMER("w comm", my_rank_world);
        // Send the ghost shells
        for (size_t s = 0; s < 4; s++) {
            if (neighs[s] != MPI_PROC_NULL) {
                // Copy corresponding side in send_buf
                double *send_buf_side = send_buf_sides[s];
                double *recv_buf_side = recv_buf_sides[s];

                // Copy the outer cells to send buf
                for (size_t i = 0; i < gh_cells_counts[s % 2]; i++) {
                    if ((s % 2) == 0) // top and bottom
                        send_buf_side[i] =
                            u0[s == 0 ? 1 : (tot_n_cells - 2)][i];
                    else // left and right
                        send_buf_side[i] =
                            u0[i][s == 1 ? 1 : (tot_n_cells - 2)];
                }

                // Send the buf
                MPI_Wait(&gh_cells_req[s], &last_status);
                MPI_Isend(send_buf_side, gh_cells_counts[s % 2], MPI_DOUBLE,
                          neighs[s], GH_CELLS_TAG, comm2d, &gh_cells_req[s]);
                // Wait and recv the corresponding side
                MPI_Recv(recv_buf_side, gh_cells_counts[s % 2], MPI_DOUBLE,
                         neighs[s], GH_CELLS_TAG, comm2d, &last_status);

                // Copy back recved gh_cells
                for (size_t i = 0; i < gh_cells_counts[s % 2]; i++) {
                    if ((s % 2) == 0) // top and bottom
                        u0[s == 0 ? 0 : (tot_n_cells - 1)][i] =
                            recv_buf_side[i];
                    else // left and right
                        u0[i][s == 1 ? 0 : (tot_n_cells - 1)] =
                            recv_buf_side[i];
                }
            }
        }
        END_TIMER;

        // move the queue
        u_tmp = u2;
        u2 = u1;
        u1 = u0;
        u0 = u_tmp;

        START_TIMER("w send m", my_rank_world);
        // Send frame (u1) to master
        if (step % c.save_period == 0) {
            MPI_Wait(&send_to_mastert_req, &last_status);
            MPI_Isend(u1, tot_n_cells * tot_n_cells, MPI_DOUBLE, master_rank,
                      SEND_MASTER_TAG, MPI_COMM_WORLD, &send_to_mastert_req);
        }
        END_TIMER;
    }

    LOGF("W[%d, CT: %d]: finished.", my_rank, my_cart_rank);
    // To be sure that every communication has ended
    MPI_Wait(&send_to_mastert_req, &last_status);
    MPI_Barrier(comm2d);

    area_deinit(&arena);

    MPI_Comm_free(&comm2d);

    return 0;
}
