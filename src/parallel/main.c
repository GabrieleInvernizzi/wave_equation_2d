#include <math.h>
#include <stdio.h>

#include <mpi.h>

#include "master_worker.h"

int main(int argc, char *argv[]) {
    int n_procs_world, my_rank_world;
    int side_len;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs_world);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_world);

    side_len = (int)sqrt(n_procs_world - 1);

    if ((n_procs_world - 1) != (side_len * side_len)) {
        if (my_rank_world == 0)
            printf("M[]: the number of processors must be a perfect square + 1 "
                   "but is %d.\n",
                   n_procs_world);
        MPI_Finalize();
        return 1;
    }

    // Split master and workers
    MPI_Comm worker_comm = MPI_COMM_NULL;
    MPI_Comm_split(
        MPI_COMM_WORLD,
        (my_rank_world == MASTER_RANK ? MPI_UNDEFINED : WORKER_COLOR),
        my_rank_world, &worker_comm);

    int ret;
    if (my_rank_world == MASTER_RANK) {
        ret = master(n_procs_world);
    } else {
        ret = worker(worker_comm, my_rank_world, MASTER_RANK);
    }

    if (worker_comm != MPI_COMM_NULL)
        MPI_Comm_free(&worker_comm);

    MPI_Finalize();

    return ret;
}
