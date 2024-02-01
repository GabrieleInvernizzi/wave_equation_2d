#pragma once

#include <mpi.h>

#include "sim_conf.h"

#define MASTER_RANK (0)
#define N_DIMS (2)
#define GH_CELLS_TAG (1)
#define SEND_MASTER_TAG (2)
#define WORKER_COLOR (1)

int master(SimConf c, int n_procs_world);
int worker(SimConf c ,MPI_Comm comm, int my_rank_world, int master_rank);
