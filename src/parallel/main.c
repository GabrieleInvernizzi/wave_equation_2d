#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "sim_conf.h"

#define MASTER_RANK (0)
#define N_DIMS (2)
#define GH_CELLS_TAG (1)
#define WORKERS_COLOR (1)


int master(int n_procs_world) {
	FILE*f = NULL;

	SimConf c = get_sim_conf();

	printf("M[]: started.\n");

	f = fopen(c.filepath, "wb");
	if (f == NULL) {
		printf("Could not open the file \"%s\".", c.filepath);
		MPI_Abort(MPI_COMM_WORLD, 1);
		return 1;
	}
	
	fprintf(f, "%zu-%zu-%zu-%zu-%u\n", sizeof(double), c.tot_rows, c.tot_cols, (size_t) (c.n_steps / c.save_period), c.framerate);

	fclose(f);

	printf("M[]: finished.\n");

	return 0;
}



int worker(MPI_Comm comm, int master_rank) {
	int my_rank;
	int n_procs;

	SimConf c = get_sim_conf();

	MPI_Comm_size(comm, &n_procs);
	MPI_Comm_rank(comm, &my_rank);

	if (n_procs != 4) {
		if (my_rank == 0) {
			printf("W[0]: worker processes must be 4 but are %d.", n_procs);
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
	}

	printf("W[%d]: started.\n", my_rank);

	// Create cartesian topology
	int dims[N_DIMS] = { 0 };
	MPI_Dims_create(n_procs, N_DIMS, dims);
	if (my_rank == 0)
		printf("W[0]: worker grid = [%d, %d]\n", dims[0], dims[1]);

	// Create cartesian mapping
	MPI_Comm comm2d;
	int wrap_around[N_DIMS] = { 0 };
	int reorder = 1;

	int err = 0;
	err = MPI_Cart_create(comm, N_DIMS, dims, wrap_around, reorder, &comm2d);
	if (err != 0) printf("W[%d]: Error creating the cart (%d).\n", my_rank, err);

	int coord[N_DIMS];
	MPI_Cart_coords(comm2d, my_rank, N_DIMS, coord);
	int my_cart_rank;
	MPI_Cart_rank(comm2d, coord, &my_cart_rank);

	printf("W[%d]: my_cart_rank WCT[%d], my coords = (%d,%d).\n",
	my_rank, my_cart_rank, coord[0], coord[1]);


	// neigh = {top, left, down, right }
	int neigh[4] = { MPI_PROC_NULL };
	// Get index of neighbors
	MPI_Cart_shift(comm2d, 0, 1, &neigh[0], &neigh[2]);	
	MPI_Cart_shift(comm2d, 1, 1, &neigh[1], &neigh[3]);	
	printf("W[%d, CT: %d] - neighs: [%d, %d, %d, %d].\n",
	my_rank, my_cart_rank, neigh[0], neigh[1], neigh[2], neigh[3]);


	// Request for async
	MPI_Request gh_cells_req[4] = { MPI_REQUEST_NULL };
	MPI_Status last_status;

	printf("W[%d, CT: %d]: CT created. Starting computations.\n", my_rank, my_cart_rank);


	// Sim
	size_t cols = c.cols / 4;
	size_t rows = c.rows / 4;
	size_t tot_cols = cols + 2;
	size_t tot_rows = rows + 2;

	// Init arrays
	double (*u_tmp)[tot_cols] = NULL;
	double (*u0)[tot_cols] = malloc(tot_rows * sizeof(double [tot_cols]));	// u(k)
	double (*u1)[tot_cols] = malloc(tot_rows * sizeof(double [tot_cols]));	// u(k-1)
	double (*u2)[tot_cols] = malloc(tot_rows * sizeof(double [tot_cols]));	// u(k-2)

	// Init ghost cells
	size_t tot_gh_cells = 2*(tot_rows+tot_cols);

	double* send_buf = malloc(tot_gh_cells*sizeof(double));
	double* send_buf_sides[4] = {
		send_buf,
		send_buf + tot_cols,
		send_buf + tot_cols + tot_rows,
		send_buf + 2*tot_cols + tot_rows
	};
	double* recv_buf = malloc(tot_gh_cells*sizeof(double));
	// {top, left, down, right} as before
	double* recv_buf_sides[4] = {
		recv_buf,
		recv_buf + tot_cols,
		recv_buf + tot_cols + tot_rows,
		recv_buf + 2*tot_cols + tot_rows
	};
	double gh_cells_counts[2] = {tot_cols, tot_rows};


	// Courant numbers
	double Cx = c.c*(c.dt/c.dx);
	double Cx_sq = Cx*Cx;
	double Cy = c.c*(c.dt/c.dy);
	double Cy_sq = Cy*Cy;



	// u1 and u2 initial conditions
	for (size_t j = 0; j < cols; j++) {
		for (size_t i = 0; i < rows; i++) {
			u2[i][j] = 0.0;
			u1[i][j] = 0.0;
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
				u0[i][j] =   2*u1[i][j]-u2[i][j]
					+ (0.5*Cx_sq)*(u1[i][j-1]-2*u1[i][j]+u1[i][j+1])
					+ (0.5*Cy_sq)*(u1[i-1][j]-2*u1[i][j]+u1[i+1][j]);
			}
		}


		// Forcing term (only in the region where it is applied)
		if (coord[0] == 1 && coord[1] == 1) {
			u0[0][0] +=
			c.dt*c.dt*200*sin(2*M_PI*2*t);
		}


		// Send the ghost shells
		for (size_t i = 0; i < 4; i++) {
			if (neigh[i] != MPI_PROC_NULL) {
				// Copy corresponding side in send_buf
				double* send_buf_side = send_buf_sides[i];
				double* recv_buf_side = recv_buf_sides[i];
				for (size_t j = 0; j < gh_cells_counts[i%2]; j++) {
					if ((i % 2) == 0)	// top and bottom
						send_buf_side[j] = u1[i == 0 ? 0 : (tot_rows-1)][j];
					else			// left and right
						send_buf_side[j] = u1[j][i == 1 ? 0 : (tot_cols-1)];
				}
				// Send the buf
				MPI_Isend(
					send_buf_side, 
					gh_cells_counts[i%2], 
					MPI_DOUBLE, 
					neigh[i], 
					GH_CELLS_TAG, 
					comm2d, 
				&gh_cells_req[i]);
				// Wait and recv the corresponding side
				MPI_Recv(
					recv_buf_side, 
					gh_cells_counts[i%2], 
					MPI_DOUBLE,
					neigh[i], 
					GH_CELLS_TAG,
					comm2d, 
				&last_status);
				// Copy back recved gh_cells
				for (size_t j = 0; j < gh_cells_counts[i%2]; j++) {
					if ((i % 2) == 0)	// top and bottom
						u1[i == 0 ? 0 : (tot_rows-1)][j] = recv_buf_side[j];
					else			// left and right
						u1[j][i == 1 ? 0 : (tot_cols-1)] = recv_buf_side[j];
				}
			}
		}

		// calc top and down out cells
		for (size_t s = 0; s < 4; s+=2) {
			double* gh   = recv_buf_sides[s];
			size_t i = (s == 0 ? 0 : (tot_cols - 1));
			if (neigh[s] != MPI_PROC_NULL) {		// There is a neighbor so we use it's ghost cells
				for (size_t j = 0; j < tot_cols; j++) {
					u0[i][j] =   2*u1[i][j]-u2[i][j]
						+ (0.5*Cx_sq)*(u1[i][j-1]-2*u1[i][j]+u1[i][j+1])
						+ (0.5*Cy_sq)*(u1[i-1][j]-2*u1[i][j]+u1[i+1][j]);
				}
			} else {					// There is no neighbor so we enforce boundary conds
				for (size_t j = 0; j < tot_cols; j++)
					u0[i][j] = (s == 0 ? u0[2][j] : u0[i - 2][j]);
			}
		}


		// calc left and right out cells
		for (size_t s = 1; s < 4; s+=2) {
			double* gh   = recv_buf_sides[s];
			size_t i = (s == 1 ? 0 : (tot_rows - 1));
			if (neigh[s] != MPI_PROC_NULL) {		// There is a neighbor so we use it's ghost cells
				for (size_t j = 0; j < tot_cols; j++) {
					u0[i][j] =   2*u1[i][j]-u2[i][j]
						+ (0.5*Cx_sq)*(u1[i][j-1]-2*u1[i][j]+u1[i][j+1])
						+ (0.5*Cy_sq)*(u1[i-1][j]-2*u1[i][j]+u1[i+1][j]);
				}
			} else {					// There is no neighbor so we enforce boundary conds
				for (size_t j = 0; j < tot_cols; j++)
					u0[j][i] = (s == 0 ? u0[j][2] : u0[j][i - 2]);
			}
		}


		// move the queue
		u_tmp = u2;
		u2 = u1;
		u1 = u0;
		u0 = u_tmp;


		// Send frame(u1) to master
		// if (step % c.save_period == 0)  
		// 	fwrite((void*)u1, sizeof(double), c.tot_cols * c.tot_rows, f);
	}

	printf("W[%d, CT: %d]: finished.\n", my_rank, my_cart_rank);
	// To be sure that every communication has ended
	MPI_Barrier(comm2d);

	free(u2);
	free(u1);
	free(u0);
	free(send_buf);
	free(recv_buf);



	MPI_Comm_free(&comm2d);


	return 0;
}


int main(int argc, char *argv[]) {
	int n_procs_world, my_rank_world;
	int n_rows, n_cols;
	n_rows = 2;
	n_cols = n_rows;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs_world);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_world);

	if (n_procs_world  != ((n_cols*n_rows) + 1)) {
		if (my_rank_world == 0)
			printf("M[]: the number of processors must be 5 but is %d.\n", n_procs_world);
		MPI_Finalize();
		return 1;
	}
	

	// Split master and workers
	MPI_Comm worker_comm = MPI_COMM_NULL;
	MPI_Comm_split(
		MPI_COMM_WORLD, 
		(my_rank_world == MASTER_RANK ? MPI_UNDEFINED : WORKERS_COLOR), 
		my_rank_world, 
		&worker_comm);
	
	int ret;
	if (my_rank_world == MASTER_RANK) {	
		ret = master(n_procs_world);	
	} else {
		ret = worker(worker_comm, MASTER_RANK);	
	}

	if (worker_comm != MPI_COMM_NULL)
		MPI_Comm_free(&worker_comm);

	MPI_Finalize();

	return ret;
}
