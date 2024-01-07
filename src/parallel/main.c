#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "sim_conf.h"


#define N_DIMS (2)
#define GH_CELLS_TAG (1)


int main(int argc, char *argv[]) {
	int n_procs, my_rank;
	int n_rows, n_cols;
	n_rows = 2;
	n_cols = n_rows;

	SimConf c = get_sim_conf();

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (n_procs != (n_cols*n_rows)) {
		if (my_rank == 0)
			printf("The number of processors is different from 4.\n");
		MPI_Finalize();
		return 1;
	}

	// Create cartesian topology
	int dims[N_DIMS] = { 0 };
	MPI_Dims_create(n_procs, N_DIMS, dims);
	if (my_rank == 0)
		printf("P x dim = [%d, %d]\n", dims[0], dims[1]);

	// Create cartesian mapping
	MPI_Comm comm2d;
	int wrap_around[N_DIMS] = { 0 };
	int reorder = 1;

	int err = 0;
	err = MPI_Cart_create(MPI_COMM_WORLD, N_DIMS, dims, wrap_around, reorder, &comm2d);
	if (err != 0) printf("Error creating the cart (%d).\n", err);
	
	int coord[N_DIMS];
	MPI_Cart_coords(comm2d, my_rank, N_DIMS, coord);
	int my_cart_rank;
	MPI_Cart_rank(comm2d, coord, &my_cart_rank);

	printf("P[%d]: my_cart_rank PCM[%d], my coords = (%d,%d)\n",
		my_rank, my_cart_rank, coord[0], coord[1]);


	// neigh = {top, left, down, right }
	int neigh[4] = { MPI_PROC_NULL };
	// Get index of neighbors
	MPI_Cart_shift(comm2d, 0, 1, &neigh[0], &neigh[2]);	
	MPI_Cart_shift(comm2d, 1, 1, &neigh[1], &neigh[3]);	
	printf("P[%d, cart: %d] - neighs: %d, %d, %d, %d\n",
		my_rank, my_cart_rank, neigh[0], neigh[1], neigh[2], neigh[3]);


	// Request for async
	MPI_Request gh_cells_req[4] = { MPI_REQUEST_NULL };
	MPI_Status last_status;


	// Sim
	
	size_t cols = c.cols / 4;
	size_t rows = c.rows / 4;
	size_t inn_cols = cols - 2;
	size_t inn_rows = rows - 2;

	// Init arrays
	double (*u_tmp)[inn_cols] = NULL;
	double (*u0)[inn_cols] = malloc(inn_rows * sizeof(double [inn_cols]));	// u(k)
	double (*u1)[inn_cols] = malloc(inn_rows * sizeof(double [inn_cols]));	// u(k-1)
	double (*u2)[inn_cols] = malloc(inn_rows * sizeof(double [inn_cols]));	// u(k-2)
	
	// Init ghost cells
	size_t tot_gh_cells = 2*(rows+cols);
	// Outer cells
	double* out_cells = malloc(tot_gh_cells*sizeof(double));
	double* out_cells_sides[4] = {
		out_cells,
		out_cells + cols,
		out_cells + cols + rows,
		out_cells + 2*cols + rows
	};

	double* gh_cells = malloc(tot_gh_cells*sizeof(double));
	// {top, left, down, right} as before
	double* gh_cells_sides[4] = {
		gh_cells,
		gh_cells + cols,
		gh_cells + cols + rows,
		gh_cells + 2*cols + rows
	};
	double gh_cells_counts[2] = {cols, rows};


	// Courant numbers
	double Cx = c.c*(c.dt/c.dx);
	double Cx_sq = Cx*Cx;
	double Cy = c.c*(c.dt/c.dy);
	double Cy_sq = Cy*Cy;



	// u2 initial conditions
	for (size_t j = 0; j < cols; j++)
		for (size_t i = 0; i < rows; i++)
			u2[i][j] = 0.0;

	// u1 initial conditions
	for (size_t j = 1; j < cols; j++)
		for (size_t i = 1; i < rows; i++)
			u1[i][j] = 0.0;

	// init gh_cells and outer cells
	for (size_t i = 0; i < tot_gh_cells; i++) {
		gh_cells[i] = 0.0;
		out_cells[i] = 1.0;
	}


	double t = c.dt;
	for (size_t step = 0; step < c.n_steps; step++) {
		t += c.dt;


		// Inner cells calc
		for (size_t j = 1; j < inn_cols - 1; j++) {
			for (size_t i = 1; i < inn_rows - 1; i++) {
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
				MPI_Isend(
					out_cells_sides[i], 
					gh_cells_counts[i%2], 
					MPI_DOUBLE, 
					neigh[i], 
					GH_CELLS_TAG, 
					comm2d, 
					&gh_cells_req[i]);
				MPI_Recv(
					gh_cells_sides[i], 
					gh_cells_counts[i%2], 
					MPI_DOUBLE,
					neigh[i], 
					GH_CELLS_TAG,
					comm2d, 
					&last_status);
			}
		}

		// calc out cells
		// for (size_t i = 0; i < rows; i++) {
		// 	out_cells_sides[0][i] =   2*u1[i][j]-u2[i][j]
		// 			+ (0.5*Cx_sq)*(u1[i][j-1]-2*u1[i][j]+u1[i][j+1])
		// 		        + (0.5*Cy_sq)*(u1[i-1][j]-2*u1[i][j]+u1[i+1][j]);
		// }
		

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
		// if (step % c.save_period == 0) 
		// 	fwrite((void*)u1, sizeof(double), c.tot_cols * c.tot_rows, f);
	}

	free(u2);
	free(u1);
	free(u0);
	free(gh_cells);



	MPI_Comm_free(&comm2d);
	MPI_Finalize();

	return 0;
}
