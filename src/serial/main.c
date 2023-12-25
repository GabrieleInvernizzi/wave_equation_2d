#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Constants
const char* filepath = "test.sim";

const double SIM_DURATION = 10;
const double dt = 0.01;
const size_t TIMESTEPS = (size_t)(SIM_DURATION / dt);

const size_t save_period = 2;
const unsigned int framerate = 1.0/(dt*save_period);

const double c = 1;
const double domain_width = 10;
const double domain_height = 10;
const double dx = 0.01;
const double dy = dx;
const size_t COLS = (size_t)(domain_width / dx);
const size_t ROWS = (size_t)(domain_height / dy);
const size_t TOT_COLS = COLS + 2;
const size_t TOT_ROWS = ROWS + 2;




int main() {
	// Create the file to save all the frames
	FILE* f = fopen(filepath, "wb");
	if (!f) {
		printf("Can't create the file \"%s\". Exiting.\n", filepath);
		return 1;
	}
	fprintf(f, "%zu-%zu-%zu-%u\n", TOT_ROWS, TOT_COLS, (size_t) (TIMESTEPS / save_period), framerate);

	// Init arrays
	double (*u_tmp)[TOT_COLS] = NULL;
	double (*u0)[TOT_COLS] = malloc(TOT_ROWS * sizeof(double [TOT_COLS]));	// u(k)
	double (*u1)[TOT_COLS] = malloc(TOT_ROWS * sizeof(double [TOT_COLS]));	// u(k-1)
	double (*u2)[TOT_COLS] = malloc(TOT_ROWS * sizeof(double [TOT_COLS]));	// u(k-2)
	

	// Courant numbers
	double Cx = c*(dt/dx);
	double Cx_sq = Cx*Cx;
	double Cy = c*(dt/dy);
	double Cy_sq = Cy*Cy;
	

	// u2 initial conditions
	for (size_t j = 0; j < TOT_COLS - 1; j++) {
		double y = j * dy;
		for (size_t i = 0; i < TOT_ROWS; i++) {
			double x = j * dx;
			u2[i][j] = sin(2*M_PI*(x+y));
		}
	}

	// u1 initial conditions
	for (size_t j = 1; j < TOT_COLS - 1; j++) {
		for (size_t i = 1; i < TOT_ROWS - 1; i++) {
			u1[i][j] = u2[i][j]+(0.25*Cx_sq)*(u2[i-1][j]-2*u2[i][j]+u2[i+1][j])
					   +(0.25*Cy_sq)*(u2[i][j-1]-2*u2[i][j]+u2[i][j+1]);
		}
	}

	double t = dt;
	for (size_t step = 0; step < TIMESTEPS; step++) {
		t += dt;


		for (size_t j = 1; j < TOT_COLS - 1; j++) {
			for (size_t i = 1; i < TOT_ROWS - 1; i++) {
				u0[i][j] =   2*u1[i][j]-u2[i][j]
					   + (0.5*Cx_sq)*(u1[i-1][j]-2*u1[i][j]+u1[i+1][j])
					   + (0.5*Cy_sq)*(u1[i][j-1]-2*u1[i][j]+u1[i][j+1]);
			}
		}

		// u0 boundary conditions
		for (size_t j = 0; j < TOT_COLS; j++) {
			double y = j * dy;
			double x = (TOT_ROWS - 1) * dx;
			u0[0][j] = sin(2*M_PI*(y+dt));
			u0[TOT_ROWS - 1][j] = sin(2*M_PI*(x+y+dt));
		}

		for (size_t i = 0; i < TOT_ROWS; i++) {
			double x = i * dx;
			double y = (TOT_COLS - 1) * dy;
			u0[i][0] = sin(2*M_PI*(x+t));
			u0[i][TOT_COLS - 1] = sin(2*M_PI*(x+y+t));
		}


		// move the queue
		u_tmp = u2;
		u2 = u1;
		u1 = u0;
		u0 = u_tmp;


		// Save the frame
		if (step % save_period == 0)
			fwrite((void*)u1, TOT_ROWS * TOT_COLS, sizeof(double), f);
	}

	fclose(f);

	free(u2);
	free(u1);
	free(u0);

	return 0;
}
