#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "config.h"  
#include "utils.h"

double *u_k   = NULL;
double *u_kp1 = NULL;
double *temp  = NULL;

void initialise(double * u_k, double * u_kp1, int nx, int ny) {
	int i,j;
	// We are setting the boundary (left and right) values here, in the parallel version this should be exactly the same and no changed required
	for (i=0;i<nx+1;i++) {
		u_k[i*(ny+2)]=LEFT_VALUE;
		u_k[(ny+1)+(i*(ny+2))]=RIGHT_VALUE;
	}
	for (j=0;j<=nx+1;j++) {
		for (i=1;i<=ny;i++) {
			u_k[i+(j*(ny+2))]=0.0;
		}
	}
	if (u_kp1 != NULL) {
		for (j=0;j<=nx+1;j++) {
			for (i=0;i<=ny+1;i++) {
				u_kp1[i+(j*(ny+2))]=u_k[i+(j*(ny+2))];
			}
		}
	}
}


void run_solver(int local_nx, int ny, int myrank, int size, double convergence_accuracy, void (*init_solver)(int, int), void (*go_solve)(int, int, int, double*, double*), void (*finalise_solver)()) {	
	int mem_size_x=local_nx+2;
	int mem_size_y=ny+2;

	init_solver(mem_size_x, mem_size_y);

	double start_time;

	initialise(u_k, u_kp1, local_nx, ny);

	double rnorm=0.0, bnorm=0.0, norm, tmpnorm=0.0;
	bnorm=get_residual(u_k, local_nx, ny, mem_size_y);

	int k;
	start_time=MPI_Wtime();
	for (k=0;k<MAX_ITERATIONS;k++) {
		perform_halo_swap(myrank, size, local_nx, ny, mem_size_y, u_k);
		rnorm=get_residual(u_k, local_nx, ny, mem_size_y);
		go_solve(local_nx, ny, mem_size_y, u_k, u_kp1);
		if (u_kp1 != NULL) {
			temp=u_kp1;
			u_kp1=u_k;
			u_k=temp;
		}
		norm=rnorm/bnorm;
		if (norm < convergence_accuracy) break;
		if (norm > 1) break;
		if (k % REPORT_NORM_PERIOD == 0 && myrank==0) printf("Iteration= %d Relative Norm=%e\n", k, norm);
	}
	if (myrank==0) {
		if (norm > 1) printf("Aborted due to divergence\n");
		printf("\nTerminated on %d iterations, Relative Norm=%e, Total time=%e seconds\n", k, norm,
				MPI_Wtime() - start_time);			
	}
	finalise_solver();
}


// Performs the halo swap in one dimension

void perform_halo_swap(int myrank, int size, int local_nx, int ny, int mem_size_y, double * data) {
	MPI_Request requests[]={MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};	

	if (myrank > 0) {
		MPI_Isend(&data[1+mem_size_y], ny, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD, &requests[0]);
		MPI_Irecv(&data[1], ny, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD, &requests[1]);
    	}
   	if (myrank < size-1) {
      		MPI_Isend(&data[(local_nx*mem_size_y)+1], ny, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, &requests[2]);
		MPI_Irecv(&data[((local_nx+1)*mem_size_y)+1], ny, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, &requests[3]);
   	}

	MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
}


// Determines the global residual of some data provided

double get_residual(double * data, int local_nx, int ny, int mem_size_y) {
	double tmpnorm=0.0, global_norm;
	int i,j;
	// Calculate the initial residual norm
	for (j=1;j<=local_nx;j++) {
		for (i=1;i<=ny;i++) {
			tmpnorm=tmpnorm+pow(data[i+(j*mem_size_y)]*4-data[(i-1)+(j*mem_size_y)]-
					data[(i+1)+(j*mem_size_y)]-data[i+((j-1)*mem_size_y)]-data[i+((j+1)*mem_size_y)], 2);
		}
	}
	MPI_Allreduce(&tmpnorm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	global_norm=sqrt(global_norm);
	return global_norm;
