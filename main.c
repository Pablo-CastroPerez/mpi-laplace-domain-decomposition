#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "config.h"   // constants, LEFT_VALUE, RIGHT_VALUE, etc.
#include "solvers.h"  // init/finalise + solver kernels
#include "utils.h"    // run_solver, initialise, halo swap, residual

int main(int argc, char * argv[]) {
	int size, myrank, nx, ny, max_its;
	double convergence_accuracy;
  
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (argc != 5) {
		if (myrank==0) {
			printf("You should provide four command line arguments, the global size in X, the global size in Y, convergence accuracy and max number iterations\n");		
			printf("In the absence of this defaulting to x=128, y=1024, convergence=3e-3, no max number of iterations\n");
		}
		nx=128;
		ny=1024;
		convergence_accuracy=3e-3;
		max_its=0;
	} else {
		nx=atoi(argv[1]);
		ny=atoi(argv[2]);
		convergence_accuracy=atof(argv[3]);
		max_its=atoi(argv[4]);
	}

	if (myrank==0) {
    		printf("Number Processes in X=%d\n", size);
    		printf("Global size in X=%d, Global size in Y=%d\n\n", nx, ny);
  	}

	int local_nx=nx/size;
  	if (local_nx * size < nx) {
    		if (myrank < nx - local_nx * size) local_nx++;
  	}


// Choose which solver to use, SOLVER_TO_USE is defined in config.h
  
#if SOLVER_TO_USE == 0
	run_solver(local_nx, ny, myrank, size, convergence_accuracy, init_jacobi_solvers, jacobi_solver, finalise_jacobi_solvers);
#elif SOLVER_TO_USE == 1
	run_solver(local_nx, ny, myrank, size, convergence_accuracy, init_jacobi_solvers, jacobi_sor_solver, finalise_jacobi_solvers);
#elif SOLVER_TO_USE == 2
	run_solver(local_nx, ny, myrank, size, convergence_accuracy, init_gauss_seidel_solvers, gauss_seidel_solver, finalise_gauss_seidel_solvers);
#elif SOLVER_TO_USE == 3
	run_solver(local_nx, ny, myrank, size, convergence_accuracy, init_gauss_seidel_solvers, gauss_seidel_sor_solver, finalise_gauss_seidel_solvers);
#endif

	MPI_Finalize();
	return 0;
}

