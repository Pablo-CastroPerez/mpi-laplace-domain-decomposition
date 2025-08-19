#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

// Expose the global solver arrays (allocated in solvers.c, used in utils.c and main.c)
extern double *u_k;
extern double *u_kp1;
extern double *temp;

// Initialise the arrays with boundary conditions
void initialise(double *u_k, double *u_kp1, int nx, int ny);

// Run the iterative solver with function pointers to the solver routines
void run_solver(int local_nx, int ny, int myrank, int size,
                double convergence_accuracy,
                void (*init_solver)(int, int),
                void (*go_solve)(int, int, int, double*, double*),
                void (*finalise_solver)());

// Exchange halo regions between neighbouring MPI ranks
void perform_halo_swap(int myrank, int size, int local_nx,
                       int ny, int mem_size_y, double *data);

// Compute the global residual
double get_residual(double *data, int local_nx, int ny, int mem_size_y);

#endif // UTILS_H
