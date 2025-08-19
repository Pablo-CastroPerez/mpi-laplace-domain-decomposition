#ifndef SOLVERS_H
#define SOLVERS_H


/* Jacobi */
void init_jacobi_solvers(int mem_size_x, int mem_size_y);
void finalise_jacobi_solvers();
void jacobi_solver(int local_nx, int ny, int mem_size_y,
                   double *u_k, double *u_kp1);
void jacobi_sor_solver(int local_nx, int ny, int mem_size_y,
                       double *u_k, double *u_kp1);

/* Gauss–Seidel */
void init_gauss_seidel_solvers(int mem_size_x, int mem_size_y);
void finalise_gauss_seidel_solvers();
void gauss_seidel_solver(int local_nx, int ny, int mem_size_y,
                         double *u_k, double *u_kp1);
void gauss_seidel_sor_solver(int local_nx, int ny, int mem_size_y,
                             double *u_k, double *u_kp1);

#endif /* SOLVERS_H */
