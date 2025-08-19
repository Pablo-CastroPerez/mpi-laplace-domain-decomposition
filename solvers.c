#include <stdlib.h>
#include "config.h"   
#include "solvers.h"  


extern double *u_k;
extern double *u_kp1;


void init_jacobi_solvers(int mem_size_x, int mem_size_y) {
	u_k = malloc(sizeof(double) * mem_size_x * mem_size_y);
	u_kp1 = malloc(sizeof(double) * mem_size_x * mem_size_y);	
}

void finalise_jacobi_solvers() {
	free(u_k);
	free(u_kp1);
}


// Jacobi solver with over relaxation
 
void jacobi_sor_solver(int local_nx, int ny,  int mem_size_y, double * u_k, double * u_kp1) {
	int i, j;
	for (j=1;j<=local_nx;j++) {
		for (i=1;i<=ny;i++) {
			u_kp1[i+(j*mem_size_y)]= ((1-W) * u_kp1[i+(j*mem_size_y)]) + W * 0.25 * (u_k[(i-1)+(j*mem_size_y)]+
				u_k[(i+1)+(j*mem_size_y)]+u_k[i+((j-1)*mem_size_y)]+u_k[i+((j+1)*mem_size_y)]);
		}
	}
}

 // Simple Jacobi solver

void jacobi_solver(int local_nx, int ny,  int mem_size_y, double * u_k, double * u_kp1) {
	int i, j;
	for (j=1;j<=local_nx;j++) {
		for (i=1;i<=ny;i++) {
			u_kp1[i+(j*mem_size_y)]=0.25 * (u_k[(i-1)+(j*mem_size_y)]+u_k[(i+1)+(j*mem_size_y)]+
				u_k[i+((j-1)*mem_size_y)]+u_k[i+((j+1)*mem_size_y)]);
		}
	}
}

void init_gauss_seidel_solvers(int mem_size_x, int mem_size_y) {
	u_k = malloc(sizeof(double) * mem_size_x * mem_size_y);	
}

void finalise_gauss_seidel_solvers() {
	free(u_k);
}


 // Gauss Seidel solver solver with over relaxation (SOR)

void gauss_seidel_sor_solver(int local_nx, int ny,  int mem_size_y, double * u_k, double * u_kp1) {
	int i, j;
	for (j=1;j<=local_nx;j++) {
		for (i=1;i<=ny;i++) {
			u_k[i+(j*mem_size_y)]= ((1-W) * u_k[i+(j*mem_size_y)]) + W * 0.25 * (u_k[(i-1)+(j*mem_size_y)]+
				u_k[(i+1)+(j*mem_size_y)]+u_k[i+((j-1)*mem_size_y)]+u_k[i+((j+1)*mem_size_y)]);
		}
	}
}


// Simple Gauss Seidel solver
 
void gauss_seidel_solver(int local_nx, int ny,  int mem_size_y, double * u_k, double * u_kp1) {
	int i, j;
	for (j=1;j<=local_nx;j++) {
		for (i=1;i<=ny;i++) {
			u_k[i+(j*mem_size_y)]=0.25 * (u_k[(i-1)+(j*mem_size_y)]+u_k[(i+1)+(j*mem_size_y)]+
				u_k[i+((j-1)*mem_size_y)]+u_k[i+((j+1)*mem_size_y)]);
		}
	}
}
