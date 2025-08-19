#ifndef CONFIG_H
#define CONFIG_H

// Boundary values 
#define LEFT_VALUE 1.0
#define RIGHT_VALUE 10.0

// The maximum number of iterations
#define MAX_ITERATIONS 1000000

// How often to report the norm
#define REPORT_NORM_PERIOD 1000

// The solver to use, 0=jacobi, 1=jacobi relaxed, 2=gauss seidel, 3=gauss seidel relaxed (SOR)
#define SOLVER_TO_USE 0

// W is the amount to overrelax in the SOR, it should be larger than 1 and smaller than 2 (but might diverge before this)
#define W 1.4
